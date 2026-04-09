//==============================================================================
//
//	Copyright (c) 2026-
//	Authors:
//	* Qiyi Tang
//	* Franck van Breugel
//	* Zainab Fatmi
//
//------------------------------------------------------------------------------
//
//	This file is part of PRISM.
//
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//==============================================================================

package explicit.bisim.distances;

import explicit.DTMC;
import explicit.rewards.Rewards;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.linearsolver.MPVariable;

import prism.PrismComponent;

/**
 * Computes the proximity for all state pairs of a labelled Markov chain by
 * means of simple policy iteration. This is the maximum probability that a
 * state pair reaches the diagonal under any policy.
 */
public class Proximity<Value> extends Distances<Value> {

    /** {@code proximity[s * numStates + t]} is the proximity of s and t. */
    private double[] proximity;

    /**
     * {@code canReachDiagonal[s * numStates + t]} is true iff (s, t) can reach
     * any (u, u) for any policy.
     */
    private boolean[] canReachDiagonal;

    public Proximity(PrismComponent parent) {
        super(parent);
    }

    /**
     * Initializes this computation of proximity for the given labelled Markov
     * chain.
     *
     * @param chain   a DTMC
     * @param propBSs the propositions used to determine labels
     * @param rewards the reward structure used to determine labels
     */
    protected void initialize(DTMC<Value> chain, List<BitSet> propBSs, Rewards<Value> rewards) {
        super.initialize(chain, propBSs, rewards);
        this.proximity = new double[this.numIndices];
        // compute the set of pairs of states that can reach a diagonal pair
        this.canReachDiagonal = new boolean[this.numIndices];
        for (int s = 0; s < this.numStates; s++) {
            this.canReachDiagonal[s * this.numStates + s] = true;
            this.proximity[s * this.numStates + s] = 1;
        }
        boolean done;
        do {
            done = true;
            for (int s = 0; s < this.numStates; s++) {
                for (int t = s + 1; t < this.numStates; t++) {
                    int index = s * this.numStates + t;
                    if (!this.canReachDiagonal[index] && !this.differentLabels[index]) {
                        checkSuccessors:
						for (int u = 0; u < this.numStates; u++) {
                            if (this.probabilities[s][u] > 0) {
                                for (int v = 0; v < this.numStates; v++) {
                                    if (this.probabilities[t][v] > 0
                                            && this.canReachDiagonal[u * this.numStates + v]) {
                                        this.canReachDiagonal[s * this.numStates + t] = true;
                                        this.toCompute[s * this.numStates + t] = true;
                                        this.canReachDiagonal[t * this.numStates + s] = true;
                                        this.toCompute[t * this.numStates + s] = true;
                                        this.setCoupling(s, t); // symmetric
                                        done = false;
                                        break checkSuccessors;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } while (!done);
    }

    /**
     * Sets the probability of reaching the diagonal based on the current policy.
     */
    private void setProximity() {
        boolean[] reachable = getReachable();
        List<Integer> activePairs = new ArrayList<>();
        // Map from flat index -> row in the linear system
        int[] indexToRow = new int[this.numIndices];
        Arrays.fill(indexToRow, -1);
        int rows = 0;
        for (int s = 0; s < this.numStates; s++) {
            for (int t = s + 1; t < this.numStates; t++) {
                int index = s * this.numStates + t;
                if (toCompute[index] && reachable[index]) {
                    activePairs.add(index);
                    indexToRow[index] = rows;
                    indexToRow[t * this.numStates + s] = rows;
                    rows++;
                }
            }
        }
        if (rows == 0) return;

        // Solve the linear system using OR-Tools
        // (LP with equality constraints and arbitrary objective)
        MPSolver solver = MPSolver.createSolver("GLOP");
        MPVariable[] vars = new MPVariable[rows];

        double[][] matrix = new double[rows][rows];
        double[] vector = new double[rows];
        for (int i = 0; i < rows; i++) {
            vars[i] = solver.makeNumVar(0.0, 1.0, "x_" + i);
            matrix[i][i] = 1.0; // identity diagonal
            int st = activePairs.get(i);
            for (int u = 0; u < this.numStates; u++) {
                for (int v = 0; v < this.numStates; v++) {
                    int uv = u * this.numStates + v;
                    if (this.canReachDiagonal[uv]) {
                        double value = this.policy[st][uv];
                        if (u == v) {
                            vector[i] += value;
                        } else if (indexToRow[uv] >= 0) {
                            matrix[i][indexToRow[uv]] -= value; // unknown pair: subtract from LHS
                        }
                    }
                }
            }
        }
        for (int i = 0; i < rows; i++) {
            // Set up the equality: matrix[i] * vars = vector[i]
            MPConstraint constraint = solver.makeConstraint(vector[i], vector[i], "eq_" + i);
            for (int j = 0; j < rows; j++) {
                if (matrix[i][j] != 0.0) {
                    constraint.setCoefficient(vars[j], matrix[i][j]);
                }
            }
        }
        solver.solve();
        for (int i = 0; i < rows; i++) {
            int st = activePairs.get(i);
            int s = st / this.numStates;
            int t = st % this.numStates;
            this.proximity[st] = vars[i].solutionValue();
            this.proximity[t * this.numStates + s] = vars[i].solutionValue();
        }
    }

    /**
     * Determine which state pairs reach the diagonal under the current policy.
     *
     * @return an array that is true at index {@code s * numStates + t} iff
     * (s, t) can reach any (u, u) under the current policy.
     */
    private boolean[] getReachable() {
        boolean[] reachable = new boolean[this.numIndices];
        for (int u = 0; u < this.numStates; u++) {
            reachable[u * this.numStates + u] = true;
        }
        boolean reachChanged;
        do {
            reachChanged = false;
            for (int s = 0; s < this.numStates; s++) {
                for (int t = s + 1; t < this.numStates; t++) {
                    int index = s * this.numStates + t;
                    if (toCompute[index] && !reachable[index]) {
                        checkReach:
                        for (int u = 0; u < this.numStates; u++) {
                            for (int v = 0; v < this.numStates; v++) {
                                if (this.policy[index][u * this.numStates + v] > 0
                                        && reachable[u * this.numStates + v]) {
                                    reachable[index] = true;
                                    reachable[t * this.numStates + s] = true; // symmetric
                                    reachChanged = true;
                                    break checkReach;
                                }
                            }
                        }
                    }
                }
            }
        } while (reachChanged);
        return reachable;
    }

    /**
     * Tests whether the coupling for the given states is optimal.
     * If not, updates the couplings for the given states.
     * 
     * @param s a state
     * @param t a state
     * @return {@code true} iff the coupling for the given states is optimal
     */
    private boolean isOptimalCoupling(int s, int t) {
        SolutionPair solution = this.getOptimalSolution(s, t, this.proximity, false);
        double value = solution.getValue();
        if (value - ACCURACY <= this.proximity[s * this.numStates + t]) {
            return true;
        } else {
            this.policy[s * this.numStates + t] = solution.getPoint();
            setSymmetricCoupling(s, t);
            return false;
        }
    }

    /**
     * Returns the proximity for the given labelled Markov chain.
     * 
     * @return the proximity for the given labelled Markov chain
     */
    public double[] compute(DTMC<Value> chain, List<BitSet> propBSs, Rewards<Value> rewards) {
        long timer = System.currentTimeMillis();
        this.initialize(chain, propBSs, rewards);
        this.setProximity();
        int iterations = 1;
        boolean allOptimal;
        do {
            allOptimal = true;
            for (int s = 0; s < this.numStates && allOptimal; s++) {
                for (int t = s + 1; t < this.numStates && allOptimal; t++) {
                    if (this.toCompute[s * this.numStates + t] && !this.isOptimalCoupling(s, t)) {
                        allOptimal = false;
                        this.setProximity();
                    }
                }
            }
            iterations++;
        } while (!allOptimal);
        timer = System.currentTimeMillis() - timer;
        mainLog.println("Iterations: " + iterations);
        mainLog.println("Time for proximity computation: " + timer / 1000.0 + " seconds.");
        printArray(this.proximity);
        return this.proximity;
    }
}
