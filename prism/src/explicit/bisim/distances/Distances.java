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
import explicit.bisim.RobustBisimulation;
import explicit.rewards.Rewards;

import prism.PrismComponent;

import java.util.Arrays;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import com.google.ortools.Loader;
import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPObjective;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.linearsolver.MPVariable;

/**
 * Parent class for bisimilarity distance computations.
 */
public abstract class Distances<Value> extends RobustBisimulation<Value> {

    // Load the native libraries for Google's OR-Tools.
    static {
        Loader.loadNativeLibraries();
    }

    /** Desired numerical accuracy for comparisons. */
    protected static final double ACCURACY = 1E-6;

    /** Transition probabilities of the labelled Markov chain. */
    protected double[][] probabilities;

    /** n^2, where n is the number of states in the labelled Markov chain. */
    protected int numIndices;

    /**
     * {@code toCompute[s * numStates + t]} stores whether the distance of
     * s and t still needs to be computed.
     */
    protected boolean[] toCompute;

    /**
     * For states s, t, with s > t, if {@code toCompute[s * numStates + t]}
     * then {@code policy[s * numStates + t]} is a coupling of
     * {@code probabilities[s]} and {@code probabilities[t]}.
     */
    protected double[][] policy;

    /**
     * {@code distanceOne[s * numStates + t]} is true iff states s and t
     * have distance one.
     */
    protected boolean[] distanceOne;

    /**
     * Represents a solution returned by the LP solver.
     */
    public static class SolutionPair {
        private final double[] point;
        private final double value;

        public SolutionPair(double[] point, double value) {
            this.point = point;
            this.value = value;
        }

        public double[] getPoint() {
            return point;
        }

        public double getValue() {
            return value;
        }
    }

    public Distances(PrismComponent parent) {
        super(parent);
    }

    /**
     * Initializes the properties of the given labelled Markov chain.
     *
     * @param chain a DTMC
     */
    protected void initialize(DTMC<Value> chain, List<BitSet> propBSs, Rewards<Value> rewards) {
        this.numStates = chain.getNumStates();
        this.probabilities = new double[this.numStates][this.numStates];
        for (int s = 0; s < this.numStates; s++) {
            Iterator<Map.Entry<Integer, Value>> iter = chain.getTransitionsIterator(s);
            while (iter.hasNext()) {
                Map.Entry<Integer, Value> transition = iter.next();
                this.probabilities[s][transition.getKey()] = chain.getEvaluator().toDouble(transition.getValue());
            }
        }
        this.numIndices = this.numStates * this.numStates;
        this.toCompute = new boolean[this.numIndices];
        this.policy = new double[this.numIndices][this.numIndices];
        // use initial partition to determine if states have the same label
        this.distanceOne = new boolean[this.numIndices];
        initialisePartitionInfo(chain, propBSs, rewards);
        for (int s = 0; s < this.numStates; s++) {
            for (int t = s + 1; t < this.numStates; t++) {
                if (this.partition[s] != this.partition[t]) {
                    this.distanceOne[s * this.numStates + t] = true;
                    this.distanceOne[t * this.numStates + s] = true; // symmetric
                }
            }
        }
    }

    /**
     * Determines whether (s, t) is a predecessor of R, that is, there exist
     * states u and v such that u is a successor of s, v is a successor of t,
     * and (u, v) is an element of R (i.e. {@code R[u * numStates + v] == true}).
     *
     * @param s a state
     * @param t a state
     * @param R a set of pairs of states, represented as an array of size numIndices
     * @return {@code true} iff (s, t) is a predecessor of some pair in R.
     */
    protected boolean isPredecessor(int s, int t, boolean[] R) {
        for (int u = 0; u < this.numStates; u++) {
            if (this.probabilities[s][u] > 0) {
                for (int v = 0; v < this.numStates; v++) {
                    if (this.probabilities[t][v] > 0 && R[u * this.numStates + v]) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    /**
     * Sets a coupling for the probability distributions of the given states in
     * {@code this.policy[s * numStates + t]} and {@code this.policy[t * numStates + s]}.
     * The coupling is computed using the North-West corner method.
     *
     * @param s a state
     * @param t a state
     */
    protected void setCoupling(int s, int t) {
        double[] source = Arrays.copyOf(this.probabilities[s], this.numStates);
        double[] target = Arrays.copyOf(this.probabilities[t], this.numStates);
        for (int u = 0; u < this.numStates; u++) {
            for (int v = 0; v < this.numStates; v++) {
                double minimum = Math.min(source[u], target[v]);
                this.policy[s * this.numStates + t][u * this.numStates + v] = minimum;
                this.policy[t * this.numStates + s][v * this.numStates + u] = minimum; // symmetric
                source[u] -= minimum;
                target[v] -= minimum;
            }
        }
    }

    /**
     * Sets the symmetric coupling for the probability distributions of the
     * given states in {@code this.policy[t * numStates + s]} from
     * {@code this.policy[s * numStates + t]}.
     *
     * @param s a state
     * @param t a state
     */
    protected void setSymmetricCoupling(int s, int t) {
        for (int u = 0; u < this.numStates; u++) {
            for (int v = 0; v < this.numStates; v++) {
                this.policy[t * this.numStates + s][v * this.numStates + u] =
                        this.policy[s * this.numStates + t][u * this.numStates + v];
            }
        }
    }

    /**
     * Solves the LP for an optimal coupling of the transition distributions of
     * the given states s and t, minimising the expected distance.
     *
     * @param s        a state
     * @param t        a state
     * @param distance the distance function
     * @return an optimal solution
     */
    protected SolutionPair getOptimalSolution(int s, int t, double[] distance) {
        return getOptimalSolution(s, t, distance, true, null);
    }

    /**
     * Solves the LP for an optimal coupling of the transition distributions of
     * the given states s and t.
     *
     * @param s        a state
     * @param t        a state
     * @param distance the distance function
     * @param minimize whether to minimise the expected distance
     * @return an optimal solution
     */
    protected SolutionPair getOptimalSolution(int s, int t, double[] distance, boolean minimize) {
        return getOptimalSolution(s, t, distance, minimize, null);
    }

    /**
     * Solves the LP for an optimal coupling of the transition distributions of
     * the given states s and t.
     * If {@code avoid} is non-null, additional equality constraints are added
     * to restrict the coupling to pairs of states outside {@code avoid}.
     * Returns {@code null} if the restricted LP is infeasible.
     *
     * @param s        a state
     * @param t        a state
     * @param distance the distance function
     * @param minimize whether to minimise the expected distance
     * @param avoid    a boolean array of size n^2 indicating pairs of
     *                 states on which the coupling mass must be zero
     * @return an optimal solution or {@code null} if infeasible
     */
    protected SolutionPair getOptimalSolution(int s, int t, double[] distance, boolean minimize, boolean[] avoid) {
        MPSolver solver = MPSolver.createSolver("GLOP");
        MPVariable[] vars = new MPVariable[this.numIndices];
        for (int i = 0; i < this.numIndices; i++) {
            vars[i] = solver.makeNumVar(0.0, 1.0, "x_" + i);
        }
        // objective function
        MPObjective objective = solver.objective();
        for (int i = 0; i < this.numIndices; i++) {
            objective.setCoefficient(vars[i], distance[i]);
        }
        if (minimize) {
            objective.setMinimization();
        } else {
            objective.setMaximization();
        }
        // marginal constraints
        for (int u = 0; u < this.numStates; u++) {
            MPConstraint sConstraint = solver.makeConstraint(this.probabilities[s][u],
                    this.probabilities[s][u], "s_marginal_" + u);
            MPConstraint tConstraint = solver.makeConstraint(this.probabilities[t][u],
                    this.probabilities[t][u], "t_marginal_" + u);
            for (int v = 0; v < this.numStates; v++) {
                sConstraint.setCoefficient(vars[u * this.numStates + v], 1.0);
                tConstraint.setCoefficient(vars[v * this.numStates + u], 1.0);
            }
        }
        // forbidden pairs
        if (avoid != null) {
            MPConstraint blockConstraint = solver.makeConstraint(0, 0, "avoid");
            for (int u = 0; u < this.numStates; u++) {
                for (int v = u + 1; v < this.numStates; v++) {
                    if (avoid[u * this.numStates + v]) {
                        blockConstraint.setCoefficient(vars[u * this.numStates + v], 1.0);
                        blockConstraint.setCoefficient(vars[v * this.numStates + u], 1.0);
                    }
                }
            }
        }
        // attempt to solve the linear program
        final MPSolver.ResultStatus resultStatus = solver.solve();
        if (resultStatus == MPSolver.ResultStatus.OPTIMAL || resultStatus == MPSolver.ResultStatus.FEASIBLE) {
            double[] point = new double[this.numIndices];
            for (int i = 0; i < this.numIndices; i++) {
                point[i] = vars[i].solutionValue();
            }
            return new SolutionPair(point, objective.value());
        } else {
            return null;
        }
    }

    /**
     * Prints the given array in a matrix format.
     * 
     * @param array the array of size n^2 to print
     */
    protected void printArray(double[] array) {
        for (int s = 0; s < this.numStates; s++) {
            for (int t = 0; t < this.numStates; t++) {
                mainLog.print(String.format("%.3f\t", array[s * this.numStates + t]));
            }
            mainLog.println();
        }
    }

    /**
     * Prints the complement of given array in a matrix format.
     *
     * @param array a boolean array of size n^2
     */
    protected void printComplementArray(boolean[] array) {
        for (int s = 0; s < this.numStates; s++) {
            for (int t = 0; t < this.numStates; t++) {
                mainLog.print(!array[s * this.numStates + t] + "\t");
            }
            mainLog.println();
        }
    }
}
