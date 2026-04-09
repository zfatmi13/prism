//==============================================================================
//
//	Copyright (c) 2026-
//	Authors:
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

import java.util.BitSet;
import java.util.List;

import prism.PrismComponent;

/**
 * Decides the continuity of the probabilistic bisimilarity distance function.
 */
public class Continuity<Value> extends SimplePolicyIteration<Value> {

    /** Stores pairs of states for which the distance is discontinuous */
    private boolean[] discontinuous;

    public Continuity(PrismComponent parent) {
        super(parent);
    }

    /**
     * Decides the continuity of the distance function for all pairs of states.
     *
     * @param chain   a labelled Markov chain
     * @param propBSs the propositions used to determine labels
     * @param rewards the reward structure used to determine labels
     * @return an array indicating pairs of states for which the distance is
     *         discontinuous
     */
    public boolean[] decide(DTMC<Value> chain, List<BitSet> propBSs, Rewards<Value> rewards) {
        // compute the distances
        compute(chain, propBSs, rewards);
        long timer = System.currentTimeMillis();
        // compute robustly bisimilar pairs of states
        robust(chain, chain.getEvaluator(), 1);
        // initialize discontinuous pairs of states to those that are
        // bisimilar but not robustly bisimilar
        this.discontinuous = new boolean[this.numIndices];
        for (int s = 0; s < this.numStates; s++) {
            for (int t = s + 1; t < this.numStates; t++) {
                if (this.bisimilar[s * this.numStates + t] && this.partition[s] != this.partition[t]) {
                    this.discontinuous[s * this.numStates + t] = true;
                    this.discontinuous[t * this.numStates + s] = true; // symmetric
                }
            }
        }
        // decide continuity
        boolean changed;
        do {
            changed = false;
            for (int s = 0; s < this.numStates; s++) {
                for (int t = s + 1; t < this.numStates; t++) {
                    int index = s * this.numStates + t;
                    if (!discontinuous[index] && !bisimilar[index] && !differentLabels[index] &&
                            isPredecessor(s, t, discontinuous)) {
                        if (isDiscontinuous(s, t)) {
                            discontinuous[index] = true;
                            discontinuous[t * this.numStates + s] = true; // symmetric
                            changed = true;
                        }
                    }
                }
            }
        } while (changed);
        timer = System.currentTimeMillis() - timer;
        mainLog.println("Time to decide continuity: " + timer / 1000.0 + " seconds.");
        printComplementArray(this.discontinuous);
        return this.discontinuous;
    }

    /**
     * Returns {@code true} iff every optimal coupling of s and t has
     * discontinuous pairs in its support, i.e. no there is no optimal
     * coupling whose support lies entirely outside discontinuous pairs.
     *
     * @param s a state
     * @param t a state
     * @return {@code true} iff all optimal couplings have support
     *         intersecting known discontinuous pairs
     */
    private boolean isDiscontinuous(int s, int t) {
        SolutionPair restrictedSolution = getOptimalSolution(s, t, this.distance, true, this.discontinuous);
        return restrictedSolution == null ||
                !((Math.abs(restrictedSolution.getValue() - this.distance[s * this.numStates + t])) < ACCURACY);
    }
}
