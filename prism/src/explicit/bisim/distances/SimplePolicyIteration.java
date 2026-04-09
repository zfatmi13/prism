//==============================================================================
//
//	Copyright (c) 2025-
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
import java.util.BitSet;
import java.util.List;

import prism.PrismComponent;

/**
 * Computes the probabilistic bisimilarity distances by
 * means of simple policy iteration.
 */
public class SimplePolicyIteration<Value> extends Distances<Value> {

	/** {@code distance[s * numStates + t]} is the distance between s and t. */
	protected double[] distance;

	/**
	 * {@code bisimilar[s * numStates + t]} is true iff states s and t are
	 * bisimilar.
	 */
	protected boolean[] bisimilar;

	public SimplePolicyIteration(PrismComponent parent) {
		super(parent);
	}

	/**
	 * Initializes this computation of probabilistic bisimilarity
	 * distances for the given labelled Markov chain.
	 *
	 * @param chain   a labelled Markov chain
	 * @param propBSs the propositions used to determine labels
	 * @param rewards the reward structure used to determine labels
	 */
	protected void initialize(DTMC<Value> chain, List<BitSet> propBSs, Rewards<Value> rewards) {
		super.initialize(chain, propBSs, rewards);
		this.distance = new double[this.numIndices];
		this.bisimilar = new boolean[this.numIndices];
		// compute bisimilarity
		bisimilarity(chain, chain.getEvaluator(), 1);
		for (int s = 0; s < this.numStates; s++) {
			this.bisimilar[s * this.numStates + s] = true;
			for (int t = s + 1; t < this.numStates; t++) {
				if (this.partition[s] == this.partition[t]) {
					this.bisimilar[s * this.numStates + t] = true;
					this.bisimilar[t * this.numStates + s] = true;
				} else if (this.differentLabels[s * this.numStates + t]) {
					this.distance[s * this.numStates + t] = 1;
					this.distance[t * this.numStates + s] = 1;
				} else {
					this.toCompute[s * this.numStates + t] = true;
					this.toCompute[t * this.numStates + s] = true;
					this.setCoupling(s, t);
				}
			}
		}
	}

	/**
	 * Returns the probabilistic bisimilarity distances.
	 *
	 * @return the probabilistic bisimilarity distances.
	 */
	public double[] getDistance() {
		return this.distance;
	}

	/**
	 * Sets the distances based on the current couplings.
	 */
	private void setDistance() {
		double[] lower = new double[this.numIndices];
		double[] upper = new double[this.numIndices];

		for (int u = 0; u < this.numStates; u++) {
			for (int v = 0; v < this.numStates; v++) {
				if (this.bisimilar[u * this.numStates + v]) {
					// lower[u * this.numStates + v] = 0;
					// upper[u * this.numStates + v] = 0;
				} else if (differentLabels[u * this.numStates + v]) {
					lower[u * this.numStates + v] = 1;
					upper[u * this.numStates + v] = 1;
				} else {
					// lower[u * this.numStates + v] = 0;
					upper[u * this.numStates + v] = 1;
				}
			}
		}

		boolean[] farApart = Arrays.copyOf(this.toCompute, this.toCompute.length);
		boolean done;
		do {
			done = true;
			double[] previousLower = Arrays.copyOf(lower, lower.length);
			double[] previousUpper = Arrays.copyOf(upper, upper.length);
			for (int s = 0; s < this.numStates; s++) {
				for (int t = 0; t < this.numStates; t++) {
					if (farApart[s * this.numStates + t]) {
						lower[s * this.numStates + t] = 0;
						upper[s * this.numStates + t] = 0;
						for (int u = 0; u < this.numStates; u++) {
							for (int v = 0; v < this.numStates; v++) {
								if (this.bisimilar[u * this.numStates + v]) {
									// do nothing
								} else {
									double value = this.policy[s * this.numStates + t][u * this.numStates + v];
									if (differentLabels[u * this.numStates + v]) {
										lower[s * this.numStates + t] += value;
										upper[s * this.numStates + t] += value;
									} else {
										lower[s * this.numStates + t] += previousLower[u * this.numStates + v] * value;
										upper[s * this.numStates + t] += previousUpper[u * this.numStates + v] * value;
									}
								}
							}
						}

						if (upper[s * this.numStates + t] - lower[s * this.numStates + t] < ACCURACY) {
							farApart[s * this.numStates + t] = false;
						} else {
							done = false;
						}
					}
				}
			}
		} while (!done);

		for (int s = 0; s < this.numStates; s++) {
			for (int t = 0; t < this.numStates; t++) {
				if (this.toCompute[s * this.numStates + t]) {
					this.distance[s * this.numStates + t] =
							(lower[s * this.numStates + t] + upper[s * this.numStates + t]) / 2;
				}
			}
		}
	}

	/**
	 * Tests whether the coupling for the given states is optimal.
	 * If not, updates the coupling for the given states.
	 *
	 * @param s a state
	 * @param t a state
	 * @pre. s > t
	 * @return true if the coupling for the given states is optimal,
	 *         false otherwise
	 */
	private boolean isOptimalCoupling(int s, int t) {
		SolutionPair solution = this.getOptimalSolution(s, t, this.distance);
		double value = solution.getValue();
		if (value + ACCURACY >= this.distance[s * this.numStates + t]) {
			return true;
		} else {
			this.policy[s * this.numStates + t] = solution.getPoint();
			setSymmetricCoupling(s, t);
			return false;
		}
	}

	/**
	 * Returns the probabilistic bisimilarity distances for the
	 * given labelled Markov chain.
	 */
	public double[] compute(DTMC<Value> chain, List<BitSet> propBSs, Rewards<Value> rewards) {
		long timer = System.currentTimeMillis();
		this.initialize(chain, propBSs, rewards);
		this.setDistance();
		int iterations = 1;
		boolean allOptimal;
		do {
			allOptimal = true;
			for (int s = 0; s < this.numStates && allOptimal; s++) {
				for (int t = s + 1; t < this.numStates && allOptimal; t++) {
					if (this.toCompute[s * this.numStates + t] && !this.isOptimalCoupling(s, t)) {
						allOptimal = false;
						this.setDistance();
					}
				}
			}
			iterations++;
		} while (!allOptimal);
		timer = System.currentTimeMillis() - timer;
		mainLog.println("Iterations: " + iterations);
		mainLog.println("Time for distance computation: " + timer / 1000.0 + " seconds.");
		printArray(this.distance);
		return this.distance;
	}
}
