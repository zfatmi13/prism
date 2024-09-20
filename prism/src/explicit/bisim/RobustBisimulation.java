//==============================================================================
//
//	Copyright (c) 2024-
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

package explicit.bisim;

import java.util.*;
import java.util.Map.Entry;

import explicit.DTMC;
import prism.Evaluator;
import prism.PrismComponent;
import prism.PrismException;

/**
 * Performs robust bisimulation minimisation for explicit-state models.
 */
public class RobustBisimulation<Value> extends BisimulationMethodNew<Value> {
	// Local storage of partition info
	ArrayList<Map<Integer, Value>> transitions;
	List<Set<Integer>> R;
	List<Set<Integer>> Q;
	boolean changedR;
	boolean changedQ;
	boolean firstIteration = true;

	/**
	 * Construct a new RobustBisimulation object.
	 */
	public RobustBisimulation(PrismComponent parent) throws PrismException {
		super(parent);
	}

	@Override
	protected boolean minimiseDTMC(DTMC<Value> dtmc) {
		Evaluator<Value> eval = dtmc.getEvaluator();
		bisimilarity(dtmc, eval, 1);
		if (numStates == numBlocks) {
			return false;
		}
		initialize(dtmc);
		// printStuff("after init");
		do {
			// mainLog.println("main iteration");
			refine(eval);
			// printStuff("after refine");
			changedR = false;
			for (int s = 0; s < numStates; s++) {
				if (!R.get(s).isEmpty()) {
					changedR = true;
					break;
				}
			}
			R = Q;
			for (int s = 0; s < numStates; s++) {
				R.get(s).remove(s);
			} // make it a set and check if not in Q
		} while (changedR);
		int block = 0; // compute numBlocks
		Arrays.fill(partition, -1);
		for (int s = 0; s < numStates; s++) {
			if (partition[s] == -1) {
				partition[s] = block;
				for (int t : R.get(s)) {
					partition[t] = block;
				}
				block++;
			}
			// System.out.print(partition[s] + " ");
		}
		// System.out.println("blocks = " + block);
		numBlocks = block;
		if (numStates == numBlocks) {
			return false;
		}
		lifting(dtmc, eval);
		// printProbabilities(dtmc, eval);
		return true;
	}

	protected void initialize(DTMC<Value> dtmc) {
		/* transitions.get(s).get(t) is the probability of transitioning from state s to state t */
		transitions = new ArrayList<Map<Integer, Value>>(numStates);
		/* construct R = bisimilar pairs */
		R = new ArrayList<Set<Integer>>(numStates);
		for (int source = 0; source < numStates; source++) {
			R.add(new HashSet<Integer>());
		}
		for (int source = 0; source < numStates; source++) {
			Map<Integer, Value> successors = new HashMap<Integer, Value>();
			Iterator<Entry<Integer, Value>> iter = dtmc.getTransitionsIterator(source);
			while (iter.hasNext()) {
				Entry<Integer, Value> transition = iter.next();
				successors.put(transition.getKey(), transition.getValue());
			}
			transitions.add(successors);

			if (R.get(source).isEmpty()) {
				for (int target = source + 1; target < numStates; target++) {
					if (partition[target] == partition[source]) {
						for (Integer state : R.get(source)) {
							R.get(state).add(target);
						}
						R.get(source).add(target);
					}
				}
			}
		}
	}

	protected void refine(Evaluator<Value> eval) {
		/* construct Q = diagonal pairs */
		Q = new ArrayList<Set<Integer>>(numStates);
		for (int source = 0; source < numStates; source++) {
			Q.add(new HashSet<Integer>());
			Q.get(source).add(source);
		}
		// printStuff("before refine");
		do {
			// mainLog.println("  refine iteration");
			changedQ = false;
			for (int s = 0; s < numStates; s++) {
				for (Iterator<Integer> it = R.get(s).iterator(); it.hasNext(); ) {
					int t = it.next();
					if (firstIteration) {
						if (preQ(s, t)) {
							Q.get(s).add(t);
							it.remove();
							changedQ = true;
						}
					} else if (coupling(s, t, eval)) {
						Q.get(s).add(t);
						it.remove();
						changedQ = true;
						// printStuff("after change");
					}
				}
			}
		} while (changedQ);
		firstIteration = false;
	}

	protected boolean preQ(int s, int t) {
		for (Entry<Integer, Value> successorS : transitions.get(s).entrySet()) {
			for (Entry<Integer, Value> successorT : transitions.get(t).entrySet()) {
				if (Q.get(successorS.getKey()).contains(successorT.getKey()) || Q.get(successorT.getKey()).contains(successorS.getKey())) {
					return true;
				}
			}
		}
		// do for pairs in Q, if transitions.get(s).contains() && transitions.get(t).contains() || vice versa
		return false;
	}

	/**
	 * Checks whether there exists a coupling ω ∈ ΩR(τ(s),τ(t)) such that
	 * support(ω) ⊆ R and support(ω) ∩ Q  ̸= ∅.
	 *
	 * @param s a state.
	 * @param t another state.
	 * @return true if such a coupling exists, false otherwise.
	 */
	private boolean coupling(int s, int t, Evaluator<Value> eval) {
		boolean preQ = false;
		int n = (2 * numStates) + 2;
		int source = n - 1;
		int target = n - 2;
		MaximumFlow<Value> dinic = new CycleCancelling<Value>(n, source, target, eval);
		for (Entry<Integer, Value> successorS : transitions.get(s).entrySet()) {
			dinic.addEdge(source, successorS.getKey(), successorS.getValue());
			for (Entry<Integer, Value> successorT : transitions.get(t).entrySet()) {
				if (Q.get(successorS.getKey()).contains(successorT.getKey()) || Q.get(successorT.getKey()).contains(successorS.getKey())) {
					preQ = true;
					// System.out.println("preQ " + successorS.getKey() + " " + successorT.getKey());
					dinic.addEdge(successorS.getKey(), successorT.getKey() + numStates, eval.one());
				} else if (R.get(successorS.getKey()).contains(successorT.getKey()) || R.get(successorT.getKey()).contains(successorS.getKey())) {
					dinic.addEdge(successorS.getKey(), successorT.getKey() + numStates, eval.one(), eval.one());
				}
			}
		}
		if (!preQ) {
			return false;
		}
		for (Entry<Integer, Value> successorT : transitions.get(t).entrySet()) {
			dinic.addEdge(successorT.getKey() + numStates, target, successorT.getValue());
		}
		// System.out.println("coupling " + s + " " + t + " " + eval.toDouble(dinic.getMaxFlow()));
		return eval.isOne(dinic.getMaxFlow()) && eval.gt(eval.one(), dinic.getMinCost()); // or !eval.isOne
	}

	private void printStuff(String where) {
		System.out.println(where);
		System.out.println("R");
		for (int s = 0; s < numStates; s++) {
			System.out.print("s = " + s + ": ");
			for (int t :R.get(s)) {
				System.out.print(t + " ");
			}
			System.out.println();
		}
		if (Q != null) {
			System.out.println("Q");
			for (int s = 0; s < numStates; s++) {
				System.out.print("s = " + s + ": ");
				for (int t :Q.get(s)) {
					System.out.print(t + " ");
				}
				System.out.println();
			}
		}
	}
}
