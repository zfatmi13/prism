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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import explicit.DTMC;
import prism.Evaluator;
import prism.PrismComponent;
import prism.PrismException;

/**
 * Performs robust bisimulation minimisation for explicit-state models.
 */
public class Robust1<Value> extends BisimulationMethodNew<Value> {
    // Local storage of partition info
    Map<Integer, List<Integer>> successors;
    List<List<Integer>> R;
    List<Map<Integer, Boolean>> Q;

    /**
     * Construct a new RobustBisimulation object.
     */
    public Robust1(PrismComponent parent) throws PrismException {
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
        // printPartition("before");
        while (filter()) {
            // printStuff("after refine");
            prune();
            bisimilarity(dtmc, eval, 1);
            // printPartition("after stabilize");
        }
        if (numStates == numBlocks) {
            return false;
        }
        lifting(dtmc, eval);
        // printProbabilities(dtmc, eval);
        return true;
    }

    /* computes the list of successors */
    protected void initialize(DTMC<Value> dtmc) {
        successors = new HashMap<Integer, List<Integer>>(numStates);
        List<Integer> transitions;
        for (int source = 0; source < numStates; source++) {
            transitions = new ArrayList<Integer>();
            Iterator<Entry<Integer, Value>> iter = dtmc.getTransitionsIterator(source);
            while (iter.hasNext()) {
                transitions.add(iter.next().getKey());
            }
            successors.put(source, transitions);
        }
    }

    protected boolean filter() {
        /* construct R from the partition */
        R = new ArrayList<List<Integer>>();
        for (int block = 0; block < numBlocks; block++) {
            R.add(new ArrayList<Integer>());
        }
        int sizeR = 0;
        for (int state = 0; state < numStates; state++) {
            R.get(partition[state]).add(state);
            sizeR += R.get(partition[state]).size();
        }
        /* initialize Q to the diagonal pairs */
        Q = new ArrayList<Map<Integer, Boolean>>(numStates);
        int sizeQ = 0;
        for (int state = 0; state < numStates; state++) {
            Q.add(new HashMap<Integer, Boolean>());
            Q.get(state).put(state, Boolean.TRUE);
            sizeQ++;
        }
        // printStuff("before refine");
        boolean changedQ;
        do {
            changedQ = false;
            for (int block = 0; block < numBlocks; block++) {
                for (int i = 0; i < R.get(block).size(); i++) {
                    int s = R.get(block).get(i);
                    for (int j = i + 1; j < R.get(block).size(); j++) {
                        int t = R.get(block).get(j);
                        if (!Q.get(s).containsKey(t) && preQ(s, t)) {
                            Q.get(t).put(s, Boolean.TRUE);
                            Q.get(s).put(t, Boolean.TRUE);
                            sizeQ++;
                            changedQ = true;
                        }
                    }
                }
            }
        } while (changedQ);
        return sizeQ != sizeR; /* whether R has changed */
    }

    protected boolean preQ(int s, int t) {
        for (int successorS : successors.get(s)) {
            for (int successorT : successors.get(t)) {
                if (Q.get(successorS).containsKey(successorT)) {
                    return true;
                }
            }
        }
        return false;
    }

    protected void prune() {
        for (int s = 0; s < numStates; s++) {
            for (Map.Entry<Integer, Boolean> t : Q.get(s).entrySet()) {
                /* skip the diagonal pair */
                if (s != t.getKey() && markRemoved(s, t.getKey())) {
                    t.setValue(Boolean.FALSE);
                }
            }
        }
        /* construct the new partition */
        numBlocks = 0;
        Arrays.fill(partition, -1);
        for (int state = 0; state < numStates; state++) {
            if (partition[state] == -1) {
                for (Map.Entry<Integer, Boolean> entry : Q.get(state).entrySet()) {
                    if (entry.getValue()) {
                        partition[entry.getKey()] = numBlocks;
                    }
                }
                numBlocks++;
            } /* else {
                // test stuff
                int p = partition[state];
                for (Map.Entry<Integer, Boolean> entry : Q.get(state).entrySet()) {
                    if (entry.getValue()) {
                        if (p != partition[entry.getKey()]) {
                            mainLog.println("ERROR: not an equivalence?");
                        }
                    }
                }
                for (int index = 0; index < numStates; index++) {
                    if (p == partition[index]) {
                        if (!Q.get(state).containsKey(index) || !Q.get(state).get(index)) {
                            mainLog.println("ERROR: not an equivalence?");
                        }
                    }
                }
            } */
        }
        // printPartition("after prune");
    }

    protected boolean markRemoved(int s, int t) {
        boolean removeT = false;
        for (Map.Entry<Integer, Boolean> u : Q.get(t).entrySet()) {
            if (!Q.get(s).containsKey(u.getKey())) {
                u.setValue(Boolean.FALSE);
                removeT = true;
            }
        }
        return removeT;
    }

    private void printStuff(String where) {
        System.out.println(where);
        System.out.println("R");
        for (int s = 0; s < numBlocks; s++) {
            System.out.print("s = " + s + ": ");
            for (int t : R.get(s)) {
                System.out.print(t + " ");
            }
            System.out.println();
        }
        if (Q != null) {
            System.out.println("Q");
            for (int s = 0; s < numStates; s++) {
                System.out.print("s = " + s + ": ");
                for (int t :Q.get(s).keySet()) {
                    System.out.print(t + " ");
                }
                System.out.println();
            }
        }
    }

    private void printPartition(String where) {
        System.out.println(where);
        System.out.println("partition");
        for (int s = 0; s < numStates; s++) {
            System.out.println("s" + s + " in " + partition[s]);
        }
    }
}
