//==============================================================================
//	
//	Copyright (c) 2020-
//	Authors:
//	* Zainab Fatmi
//	* Franck van Breugel
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

package explicit;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import prism.Evaluator;
import prism.PrismComponent;
import prism.PrismException;

/**
 * Method to perform bisimulation minimisation for explicit-state models.
 */
public class BisimulationMethodNew<Value> extends Bisimulation<Value>
{
	/**
	 * The signatures of each block in the final partition.
	 */
	private HashMap<Integer, Distribution<Value>> probabilities;

	/**
	 * Construct a new BisimulationMethodNew object.
	 */
	public BisimulationMethodNew(PrismComponent parent) throws PrismException
	{
		super(parent);
	}

	@Override
	protected boolean minimiseDTMC(DTMC<Value> dtmc)
	{
		Evaluator<Value> eval = dtmc.getEvaluator();

		/* the blocks used as splitters, exclude the last block which contains all states without a label */
		BitSet splitters = new BitSet(numBlocks);
		splitters.flip(0, numBlocks - 1);

		int numberOfBlocksOld;
		int[] partitionOld;
		BitSet splittersOld;
		do {
			numberOfBlocksOld = numBlocks;
			numBlocks = 0;
			partitionOld = partition;
			partition = new int[numStates];
			splittersOld = splitters;
			splitters = new BitSet(numberOfBlocksOld);

			/* blocks.get(b) is the list of sub-blocks that partition block b */
			List<List<SubBlock<Value>>> blocks = new ArrayList<List<SubBlock<Value>>>(numberOfBlocksOld);
			for (int b = 0; b < numberOfBlocksOld; b++) {
				List<SubBlock<Value>> empty = new ArrayList<SubBlock<Value>>();
				blocks.add(empty);
			}

			for (int source = 0; source < numStates; source++) {
				SubBlock<Value> subBlock = null; /* sub-block that contains the source state */
				Iterator<Entry<Integer, Value>> iter = dtmc.getTransitionsIterator(source);
				while (iter.hasNext()) {
					Entry<Integer, Value> transition = iter.next();
					int blockOfTarget = partitionOld[transition.getKey()];
					Value probabilityOld;
					if (splittersOld.get(blockOfTarget)) {
						if (subBlock == null) {
							subBlock = new NonEmptySubBlock<Value>(eval);
							probabilityOld = eval.zero();
						} else {
							probabilityOld = subBlock.get(blockOfTarget);
							if (probabilityOld == null) {
								probabilityOld = eval.zero();
							}
						}
						subBlock.put(blockOfTarget, eval.add(probabilityOld, transition.getValue()));
					}
				}
				if (subBlock == null) {
					subBlock = new EmptySubBlock<Value>();
				}

				int blockOfSource = partitionOld[source];
				List<SubBlock<Value>> subBlocks = blocks.get(blockOfSource);
				int index = subBlocks.indexOf(subBlock);
				if (index == -1) { /* there is no sub-block with the same lifting as the source state */
					subBlock.setID(numBlocks);
					subBlocks.add(subBlock);
					partition[source] = numBlocks;
					/* the first sub-block is not used as a splitter in the next refinement step */
					if (subBlocks.size() > 1) {
						splitters.set(numBlocks);
					}
					numBlocks++;
				} else {
					partition[source] = subBlocks.get(index).getID();
				}
			}
		} while (numBlocks != numberOfBlocksOld);

		if (numStates == numBlocks) {
			return false;
		}
		/* probabilities.get(s) is the final lifting of the probability distribution of states in block s */
		probabilities = new HashMap<Integer, Distribution<Value>>(numBlocks);
		for (int source = 0; source < numStates; source++) {
			int blockOfSource = partition[source];
			if (!probabilities.containsKey(blockOfSource)) {
				Distribution<Value> distrNew = new Distribution<Value>(eval);
				Iterator<Entry<Integer, Value>> iter = dtmc.getTransitionsIterator(source);
				while (iter.hasNext()) {
					Entry<Integer, Value> transition = iter.next();
					distrNew.add(partition[transition.getKey()], transition.getValue());
				}
				probabilities.put(blockOfSource, distrNew);
			}
		}
		return true;
	}

	@Override
	protected DTMCSimple<Value> buildReducedDTMC()
	{
		DTMCSimple<Value> dtmcNew = new DTMCSimple<>(numBlocks);
		for (int i = 0; i < numBlocks; i++) {
			for (Map.Entry<Integer, Value> e : probabilities.get(i)) {
				dtmcNew.setProbability(i, e.getKey(), e.getValue());
			}
		}
		return dtmcNew;
	}

	@Override
	protected CTMCSimple<Value> buildReducedCTMC()
	{
		CTMCSimple<Value> ctmcNew = new CTMCSimple<>(numBlocks);
		for (int i = 0; i < numBlocks; i++) {
			for (Map.Entry<Integer, Value> e : probabilities.get(i)) {
				ctmcNew.setProbability(i, e.getKey(), e.getValue());
			}
		}
		return ctmcNew;
	}
}
