//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
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
import java.util.BitSet;
import java.util.List;

import explicit.CTMC;
import explicit.CTMCSimple;
import explicit.DTMC;
import explicit.DTMCSimple;
import explicit.Model;
import explicit.ModelExplicit;
import parser.State;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismNotSupportedException;

/**
 * Class to perform bisimulation minimisation for explicit-state models.
 */
public abstract class Bisimulation<Value> extends PrismComponent
{
	// Local storage of partition info
	protected int numStates;
	protected int[] partition;
	protected int numBlocks;
	protected boolean minimised;

	/**
	 * Construct a new Bisimulation object.
	 */
	public Bisimulation(PrismComponent parent) throws PrismException
	{
		super(parent);
	}

	/**
	 * Perform bisimulation minimisation on a model.
	 * @param model     The model
	 * @param propNames Names of the propositions in {@code propBSs}
	 * @param propBSs   Propositions (satisfying sets of states) to be preserved by bisimulation.
	 */
	public Model<Value> minimise(Model<Value> model, List<String> propNames, List<BitSet> propBSs) throws PrismException
	{
		switch (model.getModelType()) {
			case DTMC:
				return minimiseDTMC((DTMC<Value>) model, propNames, propBSs);
			case CTMC:
				return minimiseCTMC((CTMC<Value>) model, propNames, propBSs);
			default:
				throw new PrismNotSupportedException(
						"Bisimulation minimisation not yet supported for " + model.getModelType() + "s");
		}
	}

	/**
	 * Perform bisimulation minimisation on a DTMC.
	 * @param dtmc      The DTMC
	 * @param propNames Names of the propositions in {@code propBSs}
	 * @param propBSs   Propositions (satisfying sets of states) to be preserved by bisimulation.
	 */
	protected DTMC<Value> minimiseDTMC(DTMC<Value> dtmc, List<String> propNames, List<BitSet> propBSs)
	{
		long timer;
		timer = System.currentTimeMillis();
		initialisePartitionInfo(dtmc, propBSs); /* Create initial partition based on propositions */
		minimised = minimiseDTMC(dtmc);
		// printProbabilities();
		if (minimised) {
			DTMCSimple<Value> dtmcNew = buildReducedDTMC();
			attachStatesAndLabels(dtmc, dtmcNew, propNames, propBSs);
			timer = System.currentTimeMillis() - timer;
			mainLog.println("Minimisation: " + numStates + " to " + numBlocks + " States");
			mainLog.println("Time for bisimulation computation: " + timer / 1000.0 + " seconds.");
			return dtmcNew;
		} /* If the state space was not minimised, do not create a reduced model */
		timer = System.currentTimeMillis() - timer;
		mainLog.println("Minimisation: " + numStates + " to " + numBlocks + " States");
		mainLog.println("Time for bisimulation computation: " + timer / 1000.0 + " seconds.");
		return dtmc;
	}

	/**
	 * Perform bisimulation minimisation on a CTMC.
	 * @param ctmc      The CTMC
	 * @param propNames Names of the propositions in {@code propBSs}
	 * @param propBSs   Propositions (satisfying sets of states) to be preserved by bisimulation.
	 */
	protected CTMC<Value> minimiseCTMC(CTMC<Value> ctmc, List<String> propNames, List<BitSet> propBSs)
	{
		long timer;
		timer = System.currentTimeMillis();
		initialisePartitionInfo(ctmc, propBSs); /* Create initial partition based on propositions */
		// printProbabilities();
		minimised = minimiseCTMC(ctmc);
		if (minimised) {
			CTMCSimple<Value> ctmcNew = buildReducedCTMC();
			attachStatesAndLabels(ctmc, ctmcNew, propNames, propBSs);
			timer = System.currentTimeMillis() - timer;
			mainLog.println("Minimisation: " + numStates + " to " + numBlocks + " States");
			mainLog.println("Time for bisimulation computation: " + timer / 1000.0 + " seconds.");
			return ctmcNew;
		} /* If the state space was not minimised, do not create a reduced model */
		timer = System.currentTimeMillis() - timer;
		mainLog.println("Minimisation: " + numStates + " to " + numBlocks + " States");
		mainLog.println("Time for bisimulation computation: " + timer / 1000.0 + " seconds.");
		return ctmc;
	}

	/**
	 * Partitions the states of the specified labelled Markov chain, given the
	 * initial partition, updating {@code numBlocks} and {@code partition}.
	 * States are in the same set, that is, are mapped to the same integer, if
	 * and only if they are probabilistic bisimilar.
	 * @param dtmc The DTMC.
	 * @return True if the state space is minimised, false otherwise.
	 */
	protected abstract boolean minimiseDTMC(DTMC<Value> dtmc);

	/**
	 * Partitions the states of the specified labelled Markov chain, given the
	 * initial partition, updating {@code numBlocks} and {@code partition}.
	 * States are in the same set, that is, are mapped to the same integer, if
	 * and only if they are probabilistic bisimilar.
	 * @param ctmc The CTMC.
	 * @return True if the state space is minimised, false otherwise.
	 */
	protected abstract boolean minimiseCTMC(CTMC<Value> ctmc);

	/**
	 * Build the reduced model.
	 * @return The reduced DTMC.
	 */
	protected abstract DTMCSimple<Value> buildReducedDTMC();

	/**
	 * Build the reduced model.
	 * @return The reduced CTMC.
	 */
	protected abstract CTMCSimple<Value> buildReducedCTMC();

	/**
	 * Construct the initial partition based on a set of proposition bitsets.
	 * Store info in {@code numStates}, {@code numBlocks} and {@code partition}.
	 */
	protected void initialisePartitionInfo(Model<Value> model, List<BitSet> propBSs)
	{
		BitSet bs1, bs0;
		numStates = model.getNumStates();
		partition = new int[numStates];

		// Compute all non-empty combinations of propositions
		List<BitSet> all = new ArrayList<BitSet>();
		bs1 = (BitSet) propBSs.get(0).clone();
		bs0 = (BitSet) bs1.clone();
		bs0.flip(0, numStates);
		all.add(bs1);
		all.add(bs0);
		int n = propBSs.size();
		for (int i = 1; i < n; i++) {
			BitSet bs = propBSs.get(i);
			int m = all.size();
			for (int j = 0; j < m; j++) {
				bs1 = all.get(j);
				bs0 = (BitSet) bs1.clone();
				bs0.andNot(bs);
				bs1.and(bs);
				if (bs1.isEmpty()) {
					all.set(j, bs0);
				} else {
					if (!bs0.isEmpty())
						all.add(bs0);
				}
			}
		}

		// Construct initial partition
		numBlocks = all.size();
		for (int j = 0; j < numBlocks; j++) {
			BitSet bs = all.get(j);
			for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i + 1)) {
				partition[i] = j;
			}
		}
	}

	/**
	 * Attach a list of states to the minimised model by adding a representative state
	 * from the original model.
	 * Also attach information about the propositions (used for bisimulation
	 * minimisation) to the minimised model, in the form of labels (stored as BitSets).
	 * @param model     The original model
	 * @param modelNew  The minimised model
	 * @param propNames The names of the propositions
	 * @param propBSs   Satisfying states (of the minimised model) for the propositions
	 */
	protected void attachStatesAndLabels(Model<Value> model, ModelExplicit<Value> modelNew, List<String> propNames, List<BitSet> propBSs)
	{
		// Attach states
		if (model.getStatesList() != null) {
			List<State> statesList = model.getStatesList();
			List<State> statesListNew = new ArrayList<State>(numBlocks);
			for (int i = 0; i < numBlocks; i++) {
				statesListNew.add(null);
			}
			for (int i = 0; i < numStates; i++) {
				if (statesListNew.get(partition[i]) == null)
					statesListNew.set(partition[i], statesList.get(i));
			}
			modelNew.setStatesList(statesListNew);
		}

		// Build/attach new labels
		int numProps = propBSs.size();
		for (int i = 0; i < numProps; i++) {
			String propName = propNames.get(i);
			BitSet propBS = propBSs.get(i);
			BitSet propBSnew = new BitSet();
			for (int j = propBS.nextSetBit(0); j >= 0; j = propBS.nextSetBit(j + 1))
				propBSnew.set(partition[j]);
			modelNew.addLabel(propName, propBSnew);
		}
	}

	/**
	 * Did bisimulation minimisation reduce the state space of the model?
	 */
	public boolean minimised()
	{
		return minimised;
	}

	/**
	 * Display the current partition, showing the states in each block.
	 */
	@SuppressWarnings("unused")
	private void printPartition(Model<Value> model)
	{
		for (int i = 0; i < numBlocks; i++) {
			mainLog.print(i + ":");
			for (int j = 0; j < numStates; j++)
				if (partition[j] == i)
					if (model.getStatesList() != null)
						mainLog.print(" " + model.getStatesList().get(j));
					else
						mainLog.print(" " + j);
			mainLog.println();
		}
	}

	/**
	 * Display the current partition, showing the states in each block.
	 */
	@SuppressWarnings("unused")
	private void printProbabilities()
	{
		for (int i = 0; i < numBlocks; i++) {
			mainLog.print(i + ":");
			int count = 0;
			for (int j = 0; j < numStates; j++)
				if (partition[j] == i)
					mainLog.print(" " + j);
					// count++;
			// mainLog.println(" " + count);
			mainLog.println();
		}
	}
}
