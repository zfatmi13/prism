//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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

import java.util.Arrays;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import parser.State;
import parser.Values;
import prism.PrismException;
import prism.PrismNotSupportedException;
import explicit.rewards.MCRewards;

/**
* Simple explicit-state representation of a DTMC, constructed (implicitly) as the uniformised DTMC of a CTMC.
* This class is read-only: most of data is pointers to other model info.
* <br>
* Note: This implicitly constructed DTMC does not provide implementations for
* all methods of a full DTMCExplicit model. See {@link CTMCSimple#buildUniformisedDTMC} for
* a method to obtain an explicit uniformised DTMC.
*/
public class DTMCUniformisedSimple<Value> extends DTMCExplicit<Value>
{
	// Parent CTMC
	protected CTMCSimple<Value> ctmc;
	// Uniformisation rate
	protected Value q;
	// Number of extra transitions added (just for stats)
	protected int numExtraTransitions;

	/**
	 * Constructor: create from CTMC and uniformisation rate q.
	 */
	public DTMCUniformisedSimple(CTMCSimple<Value> ctmc, Value q)
	{
		this.ctmc = ctmc;
		this.numStates = ctmc.getNumStates();
		this.q = q;
		numExtraTransitions = 0;
		for (int i = 0; i < numStates; i++) {
			if (!ctmc.getTransitions(i).contains(i) && !getEvaluator().geq(ctmc.getTransitions(i).sumAllBut(i), q)) {
				numExtraTransitions++;
			}
		}
	}

	/**
	 * Constructor: create from CTMC and its default uniformisation rate.
	 */
	public DTMCUniformisedSimple(CTMCSimple<Value> ctmc)
	{
		this(ctmc, ctmc.getDefaultUniformisationRate());
	}

	@Override
	public void buildFromPrismExplicit(String filename) throws PrismException
	{
		throw new PrismNotSupportedException("Not supported");
	}
	
	// Accessors (for Model)

	public int getNumStates()
	{
		return ctmc.getNumStates();
	}

	public int getNumInitialStates()
	{
		return ctmc.getNumInitialStates();
	}

	public Iterable<Integer> getInitialStates()
	{
		return ctmc.getInitialStates();
	}

	public int getFirstInitialState()
	{
		return ctmc.getFirstInitialState();
	}

	public boolean isInitialState(int i)
	{
		return ctmc.isInitialState(i);
	}

	@Override
	public int getNumDeadlockStates()
	{
		return ctmc.getNumDeadlockStates();
	}

	@Override
	public Iterable<Integer> getDeadlockStates()
	{
		return ctmc.getDeadlockStates();
	}

	@Override
	public StateValues getDeadlockStatesList()
	{
		return ctmc.getDeadlockStatesList();
	}

	@Override
	public int getFirstDeadlockState()
	{
		return ctmc.getFirstDeadlockState();
	}

	@Override
	public boolean isDeadlockState(int i)
	{
		return ctmc.isDeadlockState(i);
	}

	public List<State> getStatesList()
	{
		return ctmc.getStatesList();
	}
	
	public Values getConstantValues()
	{
		return ctmc.getConstantValues();
	}
	
	public int getNumTransitions()
	{
		return ctmc.getNumTransitions() + numExtraTransitions;
	}

	public int getNumTransitions(int s)
	{
		// TODO
		throw new RuntimeException("Not implemented yet");
	}

	public SuccessorsIterator getSuccessors(final int s)
	{
		// TODO
		throw new Error("Not yet supported");
	}

	public int getNumChoices(int s)
	{
		// Always 1 for a DTMC
		return 1;
	}

	public void findDeadlocks(boolean fix) throws PrismException
	{
		// No deadlocks by definition
	}

	public void checkForDeadlocks() throws PrismException
	{
		// No deadlocks by definition
	}

	public void checkForDeadlocks(BitSet except) throws PrismException
	{
		// No deadlocks by definition
	}

	@Override
	public String infoString()
	{
		String s = "";
		s += getNumStates() + " states (" + getNumInitialStates() + " initial)";
		s += ", " + getNumTransitions() + " transitions (incl. " + numExtraTransitions + " self-loops)";
		return s;
	}

	// Accessors (for DTMC)

	public Iterator<Entry<Integer,Value>> getTransitionsIterator(int s)
	{
		// TODO
		throw new RuntimeException("Not implemented yet");
	}

	@Override
	public double mvMultSingle(int s, double vect[])
	{
		double qDouble = getEvaluator().toDouble(q);
		Distribution<Value> distr = ctmc.getTransitions(s);
		double sum = 0.0, d = 0.0;
		for (Map.Entry<Integer, Value> e : distr) {
			int k = e.getKey();
			double prob = getEvaluator().toDouble(e.getValue());
			// Non-diagonal entries
			if (k != s) {
				sum += prob;
				d += (prob / qDouble) * vect[k];
			}
		}
		// Diagonal entry
		if (sum < qDouble) {
			d += (1 - sum/qDouble) * vect[s];
		}

		return d;
	}

	@Override
	public double mvMultJacSingle(int s, double vect[])
	{
		double qDouble = getEvaluator().toDouble(q);
		Distribution<Value> distr = ctmc.getTransitions(s);
		double sum = 0.0, d = 0.0;
		for (Map.Entry<Integer, Value> e : distr) {
			int k = e.getKey();
			double prob = (Double) e.getValue();
			// Non-diagonal entries only
			if (k != s) {
				sum += prob;
				d += (prob / qDouble) * vect[k];
			}
		}
		// Diagonal entry is 1 - sum/q
		d /= (sum / qDouble);

		return d;
	}

	public double mvMultRewSingle(int s, double vect[], MCRewards<Double> mcRewards)
	{
		// TODO
		throw new Error("Not yet supported");
	}

	@Override
	public void vmMult(double vect[], double result[])
	{
		double qDouble = getEvaluator().toDouble(q);
		// Initialise result to 0
		Arrays.fill(result, 0);
		// Go through matrix elements (by row)
		for (int state = 0; state < numStates; state++) {
			double sum = 0.0;
			for (Iterator<Entry<Integer, Value>> transitions = ctmc.getTransitionsIterator(state); transitions.hasNext();) {
				Entry<Integer, Value> trans = transitions.next();
				int target  = trans.getKey();
				double prob = getEvaluator().toDouble(trans.getValue()) / qDouble;
				// Non-diagonal entries only
				if (target != state) {
					sum += prob;
					result[target] += prob * vect[state];
				}
			}
			// Diagonal entry is 1 - sum
			result[state] += (1 - sum) * vect[state];
		}
	}

	@Override
	public String toString()
	{
		String s = "";
		s += "ctmc: " + ctmc;
		s = ", q: " + q;
		return s;
	}

	@Override
	public boolean equals(Object o)
	{
		if (o == null || !(o instanceof DTMCUniformisedSimple))
			return false;
		DTMCUniformisedSimple<?> dtmc = (DTMCUniformisedSimple<?>) o;
		if (!ctmc.equals(dtmc.ctmc))
			return false;
		if (q != dtmc.q)
			return false;
		if (numExtraTransitions != dtmc.numExtraTransitions)
			return false;
		return true;
	}
}
