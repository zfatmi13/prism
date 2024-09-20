//==============================================================================
//
//	Copyright (c) 2025-
//	Authors:
//	* Zainab Fatmi
//	* Dave Parker <david.parker@cs.ox.ac.uk> (University of Oxford)
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

import explicit.Model;
import prism.PrismComponent;
import prism.PrismException;

import java.util.BitSet;
import java.util.List;

/**
 * Handles bisimulation minimisation for explicit-state models. Instantiates the
 * bisimulation minimisation method specified by the {@code bisimmethod} option.
 * If the class has a registered alias, that can be used. Otherwise, reflection
 * is used to instantiate the class. The default is BisimulationMethodNew.
 */
public class BisimulationTools<Value> {

	protected Bisimulation<Value> bisim;

	// aliases
	protected final String PACKAGE = "explicit.bisim.";
	protected final String DEFAULT = "new";
	protected final String OLD = "old";
	protected final String ROBUST = "robust";

	/**
	 * Perform bisimulation minimisation on a model.
	 *
	 * @param parent      The parent component
	 * @param bisimMethod The bisimulation minimization method
	 * @param model       The model
	 * @param propNames   Names of the propositions in {@code propBSs}
	 * @param propBSs     Propositions (satisfying sets of states) to be preserved by bisimulation.
	 */
	@SuppressWarnings("unchecked")
	public Model<Value> minimise(PrismComponent parent, String bisimMethod, Model<Value> model, List<String> propNames, List<BitSet> propBSs) throws PrismException {
		if (bisimMethod == null) {
			/* use default */
			bisim = new BisimulationMethodNew<Value>(parent);
		} else {
			/* use alias */
			switch (bisimMethod) {
				case (OLD):
					bisim = new BisimulationMethodOld<Value>(parent);
					break;
				case (DEFAULT):
					bisim = new BisimulationMethodNew<Value>(parent);
					break;
				case (ROBUST):
					bisim = new RobustBisimulation<Value>(parent);
					break;
				default:
					/* use reflection */
					try {
						/* class is in the bisimulation package */
						Class<?> bisimClass = Class.forName(PACKAGE + bisimMethod);
						bisim = (Bisimulation<Value>) bisimClass.getConstructor(PrismComponent.class).newInstance(parent);
					} catch (ClassNotFoundException e) {
						/* class may be elsewhere */
						try {
							Class<?> bisimClass = Class.forName(bisimMethod);
							bisim = (Bisimulation<Value>) bisimClass.getConstructor(PrismComponent.class).newInstance(parent);
						} catch (ClassNotFoundException exception) {
							throw new PrismException("Unknown bisimulation minimisation method");
						} catch (ReflectiveOperationException exception) {
							throw new PrismException("Could not instantiate " + bisimMethod + " using the required constructor");
						}
					} catch (ReflectiveOperationException e) {
						throw new PrismException("Could not instantiate " + bisimMethod + " using the required constructor");
					}
			}
		}
		return bisim.minimise(model, propNames, propBSs);
	}

	/**
	 * Did bisimulation minimisation reduce the state space of the model?
	 */
	public boolean minimised() {
		if (bisim == null) {
			return false;
		}
		return bisim.minimised();
	}
}
