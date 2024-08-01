//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@cs.ox.ac.uk> (University of Oxford)
//	* Ludwig Pauly <ludwigpauly@gmail.com> (TU Dresden)
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

package symbolic.build;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.BitSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Vector;
import java.util.Map.Entry;

import common.IterableStateSet;
import jdd.JDD;
import jdd.JDDNode;
import jdd.JDDVars;
import parser.Values;
import parser.VarList;
import prism.ExplicitFilesRewardGenerator;
import prism.ModelInfo;
import prism.ModelType;
import prism.Prism;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismNotSupportedException;
import symbolic.model.Model;
import symbolic.model.ModelSymbolic;
import symbolic.model.ModelVariablesDD;
import symbolic.model.NondetModel;
import symbolic.model.ProbModel;
import symbolic.model.StochModel;


/**
 * Class to convert explicit-state file storage of a model to symbolic representation.
 */
public class ExplicitFiles2MTBDD
{
	// Prism stuff
	private Prism prism;
	private PrismLog mainLog;

	// Files to read in from
	private File statesFile;
	private File transFile;
	private File labelsFile;

	// Model info
	private ModelInfo modelInfo;
	private ModelType modelType;
	private VarList varList;
	private int numVars;
	private int numStates;

	// Explicit storage of states
	private int statesArray[][] = null;

	// mtbdd stuff

	// dds/dd vars - whole system
	private JDDNode trans; // transition matrix dd
	private JDDNode start; // dd for start state
	private JDDNode stateRewards[]; // array of dds for state rewards
	private JDDNode transRewards[]; // array of dds for transition rewards
	private JDDVars allDDRowVars; // all dd vars (rows)
	private JDDVars allDDColVars; // all dd vars (cols)
	private JDDVars allDDNondetVars; // all dd vars (all non-det.)
	// dds/dd vars - variables
	private JDDVars[] varDDRowVars; // dd vars (row/col) for each module variable
	private JDDVars[] varDDColVars;
	// dds/dd vars - nondeterminism

	private ModelVariablesDD modelVariables;
	
	// action info
	private Vector<String> synchs; // list of action names
	private JDDNode transActions; // dd for transition action labels (MDPs)
	private Vector<JDDNode> transPerAction; // dds for transition action labels (D/CTMCs)

	private int maxNumChoices = 0;
	private LinkedHashMap<String, JDDNode> labelsDD;

	protected ExplicitFilesRewardGenerator efrg; // reward generator

	public ExplicitFiles2MTBDD(Prism prism)
	{
		this.prism = prism;
		mainLog = prism.getMainLog();
	}

	/**
	 * Build a Model corresponding to the passed in states/transitions/labels files.
	 * Variable info and model type is taken from a {@code ModelInfo} object.
	 * The number of states should also be passed in as {@code numStates}.
	 */
	public Model build(File statesFile, File transFile, File labelsFile, ModelInfo modelInfo, int numStates,
					   ExplicitFilesRewardGenerator efrg) throws PrismException
	{
		this.statesFile = statesFile;
		this.transFile = transFile;
		this.labelsFile = labelsFile;
		this.modelInfo = modelInfo;
		modelType = modelInfo.getModelType();
		varList = modelInfo.createVarList();
		numVars = varList.getNumVars();
		this.numStates = numStates;
		modelVariables = new ModelVariablesDD();

		this.efrg = efrg;
		
		// Build states list, if info is available
		if (statesFile != null) {
			readStatesFromFile();
		}

		return buildModel();
	}

	/** read info about reachable state space from file and store explicitly */
	private void readStatesFromFile() throws PrismException
	{
		String s, ss[];
		int i, j, lineNum = 0;

		// create arrays for explicit state storage
		statesArray = new int[numStates][];

		// open file for reading, automatic close when done
		try (BufferedReader in = new BufferedReader(new FileReader(statesFile))) {
			// skip first line
			in.readLine();
			lineNum = 1;
			// read remaining lines
			s = in.readLine();
			lineNum++;
			while (s != null) {
				// skip blank lines
				s = s.trim();
				if (s.length() > 0) {
					// split into two parts
					ss = s.split(":");
					// determine which state this line describes
					i = Integer.parseInt(ss[0]);
					// now split up middle bit and extract var info
					ss = ss[1].substring(ss[1].indexOf('(') + 1, ss[1].indexOf(')')).split(",");
					if (ss.length != numVars)
						throw new PrismException("(wrong number of variable values) ");
					if (statesArray[i] != null)
						throw new PrismException("(duplicated state) ");
					statesArray[i] = new int[numVars];
					for (j = 0; j < numVars; j++) {
						statesArray[i][j] = varList.encodeToIntFromString(j, ss[j]);
					}
				}
				// read next line
				s = in.readLine();
				lineNum++;
			}
		} catch (IOException e) {
			throw new PrismException("File I/O error reading from \"" + statesFile + "\"");
		} catch (PrismException e) {
			throw new PrismException("Error detected " + e.getMessage() + "at line " + lineNum + " of states file \"" + statesFile + "\"");
		}
	}

	/** build model */
	private Model buildModel() throws PrismException
	{
		ModelSymbolic model = null;
		JDDNode tmp, tmp2;
		JDDVars ddv;
		int i;

		// for an mdp, compute the max number of choices in a state
		if (modelType == ModelType.MDP)
			computeMaxChoicesFromFile();

		// allocate dd variables
		allocateDDVars();

		// construct transition matrix from file
		buildTrans();

		// get rid of any nondet dd variables not needed
		if (modelType == ModelType.MDP) {
			tmp = JDD.GetSupport(trans);
			tmp = JDD.ThereExists(tmp, allDDRowVars);
			tmp = JDD.ThereExists(tmp, allDDColVars);
			tmp2 = tmp;
			ddv = new JDDVars();
			while (!tmp2.equals(JDD.ONE)) {
				ddv.addVar(JDD.Var(tmp2.getIndex()));
				tmp2 = tmp2.getThen();
			}
			JDD.Deref(tmp);
			allDDNondetVars.derefAll();
			allDDNondetVars = ddv;
		}

		// 		// print dd variables actually used (support of trans)
		// 		mainLog.print("\nMTBDD variables used (" + allDDRowVars.n() + "r, " + allDDRowVars.n() + "c");
		// 		if (modelType == ModelType.MDP) mainLog.print(", " + allDDNondetVars.n() + "nd");
		// 		mainLog.print("):");
		// 		tmp = JDD.GetSupport(trans);
		// 		tmp2 = tmp;
		// 		while (!tmp2.isConstant()) {
		// 			//mainLog.print(" " + tmp2.getIndex() + ":" + ddVarNames.elementAt(tmp2.getIndex()));
		// 			mainLog.print(" " + ddVarNames.elementAt(tmp2.getIndex()));
		// 			tmp2 = tmp2.getThen();
		// 		}
		// 		mainLog.println();
		// 		JDD.Deref(tmp);

		// load labels
		loadLabels();

		// calculate dd for initial state
		buildInit();

		// compute state rewards
		buildStateRewards();

		Values constantValues = new Values(); // no constants

		// transition reward not yet supported - just make 0
		JDDNode transRewardsArray[] = new JDDNode[stateRewards.length];
		for (int k = 0; k < stateRewards.length; k++) {
			transRewardsArray[k] = JDD.Constant(0);
		}

		// create new Model object to be returned
		// they need a module name list, so we fake that
		int numModules = 1;
		String moduleNames[] = new String[] { "M" };
		String rewardStructNames[] = efrg.getRewardStructNames().toArray(new String[0]);

		if (modelType == ModelType.DTMC) {
			model = new ProbModel(trans, start, allDDRowVars, allDDColVars, modelVariables,
					varList, varDDRowVars, varDDColVars);
		} else if (modelType == ModelType.MDP) {
			model = new NondetModel(trans, start, allDDRowVars, allDDColVars, allDDNondetVars, modelVariables,
					varList, varDDRowVars, varDDColVars);
		} else if (modelType == ModelType.CTMC) {
			model = new StochModel(trans, start, allDDRowVars, allDDColVars, modelVariables,
					varList, varDDRowVars, varDDColVars);
		}
		model.setRewards(stateRewards, transRewardsArray, rewardStructNames);
		model.setConstantValues(constantValues);
		// set action info
		// TODO: disable if not required?
		model.setSynchs(synchs);
		if (modelType != ModelType.MDP) {
			((ProbModel) model).setTransPerAction((JDDNode[]) transPerAction.toArray(new JDDNode[0]));
		} else {
			((NondetModel) model).setTransActions(transActions);
		}

		// do reachability (or not)
		if (prism.getDoReach()) {
			mainLog.print("\nComputing reachable states...\n");
			model.doReachability();
			model.filterReachableStates();
		} else {
			mainLog.print("\nSkipping reachable state computation.\n");
			model.skipReachability();
			model.filterReachableStates();
		}

		// Print some info (if extraddinfo flag on)
		if (prism.getExtraDDInfo()) {
			mainLog.print("Reach: " + JDD.GetNumNodes(model.getReach()) + " nodes\n");
		}

		// find any deadlocks
		model.findDeadlocks(prism.getFixDeadlocks());

		// attach labels
		attachLabels(model);

		// deref spare dds
		if (labelsDD != null) {
			for (JDDNode d : labelsDD.values()) {
				JDD.Deref(d);
			}
		}

		return model;
	}

	/** for an mdp, compute max number of choices in a state (from transitions file) */
	private void computeMaxChoicesFromFile() throws PrismException
	{
		String s, ss[];
		int j, lineNum = 0;

		// open file for reading, automatic close when done
		try (BufferedReader in = new BufferedReader(new FileReader(transFile))) {
			// skip first line
			in.readLine();
			lineNum = 1;
			// read remaining lines
			s = in.readLine();
			lineNum++;
			maxNumChoices = 0;
			while (s != null) {
				s = s.trim();
				if (s.length() > 0) {
					ss = s.split(" ");
					if (ss.length < 4 || ss.length > 5)
						throw new PrismException("");
					j = Integer.parseInt(ss[1]);
					if (j + 1 > maxNumChoices)
						maxNumChoices = j + 1;
				}
				s = in.readLine();
				lineNum++;
			}
		} catch (IOException e) {
			throw new PrismException("File I/O error reading from \"" + transFile + "\": " + e.getMessage());
		} catch (NumberFormatException e) {
			throw new PrismException("Error detected at line " + lineNum + " of transition matrix file \"" + transFile + "\"");
		} catch (PrismException e) {
			throw new PrismException("Error detected " + e.getMessage() + "at line " + lineNum + " of transition matrix file \"" + transFile + "\"");
		}
	}

	/**
	 * allocate DD vars for system
	 * i.e. decide on variable ordering and request variables from CUDD
	 */
	private void allocateDDVars()
	{
		int i, j, n;

		modelVariables = new ModelVariablesDD();

		// create arrays/etc. first

		// module variable (row/col) vars
		varDDRowVars = new JDDVars[numVars];
		varDDColVars = new JDDVars[numVars];
		for (i = 0; i < numVars; i++) {
			varDDRowVars[i] = new JDDVars();
			varDDColVars[i] = new JDDVars();
		}

		// now allocate variables

		// allocate nondeterministic variables
		if (modelType == ModelType.MDP) {
			allDDNondetVars = new JDDVars();
			for (i = 0; i < maxNumChoices; i++) {
				allDDNondetVars.addVar(modelVariables.allocateVariable("l" + i));
			}
		}

		// allocate dd variables for module variables (i.e. rows/cols)
		// go through all vars in order (incl. global variables)
		// so overall ordering can be specified by ordering in the input file
		allDDRowVars = new JDDVars();
		allDDColVars = new JDDVars();
		for (i = 0; i < numVars; i++) {
			// get number of dd variables needed
			// (ceiling of log2 of range of variable)
			n = varList.getRangeLogTwo(i);
			// add pairs of variables (row/col)
			for (j = 0; j < n; j++) {
				// new dd row variable
				varDDRowVars[i].addVar(modelVariables.allocateVariable(varList.getName(i) + "." + j));
				// new dd col variable
				varDDColVars[i].addVar(modelVariables.allocateVariable(varList.getName(i) + "'." + j));
			}
			allDDRowVars.copyVarsFrom(varDDRowVars[i]);
			allDDColVars.copyVarsFrom(varDDColVars[i]);
		}
	}

	/** Construct transition matrix from file */
	private void buildTrans() throws PrismException
	{
		String s, ss[], a;
		int i, j, r, c, k = 0, lineNum = 0;
		double d;
		JDDNode elem, tmp;

		// initialise action list
		synchs = new Vector<String>();

		// initialise mtbdds
		trans = JDD.Constant(0);
		if (modelType != ModelType.MDP) {
			transPerAction = new Vector<JDDNode>();
			transPerAction.add(JDD.Constant(0));
		} else {
			transActions = JDD.Constant(0);
		}

		// open file for reading, automatic close when done
		try (BufferedReader in = new BufferedReader(new FileReader(transFile))) {
			// skip first line
			in.readLine();
			lineNum = 1;
			// read remaining lines
			s = in.readLine();
			lineNum++;
			while (s != null) {
				// skip blank lines
				s = s.trim();
				if (s.length() > 0) {
					a = "";
					// parse line, split into parts
					ss = s.split(" ");
					// case for dtmcs/ctmcs...
					if (modelType != ModelType.MDP) {
						if (ss.length < 3 || ss.length > 4)
							throw new PrismException("");
						r = Integer.parseInt(ss[0]);
						c = Integer.parseInt(ss[1]);
						d = Double.parseDouble(ss[2]);
						if (ss.length == 4) {
							a = ss[3];
						}
					}
					// case for mdps...
					else {
						if (ss.length < 4 || ss.length > 5)
							throw new PrismException("");
						r = Integer.parseInt(ss[0]);
						k = Integer.parseInt(ss[1]);
						c = Integer.parseInt(ss[2]);
						d = Double.parseDouble(ss[3]);
						if (ss.length == 5) {
							a = ss[4];
						}
					}
					// construct element of matrix mtbdd
					// case where we don't have a state list...
					if (statesFile == null) {
						/// ...for dtmcs/ctmcs...
						if (modelType != ModelType.MDP) {
							elem = JDD.SetMatrixElement(JDD.Constant(0), varDDRowVars[0], varDDColVars[0], r, c, 1.0);
						}
						/// ...for mdps...
						else {
							elem = JDD.Set3DMatrixElement(JDD.Constant(0), varDDRowVars[0], varDDColVars[0], allDDNondetVars, r, c, k, 1.0);
						}
					}
					// case where we do have a state list...
					else {
						elem = JDD.Constant(1);
						for (i = 0; i < numVars; i++) {
							elem = JDD.Apply(JDD.TIMES, elem, JDD.SetVectorElement(JDD.Constant(0), varDDRowVars[i], statesArray[r][i], 1));
							elem = JDD.Apply(JDD.TIMES, elem, JDD.SetVectorElement(JDD.Constant(0), varDDColVars[i], statesArray[c][i], 1));
						}
						if (modelType == ModelType.MDP) {
							elem = JDD.Apply(JDD.TIMES, elem, JDD.SetVectorElement(JDD.Constant(0), allDDNondetVars, k, 1));
						}
					}
					// add it into mtbdds for transition matrix and transition rewards
					JDD.Ref(elem);
					trans = JDD.Apply(JDD.PLUS, trans, JDD.Apply(JDD.TIMES, JDD.Constant(d), elem));
					// look up action name
					if (!("".equals(a))) {
						j = synchs.indexOf(a);
						// add to list if first time seen 
						if (j == -1) {
							synchs.add(a);
							j = synchs.size() - 1;
						}
						j++;
					} else {
						j = 0;
					}
					/// ...for dtmcs/ctmcs...
					if (modelType != ModelType.MDP) {
						// get (or create) dd for action j
						if (j < transPerAction.size()) {
							tmp = transPerAction.get(j);
						} else {
							tmp = JDD.Constant(0);
							transPerAction.add(tmp);
						}
						// add element to matrix
						JDD.Ref(elem);
						tmp = JDD.Apply(JDD.PLUS, tmp, JDD.Apply(JDD.TIMES, JDD.Constant(d), elem));
						transPerAction.set(j, tmp);
					}
					/// ...for mdps...
					else {
						JDD.Ref(elem);
						tmp = JDD.ThereExists(elem, allDDColVars);
						// use max here because we see multiple transitions for a sinlge choice
						transActions = JDD.Apply(JDD.MAX, transActions, JDD.Apply(JDD.TIMES, JDD.Constant(j), tmp));
					}
					// deref element dd
					JDD.Deref(elem);
				}
				// read next line
				s = in.readLine();
				lineNum++;
			}
		} catch (IOException e) {
			throw new PrismException("File I/O error reading from \"" + transFile + "\": " + e.getMessage());
		} catch (NumberFormatException e) {
			throw new PrismException("Error detected at line " + lineNum + " of transition matrix file \"" + transFile + "\"");
		} catch (PrismException e) {
			throw new PrismException("Error detected " + e.getMessage() + "at line " + lineNum + " of transition matrix file \"" + transFile + "\"");
		}
	}

	/** calculate dd for initial state */
	private void buildInit() throws PrismException
	{
		// If no labels are provided, just use state 0 (i.e. min value for each var) 
		if (labelsDD == null) {
			start = JDD.Constant(1);
			for (int i = 0; i < numVars; i++) {
				JDDNode tmp = JDD.SetVectorElement(JDD.Constant(0), varDDRowVars[i], 0, 1);
				start = JDD.And(start, tmp);
			}
		}
		// Otherwise, construct from the labels
		else {
			start = labelsDD.get("init").copy();
			if (start == null || start.equals(JDD.ZERO)) {
				throw new PrismException("No initial states found in labels file");
			}
		}
	}

	/** Load label information form label file (if exists) */
	private void loadLabels() throws PrismException
	{
		if (labelsFile == null) {
			return;
		}

		// we construct BitSets for the labels
		Map<String, BitSet> labels = explicit.StateModelChecker.loadLabelsFile(labelsFile.getAbsolutePath());
		labelsDD = new LinkedHashMap<>();

		for (Entry<String, BitSet> e : labels.entrySet()) {
			// compute dd

			JDDNode labelStatesDD = JDD.Constant(0);

			for (int state : new IterableStateSet(e.getValue(), numStates)) {
				JDDNode tmp;
				// case where we don't have a state list...
				if (statesFile == null) {
					tmp = JDD.SetVectorElement(JDD.Constant(0), varDDRowVars[0], state, 1);
				}
				// case where we do have a state list...
				else {
					tmp = JDD.Constant(1);
					for (int i = 0; i < numVars; i++) {
						tmp = JDD.Apply(JDD.TIMES, tmp, JDD.SetVectorElement(JDD.Constant(0), varDDRowVars[i], statesArray[state][i], 1));
					}
				}

				// add it into bdd
				labelStatesDD = JDD.Or(labelStatesDD, tmp);
			}

			// store the dd
			labelsDD.put(e.getKey(), labelStatesDD);
		}
	}

	/** Attach the computed label information to the model */
	private void attachLabels(ModelSymbolic model) throws PrismNotSupportedException
	{
		if (labelsDD == null)
			return;

		for (Entry<String, JDDNode> e : labelsDD.entrySet()) {
			if (e.equals("init")) {
				// handled in buildInit()
				continue;
			}
			if (e.equals("deadlock")) {
				// ignored (info is recomputed)
				continue;
			}
			model.addLabelDD(e.getKey(), e.getValue().copy());
		}
	}

	/**
	 * Load state rewards from the ExplicitFilesRewardGenerator.
	 */
	private void buildStateRewards() throws PrismException
	{
		int numRewardStructs = efrg.getNumRewardStructs();
		stateRewards = new JDDNode[numRewardStructs];
		for (int r = 0; r < numRewardStructs; r++) {
			stateRewards[r] = JDD.Constant(0);
			int finalR = r;
			efrg.extractStateRewards(r, (i, d) -> storeReward(finalR, i, d));
		}
	}

	/**
	 * Stores stateRewards in the required format for mtbdd.
	 *
	 * @param rewardStructIndex reward structure index
	 * @param i state index
	 * @param d reward value
	 */
	protected void storeReward(int rewardStructIndex, int i, double d)
	{
		// construct element of vector mtbdd
		// case where we don't have a state list...
		JDDNode tmp;
		if (statesArray == null) {
			tmp = JDD.SetVectorElement(JDD.Constant(0), varDDRowVars[0], i, 1.0);
		}
		// case where we do have a state list...
		else {
			tmp = JDD.Constant(1);
			for (int j = 0; j < numVars; j++) {
				tmp = JDD.Times(tmp, JDD.SetVectorElement(JDD.Constant(0), varDDRowVars[j], statesArray[i][j], 1));
			}
		}
		// add it into mtbdd for state rewards
		stateRewards[rewardStructIndex] = JDD.Plus(stateRewards[rewardStructIndex], JDD.Times(JDD.Constant(d), tmp));
	}
}

//------------------------------------------------------------------------------
