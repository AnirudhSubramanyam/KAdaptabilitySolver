/***************************************************************************************/
/*                                                                                     */
/*  Copyright 2018 by Anirudh Subramanyam, Chrysanthos Gounaris and Wolfram Wiesemann  */
/*                                                                                     */
/*  Licensed under the FreeBSD License (the "License").                                */
/*  You may not use this file except in compliance with the License.                   */
/*  You may obtain a copy of the License at                                            */
/*                                                                                     */
/*  https://www.freebsd.org/copyright/freebsd-license.html                             */
/*                                                                                     */
/***************************************************************************************/


#include "problemInfo_psp.hpp"
#include <cassert>

//-----------------------------------------------------------------------------------

#define FACTOR_ALPHA 0.5
#define FACTOR_BETA 1.0
#define DO_MAX 1
#define SPECIAL_RESTRICTIVE_LDR 0

//-----------------------------------------------------------------------------------

void KAdaptableInfo_PSP::makeUncSet() {
	U.clear();

	if (data.special) {
		// compute k
		int k = (data.N - 1)/3;

		// define primary uncertain parameters
		for (int i = 1; i <= k; ++i) {
			U.addParam(0.5, 0, 1);
		}

		// define secondary parameters representing absolute deviation from 0.5
		for (int i = 1; i <= k; ++i) {
			U.addParam(0, 0, 0.5);
		}

		// define budget constraint
		std::vector<std::pair<int, double> > constraint;
		for (int i = k + 1; i <= (2 * k); ++i) {
			constraint.emplace_back(i, 1);
		}
		U.addFacet(constraint, 'L', 0.5);

		// enforce absolute value definition
		for (int i = 1; i <= k; i++) {
			constraint.clear();
			constraint.emplace_back(i + k, 1);
			constraint.emplace_back(i, -1);
			U.addFacet(constraint, 'G', -0.5);

			constraint.clear();
			constraint.emplace_back(i + k, 1);
			constraint.emplace_back(i, +1);
			U.addFacet(constraint, 'G', +0.5);
		}
	}
	else {
		// define uncertain risk factors
		for (unsigned int f = 1; f < data.phi[0].size(); ++f) {
			U.addParam(0, -1, 1);
		}

		// define beta-net-alpha
		std::vector<std::pair<int, double> > constraint;
		for (unsigned int f = 1; f < data.phi[0].size(); ++f) {
			constraint.emplace_back(f, 1);
		}

		U.addFacet(constraint, 'L', FACTOR_BETA * ((int)data.phi[0].size() - 1) );
		U.addFacet(constraint, 'G', -FACTOR_BETA * ((int)data.phi[0].size() - 1));
	}
}

//-----------------------------------------------------------------------------------

void KAdaptableInfo_PSP::makeVars() {
	X.clear();
	Y.clear();

	// objective variable (considered to be 1st-stage)
	X.addVarType("O", 'C', -CPX_INFBOUND, +CPX_INFBOUND, 1);

	// x(i) : resource allocated to task i
	if (!data.special) X.addVarType("x", 'C', 0, 0.5, 1 + data.N);

	// y(i) : start time of task i
	Y.addVarType("y", 'C', (data.affine ? (data.special ? (SPECIAL_RESTRICTIVE_LDR ? 0 : -CPX_INFBOUND) : -CPX_INFBOUND) : 0), +CPX_INFBOUND, 1 + data.N);

	// gamma(i,j) : slopes of affine decision rules
	if (data.affine) {
		if (!data.special) {
			Y.addVarType("g", 'C', -CPX_INFBOUND, +CPX_INFBOUND, data.N, data.phi[0].size());
		} else {
			Y.addVarType("g", 'C', (SPECIAL_RESTRICTIVE_LDR ? 0 : -CPX_INFBOUND), +CPX_INFBOUND, data.N, 1 + ((data.N - 1) / 3));
		}
	}

	// All arrays/atrices are 1-indexed
	if (!data.special) X.setUndefinedVar("x", 0);
	Y.setUndefinedVar("y", 0);
	if (data.affine) {
		for (int i = 0; i < data.N; i++) Y.setUndefinedVar("g", i, 0);
		for (int f = 0, flim = (data.special ? (1 + ((data.N - 1)/3)) : data.phi[0].size()); f < flim; f++) Y.setUndefinedVar("g", 0, f);
	}

	// cannot adjust uncertain durations of dummy nodes
	if (data.affine) for (int f = 0, flim = (data.special ? (1 + ((data.N - 1)/3)) : data.phi[0].size()); f < flim; f++) {
		Y.setUndefinedVar("g", 1, f);
	}

	// define objective variable
	X.setVarObjCoeff(1.0, "O", 0);

	// fix allocations of dummy nodes
	if (!data.special) X.setVarUB(0, "x", 1);
	if (!data.special) X.setVarUB(0, "x", data.N);

	// set early-start schedule
	Y.setVarLB(0, "y", 1);
	Y.setVarUB(0, "y", 1);
	if (data.affine) {
		Y.setVarLB(0, "y", data.N);
	}

	// restrict slopes of affine decision rules
	if (data.special && data.affine) {
		int k = (data.N - 1)/3;
		for (int i = 1; i <= k; i++) for (int f = i + 1; f <= k; f++) {
			Y.setVarUB(0, "g", (3*i) - 1, f);
			Y.setVarUB(0, "g", (3*i), f);
			if (i != k) Y.setVarUB(0, "g", (3*i) + 1, f);
		}
	}
}

//-----------------------------------------------------------------------------------

void KAdaptableInfo_PSP::makeConsX() {
	ConstraintExpression temp;

	/////////
	// B_X //
	/////////
	B_X.clear();

	// bounds on x(i)
	if (!data.special) for (int i = 1; i <= data.N; ++i) {
		temp.clear();
		temp.addTermX(getVarIndex_1("x", i), 1);

		temp.rowname("LB_x(" + std::to_string(i) + ")");
		temp.sign('G');
		temp.RHS(0);
		B_X.emplace_back(temp);

		temp.rowname("UB_x(" + std::to_string(i) + ")");
		temp.sign('L');
		temp.RHS(0.5 * ((i > 1) && (i < data.N)));
		B_X.emplace_back(temp);
	}



	/////////
	// C_X //
	/////////
	C_X.clear();
	
	// Sum(x) <= B
	if (!data.special) {
		temp.clear();
		temp.rowname("BUDGET");
		temp.sign('L');
		temp.RHS(data.B);
		for (int i = 1; i <= data.N; ++i) {
			temp.addTermX(getVarIndex_1("x", i), 1);
		}
		C_X.emplace_back(temp);
	}



	//////////
	// C_XQ //
	//////////
	C_XQ.clear();

	// nothing to do
}

//-----------------------------------------------------------------------------------

void KAdaptableInfo_PSP::makeConsY(unsigned int l) {
	assert(C_XY.size() == B_Y.size());
	assert(C_XY.size() == C_XYQ.size());
	assert(numPolicies >= l);
	if (l == 0) {
		B_Y.clear();
		C_XY.clear();
		C_XYQ.clear();
	}
	ConstraintExpression temp;

	/////////
	// B_Y //
	/////////
	for (unsigned int k = B_Y.size(); k <= l; k++) {
		if (B_Y.size() < k + 1) B_Y.resize(k + 1);

		B_Y[k].clear();

		// early-start schedule
		temp.clear();
		temp.addTermX(getVarIndex_2(k, "y", 1), 1);
		
		temp.rowname("LB_y(1," + std::to_string(k) + ")");
		temp.sign('G');
		temp.RHS(0);
		B_Y[k].emplace_back(temp);

		temp.rowname("UB_y(1," + std::to_string(k) + ")");
		temp.sign('L');
		temp.RHS(0);
		B_Y[k].emplace_back(temp);

		// non-negative makespan
		temp.clear();
		temp.addTermX(getVarIndex_2(k, "y", data.N), 1);
		temp.rowname("LB_y(" + std::to_string(data.N) + "," + std::to_string(k) + ")");
		temp.sign('G');
		temp.RHS(0);
		B_Y[k].emplace_back(temp);

		// bounds on y(i) start times
		if (!data.affine || (data.special && SPECIAL_RESTRICTIVE_LDR)) for (int i = 2; i < data.N; i++) {
			temp.clear();
			temp.addTermX(getVarIndex_2(k, "y", i), 1);
			temp.rowname("LB_y(" + std::to_string(i) + "," + std::to_string(k) + ")");
			temp.sign('G');
			temp.RHS(0);
			B_Y[k].emplace_back(temp);

			// bound slopes of decision rules g(i,f) to be non-negative
			if (data.affine && data.special && SPECIAL_RESTRICTIVE_LDR) for (int f = 1; f <= ((data.N - 1)/3); f++) {
				temp.clear();
				temp.addTermX(getVarIndex_2(k, "g", i, f), 1);
				temp.rowname("LB_g(" + std::to_string(i) + "," + std::to_string(f) + "," + std::to_string(k) + ")");
				temp.sign('G');
				temp.RHS(0);
				B_Y[k].emplace_back(temp);
			}
		}
	}


	//////////
	// C_XY //
	//////////
	for (unsigned int k = C_XY.size(); k <= l; k++) {
		if (C_XY.size() < k + 1) C_XY.resize(k + 1);

		C_XY[k].clear();

		// objective function
		temp.clear();
		temp.rowname("OBJ_CONSTRAINT(" + std::to_string(k) + ")");
		temp.sign('G');
		temp.RHS(0);
		temp.addTermX(getVarIndex_1("O", 0), 1);
		temp.addTermX(getVarIndex_2(k, "y", data.N), -1);
		C_XY[k].emplace_back(temp);


		// precedence relationships
		if (!data.special && !data.affine && DO_MAX) for (int i = 2; i < data.N; i++) for (int j : data.successor[i]) {
			temp.clear();
			temp.rowname("PRE(" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + ")");
			temp.sign('G');
			temp.RHS(data.duration[i]);
			temp.addTermX(getVarIndex_2(k, "y", j), 1);
			temp.addTermX(getVarIndex_2(k, "y", i), -1);
			temp.addTermX(getVarIndex_1("x", i), data.duration[i]);
			C_XY[k].emplace_back(temp);
		}
	}	

		
	///////////
	// C_XYQ //
	///////////
	for (unsigned int k = C_XYQ.size(); k <= l; k++) {
		if (C_XYQ.size() < k + 1) C_XYQ.resize(k + 1);
		
		C_XYQ[k].clear();


		// precedence relationships
		for (int i = 2; i < data.N; i++) for (int j : data.successor[i]) {
			temp.clear();
			temp.rowname("PRE(" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + ")");
			temp.sign('G');
			temp.RHS(data.special ? 0 : data.duration[i]);
			temp.addTermX(getVarIndex_2(k, "y", j), 1);
			if (data.affine) if (j < data.N) for (int f = 1, flim = (data.special ? (1 + ((data.N - 1)/3)) : data.phi[i].size()); f < flim; f++) {
				temp.addTermProduct(getVarIndex_2(k, "g", j, f), U.getParamIndex('q', f), 1);
			}
			temp.addTermX(getVarIndex_2(k, "y", i), -1);
			if (data.affine) for (int f = 1, flim = (data.special ? (1 + ((data.N - 1)/3)) : data.phi[i].size()); f < flim; f++) {
				temp.addTermProduct(getVarIndex_2(k, "g", i, f), U.getParamIndex('q', f), -1);
			}
			if (!data.special) {
				temp.addTermX(getVarIndex_1("x", i), data.duration[i]);
				for (unsigned int f = 1; f < data.phi[i].size(); f++) {
					double coef = data.phi[i][f] * data.duration[i] * FACTOR_ALPHA;

					if (coef != 0.0) {
						temp.addTermProduct(getVarIndex_1("x", i), U.getParamIndex('q', f), coef);
						temp.addTermQ(U.getParamIndex('q', f), -coef);
					}
				}
			}
			else {
				if (data.duration[i] > 0.5) {
					int l = (i + 1)/3;
					temp.addTermQ(U.getParamIndex('q', l), -1);
				}
				else if (data.duration[i] < -0.5) {
					int l = i / 3;
					temp.RHS(1);
					temp.addTermQ(U.getParamIndex('q', l), 1);
				}
				else if (!data.affine) {
					C_XY[k].emplace_back(temp);
					continue;
				}
			}
			C_XYQ[k].emplace_back(temp);
		}


		if (!data.special && data.affine && DO_MAX) for (int i = 2; i < data.N; i++) for (int j : data.successor[i]) {
			temp.clear();
			temp.rowname("PRE(" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + ")");
			temp.sign('G');
			temp.RHS(data.duration[i]);
			temp.addTermX(getVarIndex_2(k, "y", j), 1);
			if (j < data.N) for (unsigned int f = 1; f < data.phi[i].size(); f++) {
				temp.addTermProduct(getVarIndex_2(k, "g", j, f), U.getParamIndex('q', f), 1);
			}
			temp.addTermX(getVarIndex_2(k, "y", i), -1);
			for (unsigned int f = 1; f < data.phi[i].size(); f++) {
				temp.addTermProduct(getVarIndex_2(k, "g", i, f), U.getParamIndex('q', f), -1);
			}
			temp.addTermX(getVarIndex_1("x", i), data.duration[i]);
			C_XYQ[k].emplace_back(temp);
		}



		// non-negative start times
		if (data.affine) if (data.special ? (SPECIAL_RESTRICTIVE_LDR ? 0 : 1) : 1) for (int i = 2; i < data.N; i++) {
			temp.clear();
			temp.rowname("GREATER(" + std::to_string(i) + "," + std::to_string(k) + ")");
			temp.sign('G');
			temp.RHS(0);
			temp.addTermX(getVarIndex_2(k, "y", i), 1);
			for (int f = 1, flim = (data.special ? (1 + ((data.N - 1)/3)) : data.phi[i].size()); f < flim; f++) {
				temp.addTermProduct(getVarIndex_2(k, "g", i, f), U.getParamIndex('q', f), 1);
			}
			C_XYQ[k].emplace_back(temp);
		}
	}
}

//-----------------------------------------------------------------------------------

void KAdaptableInfo_PSP::setInstance(const PSP& d) {
	data = d;
	hasInteger = 0;
	objectiveUnc = 0;
	existsFirstStage = !data.special;
	numFirstStage = 1 + (data.special ? 0 : data.N);
	numSecondStage = data.N + (data.affine * (data.N - 2) * (data.special ? ((data.N - 1)/3) : ((int)data.phi[0].size() - 1)));
	numPolicies = 1;
	solfilename = data.solfilename;

	makeUncSet();
	makeVars();
	makeConsX();
	makeConsY(0);
	assert(isConsistentWithDesign());
}

//-----------------------------------------------------------------------------------