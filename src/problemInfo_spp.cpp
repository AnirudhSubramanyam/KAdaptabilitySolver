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


#include "problemInfo_spp.hpp"
#include <cassert>

//-----------------------------------------------------------------------------------

#define SPP_UPPER_BOUND_BUDGET 3

//-----------------------------------------------------------------------------------

void KAdaptableInfo_SPP::makeUncSet() {
	U.clear();

	// define uncertain costs
	for (int i = 1; i <= data.A; ++i) {
		U.addParam(0, 0, 1);
	}

	// define budget constraint
	std::vector<std::pair<int, double> > constraint;
	for (int i = 1; i <= data.A; ++i) {
		constraint.emplace_back(i, 1);
	}
	U.addFacet(constraint, 'L', SPP_UPPER_BOUND_BUDGET);
}

//-----------------------------------------------------------------------------------

void KAdaptableInfo_SPP::makeVars() {
	X.clear();
	Y.clear();

	// objective variable (considered to be 1st-stage)
	X.addVarType("O", 'C', -CPX_INFBOUND, +CPX_INFBOUND, 1);

	// y(i,j) : travel along arc (i,j)
	Y.addVarType("y", 'B', 0, 1, 1 + data.N, 1 + data.N);

	// All arrays/atrices are 1-indexed
	Y.setUndefinedVar("y", 0, 0);
	for (int i = 1; i <= data.N; ++i) {
		Y.setUndefinedVar("y", 0, i);
		Y.setUndefinedVar("y", i, 0);
		Y.setUndefinedVar("y", i, i);
	}

	// delete non-existent arcs
	for (int i = 1; i <= data.N; ++i) for (int j = 1; j <= data.N; ++j) if (i != j) {
		assert(data.costs[i][j] == -1.0 || data.costs[i][j] >= 0);
		if (data.costs[i][j] < 0) {
			Y.setUndefinedVar("y", i, j);
		}
	}
	assert(Y.getTotalDefVarSize() == data.A);

	// define objective variable
	X.setVarObjCoeff(1.0, "O", 0);
}

//-----------------------------------------------------------------------------------

void KAdaptableInfo_SPP::makeConsX() {


	/////////
	// B_X //
	/////////
	B_X.clear();



	/////////
	// C_X //
	/////////
	C_X.clear();



	//////////
	// C_XQ //
	//////////
	C_XQ.clear();

	// nothing to do
}

//-----------------------------------------------------------------------------------

void KAdaptableInfo_SPP::makeConsY(unsigned int l) {
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

		// bounds on y(i,j)
		for (int i = 1; i <= data.N; ++i) for (int j = 1; j <= data.N; ++j) {
			if (data.costs[i][j] > 0.0) {
				temp.clear();
				temp.addTermX(getVarIndex_2(k, "y", i, j), 1);

				temp.rowname("LB_y(" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + ")");
				temp.sign('G');
				temp.RHS(0);
				B_Y[k].emplace_back(temp);

				temp.rowname("UB_y(" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + ")");
				temp.sign('L');
				temp.RHS(1);
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

		// Flow balances
		for (int j = 1; j <= data.N; ++j) {
			temp.clear();
			temp.rowname("FLOW(" + std::to_string(j) + "," + std::to_string(k) + ")");
			temp.sign('G');
			temp.RHS((j == data.src) ? 1.0 : (j == data.tgt ? -1.0 : 0));
			for (int m = 1; m <= data.N; ++m) if (data.costs[j][m] >= -0.1) {
				temp.addTermX(getVarIndex_2(k, "y", j, m), 1.0);
			}
			for (int i = 1; i <= data.N; ++i) if (data.costs[i][j] >= -0.1) {
				temp.addTermX(getVarIndex_2(k, "y", i, j), -1.0);
			}
			C_XY[k].emplace_back(temp);
		}
	}	

		
	///////////
	// C_XYQ //
	///////////
	for (unsigned int k = C_XYQ.size(); k <= l; k++) {
		if (C_XYQ.size() < k + 1) C_XYQ.resize(k + 1);
		
		C_XYQ[k].clear();

		// objective function
		temp.clear();
		temp.rowname("OBJ_CONSTRAINT(" + std::to_string(k) + ")");
		temp.sign('G');
		temp.RHS(0);
		temp.addTermX(getVarIndex_1("O", 0), 1);
		for (int i = 1; i <= data.N; ++i) for (int j = 1; j <= data.N; ++j) {
			if (data.costs[i][j] > 0.0) {
				const int qind = getVarIndex_2(0, "y", i, j);
				const int yind = getVarIndex_2(k, "y", i, j);
				temp.addTermX(yind, -data.costs[i][j]);
				temp.addTermProduct(yind, qind, -data.costs[i][j] / 2.0);
			}
		}
		C_XYQ[k].emplace_back(temp);
	}
}

//-----------------------------------------------------------------------------------

void KAdaptableInfo_SPP::setInstance(const SPP& d) {
	data = d;
	hasInteger = 1;
	objectiveUnc = 1;
	existsFirstStage = 0;
	numFirstStage = 1;
	numSecondStage = data.A;
	numPolicies = 1;
	solfilename = data.solfilename;

	makeUncSet();
	makeVars();
	makeConsX();
	makeConsY(0);
	assert(isConsistentWithDesign());
}

//-----------------------------------------------------------------------------------
