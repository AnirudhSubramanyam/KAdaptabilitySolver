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


#include "problemInfo.hpp"
#include <cassert>


//-----------------------------------------------------------------------------------

void KAdaptableInfo::resize(unsigned int K) {
	assert(isConsistentWithDesign());
	unsigned int l = numPolicies;
	if (numPolicies < K) {
		numPolicies = K;
	}
	for (unsigned int k = l; k < numPolicies; k++) {
		makeConsY(k);
	}
	assert(isConsistentWithDesign());
}

//-----------------------------------------------------------------------------------


bool KAdaptableInfo::isConsistentWithDesign() const {
	if (numPolicies < 1) return 0;
	if (C_XY.size() != numPolicies) return 0;
	if (C_XYQ.size() != numPolicies) return 0;
	if (B_Y.size() != numPolicies) return 0;
	if (X.getTotalDefVarSize() <= 0) return 0;
	if (Y.getTotalDefVarSize() < 0) return 0;
	if (X.getTotalDefVarSize() != getNumFirstStage()) return 0;
	if (Y.getTotalDefVarSize() != getNumSecondStage()) return 0;
	if (existsFirstStage) {
		if (getNumFirstStage() <= 1) return 0;
	}
	if (U.getNoOfUncertainParameters() == 0) {
		if (!C_XQ.empty()) return 0;
		for (unsigned int p = 0; p < numPolicies; p++) {
			if (!C_XYQ[p].empty()) return 0;
		}
	}
	if (objectiveUnc) {
		if (!C_XQ.empty()) return 0;
		for (unsigned int p = 0; p < numPolicies; p++) {
			if (C_XYQ[p].size() != 1) return 0;
		}
	}
	bool check = 0;
	for (int x = 0; x < X.getTotalVarSize(); x++) {
		if (X.getVarColType(x) == 'B') {
			check = 1;
			if (!hasInteger) return 0;
		}
	}
	for (int y = 0; y < Y.getTotalVarSize(); y++) {
		if (Y.getVarColType(y) == 'B') {
			check = 1;
			if (!hasInteger) return 0;
		}
	}
	if (hasInteger && !check) return 0;
	for (const auto& con: C_X) {
		if (con.isEmpty()) return 0;
		if (con.existBilinearTerms()) return 0;
		if (con.existConstQTerms()) return 0;
	}
	for (const auto& con: C_XQ) {
		if (con.isEmpty()) return 0;
		if (!con.existBilinearTerms() && !con.existConstQTerms()) return 0;
	}
	for (const auto& con_K: C_XY) {
		for (const auto& con: con_K) {
			if (con.isEmpty()) return 0;
			if (con.existBilinearTerms()) return 0;
			if (con.existConstQTerms()) return 0;
		}
	}
	for (const auto& con_K: C_XYQ) {
		for (const auto& con: con_K) {
			if (con.isEmpty()) return 0;
			if (!con.existBilinearTerms() && !con.existConstQTerms()) return 0;
		}
	}
	for (const auto& con: B_X) {
		if (con.isEmpty()) return 0;
		if (con.existBilinearTerms()) return 0;
		if (con.existConstQTerms()) return 0;
	}
	for (const auto& con_K: B_Y) {
		for (const auto& con: con_K) {
			if (con.isEmpty()) return 0;
			if (con.existBilinearTerms()) return 0;
			if (con.existConstQTerms()) return 0;
		}
	}
	return 1;
}