/***************************************************************************************/
/*                                                                                     */
/*  Copyright 2018 by Nikolaos Lappas and Chrysanthos Gounaris                         */
/*                                                                                     */
/*  Licensed under the FreeBSD License (the "License").                                */
/*  You may not use this file except in compliance with the License.                   */
/*  You may obtain a copy of the License at                                            */
/*                                                                                     */
/*  https://www.freebsd.org/copyright/freebsd-license.html                             */
/*                                                                                     */
/***************************************************************************************/
// Variable Indexing and Querying Class
// supporting up to 4 indices

#include "indexingTools.h"
#include <iostream>
#include <algorithm>
//################################
static inline void abortProgram(const std::string &customMessage = "", const std::string &doNotUseThis = ""){
#define abortProgram(...) abortProgram(std::string(__FILE__) + " -- Line " + std::to_string(__LINE__), std::string(__VA_ARGS__))
	std::cerr << " == ABORTING PROGRAM >>> " << (doNotUseThis.empty() ? customMessage : doNotUseThis) << " == " << std::endl;
	abort();
	return;
}
//################################



bool VarInfo::isVarTypeConsistent(const int & ind1, const int & ind2, const int & ind3, const int & ind4, const int & ind5){
	if(ind1 < 0) return false;
	if(ind2 > -1 && ind1 <0 ) return false;
	if((ind3 > -1 && ind2 <0) || (ind3 > -1 && ind1 <0)) return false;
	if((ind4 > -1 && ind3 < 0) || (ind4 > -1 && ind2 <0) || (ind4 > -1 && ind1 <0)) return false;
	if((ind5 > -1 && ind4 < 0) || (ind5 > -1 && ind3 <0) || (ind5 > -1 && ind2 <0) || (ind5 > -1 && ind1 <0)) return false;
	return true;
}


void VarInfo::linTo1dIndex(const int & linIndex, const std::tuple <int, int, int, int, int, int, int, int> & VarType, int & ind1) const {
	int diff = linIndex - std::get<0>(VarType); assert(diff >= 0); assert(diff <= std::get<1> (VarType));
	ind1 = diff;
}

void VarInfo::linTo2dIndex(const int & linIndex, const std::tuple <int, int, int, int, int, int, int, int> & VarType, int & ind1, int & ind2) const {
	int diff = linIndex - std::get<0>(VarType); assert(diff >= 0); assert(diff <= std::get<1> (VarType));
	ind1 = diff / (std::get<4>(VarType) ); diff -= ind1 * std::get<4>(VarType);
	ind2 = diff;
}

void VarInfo::linTo3dIndex(const int & linIndex, const std::tuple <int, int, int, int, int, int, int, int> & VarType, int & ind1, int & ind2, int & ind3) const {
	int diff = linIndex - std::get<0>(VarType); assert(diff >= 0); assert(diff <= std::get<1> (VarType));
	ind1 = diff / (std::get<4>(VarType) * std::get<5>(VarType) ); diff -= ind1 * std::get<4>(VarType) * std::get<5>(VarType);
	ind2 = diff /  std::get<5>(VarType); diff -= ind2 * std::get<5>(VarType);
	ind3 = diff;
}

void VarInfo::linTo4dIndex(const int & linIndex, const std::tuple <int, int, int, int, int, int, int, int> & VarType, int & ind1, int & ind2, int & ind3, int & ind4) const {
	int diff = linIndex - std::get<0>(VarType); assert(diff >= 0); assert(diff <= std::get<1> (VarType));
	ind1 = diff / (std::get<4>(VarType) * std::get<5>(VarType) * std::get<6>(VarType)  ); diff -= ind1 * std::get<4>(VarType) * std::get<5>(VarType) * std::get<6>(VarType);
	ind2 = diff / (std::get<5>(VarType) * std::get<6>(VarType) ); diff -= ind2 * std::get<5>(VarType) * std::get<6>(VarType);
	ind3 = diff /  std::get<6>(VarType); diff -= ind3 * std::get<6>(VarType);
	ind4 = diff;
}

void VarInfo::linTo5dIndex(const int & linIndex, const std::tuple <int, int, int, int, int, int, int, int> & VarType, int & ind1, int & ind2, int & ind3, int & ind4, int & ind5) const {
	int diff = linIndex - std::get<0>(VarType); assert(diff >= 0); assert(diff <= std::get<1> (VarType));
	ind1 = diff / (std::get<4>(VarType) * std::get<5>(VarType) * std::get<6>(VarType) * std::get<7>(VarType)  ); diff -= ind1 * std::get<4>(VarType) * std::get<5>(VarType) * std::get<6>(VarType) * std::get<7>(VarType);
	ind2 = diff / (std::get<5>(VarType) * std::get<6>(VarType) * std::get<7>(VarType) ); diff -= ind2 * std::get<5>(VarType) * std::get<6>(VarType) * std::get<7>(VarType);
	ind3 = diff / (std::get<6>(VarType) * std::get<7>(VarType) ); diff -= ind3 * std::get<6>(VarType) * std::get<7>(VarType);
	ind4 = diff /  std::get<7>(VarType); diff -= ind4 * std::get<7>(VarType);
	ind5 = diff;
}

void VarInfo::linToIndex(const int & linIndex, const std::tuple <int, int, int, int, int, int, int, int> & VarType, int & ind1, int & ind2, int & ind3, int & ind4, int & ind5) const {
	int dimensions = std::get<2> (VarType);
	if(dimensions == 1) {
		linTo1dIndex(linIndex, VarType, ind1);
		ind2 = -1;
		ind3 = -1;
		ind4 = -1;
		ind5 = -1;
	}
	else if(dimensions == 2) {
		linTo2dIndex(linIndex, VarType, ind1, ind2);
		ind3 = -1;
		ind4 = -1;
		ind5 = -1;
	}
	else if(dimensions == 3) {
	  linTo3dIndex(linIndex, VarType, ind1, ind2, ind3);
	  ind4 = -1;
	  ind5 = -1;
	}
	else if(dimensions == 4) {
		linTo4dIndex(linIndex, VarType, ind1, ind2, ind3, ind4);
		ind5 = -1;
	}
	else if(dimensions == 5) {
		linTo5dIndex(linIndex, VarType, ind1, ind2, ind3, ind4, ind5);
	}
	else {
		abortProgram("invalid number of dimensions for this Variable");
	}
}

int VarInfo::xdIndicesToLinIndex(const std::tuple <int, int, int, int, int, int, int, int> & VarType, const int & ind1, const int & ind2, const int & ind3, const int & ind4, const int & ind5) const {
	int index;
	int dimensions = std::get<2> (VarType);
	if(dimensions == 1) {
		assert(ind1 > -1);
		assert(ind2 < 0 && ind3 < 0 && ind4 < 0 && ind5 < 0);
		index = ind1;
	}
	else if(dimensions == 2) {
		assert(ind1 > -1 && ind2 > -1);
		assert(ind3 < 0 && ind4 < 0 && ind5 < 0);
		index = ind1 * std::get<4>(VarType) + ind2;
	}
	else if(dimensions == 3) {
		assert(ind1 > -1 && ind2 > -1 && ind3 > -1);
		assert(ind4 < 0 && ind5 < 0);
	  index = ind1 * std::get<4>(VarType) * std::get<5>(VarType) + ind2 * std::get<5>(VarType) + ind3;
  }
  else if(dimensions == 4) {
    assert(ind1 > -1 && ind2 > -1 && ind3 > -1 && ind4 > -1);
    assert(ind5 < 0);
    index = ind1 * std::get<4>(VarType) * std::get<5>(VarType) * std::get<6>(VarType) + ind2 * std::get<5>(VarType) * std::get<6>(VarType) + ind3 * std::get<6>(VarType) + ind4;
  }
  else if(dimensions == 5) {
    assert(ind1 > -1 && ind2 > -1 && ind3 > -1 && ind4 > -1 && ind5 > -1);
    index = ind1 * std::get<4>(VarType) * std::get<5>(VarType) * std::get<6>(VarType) * std::get<7>(VarType) + ind2 * std::get<5>(VarType) * std::get<6>(VarType) * std::get<7>(VarType) + ind3 * std::get<6>(VarType) * std::get<7>(VarType) + ind4 * std::get<7>(VarType) + ind5;
  }
  else {
    abortProgram("invalid number of dimensions for this Variable (xdIndicesToLinIndex)");
  }

	index += std::get<0>(VarType);
	assert(index <= std::get<1>(VarType));
	return index;
}

// Clear the object
void VarInfo::clear() {
	VarType.clear();
	totalVars = 0;
	UpperBound.clear();
	LowerBound.clear();
	ColumnType.clear();
	ObjCoefficient.clear();
	UndefinedVar.clear();
	UndefinedVarCount.clear();
}

// Add a new variable type
void VarInfo::addVarType(const std::string & name, const char & type, const double & lb, const double & ub, const int & ind1, const int & ind2, const int & ind3, const int & ind4, const int & ind5) {
	int first = 0;
	if(VarType.size() > 0) {
		first = totalVars;
	}
	if (type != 'B' && type != 'C') abortProgram("Column type must be 'B' or 'C'");
	// Calculate the size
	int size = 0;
	int numDimensions = 0;
	if(ind1 > -1) {
		size = ind1;
		numDimensions++;
		if(ind2 > -1 ) {
			size *= ind2;
			numDimensions++;
			if(ind3 > -1) {
				size *= ind3;
				numDimensions++;
				if(ind4 > -1 ) {
					size *= ind4;
					numDimensions++;
					if(ind5 > -1 ) {
						size *= ind5;
						numDimensions++;
					}
				}
			}
		}
	}
	totalVars = first + size;
	int last = first + size -1;
	auto temp = std::make_tuple(first, last, numDimensions, ind1, ind2, ind3, ind4, ind5);
	VarType.emplace(name, temp);
	UpperBound.resize(totalVars, ub);
	LowerBound.resize(totalVars, lb);
	ColumnType.resize(totalVars, type);
	ObjCoefficient.resize(totalVars, 0.0);
	UndefinedVar.resize(totalVars, 0);
	UndefinedVarCount.resize(totalVars, (UndefinedVarCount.empty() ? 0 : UndefinedVarCount.back()) );

	if(isVarTypeConsistent(ind1, ind2, ind3, ind4, ind5) != true) abortProgram("Var Type not valid");
}

// // set the vector of undefined vars
// void VarInfo::setUndefinedVars(const std::vector <bool> & Undefined) {
// 	UndefinedVar.clear(); UndefinedVar.resize(totalVars, 0);
// 	assert((int)Undefined.size() == totalVars);
// 	UndefinedVar = Undefined;
// 	UndefinedVarCount.resize(totalVars, 0);
// 	if (totalVars > 0) {
// 		UndefinedVarCount[0] = static_cast<int>(UndefinedVar[0]);
// 	}
// 	for (int i = 1; i < totalVars; i++) {
// 		UndefinedVarCount[i] = UndefinedVarCount[i-1] + static_cast<int>(UndefinedVar[i]);
// 		assert(std::accumulate(UndefinedVar.begin(), UndefinedVar.begin() + i + 1, int(0)) == UndefinedVarCount[i]);
// 	}
// }

// // set the vector with the LB of the vars
// void VarInfo::setVarsLB(const std::vector <double> & Lbounds) {
// 	LowerBound.clear(); LowerBound.resize(totalVars, 0.0);
// 	assert((int)Lbounds.size() == totalVars);
// 	LowerBound = Lbounds;
// }

// // set the vector with the UB of the vars
// void VarInfo::setVarsUB(const std::vector <double> & Ubounds) {
// 	UpperBound.clear(); UpperBound.resize(totalVars, 0.0);
// 	assert((int)Ubounds.size() == totalVars);
// 	UpperBound = Ubounds;
// }

// // set the vector with the column types
// void VarInfo::setVarsColType(const std::vector <char> & Ctype) {
// 	ColumnType.clear(); ColumnType.resize(totalVars, 0.0);
// 	assert((int)Ctype.size() == totalVars);
// 	ColumnType = Ctype;
// }

// // set the vector with the objective coefficients
// void VarInfo::setVarsObjCoeff(const std::vector <double> & Ocoeff) {
// 	ObjCoefficient.clear(); ObjCoefficient.resize(totalVars, 0.0);
// 	assert((int)Ocoeff.size() == totalVars);
// 	ObjCoefficient = Ocoeff;
// }

// set the given var to be undefined
void VarInfo::setUndefinedVar(const int & index) {
	assert(index > -1 && index < totalVars);
	if (!UndefinedVar[index]) {
		UndefinedVar[index] = 1;
		for (int i = index; i < totalVars; i++) UndefinedVarCount[i]++;
	}
}

// set the UB of the given var
void VarInfo::setVarUB(const double & Ub, const int & index) {
	assert(index > -1 && index < totalVars);
	UpperBound[index] = Ub;
}

// set the LB of the given var
void VarInfo::setVarLB(const double & Lb, const int & index) {
	assert(index > -1 && index < totalVars);
	LowerBound[index] = Lb;
}

// set the column type of the given var
void VarInfo::setVarColType(const char & Ctype, const int & index) {
	assert(index > -1 && index < totalVars);
	ColumnType[index] = Ctype;
}

// set the objective coefficient of the given var
void VarInfo::setVarObjCoeff(const double & Ocoeff, const int & index) {
	assert(index > -1 && index < totalVars);
	ObjCoefficient[index] = Ocoeff;
}

// set the given var to be undefined
void VarInfo::setUndefinedVar(const std::string & type, const int ind1, const int ind2, const int ind3, const int ind4, const int ind5) {
	setUndefinedVar(getVarLinIndex(type, ind1, ind2, ind3, ind4, ind5));
}

// set the UB of the given var
void VarInfo::setVarUB(const double & Ub, const std::string & type, const int ind1, const int ind2, const int ind3, const int ind4, const int ind5) {
	setVarUB(Ub, getVarLinIndex(type, ind1, ind2, ind3, ind4, ind5));
}

// set the LB of the given var
void VarInfo::setVarLB(const double & Lb, const std::string & type, const int ind1, const int ind2, const int ind3, const int ind4, const int ind5) {
	setVarLB(Lb, getVarLinIndex(type, ind1, ind2, ind3, ind4, ind5));
}

// set the column type of the given var
void VarInfo::setVarColType(const char & Ctype, const std::string & type, const int ind1, const int ind2, const int ind3, const int ind4, const int ind5) {
	setVarColType(Ctype, getVarLinIndex(type, ind1, ind2, ind3, ind4, ind5));
}

// set the objective coefficient of the given var
void VarInfo::setVarObjCoeff(const double & Ocoeff, const std::string & type, const int ind1, const int ind2, const int ind3, const int ind4, const int ind5) {
	setVarObjCoeff(Ocoeff, getVarLinIndex(type, ind1, ind2, ind3, ind4, ind5));
}


// get the number of vars that this object holds
int VarInfo::getTotalVarSize() const {
	return totalVars;
}

// get the number of defined vars that this object holds
int VarInfo::getTotalDefVarSize() const {
	return totalVars - (totalVars > 0 ? UndefinedVarCount.back() : 0);
}

// get the name and the indices, given the linear index
void VarInfo::getVarInfo(const int & linIndex, std::string &type, int & ind1, int & ind2, int & ind3, int & ind4, int & ind5) const {
	bool found = false;
	for(auto it = VarType.cbegin(); it != VarType.cend(); ++it) {
		int firstOf = std::get<0>((*it).second);
		int lastOf  = std::get<1>((*it).second);
		if(linIndex >= firstOf && linIndex <= lastOf) {
			found = true;
			type = (*it).first;
			linToIndex(linIndex, (*it).second, ind1, ind2, ind3, ind4, ind5);
		}
	}
	if(!found) abortProgram("Variable not found (getVarInfo)");
}

// get the linear index, given the name and the indices
int VarInfo::getVarLinIndex(const std::string & type, const int ind1, const int ind2, const int ind3, const int ind4, const int ind5) const {
	if(VarType.find(type) == VarType.end()) { // key not found
		abortProgram("Unknown Variable key (getVarLinIndex)");
		return 0;
	}
	else {
		return xdIndicesToLinIndex(VarType.at(type), ind1, ind2, ind3, ind4, ind5);
	}
}

// get the linear index of defined variable, given the name and the indices
int VarInfo::getDefVarLinIndex(const std::string & type, const int ind1, const int ind2, const int ind3, const int ind4, const int ind5) const {
	if(VarType.find(type) == VarType.end()) { // key not found
		abortProgram("Unknown Variable key (getDefVarLinIndex)");
	}
	const int linIndex = xdIndicesToLinIndex(VarType.at(type), ind1, ind2, ind3, ind4, ind5);
	assert(!UndefinedVar[linIndex]);
	return linIndex - UndefinedVarCount[linIndex];
}

// get the linear index of defined variable, given the original index
int VarInfo::getDefVarLinIndex(const int & index) const {
	assert(index > -1 && index < totalVars);
	return index - UndefinedVarCount[index];
}

// get the number of variables of the requested type
int VarInfo::getVarTypeSize(const std::string & name) const {
	if(VarType.find(name) == VarType.end()) {
		abortProgram("Invalid key (getVarTypeSize)");
	}

	const int elements = std::get<1>(VarType.at(name)) - std::get<0>(VarType.at(name)) + 1;
	return elements;
}

// get the number of defined variables of the requested type
int VarInfo::getDefVarTypeSize(const std::string & name) const {
	if(VarType.find(name) == VarType.end()) {
		abortProgram("Invalid key (getVarTypeSize)");
	}
	const int last  = std::get<1>(VarType.at(name));
	const int first = std::get<0>(VarType.at(name));
	int elements = last - first + 1 - UndefinedVarCount[last];
	if (first >= 1) {
		elements += UndefinedVarCount[first - 1];
	}
	return elements;
}

// get the linear index of the first var of a given var type  (iterating tool)
int VarInfo::getFirstOfVarType(const std::string & name) const {
	if(VarType.find(name) == VarType.end()) abortProgram("invalid var key (getFirstOfVarType)");
	return std::get<0> (VarType.at(name));
}

// get the linear index of the last var of a given var type (iterating tool)
int VarInfo::getLastOfVarType(const std::string & name) const {
	if(VarType.find(name) == VarType.end()) abortProgram("invalid var key (getLastOfVarType)");
	return std::get<1> (VarType.at(name));
}

// get the linear index of the first defined var of a given var type  (iterating tool)
int VarInfo::getFirstDefOfVarType(const std::string & name) const {
	if(VarType.find(name) == VarType.end()) abortProgram("invalid var key (getFirstOfDefVarType)");
	int first = std::get<0> (VarType.at(name));
	const int last = std::get<1> (VarType.at(name));
	while (UndefinedVar[first] && first <= last) first++;
	if (first > last) abortProgram("All variables of type are undefined (getFirstOfDefVarType)");
	return getDefVarLinIndex(first);
}

// get the linear index of the last defined var of a given var type (iterating tool)
int VarInfo::getLastDefOfVarType(const std::string & name) const {
	if(VarType.find(name) == VarType.end()) abortProgram("invalid var key (getLastOfLastVarType)");
	int last = std::get<1> (VarType.at(name));
	const int first = std::get<0> (VarType.at(name));
	while(UndefinedVar[last] && last >= first) last--;
	if (last < first) abortProgram("All variables of type are undefined (getLastOfDefVarType)");
	return getDefVarLinIndex(last);
}

// get the variable name, given the linear index
std::string VarInfo::getVarName(const int & index) const {
	bool found = false;
	std::string type;
	for(auto it = VarType.cbegin(); it != VarType.cend(); ++it) {
		int firstOf = std::get<0>((*it).second);
		int lastOf  = std::get<1>((*it).second);
		if(index >= firstOf && index <= lastOf) {
			found = true;
			type = (*it).first;
			return type;
		}
	}
	if(!found) abortProgram("Variable not found (getVarName)");
	return "";
}

// query if a var is undefined
bool VarInfo::isUndefVar(const int & index) const {
	assert(index > -1 && index < totalVars);
	if(UndefinedVar.size() == 0) {
		return false; // default is 0
	}
	else {
		return UndefinedVar[index];
	}
}

// get the UB of a specific var
double VarInfo::getVarUB(const int & index) const {
	assert(index > -1 && index < totalVars);
	if(UpperBound.size() == 0) {
		return 0.0; // default bound is 0
	}
	else {
		return UpperBound[index];
	}
}


// get the LB of a specific var
double VarInfo::getVarLB(const int & index) const {
	assert(index > -1 && index < totalVars);
	if(LowerBound.size() == 0) {
		return 0.0; // default bound is 0
	}
	else {
		return LowerBound[index];
	}
}


// get the type of a specific var, given the linear index as returned by getVarLinIndex()
char VarInfo::getVarColType(const int & index) const {
	assert(index > -1 && index < totalVars);
	if(ColumnType.size() == 0) {
		return 'C'; // default type is 'C'
	}
	else {
		return ColumnType[index];
	}
}

// get the LB of a specific var
double VarInfo::getVarObjCoeff(const int & index) const {
	assert(index > -1 && index < totalVars);
	if(ObjCoefficient.size() == 0) {
		return 0.0; // default coefficient is 0
	}
	else {
		return ObjCoefficient[index];
	}
}





// void addToCPX(const VarInfo& X) {

// }




// int main(int argc, char** argv) {

// 	VarInfo X;
// 	std::vector<bool> UV;
// 	X.addVarType("z", 1); UV.emplace_back(0);
// 	X.addVarType("x", 6, 6);
// 	for (int i = 0; i <= 5; i ++) for (int j = 0; j <= 5; j++) {
// 		UV.emplace_back(0);
// 		X.setVarUB(X.getVarLinIndex("x",i,j), 1.0);
// 		if (i == j || (i == 5) || (j==0) || (i==0 && j==5)) UV.back() = 1;
// 	}
// 	X.setUndefinedVars(UV);
	

// 	for (int i = 0; i <= 5; i ++) for (int j = 0; j <= 5; j++) if (!X.isUndefVar(X.getVarLinIndex("x", i, j))) {
// 		std::cout << "(" << i << "," << j << ") = " << X.getDefVarLinIndex("x", i, j) << "\n";
// 	}

// 	std::cout << X.getDefVarLinIndex("z", 0) << std::endl;

// 	for (int i = 0; i < X.getTotalVarSize(); i++) if (!X.isUndefVar(i)) {
// 		std::cout << X.getVarName(i) << "(" << X.getDefVarLinIndex(i) << ") in [" << X.getVarLB(i) << "," << X.getVarUB(i) <<"]\n";
// 	}

// 	X.addVarType("y", 0);

// 	std::cout << X.getFirstOfVarType("x") << " " << X.getFirstOfVarType("z") << "\n";
// 	std::cout << X.getFirstDefOfVarType("x") << " " << X.getFirstDefOfVarType("z") << "\n";
// 	std::cout << X.getLastOfVarType("x") << " " << X.getLastOfVarType("z") << "\n";
// 	std::cout << X.getLastDefOfVarType("x") << " " << X.getLastDefOfVarType("z") << "\n";
// 	std::cout << X.getDefVarTypeSize("x") << " " << X.getDefVarTypeSize("z") << " " << X.getDefVarTypeSize("y") << "\n";
// 	std::cout << X.getTotalVarSize() << " " << X.getTotalDefVarSize() << "\n";

// 	return 0;
// }


