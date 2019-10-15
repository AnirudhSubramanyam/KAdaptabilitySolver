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
#ifndef IND_TOOLS
#define IND_TOOLS

// Variable Indexing and Querying Class
// supporting up to 4 indices

#include <vector>
#include <string>
#include <map>
#include <tuple>
#include <cassert>

// v2: re-arranged order of indices in VarType
// v3: includes support for a fifth index
class VarInfo {

public: 
	
	// basic internal data structure. Map with key the variable name, and a tuple 
	// with fundamental sizes
	//         "VarName", < First, Last, numDimensions,  size1, size2, size3, size4, size5>
	//                  get   0      1      2              3      4      5      6      7
	std::map<std::string, std::tuple <int, int, int, int, int, int, int, int> > VarType; 
	int totalVars;  // stores the number of elements of these vars
	std::vector <double> UpperBound; // stores the upper bounds
	std::vector <double> LowerBound; // stores the lower bounds
	std::vector <double> ObjCoefficient; // stores the objective coefficients
	std::vector <char> ColumnType; // stores the type of variable: continuous 'C' or binary 'B'
	std::vector <bool> UndefinedVar;  // stores if a specific variable is undefined or not
	std::vector <int> UndefinedVarCount;  // stores the number of variables undefined up to a specific index

	bool isVarTypeConsistent(const int & ind1, const int & ind2, const int & ind3, const int & ind4, const int & ind5);

	void linTo1dIndex(const int & linIndex, const std::tuple <int, int, int, int, int, int, int, int> & VarType, int & ind1) const;

	void linTo2dIndex(const int & linIndex, const std::tuple <int, int, int, int, int, int, int, int> & VarType, int & ind1, int & ind2) const;

	void linTo3dIndex(const int & linIndex, const std::tuple <int, int, int, int, int, int, int, int> & VarType, int & ind1, int & ind2, int & ind3) const;

	void linTo4dIndex(const int & linIndex, const std::tuple <int, int, int, int, int, int, int, int> & VarType, int & ind1, int & ind2, int & ind3, int & ind4) const;

	void linTo5dIndex(const int & linIndex, const std::tuple <int, int, int, int, int, int, int, int> & VarType, int & ind1, int & ind2, int & ind3, int & ind4, int & ind5) const;
	
	void linToIndex(const int & linIndex, const std::tuple <int, int, int, int, int, int, int, int> & VarType, int & ind1, int & ind2, int & ind3, int & ind4, int & ind5) const;

	int xdIndicesToLinIndex(const std::tuple <int, int, int, int, int, int, int, int> & VarType, const int & ind1, const int & ind2, const int & ind3, const int & ind4, const int & ind5) const;

public:

	// Clear the object 
	void clear();

	// Add a new variable type
	void addVarType(const std::string & name, const char & type, const double & lb, const double & ub, const int & ind1, const int & ind2 = -1, const int & ind3 = -1, const int & ind4 = -1, const int & ind5 = -1);

	// set the given var to be undefined
	void setUndefinedVar(const int & index);

	// set the UB of the given var
	void setVarUB(const double & Ub, const int & index);

	// set the LB of the given var
	void setVarLB(const double & Lb, const int & index);

	// set the column type of the given var
	void setVarColType(const char & Ctype, const int & index);

	// set the objective coefficient of the given var
	void setVarObjCoeff(const double & Ocoeff, const int & index);

	// set the given var to be undefined
	void setUndefinedVar(const std::string & type, const int ind1 = -1, const int ind2 = -1, const int ind3 = -1, const int ind4 = -1, const int ind5 = -1);

	// set the UB of the given var
	void setVarUB(const double & Ub, const std::string & type, const int ind1 = -1, const int ind2 = -1, const int ind3 = -1, const int ind4 = -1, const int ind5 = -1);

	// set the LB of the given var
	void setVarLB(const double & Lb, const std::string & type, const int ind1 = -1, const int ind2 = -1, const int ind3 = -1, const int ind4 = -1, const int ind5 = -1);

	// set the column type of the given var
	void setVarColType(const char & Ctype, const std::string & type, const int ind1 = -1, const int ind2 = -1, const int ind3 = -1, const int ind4 = -1, const int ind5 = -1);

	// set the objective coefficient of the given var
	void setVarObjCoeff(const double & Ocoeff, const std::string & type, const int ind1 = -1, const int ind2 = -1, const int ind3 = -1, const int ind4 = -1, const int ind5 = -1);

	// set the vector of undefined vars - must be of size getTotalVarSize()
	// void setUndefinedVars(const std::vector <bool> & Undefined);

	// // set the vector with the UB of the vars - must be of size getTotalVarSize()
	// void setVarsUB(const std::vector <double> & Ubounds);

	// // set the vector with the LB of the vars - must be of size getTotalVarSize()
	// void setVarsLB(const std::vector <double> & Lbounds);

	// // set the vector with the column types - must be of size getTotalVarSize()
	// void setVarsColType(const std::vector <char> & Ctype);

	// // set the vector with the objective coefficients - must be of size getTotalVarSize()
	// void setVarsObjCoeff(const std::vector <double> & Ocoeff);






	// get the number of vars that this object holds
	int getTotalVarSize() const;

	// get the name and the indices, given the linear index
	void getVarInfo(const int & linIndex, std::string &type, int & ind1, int & ind2, int & ind3, int & ind4, int & ind5) const;

	// get the linear index, given the name and the indices
	int getVarLinIndex(const std::string & type, const int ind1 = -1, const int ind2 = -1, const int ind3 = -1, const int ind4 = -1, const int ind5 = -1) const;

	// get the number of variables of the requested type
	int getVarTypeSize(const std::string & name) const;

	// get the linear index of the first var of a given var type  (iterating tool)
	int getFirstOfVarType(const std::string & name) const;

	// get the linear index of the last var of a given var type (iterating tool)
	int getLastOfVarType(const std::string & name) const;

	// get the variable name, given the linear index as returned by getVarLinIndex()
	std::string getVarName(const int & index) const;

	// query if a var is undefined, given the linear index as returned by getVarLinIndex()
	bool isUndefVar(const int & index) const;

	// get the UB of a specific var, given the linear index as returned by getVarLinIndex()
	double getVarUB(const int & index) const;

	// get the LB of a specific var, given the linear index as returned by getVarLinIndex()
	double getVarLB(const int & index) const;

	// get the type of a specific var, given the linear index as returned by getVarLinIndex()
	char getVarColType(const int & index) const;

	// get the objective coefficient of a specific var, given the linear index as returned by getVarLinIndex()
	double getVarObjCoeff(const int & index) const;






	// get the true number of defined vars that this object holds
	int getTotalDefVarSize() const;

	// get the true linear index of defined var, given the original index as returned by getVarLinIndex()
	int getDefVarLinIndex(const int & index) const;

	// get the true linear index of defined var, given the name and the indices
	int getDefVarLinIndex(const std::string & type, const int ind1 = -1, const int ind2 = -1, const int ind3 = -1, const int ind4 = -1, const int ind5 = -1) const;

	// get the number of defined vars of the requested type
	int getDefVarTypeSize(const std::string & name) const;

	// get the true linear index of the first defined var of a given var type  (iterating tool)
	int getFirstDefOfVarType(const std::string & name) const;

	// get the true linear index of the last defined var of a given var type (iterating tool)
	int getLastDefOfVarType(const std::string & name) const;

};


#endif



