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


#include "uncertainty.hpp"
#include "Constants.h"
#include <cassert>
#include <iostream>
#include <algorithm> // std::transform
#include <numeric>

//---------------------------------------------------------------------------//

static inline void setCPXoptions(CPXENVptr& env) {
	CPXXsetdefaults(env);
	CPXXsetintparam(env, CPXPARAM_ScreenOutput,          CPX_OFF);
	CPXXsetintparam(env, CPXPARAM_Threads,               NUM_THREADS);
	CPXXsetintparam(env, CPXPARAM_Parallel,              CPX_PARALLEL_DETERMINISTIC);
	CPXXsetintparam(env, CPXPARAM_ClockType,             2);
	CPXXsetdblparam(env, CPXPARAM_TimeLimit,             TIME_LIMIT);
	CPXXsetdblparam(env, CPXPARAM_MIP_Limits_TreeMemory, MEMORY_LIMIT);
}

//---------------------------------------------------------------------------//

static inline void initCPX(CPXENVptr& env) {
	int status;

	assert(!env);

	// Try to initialize solver environment
	env = CPXXopenCPLEX(&status);
	if (!env) throw(EXCEPTION_CPXINIT);


	// Set solver options
	setCPXoptions(env);
}

//---------------------------------------------------------------------------//

UncertaintySet::UncertaintySet() {

	int status;

	// Initialize members
	N = 0;
	nominal.reserve(1000);
	nominal.emplace_back(0);
	low = nominal;
	high = nominal;

	polytope_W.emplace_back(std::vector<double>(1, 0));
	polytope_V.emplace_back(std::vector<double>(1, 0));
	polytope_h.emplace_back(0);
	polytope_sense.emplace_back('L');


	// initialize CPLEX objects
	env = NULL;
	lp  = NULL;
	initCPX(env);

	// Create problem object/model
	lp = CPXXcreateprob(env, &status, "UncertaintySet");
	if (!lp) throw(EXCEPTION_CPXINIT);

	// Create objective function variable
	const double obj    = 1;
	const double lo     = -CPX_INFBOUND;
	const double hi     = +CPX_INFBOUND;
	const char xctype   = 'C';
	const char *colname = "O";
	if (CPXXnewcols(env, lp, 1, &obj, &lo, &hi, &xctype, &colname)) {
		throw(EXCEPTION_CPXNEWCOLS);
	}

	// Make it linear
	CPXXchgprobtype(env, lp, CPXPROB_LP);

	// Make it a maximization problem
	CPXXchgobjsen(env, lp, CPX_MAX);
}


//---------------------------------------------------------------------------//


UncertaintySet::UncertaintySet(const UncertaintySet& U) :
	N(U.N),
	nominal(U.nominal),
	low(U.low),
	high(U.high),
	polytope_W(U.polytope_W),
	polytope_V(U.polytope_V),
	polytope_h(U.polytope_h),
	polytope_sense(U.polytope_sense)
{
	int status;

	// Try to initialize solver environment
	env = NULL;
	lp  = NULL;
	initCPX(env);


	// Copy solver model object
	lp = U.getLPObject(env, &status);
	if (status) {
		throw(EXCEPTION_CPXINIT);
	}
}

//---------------------------------------------------------------------------//


UncertaintySet& UncertaintySet::operator=(const UncertaintySet& U) {
	if (this == &U) {
		return *this;
	}

	N = U.N;
	nominal = U.nominal;
	low = U.low;
	high = U.high;
	polytope_W = U.polytope_W;
	polytope_V = U.polytope_V;
	polytope_h = U.polytope_h;
	polytope_sense = U.polytope_sense;

	int status;

	// env, lp are already allocated
	assert(env);
	assert(lp);

	// delete lp
	if (lp) if (CPXXfreeprob (env, &lp)) {
		throw(EXCEPTION_CPXEXIT);
	}

	// Copy solver model object
	lp = U.getLPObject(env, &status);
	if (status) {
		throw(EXCEPTION_CPXINIT);
	}

	return *this;
}




//---------------------------------------------------------------------------//


UncertaintySet::~UncertaintySet() noexcept(false) {
	this->clear();

	assert(env);
	assert(lp);

	// Free problem object memory
	if (lp) if (CPXXfreeprob(env, &lp)) {
		throw(EXCEPTION_CPXEXIT);
	}
	

	// Free solver memory
	if (env) if (CPXXcloseCPLEX(&env)) {
		char  errmsg[CPXMESSAGEBUFSIZE];
		int status = 0;
		CPXXgeterrorstring (env, status, errmsg);
		std::cerr << errmsg << std::endl;
		throw(EXCEPTION_CPXEXIT);
	}
}


//---------------------------------------------------------------------------//

void UncertaintySet::clear() {
	assert(env);
	assert(lp);

	// get # of cols
	const int cur_numcols = CPXXgetnumcols(env, lp);

	// delete all columns except the objective
	assert(cur_numcols >= 1);
	if (cur_numcols > 1) {
		CPXXdelcols(env, lp, 1, CPXXgetnumcols(env, lp) - 1);
	}
	assert(CPXXgetnumcols(env, lp) == 1);

	// set default options
	setCPXoptions(env);
	
	// Zero all members
	N    = 0;
	nominal.assign(1, 0.0);
	low  = nominal;
	high = nominal;
	polytope_W.assign(1, std::vector<double>(1, 0));
	polytope_V.assign(1, std::vector<double>(1, 0));
	polytope_h.assign(1, 0.0);
	polytope_sense.assign(1, 'L');
}

//---------------------------------------------------------------------------//




void UncertaintySet::addParam(const double nom, const double lo, const double hi) {
	if (nom < lo)
		std::cerr << "Warning. Nominal value of uncertain parameter is lower than the lower bound.\n";

	if (nom > hi)
		std::cerr << "Warning. Nominal value of uncertain parameter is greater than the upper bound.\n";

	N++;
	nominal.emplace_back(nom);
	low.emplace_back(lo);
	high.emplace_back(hi);

	assert(nominal.size() == low.size());
	assert(high.size() == low.size());
	assert(1 + N == (int)low.size());


	// update matrices
	for (unsigned i = 0; i < polytope_W.size(); ++i) {
		polytope_W[i].emplace_back(0);
		assert((int)polytope_W[i].size() == 1 + N);
	}

	// Upper bound
	polytope_W.emplace_back(std::vector<double>(1 + N, 0));
	polytope_V.emplace_back(std::vector<double>(1, 0));
	polytope_W.back().back() = 1;
	polytope_sense.emplace_back('L');
	polytope_h.emplace_back(hi);

	// Lower bound
	polytope_W.emplace_back(std::vector<double>(1 + N, 0));
	polytope_V.emplace_back(std::vector<double>(1, 0));
	polytope_W.back().back() = 1;
	polytope_sense.emplace_back('G');
	polytope_h.emplace_back(lo);

	// Matrix sizes must match
	assert(polytope_h.size() == polytope_sense.size());
	assert(polytope_W.size() == polytope_sense.size());
	assert(polytope_V.size() == polytope_sense.size());

	// Add variable to solver model object
	double obj = 0;
	const char xctype = 'C';
	const char *colname = std::string("q(" + std::to_string(N) + ")").c_str();
	if (CPXXnewcols(env, lp, 1, &obj, &lo, &hi, &xctype, &colname)) {
		throw(EXCEPTION_CPXNEWCOLS);
	}
	assert(CPXXgetnumcols(env, lp) == 1 + N);

	return;

}

//---------------------------------------------------------------------------//


void UncertaintySet::addFacet(const std::vector<std::pair<int, double> >& data, const char sense, const double rhs) {

	for (unsigned i = 0; i < data.size(); ++i) {
		if(data[i].first < 1 || data[i].first > N) {
			std::cerr << "Warning: Facet description contains invalid indices of uncertain parameters. ";
			std::cerr << "Will not add facet.\n";
			return;
		}
	}

	// update matrices
	polytope_W.emplace_back(std::vector<double>(1 + N, 0));
	polytope_V.emplace_back(std::vector<double>(1, 0));
	polytope_h.emplace_back(rhs);
	polytope_sense.emplace_back(sense);

	for (const auto& d : data)
		polytope_W.back().at(d.first) = d.second;


	// Matrix sizes must match
	assert(polytope_h.size() == polytope_sense.size());
	assert(polytope_W.size() == polytope_sense.size());
	assert(polytope_V.size() == polytope_sense.size());

	// Update solver model object
	std::vector<int> indices;
	std::vector<double> coeffs;

	for (const auto& d : data) {
		indices.emplace_back(d.first);
		coeffs.emplace_back(d.second);
	}

	const std::string rname = "row(" + std::to_string(CPXXgetnumrows(env, lp) + 1) + ")";
	const std::vector<const char*> rowname({rname.c_str()});
	CPXNNZ rmatbeg = 0;

	if (CPXXaddrows(env, lp, 0, 1, indices.size(), &rhs, &sense, &rmatbeg, &indices[0], &coeffs[0], nullptr, &rowname[0]))
		throw(EXCEPTION_CPXNEWROWS);

	return;
}


//---------------------------------------------------------------------------//


int UncertaintySet::addVariables_DualVars(CPXCENVptr env_, CPXLPptr lp_, const std::string& dualName) const {
	/* Assumes that polytope_W, polytope_V, polytope_h and polytope_sense
	 * have all been allocated and have the same number of rows.
	 *
	 * Assumes that these matrices/vectors are all 1-indexed.
	 */
	const CPXDIM ccnt = polytope_sense.size() - 1;
	const std::vector<double> obj(ccnt, 0);
	const std::vector<char> xctype(ccnt, 'C');
	std::vector<double> lb, ub;
	std::vector<std::string> cname;
	for (int s = 1; s <= ccnt; s++) {
		cname.emplace_back("UncSetDual(" + (std::string)dualName + "," + std::to_string(s) + ")");
		switch (polytope_sense.at(s)){
			case 'L': ub.push_back(+CPX_INFBOUND); lb.push_back(0);             break;
			case 'G': ub.push_back(0);             lb.push_back(-CPX_INFBOUND); break;
			case 'E': ub.push_back(+CPX_INFBOUND); lb.push_back(-CPX_INFBOUND); break;
		}
	}

	/* ADD COLUMNS TO CPLEX */
	std::vector<const char*> colname; for (size_t i = 0; i<cname.size(); i++) colname.push_back(cname[i].c_str());
	if (CPXXnewcols(env_, lp_, ccnt, &obj[0], &lb[0], &ub[0], &xctype[0], &colname[0]))
		throw(EXCEPTION_CPXNEWCOLS);

	return ccnt;
}


//---------------------------------------------------------------------------//


double UncertaintySet::max(const std::vector<std::pair<int, double> >& data, std::vector<double>& result) const {
	
	if(N == 0) return 0;

	int probtype = CPXXgetprobtype(env, lp);
	std::vector<int> indices;
	std::vector<double> coeffs;

	for (const auto& d : data) {
		indices.emplace_back(d.first);
		coeffs.emplace_back(d.second);
	}


	// return value
	double val = 0; result.assign(1 + N, 0);

	// Set objective function of LP over the uncertainty set
	std::vector<CPXDIM> CplexIndices(1 + N, 0);
	std::vector<double> values(1 + N, 0);
	std::iota(CplexIndices.begin(), CplexIndices.end(), 0);
	for (unsigned i = 0; i < indices.size(); ++i) {
		values.at(indices.at(i)) = coeffs.at(i);
	}

	// Replace current objective function
	CPXXchgobj(env, lp, (CPXDIM)CplexIndices.size(), &CplexIndices[0], &values[0]);

	// Set objective sense
	CPXXchgobjsen(env, lp, CPX_MAX);

	
	// Solve MILP/LP
	switch(probtype) {
		case CPXPROB_LP: CPXXlpopt(env, lp); assert(CPXXgetstat(env, lp) == CPX_STAT_OPTIMAL); break;
		case CPXPROB_MILP: CPXXmipopt(env, lp); assert(CPXXgetstat(env, lp) == CPXMIP_OPTIMAL || CPXXgetstat(env, lp) == CPXMIP_OPTIMAL_TOL); break;
	}


	// Get optimal value
	CPXXgetobjval(env, lp, &val);
	
	// Get solution vector
	CPXXgetx(env, lp, &result[0], 0, N);
	

	return val;
}

//---------------------------------------------------------------------------//



double UncertaintySet::min(const std::vector<std::pair<int, double> >& data, std::vector<double>& result) const {

	auto d(data);
	for (unsigned i = 0; i < data.size(); ++i) d[i].second *= -1.0;

	return (-1.0*UncertaintySet::max(d, result));
}

//---------------------------------------------------------------------------//

double UncertaintySet::getMaximumValue(const std::vector<int>& ind, const std::vector<double>& coef) const {
	std::vector<std::pair<int, double> > input;
	input.reserve(ind.size());

	std::transform(ind.begin(), ind.end(), coef.begin(), std::back_inserter(input),
		[](int i, double c) {
			return std::make_pair(i, c);
		}
	);

	std::vector<double> result;

	return UncertaintySet::max(input, result);
}


//---------------------------------------------------------------------------//











