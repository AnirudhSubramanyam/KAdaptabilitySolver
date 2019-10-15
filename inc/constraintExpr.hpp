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


#ifndef CONSTRAINT_EXPRESSION_HPP
#define CONSTRAINT_EXPRESSION_HPP

#include <string>
#include <cstring>
#include <vector>
#include <set>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <ilcplex/cplexx.h>
#include "uncertainty.hpp"


#define CONSTRAINT_EXPRESSION_OUTPUTLEVEL 0


typedef std::vector<int> indices_t;
typedef std::vector<double> coefficients_t;
typedef std::pair<int, int> pair_t;
typedef std::vector<pair_t> pairIndices_t;




/*---------------------------------------------------------
   This class supports constraints of the following type:
   
   cx + bq + xSq (sense) d
  
   x are decision variables (varIndices)
   c are const doubles with same size as x (varCoeffs)
   q are uncertain parameters (paramIndices)
   b are const doubles with same size as q (paramCoeffs)
   S is a (nx times nq) matrix (sMatrix)
   d is the right hand side double value (rhs)
   
 *---------------------------------------------------------*/
class ConstraintExpression {
	friend class KAdaptableExpression;
private:
	char sense;/* Default: 'L', i.e., <= */
	double rhs;/* Default: 0 */
	std::string name;
	indices_t varIndices;
	coefficients_t varCoeffs;
	

	/* Following are required only for Robust counterparts */
	pairIndices_t bilinearIndices; /* pairs of (indices of) bilinear terms (xi,qj), i.e., products of decision variables with uncertain parameters */
	coefficients_t bilinearCoeffs; /* S(xi,qj) := coefficient of term xi*qj */
	indices_t paramIndices;
	coefficients_t paramCoeffs;

	inline bool isConsistentWithDesign() const {
		if (std::set<int>(varIndices.begin(), varIndices.end()).size()     != varIndices.size())   return false;
		if (std::set<int>(paramIndices.begin(), paramIndices.end()).size() != paramIndices.size()) return false;
		for (unsigned i = 0; i < bilinearIndices.size(); ++i)
			for (unsigned j = i + 1; j < bilinearIndices.size(); ++j)
				if (bilinearIndices[i] == bilinearIndices[j])return false;
		if (varIndices.size()   != varCoeffs.size())         return false;
		if (paramIndices.size() != paramCoeffs.size())       return false;
		if (bilinearIndices.size() != bilinearCoeffs.size()) return false;
		if (sense != 'L' && sense != 'E' && sense != 'G')    return false;
		/* For every bilinear term xi*qj, ci*xi and bj*qj must exist (possibly with ci=0,bj=0) */
		/* Inverse need not be true */
		for (const auto& xq : bilinearIndices) {
			if (std::find(varIndices.begin(), varIndices.end(), xq.first) == varIndices.end())        return false;
			if (std::find(paramIndices.begin(), paramIndices.end(), xq.second) == paramIndices.end()) return false;
		}
		return true;
	}
	/* Add the term cx (vec1 = varIndices, coeffs = paramCoeffs)
	 *           or bq (vec1 = paramIndices, coeffs = paramCoeffs)
	 */
	inline void addTerm(const int ind, const double val, indices_t& vec1, coefficients_t& coeffs) {

		if (val == 0.0) { if(CONSTRAINT_EXPRESSION_OUTPUTLEVEL) std::cerr << "Warning: Attempting to add term with zero coefficient. Ignoring. \n"; return; }

		const int Pos = std::find(vec1.begin(), vec1.end(), ind) - vec1.begin();

		if (Pos >= (int)vec1.size()) {
			vec1.push_back(ind);
			coeffs.push_back(val);
		}
		else {
			assert(vec1.at(Pos) == ind);
			coeffs.at(Pos) += val;
		}

		assert(isConsistentWithDesign());
		return;
	}
	inline void reserve(size_t n) {
		varIndices.reserve(n);
		varCoeffs.reserve(n);
		paramIndices.reserve(n);
		paramCoeffs.reserve(n);
		bilinearIndices.reserve(n);
		bilinearCoeffs.reserve(n);
		
		return;
	}

public:
	ConstraintExpression() : sense('L'), rhs(0) { reserve(1000); }
	ConstraintExpression(const std::string& argName) : sense('L'), rhs(0), name(argName) { reserve(1000); }
	ConstraintExpression(const std::string& argName, const char argSense, const double argRhs)
		: sense(argSense), rhs(argRhs), name(argName) { reserve(1000); assert(isConsistentWithDesign()); }
	~ConstraintExpression() {/**/}
	ConstraintExpression (const ConstraintExpression&) = default;
	ConstraintExpression (ConstraintExpression&&) = default;
	ConstraintExpression& operator=(const ConstraintExpression&) = default;
	ConstraintExpression& operator=(ConstraintExpression&&) = default;
	
	inline void sign(const char argSense) { sense = argSense; assert(isConsistentWithDesign()); return; }
	inline void RHS(const double argRhs) { rhs = argRhs; return; }
	inline void rowname(const std::string& argName) { name = argName; return; }
	inline void addTermX(const int ind, const double val) { addTerm(ind, val, varIndices, varCoeffs); return; }
	inline void addTermQ(const int ind, const double val) { addTerm(ind, val, paramIndices, paramCoeffs); return; }
	inline void clear() {
		sense = 'L';
		rhs = 0;
		name.clear();
		varIndices.clear();
		varCoeffs.clear();
		bilinearIndices.clear();
		bilinearCoeffs.clear();
		paramIndices.clear();
		paramCoeffs.clear();
		return;
	}
	/* Add the term xSq */
	inline void addTermProduct(const int xind, const int qind, const double coeff = 1.0) {
		if (coeff == 0.0) { if(CONSTRAINT_EXPRESSION_OUTPUTLEVEL) std::cerr << "Warning: Attempting to add term with zero coefficient. Ignoring. \n"; return; }

		/* Add x to cx with c = 0 if it does not already exist */
		if (std::find(varIndices.begin(), varIndices.end(), xind) == varIndices.end()) {
			varIndices.push_back(xind);
			varCoeffs.push_back(0);
		}
		/* Add q to bq with b = 0 if it does not already exist */
		if (std::find(paramIndices.begin(), paramIndices.end(), qind) == paramIndices.end()) {
			paramIndices.push_back(qind);
			paramCoeffs.push_back(0);
		}
		
		const pair_t xq(xind, qind);
		/* Add xi*qj with coefficient "coeff" if this bilinear term does not already exist */
		const int curPos = std::find(bilinearIndices.begin(), bilinearIndices.end(), xq) - bilinearIndices.begin();
		if (curPos >= (int)bilinearIndices.size()) {
			bilinearIndices.emplace_back(xq);
			bilinearCoeffs.emplace_back(coeff);
		}
		/* Else, update S(xi,qj) with "coeff" */
		else {
			assert(bilinearIndices.at(curPos) == xq);
			bilinearCoeffs.at(curPos) += coeff;
		}

		assert(isConsistentWithDesign());
		return;
	}
	
	
	
	inline bool isEmpty() const {
		for (const auto& i : varCoeffs) if (i != 0.0) return false;
		for (const auto& i : bilinearCoeffs) if (i != 0.0) return false;
		return true;
	}
	inline bool existConstQTerms() const {
		for (const auto& i : paramCoeffs) if (i != 0.0) return true;
		return false;
	}
	inline bool existBilinearTerms() const {
		for (const auto& i : bilinearCoeffs) if (i != 0.0) return true;
		return false;
	}
	inline std::vector<std::string> getVarNames(CPXCENVptr env = nullptr, CPXCLPptr lp = nullptr) const {

		std::vector<std::string> varNames;
		for (unsigned i = 0; i < varIndices.size(); i++) {
			if (env == nullptr || lp == nullptr)
				varNames.emplace_back("x(" + std::to_string(varIndices[i]) + ")");
			else {
				CPXSIZE surplus = 0; CPXXgetcolname(env, lp, nullptr, nullptr, 0, &surplus, varIndices[i], varIndices[i]);
				CPXSIZE colnamespace = -surplus;
				if (colnamespace > 0) {
					char **colname     = (char **)malloc(sizeof(char *));
					char *colnamestore = (char *) malloc(colnamespace);
					CPXXgetcolname(env, lp, colname, colnamestore, colnamespace, &surplus, varIndices[i], varIndices[i]);

					varNames.emplace_back(colname[0], std::find(colname[0], colname[0] + strlen(colname[0]), '\0'));
					if (colname != NULL)      { free(colname);      colname      = nullptr; }
					if (colnamestore != NULL) { free(colnamestore); colnamestore = nullptr; }
				}
				else varNames.emplace_back("xFake(" + std::to_string(varIndices[i]) + ")");
			}
		}

		return varNames;
	}
	inline void print(CPXCENVptr env = nullptr, CPXCLPptr lp = nullptr, std::ostream& out = std::cout) const {
		if (isEmpty()) { out << "(empty expression)\n"; return; }

		/* Determine variable names */
		std::vector<std::string> varNames = getVarNames(env, lp);
		

		if (!name.empty()) out << "Printing " << name << "->    ";
		for (unsigned i = 0; i < varCoeffs.size(); i++)
			if(varCoeffs[i] != 0.0)
				out << ( (varCoeffs[i] > 0) ? "+" : "-") << std::abs(varCoeffs[i]) << " " << varNames[i] << "  " ;
		for (unsigned i = 0; i < paramCoeffs.size(); i++)
			if (paramCoeffs[i] != 0.0)
				out << ((paramCoeffs[i] > 0) ? "+" : "-") << std::abs(paramCoeffs[i]) << " q(" << paramIndices[i] << ")  ";
		for (unsigned i = 0; i < bilinearCoeffs.size(); i++)
			if (bilinearCoeffs[i] != 0.0) {
				out << ((bilinearCoeffs[i] > 0) ? "+" : "-") << std::abs(bilinearCoeffs[i]) << " ";
				out << varNames[std::find(varIndices.begin(), varIndices.end(), bilinearIndices[i].first) - varIndices.begin()] << " q(" << bilinearIndices[i].second << ")  ";
			}
		out << (sense=='E' ? "=" : (sense=='G' ? ">=" : "<=")) << " " << rhs << "\n";
		
		return;
	}
	inline void getDeterministicConstraint( const std::vector<double>& paramValues,
											CPXNNZ& nzcnt,
											double& trueRhs,
											char& trueSense,
											std::vector<CPXDIM>& rmatind,
											std::vector<double>& rmatval) const
	{
		trueSense = sense;

		/* Compute right-hand-size as supplied rhs - bq */
		trueRhs = rhs;
		for (unsigned i = 0; i < paramCoeffs.size(); i++)
			if (paramCoeffs.at(i) != 0.0)
				trueRhs -= (paramCoeffs.at(i)*paramValues.at(paramIndices.at(i)));
		
		/* Compute non-zeros, i.e., coefficients of decision variables */
		nzcnt = 0;
		rmatind.clear();
		rmatval.clear();
		for (unsigned i = 0; i < varIndices.size(); ++i) {
			bool inExpr = 0;  /* Does this decision variable appear in this expression? */
			double coef = 0.0;/* It's coefficient should be non-zero */

			/* Check terms involving x only: cx */
			if (varCoeffs[i] != 0.0) {
				coef += varCoeffs[i];
				inExpr = 1;
			}
			/* Check product terms: xSq */
			for (unsigned k = 0; k < bilinearIndices.size(); ++k) if (bilinearIndices[k].first == varIndices[i]) {
				if (bilinearCoeffs[k] != 0.0) {
					coef += paramValues.at(bilinearIndices[k].second)*bilinearCoeffs[k];
					inExpr = 1;
				}
			}
			
			/* This variable must appear somewhere */
			if (inExpr) {
				nzcnt++;
				rmatind.emplace_back(varIndices[i]);
				rmatval.emplace_back(coef);
			}
			else if(CONSTRAINT_EXPRESSION_OUTPUTLEVEL)
				std::cerr << "Warning: Expression contains empty variable. Ignoring. \n";
		}
		assert((CPXNNZ)rmatind.size() == nzcnt);
		assert(rmatval.size() == rmatind.size());

		return;
	}
	inline void getStochasticConstraint(const std::vector<double>& varValues,
										CPXNNZ& nzcnt,
										double& trueRhs,
										char& trueSense,
										std::vector<CPXDIM>& rmatind,
										std::vector<double>& rmatval) const
	{
		trueSense = sense;

		/* Compute right-hand-size as supplied rhs - cx */
		trueRhs = rhs;
		for (unsigned i = 0; i < varCoeffs.size(); i++)
			if (varCoeffs.at(i) != 0.0)
				trueRhs -= (varCoeffs.at(i)*varValues.at(varIndices.at(i)));
		
		/* Compute non-zeros, i.e., coefficients of decision variables */
		nzcnt = 0;
		rmatind.clear();
		rmatval.clear();
		for (unsigned i = 0; i < paramIndices.size(); ++i) {
			bool inExpr = 0;  /* Does this parameter appear in this expression? */
			double coef = 0.0;/* It's coefficient should be non-zero */

			/* Check terms involving q only: bq */
			if (paramCoeffs[i] != 0.0) {
				coef += paramCoeffs[i];
				inExpr = 1;
			}
			/* Check product terms: xSq */
			for (unsigned k = 0; k < bilinearIndices.size(); ++k) if (bilinearIndices[k].second == paramIndices[i]) {
				if (bilinearCoeffs[k] != 0.0) {
					coef += varValues.at(bilinearIndices[k].first)*bilinearCoeffs[k];
					inExpr = 1;
				}
			}
			
			/* This variable must appear somewhere */
			if (inExpr) {
				nzcnt++;
				rmatind.emplace_back(paramIndices[i]);
				rmatval.emplace_back(coef);
			}
			else if(CONSTRAINT_EXPRESSION_OUTPUTLEVEL)
				std::cerr << "Warning: Expression contains empty uncertain parameter. Ignoring. \n";
		}
		assert((CPXNNZ)rmatind.size() == nzcnt);
		assert(rmatval.size() == rmatind.size());

		return;
	}
	inline int addToCplex(CPXCENVptr env,
		                  CPXLPptr lp,
						  UNCSetCPtr UncSet = nullptr,
						  const bool isDualized = false,
						  const std::vector<double>& paramValues = {}) const
	{
		if (isEmpty()) {
			if(CONSTRAINT_EXPRESSION_OUTPUTLEVEL) std::cerr << " Warning: Attempting to add empty expression to CPLEX. Ignoring. \n";
			return 0;
		}

		/* Do the uncertain parameters appear as coefficients of decision variables? */
		const bool qTermsInProduct = existBilinearTerms();
		/* Do the uncertain parameters appear independently by themselves? */
		const bool qTermsConst     = existConstQTerms();
		/* Is this a "deterministic" constraint that does not require any dualization? */
		bool deterministic = !isDualized;

		if (isDualized) if (!qTermsInProduct) if (!qTermsConst) deterministic = true;
		if (deterministic) if (qTermsInProduct || qTermsConst) if (paramValues.empty()) { std::cerr << " Error: Parameter values not supplied for deterministic model. \n"; return 1; }




		//=====================================================================
		/* THIS CONSTRAINT DOES NOT REQUIRE REFORMULATION IN THE MAIN MODEL */
		//=====================================================================
		if (deterministic) {
			char trueSense = sense;
			const CPXNNZ rmatbeg = 0;
			const char *rowname = name.c_str();

			/* Compute right-hand-size as supplied rhs - bq */
			double trueRhs = rhs;
			for (unsigned i = 0; i < paramCoeffs.size(); i++)
				if (paramCoeffs.at(i) != 0.0)
					trueRhs -= (paramCoeffs.at(i)*paramValues.at(paramIndices.at(i)));
			
			/* Compute non-zeros, i.e., coefficients of decision variables */
			CPXNNZ nzcnt = 0;
			std::vector<CPXDIM> rmatind;
			std::vector<double> rmatval;
			for (unsigned i = 0; i < varIndices.size(); ++i) {
				bool inExpr = 0;  /* Does this decision variable appear in this expression? */
				double coef = 0.0;/* It's coefficient should be non-zero */

				/* Check terms involving x only: cx */
				if (varCoeffs[i] != 0.0) {
					coef += varCoeffs[i];
					inExpr = 1;
				}
				/* Check product terms: xSq */
				for (unsigned k = 0; k < bilinearIndices.size(); ++k) if (bilinearIndices[k].first == varIndices[i]) {
					if (bilinearCoeffs[k] != 0.0) {
						coef += paramValues.at(bilinearIndices[k].second)*bilinearCoeffs[k];
						inExpr = 1;
					}
				}
				
				/* This variable must appear somewhere */
				if (inExpr) {
					nzcnt++;
					rmatind.emplace_back(varIndices[i]);
					rmatval.emplace_back(coef);
				}
				else if(CONSTRAINT_EXPRESSION_OUTPUTLEVEL)
					std::cerr << "Warning: Expression contains empty variable. Ignoring. \n";
			}
			assert((CPXNNZ)rmatind.size() == nzcnt);
			assert(rmatval.size() == rmatind.size());

			return CPXXaddrows(env, lp, 0, 1, nzcnt, &trueRhs, &trueSense, &rmatbeg, &rmatind[0], &rmatval[0], nullptr, &rowname);
		}
		//=====================================================================
		/* DUALIZE CONSTRAINT BEFORE ADDING BACK TO THE (REFORMULATED) MODEL */
		//=====================================================================
		else {
			if (UncSet == nullptr) {
				std::cerr << " Error: Null pointer supplied for uncertainty set in dualized constraint. \n";
				return 1;
			}
			if (!UncSet->isUncertain()) {
				std::cerr << " Error: Attempting to dualize constraint for deterministic problem. \n";
				return 1;
			}
			if (!qTermsInProduct) if (qTermsConst) { /* "Effectively" deterministic constraint */
				if (sense == 'E') {
					std::cerr << "Error: Attempting to dualize equality constraint in which uncertain parameters do not appear as coefficients of decision variables. \n";
					return 1;
				}
				//if(CONSTRAINT_EXPRESSION_OUTPUTLEVEL) std::cerr << " Warning: Not dualizing constraint because uncertain parameters do not appear as coefficients of decision variables. \n";
				
				/* If this is a >= constraint, move the q_terms to the right hand side and then maximize */
				std::vector<double> vecCoeffs(paramCoeffs); if (sense == 'G') for (auto& qj : vecCoeffs) qj *= (-1.0);

				/* Obtain maximum value of the q_terms over the uncertainty set */
				const double max_qTermsConst = UncSet->getMaximumValue(paramIndices, vecCoeffs);
				const double NewRhs = rhs + ((sense == 'G') ? +max_qTermsConst : -max_qTermsConst);
				
				/* Create a new constraint which is a copy of *this but with modified RHS
				 * and add that constraint deterministically
				 */
				ConstraintExpression NewExpr(name, sense, NewRhs);
				for (unsigned i = 0; i < varIndices.size(); ++i)
					NewExpr.addTermX(varIndices[i], varCoeffs[i]);
				return NewExpr.addToCplex(env, lp);
			}

			/*************************************************************/
			/* Before you add the "single" dualized constraint,
			 * create as many new dual variables as there are rows in your
			 * "W matrix" of the uncertainty set
			 */
			/*************************************************************/
			assert(qTermsInProduct);

			

			int status = 0; /* return value */

			if (sense == 'E') {
				std::cerr << "Error!!!! To-do!! Coefficient matching in equality constraints!! \n";
				exit(100);
			}

			/* r = # of new (dual) variables created on top of existing model */
			const int r = UncSet->addVariables_DualVars(env, lp, name);
			if (r == 0) { std::cerr << "Error: No dual variables created while reformulating/dualizing constraints.\n"; return 1; }

			/* UncSet: W*q + V*Ksee [sense] h */
			/* Depending on sense [>=, <=, =], the dual variables have appropriate signs and bounds */
			const auto W = UncSet->getMatrixW();
			const auto V = UncSet->getMatrixV();
			const auto h = UncSet->getMatrixH();
			assert(W.size() == h.size());
			assert(V.size() == h.size());
			assert(r + 1 == (int)h.size());

			/* Last variable index of the original model
			 * before new variables were created
			 */
			const int lastIndBeforeDual = CPXXgetnumcols(env, lp) - r - 1;




			/****/
			/****/
			


			/* Add the constraint corresponding to the objective of the inner maximization problem */
			ConstraintExpression ReformObjExpr("DualObj(" + name + ")", sense, rhs);

			/* Terms involving only x (e.g., cx) carry over as-is */
			for (unsigned i = 0; i < varIndices.size(); ++i) if(varCoeffs[i] != 0.0) ReformObjExpr.addTermX(varIndices[i], varCoeffs[i]);
			
			/* Add objective function expression of the dual of the inner problem
			 * Only add those coefficients which are true non-zeros
			 * SUM[s = 1 to r, h(s)Ksee(s)]
			 * Note: if the original constraint is a >= constraint, then multiply above term by (-1)
			 */
			for (int s = 1; s <= r; s++) if (h.at(s) != 0.0) {
				const int dualVarIndex = lastIndBeforeDual + s;
				ReformObjExpr.addTermX(dualVarIndex, ((sense == 'G') ? -1 : +1)*h.at(s));
			}

			/* Add constraint */
			status += ReformObjExpr.addToCplex(env, lp);

			

			/****/
			/****/



			/* Add the constraints of the dual problem */
			/* Assuming that rows and columns of W, V and h are 1-indexed */
			for (int ll = 1; ll <= UncSet->getNoOfUncertainParameters(); ll++) {
				const int l = UncSet->getParamIndex('q', ll);

				/* To get the rhs of the constraint of the dual problem, compute the coefficient of q(l) in the original constraint */
				/* This coefficient has an x component (coming from xSq) and a constant component (coming from bq) */
				/* The "rhs" (which is a constant) comes from the constant component bq */
				double dualConstraintRhs = 0;

				/* If b(l)*q(l) is in this expression, then rhs = b(l) if <=, -b(l) if >= */
				const int lPos = std::find(paramIndices.begin(), paramIndices.end(), l) - paramIndices.begin();
				if (lPos < (int)paramIndices.size()) dualConstraintRhs += paramCoeffs[lPos];
				if (sense == 'G') dualConstraintRhs *= (-1.0);
				
				/* Constraint of dual problem corresponding to this column of the uncertainty set */
				ConstraintExpression ReformConstrExpr("DualCon(" + name + "," + std::to_string(l) + ")");
				ReformConstrExpr.sign('E');
				ReformConstrExpr.RHS(dualConstraintRhs);

				/* Add terms involving products of q(l) with x(i), i.e., xSq */
				for (unsigned i = 0; i < bilinearIndices.size(); i++) if (bilinearIndices[i].second == l) {
					assert(lPos < (int)paramIndices.size());
					if (bilinearCoeffs[i] != 0.0)
						ReformConstrExpr.addTermX(bilinearIndices[i].first, ((sense == 'G') ? +1 : -1)*bilinearCoeffs[i]);
				}

				/* Add nonzeros corresponding to constraints of the uncertainty set */
				/* Only add those coefficients which are true non-zeros */
				/* SUM[s = 1 to r, W(s,l)Ksee(s)] */
				for (int s = 1; s <= r; s++) if (W.at(s).at(l) != 0.0) {
					const int dualVarIndex = lastIndBeforeDual + s;
					ReformConstrExpr.addTermX(dualVarIndex, W.at(s).at(l));
				}

				/* Add constraint */
				status += ReformConstrExpr.addToCplex(env, lp);
			}


			/****/
			/****/



			/* If using factor models, then we need an additional constraint
			 * Note that for any other uncertainty set, V = 0 and we'll end up adding no new constraints
			 */
			for (int f = 1; f <= UncSet->getNumBudgetsFactors(); f++) {
				bool nonzero = 0;
				for (int s = 1; s <= r; s++) if (V.at(s).at(f) != 0.0) { nonzero = 1; break; }

				
				/* Attempt to add the constraint only if at least one non-zero exists */
				if (nonzero) {
					ConstraintExpression ReformFactorExpr("DualConFactor(" + name + ",F" + std::to_string(f) + ")", 'E', 0);
					for (int s = 1; s <= r; s++) if (V.at(s).at(f) != 0.0) {
						const int dualVarIndex = lastIndBeforeDual + s;
						ReformFactorExpr.addTermX(dualVarIndex, V.at(s).at(f));
					}

					/* Add constraint */
					status += ReformFactorExpr.addToCplex(env, lp);
				}
				
				/****/
				/****/
			}

			return status;
		}
	}
	inline char getSense() const { return sense; }
	inline double getRHS() const { return rhs; }
	inline std::string getName() const { return name; }
	inline double getViolation_fixedXQ(const std::vector<double>& varValues,
	                                   const std::vector<double>& paramValues = {},
	                                   const bool normalized = 0) const
	{	
		double viol = 0;

		const bool xTerms          = !isEmpty();
		const bool qTermsInProduct = existBilinearTerms();
		const bool qTermsConst     = existConstQTerms();


		if (xTerms) {
			if(varValues.empty()) {
				std::cerr << "Error: getViolation() requires values for primal variables\n";
				return 1.0;
			}
		}
		if (qTermsInProduct || qTermsConst) {
			if (paramValues.empty()) {
				std::cerr << "Error: getViolation() requires values for uncertain parameters\n";
				return 1.0;
			}
		}
		else if (!paramValues.empty()) {
			if (CONSTRAINT_EXPRESSION_OUTPUTLEVEL)
				std::cerr << "Warning: Ignoring supplied parameter values in getViolation()\n";
		}

		// normalizing factor
		double f = (normalized ? 0.0 : 1.0);

		// -d
		viol += -rhs;

		// +bq*
		if (qTermsConst) {
			for (unsigned i = 0; i < paramCoeffs.size(); i++) {
				viol += paramCoeffs[i] * paramValues.at(paramIndices.at(i));
				if (normalized) f += paramCoeffs[i] * paramCoeffs[i];
			}
		}

		// +cx*
		if (xTerms) {
			for (unsigned i = 0; i < varCoeffs.size(); i++) {
				viol += varCoeffs[i] * varValues.at(varIndices.at(i));
			}
		}

		// +x*Sq*
		if (qTermsInProduct) {
			for (unsigned i = 0; i < bilinearCoeffs.size(); i++) {
				viol += bilinearCoeffs[i] * varValues.at(bilinearIndices.at(i).first) * paramValues.at(bilinearIndices.at(i).second);
				if (normalized) f += bilinearCoeffs[i] * varValues.at(bilinearIndices.at(i).first) * bilinearCoeffs[i] * varValues.at(bilinearIndices.at(i).first);
			}
		}

		// normalizing factor
		if (normalized) {
			assert(f != 0.0);
			viol /= std::sqrt(f);
		}

		// Make changes for >= , =
		switch (sense) {
			case 'G': viol *= -1.0; break;
			case 'E': if (viol < 0) viol *= -1.0; break;
		}

		return viol;
	}
	inline ConstraintExpression flip0() const {
		ConstraintExpression N;
		switch (sense) {
			case 'G': N.sense = 'L';break;
			case 'L': N.sense = 'G';break;
			case 'E': N.sense = 'E';break;
		}
		N.rhs = rhs;
		N.name = name;
		N.varIndices = paramIndices;
		N.varCoeffs = paramCoeffs;
		N.paramIndices = varIndices;
		N.paramCoeffs = varCoeffs;
		N.bilinearCoeffs = bilinearCoeffs;
		N.bilinearIndices = bilinearIndices;
		for (unsigned int i = 0; i < bilinearIndices.size(); i++) {
			N.bilinearIndices[i].first = bilinearIndices[i].second;
			N.bilinearIndices[i].second = bilinearIndices[i].first;
		}
		assert(N.isConsistentWithDesign());

		const int curPos = std::find(N.paramIndices.begin(), N.paramIndices.end(), 0) - N.paramIndices.begin();
		assert(curPos < (int)N.paramIndices.size());
		if (curPos < (int)N.paramIndices.size()) {
			assert(N.paramCoeffs[curPos] == 1);
			N.addTermX(0, N.paramCoeffs[curPos]);
			N.paramIndices.erase(N.paramIndices.begin() + curPos);
			N.paramCoeffs.erase(N.paramCoeffs.begin() + curPos);
		}
		assert(N.isConsistentWithDesign());
		return N;
	}
};








/*----------------------------------------------------------------------
  This class supports constraints of the following type:

      max_q min_p {cx(p) + bq + x(p)Sq - d(p)}    <= rhs
  <=> max_{q,t}   t                               <= rhs
          s.t.    t <= cx(p) + bq + x(p)Sq - d(p)        \forall p
                  q \in UncSet

[ cx(p) + bq(p) + x(p)Sq <= d(p) ] are objects of class ConstraintExpression
There must be p such objects (hence the parametrization on p)

*-----------------------------------------------------------------------*/
class KAdaptableExpression {
private:
	std::vector<ConstraintExpression> Expr; /* Collection of constraints */
	std::string name;						/* Name of this "meta-constraint" */
	double rhs;                             /* right-hand side of K-adaptable expression */

	inline bool isConsistentWithDesign() const {
		if(Expr.empty())                       return false;
		for (const auto& exp : Expr)
			if (!exp.isConsistentWithDesign()) return false;
		for (const auto& exp : Expr)
			if (exp.sense != 'L')              return false;

		return true;
	}

public:

	KAdaptableExpression(const std::vector<ConstraintExpression>& arg, const std::string arg_name, const double arg_rhs = 0)
	: Expr(arg), name(arg_name), rhs(arg_rhs)
	{
		/* Convert to <= form */
		for(auto& exp:Expr) {
			if (exp.sense == 'G') {
				exp.sign('L');
				exp.rhs *= (-1.0);
				for (auto& c : exp.varCoeffs)      c *= (-1.0);
				for (auto& b : exp.paramCoeffs)    b *= (-1.0);
				for (auto& S : exp.bilinearCoeffs) S *= (-1.0);
			}
			else if (exp.sense == 'E') {
				std::cerr << "Cannot use equality constraints in K-adaptable expression.\n";
				exit(-2);
			}
		}
		assert(isConsistentWithDesign());
	}

	inline bool existConstXTerms() const {
		for (const auto& exp : Expr) if (exp.isEmpty()) return false;
		return true;
	}
	inline bool existConstQTerms() const {
		for (const auto& exp : Expr) if(exp.existConstQTerms()) return true;
		return false;
	}
	inline bool existBilinearTerms() const {
		for (const auto& exp : Expr) if(exp.existBilinearTerms()) return true;
		return false;
	}
	
	



	/*------------------------------------------
	 | This function evaluates the value of:
	 |
	 | max_q min_p {cx*(p) + bq + x*(p)Sq - d(p)}
	 |
	 *------------------------------------------*/
	inline std::vector<double> evaluate(UNCSetCPtr UncSet,
	                                    const std::vector<double>& paramValues,
	                                    const std::vector<double>& varValues = {},
	                                    const bool normalized = 0) const
	{
		const bool qTermsInProduct = existBilinearTerms(); 
		const bool qTermsConst     = existConstQTerms();
		const bool xTermsConst     = existConstXTerms();

		const bool donotSolveLP    = (UncSet == nullptr) ? 1 : (!UncSet->isUncertain() || (UncSet->isUncertain() && !qTermsConst && !qTermsInProduct));

		/* Deterministic case: set q equal to supplied value */
		if(donotSolveLP) {

			
			if(qTermsInProduct || qTermsConst) if(paramValues.empty()) { std::cerr << " Error: Parameter values not supplied for evaluation of K-adaptable expression. \n"; exit(-1); }
			if(xTermsConst)                    if(varValues.empty())   { std::cerr << " Error: Primal variable values not supplied for evaluation of K-adaptable expression. \n"; exit(-1); }


			/* Value of the {...} term in max_q min_p {...} */
			std::vector<double> ExprEval(Expr.size(), 0.0);

			/* Normalized value --> quantity dividing the {...} term in max_q min_p {...} */
			std::vector<double> Normalization(Expr.size(), (normalized ? 0.0 : 1.0));
			
			/* -d(p) */
			for(size_t p=0; p<Expr.size(); p++)
				ExprEval[p] += -Expr[p].rhs;

			/* +bq* */
			if(qTermsConst)
				for(size_t p=0; p<Expr.size(); p++)
					for(size_t i=0; i<Expr[p].paramCoeffs.size(); i++) {
						const double t = Expr[p].paramCoeffs[i];
						ExprEval[p] += t * paramValues.at(Expr[p].paramIndices.at(i));
						if (normalized) Normalization[p] += t * t;
					}

			/* +cx*(p) */
			if(xTermsConst)
				for(size_t p=0; p<Expr.size(); p++)
					for(size_t i=0; i<Expr[p].varCoeffs.size(); i++)
						ExprEval[p] += Expr[p].varCoeffs[i]*varValues.at(Expr[p].varIndices.at(i));

			/* +x*(p)Sq* */
			if(qTermsInProduct)
				for(size_t p=0; p<Expr.size(); p++)
					for(size_t i=0; i<Expr[p].bilinearCoeffs.size(); i++) {
						const double t = Expr[p].bilinearCoeffs[i]*varValues.at(Expr[p].bilinearIndices.at(i).first);
						ExprEval[p] += t * paramValues.at(Expr[p].bilinearIndices.at(i).second);
						if (normalized) Normalization[p] += t * t;
					}
			
			if (normalized) for(size_t p=0; p<Expr.size(); p++) {
				assert(Normalization[p] != 0.0);
				ExprEval[p] /= std::sqrt(Normalization[p]);
			}

			/* Prepend evaluation to supplied q vector */
			std::vector<double> eval{*std::min_element(ExprEval.begin(), ExprEval.end())}; eval.insert(eval.end(), paramValues.begin(), paramValues.end());

			return eval;
		}
		else {

			assert(UncSet->isUncertain());
			assert(qTermsInProduct || qTermsConst);
			if (qTermsInProduct || xTermsConst) if (varValues.empty()) { std::cerr << " Error: Primal variable values not supplied for evaluation of K-adaptable expression. \n"; exit(-1); }
			
			
			
			int stat;

			/* Get ENV object from Uncertainty Set */
			CPXENVptr env = UncSet->getENVObject();
			
			/* Clone LP object from Uncertainty Set */
			CPXLPptr lp = UncSet->getLPObject(&stat);
			if (stat) {std::cerr << "Failed to clone CPLEX LP object inside evaluation of K-adaptable expression. \n"; exit(-1);}

			

			
			int ind = 0; double coef = 1; char lb = 'L'; char ub = 'U'; double lbVal = -CPX_INFBOUND; double ubVal = +CPX_INFBOUND;

			/* Make it a maximization problem */
			CPXXchgobjsen(env, lp, CPX_MAX);
			/* Make it an LP */
			CPXXchgprobtype(env, lp, CPXPROB_LP);
			/* Replace current objective function */
			CPXXchgobj(env, lp, 1, &ind, &coef); coef = 0; for (int l = 1; l <= UncSet->getNoOfUncertainParameters(); ++l) CPXXchgobj(env, lp, 1, &l, &coef);
			/* Change bounds on objective function variable to make it (-Inf, +Inf) */
			CPXXchgbds(env, lp, 1, &ind, &lb, &lbVal);
			CPXXchgbds(env, lp, 1, &ind, &ub, &ubVal);

			/* Add constraints corresponding to each of the P constraints in Expr */
			for(const auto& exp: Expr) {
				ConstraintExpression TauConstr("TauCon(" + exp.name + ")");

				/* -d(p) + cx*(p) */
				double this_rhs = -exp.rhs; if (xTermsConst) for (size_t i = 0; i < exp.varCoeffs.size(); i++) this_rhs += exp.varCoeffs[i] * varValues.at(exp.varIndices.at(i));

				TauConstr.sign('L');
				TauConstr.RHS(this_rhs);

				/* Normalized value --> quantity dividing the {...} term in max_q min_p {...} */
				double Normalization = (normalized ? 0.0 : 1.0);

				/* -bq */
				for (size_t i = 0; i < exp.paramCoeffs.size(); i++) if (exp.paramCoeffs[i] != 0.0) {
					TauConstr.addTermX( exp.paramIndices[i], -exp.paramCoeffs[i]);
					if (normalized) Normalization += exp.paramCoeffs[i] * exp.paramCoeffs[i];
				}

				/* -x*Sq */
				if (qTermsInProduct) for (size_t i = 0; i < exp.bilinearCoeffs.size(); i++) if (exp.bilinearCoeffs[i] != 0.0) if (varValues.at(exp.bilinearIndices[i].first) != 0.0) {
					const double t = (exp.bilinearCoeffs[i] * varValues.at(exp.bilinearIndices[i].first));
					TauConstr.addTermX( exp.bilinearIndices[i].second, -t );
					if (normalized) Normalization += t * t;
				}

				/* \tau */
				TauConstr.addTermX(0, std::sqrt(Normalization));

				TauConstr.addToCplex(env, lp);
			}


			/* Solve LP */
			stat = CPXXlpopt(env, lp);
			if (stat) std::cerr << "Could not solve LP inside evaluation of K-adaptable expression.\n";


			/* Assert that LP was solved to optimality */
			if (CPXXgetstat(env, lp) != CPX_STAT_OPTIMAL) {
				std::cerr << "\n\n Could not solve LP to optimality inside the evaluation problem of K-adaptable expression.\n\n";
				exit(-2);
			}

			/* Get optimal RHS value */
			double objval; CPXXgetobjval(env, lp, &objval);

			/* Get primal solution vector */
			std::vector<double> primalX(1 + UncSet->getNoOfUncertainParameters(), objval); CPXXgetx(env, lp, &primalX[1], UncSet->getParamIndex('q', 1), UncSet->getParamIndex('q', UncSet->getNoOfUncertainParameters()));

			/* Free memory */
			if (lp != nullptr) {
				stat = CPXXfreeprob(env, &lp);
				if (stat) std::cerr << "CPXfreeprob failed inside K-adaptable expression, error code " << stat << "\n";
			}


			return primalX;
		}
	}
};








static inline double getViolation(const ConstraintExpression& con,
                                  const std::vector<double>& varValues,
                                  const std::vector<double>& paramValues = {},
                                  const bool normalized = 0)
{
	if (con.getSense() == 'E') {
		return con.getViolation_fixedXQ(varValues, paramValues, normalized);
	}

	// Construct K-adaptable expression with single policy
	std::vector<ConstraintExpression> exp{con};
	KAdaptableExpression CExpr(exp, con.getName());

	// evaluate worst-case value
	const auto q = CExpr.evaluate(nullptr, paramValues, varValues, normalized);
	assert(!q.empty());

	return q[0];
}

static inline std::vector<double> getViolation(const ConstraintExpression& con,
                                               UNCSetCPtr UncSet,
                                               const std::vector<double>& varValues = {},
                                               const bool normalized = 0)
{
	// temporary var
	std::vector<double> paramValues;

	if (UncSet == nullptr && con.getSense() == 'E') {
		const double viol = con.getViolation_fixedXQ(varValues, paramValues, normalized);
		std::vector<double> q{viol};
		q.insert(q.end(), paramValues.begin(), paramValues.end());
		return q;
	}

	// Construct K-adaptable expression with single policy
	std::vector<ConstraintExpression> exp{con};
	KAdaptableExpression CExpr(exp, con.getName());

	// evaluate worst-case value
	const auto q = CExpr.evaluate(UncSet, paramValues, varValues, normalized);
	assert((int)q.size() == 1 + UncSet->getNoOfUncertainParameters());

	return q;
}

#endif