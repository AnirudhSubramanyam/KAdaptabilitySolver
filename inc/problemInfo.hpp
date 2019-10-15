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


#ifndef KADAPTABLEINFO_HPP
#define KADAPTABLEINFO_HPP

#include "indexingTools.h"
#include "constraintExpr.hpp"
#include "uncertainty.hpp"
#include <vector>
#include <string>

/**
 * Abstract class representing information that any
 * two-/multi-stage RO instance must provide.
 * All problems must be minimization problems.
 */
class KAdaptableInfo {

protected:

	/** Continuous or (mixed) integer problem? */
	bool hasInteger;

	/** Uncertainty in objective only */
	bool objectiveUnc;

	/** Any 1st-stage decisions to be made? */
	bool existsFirstStage;

	/** # of 1st-stage variables including worst-case objective */
	unsigned int numFirstStage;

	/** # of 2nd-stage variables */
	unsigned int numSecondStage;

	/** # of 2nd-stage policies */
	unsigned int numPolicies;

	/** Solution file name */
	std::string solfilename;

	/** Resident uncertainty set */
	UncertaintySet U;

	/** 1st-stage variables only */
	VarInfo X;

	/** 2nd-stage variables only */
	VarInfo Y;

	/** Constraints with 1st-stage variables only (excluding bounds) */
	std::vector<ConstraintExpression> C_X;

	/** Constraints with 2nd-(and, if applicable, 1st-) stage variables BUT NO uncertain parameters (excluding bounds) */
	std::vector<std::vector<ConstraintExpression> > C_XY;

	/** Constraints with 1st-stage variables AND uncertain parameters */
	std::vector<ConstraintExpression> C_XQ;

	/** Constraints with 2nd-(and, if applicable, 1st-) stage variables AND uncertain parameters */
	std::vector<std::vector<ConstraintExpression> > C_XYQ;

	/** Bound constraints on 1st-stage variables only */
	std::vector<ConstraintExpression> B_X;

	/** Bound constraints on 2nd-stage variables only */
	std::vector<std::vector<ConstraintExpression> > B_Y;

	/**
	 * Define 1st-stage and 2nd-stage variables
	 */
	virtual void makeVars() = 0;

	/**
	 * Define uncertainty set
	 */
	virtual void makeUncSet() = 0;

	/**
	 * Define constraints involving 1st-stage variables only
	 */
	virtual void makeConsX() = 0;
	
	/**
	 * Define constraints involving 2nd-stage variables
	 * @param k policy number to make constraints for
	 */
	virtual void makeConsY(unsigned int k = 0) = 0;

public:
	/**
	 * Virtual destructor does nothing
	 */
	virtual ~KAdaptableInfo() {/**/}

	/**
	 * Virtual (default) constructor
	 */
	virtual KAdaptableInfo* create() const = 0;

	/**
	 * Virtual (copy) constructor
	 */
	virtual KAdaptableInfo* clone() const = 0;

	/**
	 * Resizes data structures to contain at least K 2nd-stage policies
	 * @param K Number of 2nd-stage policies to be supported
	 */
	void resize(unsigned int K);

	/**
	 * Check consistency with class design
	 * @return true if object is conistent with class design
	 */
	bool isConsistentWithDesign() const;

	/**
	 * Linear variable index of 1st-stage variable
	 * [All indices start from 0]
	 * 
	 * @param  type name of variable
	 * @param  ind1 first index
	 * @param  ind2 second index
	 * @param  ind3 third index
	 * @param  ind4 fourth index
	 * @return      the linear index of the variable
	 */
	inline int getVarIndex_1(const std::string type, const int ind1, const int ind2 = -1, const int ind3 = -1, const int ind4 = -1) const {
		return X.getDefVarLinIndex(type, ind1, ind2, ind3, ind4);
	}

	/**
	 * Linear variable index of 2nd-stage variable
	 * [All indices start from 0]
	 * [Policy numbering starts from 0]
	 *
	 * @param  k    policy number in which you want linear index
	 * @param  type name of variable
	 * @param  ind1 first index
	 * @param  ind2 second index
	 * @param  ind3 third index
	 * @param  ind4 fourth index
	 * @return      the linear index of the variable
	 */
	inline int getVarIndex_2(const unsigned int k, const std::string type, const int ind1, const int ind2 = -1, const int ind3 = -1, const int ind4 = -1) const {
		return numFirstStage + (k * numSecondStage) + Y.getDefVarLinIndex(type, ind1, ind2, ind3, ind4);
	}

	/**
	 * Get total # of 1st- and 2nd-stage decisions
	 * @param  K Total number of 2nd-stage policies to consider
	 * @return   Total # of 1st- and 2nd-stage decisions for a K-policy problem
	 */
	inline int getNumVars(const unsigned int K = 1) const {
		return numFirstStage + (K * numSecondStage);
	}

	/**
	 * Get # of 1st-stage decisions
	 * @return # of 1st-stage decisions
	 */
	inline int getNumFirstStage() const {
		return numFirstStage;
	}

	/**
	 * Get # of 2nd-stage decisions
	 * @return # of 2nd-stage decisions
	 */
	inline int getNumSecondStage() const {
		return numSecondStage;
	}

	/**
	 * Get current upper bound on # of 2nd-stage policies
	 * @return current upper bound on # of 2nd-stage policies
	 */
	inline int getNumPolicies() const {
		return numPolicies;
	}

	/**
	 * Get solution file name
	 * @return solution file name
	 */
	inline std::string getSolFileName() const {
		return solfilename;
	}

	/**
	 * Get nominal parameter values
	 * @return nominal parameter values
	 */
	inline std::vector<double> getNominal() const {
		return U.getNominal();
	}

	/**
	 * Get # of uncertain parameters
	 * @return # of uncertain parameters
	 */
	inline int getNoOfUncertainParameters() const {
		return U.getNoOfUncertainParameters();
	}

	/**
	 * Indicates if problem has continuous variables only
	 * @return true if problem has continuous variables only
	 */
	inline bool isContinuous() const {
		return !hasInteger;
	}

	/**
	 * Indicates if problem has objective uncertainty only
	 * @return true if problem has objective uncertainty only
	 */
	inline bool hasObjectiveUncOnly() const {
		return objectiveUnc;
	}

	/**
	 * Indicates if problem has 2nd-stage decisions only
	 * @return true if problem has 2nd-stage decisions only
	 */
	inline bool isSecondStageOnly() const {
		return !existsFirstStage;
	}

	/**
	 * Get resident uncertainty set
	 * @return return the uncertainty set
	 */
	inline const decltype(U)& getUncSet() {
		return U;
	}

	/**
	 * Get 1st-stage variables only
	 * @return the 1st-stage variables
	 */
	inline const decltype(X)& getVarsX() {
		return X;
	}

	/**
	 * Get 2nd-stage variables only
	 * @return the 2nd-stage variables
	 */
	inline const decltype(Y)& getVarsY() {
		return Y;
	}

	/**
	 * Get constraints with 1st-stage variables only (excluding bounds)
	 * @return constraints with 1st-stage variables only (excluding bounds)
	 */
	inline const decltype(C_X)& getConstraintsX() {
		return C_X;
	}

	/**
	 * Get constraints with 2nd-(and, if applicable, 1st-) stage variables BUT NO uncertain parameters (excluding bounds)
	 * @return constraints with 2nd-(and, if applicable, 1st-) stage variables BUT NO uncertain parameters (excluding bounds)
	 */
	inline const decltype(C_XY)& getConstraintsXY() {
		return C_XY;
	}

	/**
	 * Get constraints with 1st-stage variables AND uncertain parameters
	 * @return constraints with 1st-stage variables AND uncertain parameters
	 */
	inline const decltype(C_XQ)& getConstraintsXQ() {
		return C_XQ;
	}

	/**
	 * Get constraints with 2nd-(and, if applicable, 1st-) stage variables AND uncertain parameters
	 * @return constraints with 2nd-(and, if applicable, 1st-) stage variables AND uncertain parameters
	 */
	inline const decltype(C_XYQ)& getConstraintsXYQ() {
		return C_XYQ;
	}

	/**
	 * Get bound constraints on 1st-stage variables only
	 * @return bound constraints on 1st-stage variables only
	 */
	inline const decltype(B_X)& getBoundsX() {
		return B_X;
	}

	/**
	 * Get bound constraints on 2nd-stage variables only
	 * @return bound constraints on 2nd-stage variables only
	 */
	inline const decltype(B_Y)& getBoundsY() {
		return B_Y;
	}

};


#endif