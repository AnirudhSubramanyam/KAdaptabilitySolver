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



////////////////////////////////////////////
// File: Uncertainty.hpp                  //
////////////////////////////////////////////
// Class that represents a generic        //
// uncertainty set                        //
////////////////////////////////////////////



#ifndef UNCERTAINTY_HPP
#define UNCERTAINTY_HPP


#include <ilcplex/cplexx.h>
#include <vector>
#include <string>
#include <utility>




/**
 * Generic uncertainty set 
 * All arrays are 1-indexed (indexing starts from 1, not 0)
 */
class UncertaintySet {
protected:

	/** # of uncertain parameters */
	int N;

	/** Nominal realization */
	std::vector<double> nominal;

	/** Lower bounds on uncertain parameters */
	std::vector<double> low;

	/** Upper bounds on uncertain parameters */
	std::vector<double> high;
	
	/** Constraint matrix defining uncertainty set (including bounds) */
	std::vector<std::vector<double> > polytope_W;

	/** Matrix of zeros (backward compatibility) */
	std::vector<std::vector<double> > polytope_V;
	
	/** Corresponding right hand side vector */
	std::vector<double> polytope_h;

	/* Corresponding sense of each constraint */
	std::vector<char> polytope_sense;

	/** Solver environment to carry out optimizations on uncertainty set */
	CPXENVptr env;

	/** Solver problem object to carry out above optimizations */
	CPXLPptr lp;



public:
	
	/** 
	 * Constructs an UncertaintySet object
	 */
	UncertaintySet();

	/**
	 * Copy constructor
	 */
	UncertaintySet(const UncertaintySet& U);

	/**
	 * Copy assignment operator
	 */
	UncertaintySet& operator=(const UncertaintySet& U);

	/** 
	 * Free memory allocated 
	 */
	~UncertaintySet() noexcept(false);

	/**
	 * Clear all data structures that this uncertainty set holds
	 */
	void clear();

	/**
	 * Add an uncertain parameter to the set
	 * @param nom the nominal value of the parameter
	 * @param lo  the lower bound of the parameter
	 * @param hi  the upper bound of the parameter
	 */
	void addParam(const double nom, const double lo, const double hi);

	/**
	 * Add a single facet to the uncertainty set
	 * @param input the list of uncertain parameters and their coefficients in the facet
	 * @param sense the sense of the constraint ('L', 'G', 'E')
	 * @param rhs   the rhs of the facet
	 */
	void addFacet(const std::vector<std::pair<int, double> >& input, const char sense, const double rhs);

	/**
	 * Compute the maximum of a linear function of uncertain parameters over the uncertainty set
	 * @param  input  the list of uncertain parameters and their coefficients in the linear function
	 * @param  result argmax is returned here (1-indexed)
	 * @return        the result of the optimization (objective)
	 */
	double max(const std::vector<std::pair<int, double> >& input, std::vector<double>& result) const;
	
	/**
	 * Compute the minimum of a linear function of uncertain parameters over the uncertainty set
	 * @param  input  the list of uncertain parameters and their coefficients in the linear function
	 * @param  result argmin is returned here (1-indexed)
	 * @return        the result of the optimization (objective)
	 */
	double min(const std::vector<std::pair<int, double> >& input, std::vector<double>& result) const;

	/**
	 * Overload of previous function (to maintain backward compatibility)
	 * @param  ind  the list of indices of the uncertain parameters
	 * @param  coef the list of coefficients of the above uncertain parameters
	 * @return      the result of the optimization (objective)
	 */
	double getMaximumValue(const std::vector<int>& ind, const std::vector<double>& coef) const;

	/**
	 * This function adds an entire set of (empty) columns to the given lp object.
	 * These columns correspond to dual variables associated with
	 * each of the constraints of the inner optimization problem (to be dualized),
	 * i.e., constraints of the uncertainty set.
	 *
	 * THIS FUNCTION MUST BE CALLED WHENEVER A CONSTRAINT (WITH UNCERTAIN PARAMETERS)
	 * IS TO BE DUALIZED AND ADDED TO THE MAIN MILP.
	 *
	 * @param  env      solver environment in which master MILP resides
	 * @param  lp       solver model object in which master MILP resides
	 * @param  dualName name of constraint to be attached to new columns
	 * @return          number of new dual variables added
	 */
	int addVariables_DualVars(CPXCENVptr env, CPXLPptr lp, const std::string& dualName = "") const;

	/** 
	 * Get nominal realization
	 * @return the nominal realization
	 */
	inline std::vector<double> getNominal() const { return nominal; }

	/** 
	 * Get lower bounds of uncertain parameters
	 * @return the lower bounds
	 */
	inline std::vector<double> getLowerBounds() const { return low; }

	/** 
	 * Get upper bounds of uncertain parameters
	 * @return the upper bounds
	 */
	inline std::vector<double> getUpperBounds() const { return high; }

	/**
	 * Get clone of solver model object representing uncertainty set
	 * @param  env_ pointer to environment to be used in conjunction with cloned solver model object
	 * @param  stat pointer to solve status of clone operation
	 * @return      pointer to clone of solver model object
	 */
	inline CPXLPptr getLPObject(CPXCENVptr env_, int *stat) const { return CPXXcloneprob(env_, lp, stat); }
	
	/**
	 * Get clone of solver model object representing uncertainty set
	 * @param  stat pointer to solve status of clone operation
	 * @return      pointer to clone of solver model object
	 */
	inline CPXLPptr getLPObject(int *stat) const { return getLPObject(env, stat); }

	/**
	 * Get pointer to solver environment to carry out optimizations on uncertainty set
	 * @return pointer to solver environment
	 */
	inline CPXENVptr getENVObject() const { return env; }

	/**
	 * Is the uncertainty set empty?
	 * @return bool indicating if uncertainty set is empty
	 */
	inline bool isUncertain() const { return (N != 0); }

	/**
	 * Get number of constraints defining uncertainty set (other than bound constraints)
	 * @return number of facets of uncertainty set
	 */
	inline int getNoOfFacets() const { return polytope_h.size(); }

	/**
	 * Get number of uncertain parameters
	 * @return number of uncertain parameters
	 */
	inline int getNoOfUncertainParameters() const { return N; }

	/**
	 * Get index of uncertain parameter inside uncertainty set
	 * @param  char character specifying uncertain parameter
	 * @param  ind  index of uncertain parameter
	 * @return      index of uncertain parameter
	 */
	inline int getParamIndex(const char, const int ind) const { return ind; }

	/**
	 * Get number of factors (backward compatibility)
	 * @return number of factors
	 */
	inline int getNumBudgetsFactors() const { return 0; }

	/**
	 * Get constraint matrix W in Wq [sense] h
	 * @return constraint matrix W
	 */
	inline std::vector<std::vector<double> > getMatrixW() const { return polytope_W; }

	/**
	 * Get constraint matrix V
	 * @return constraint matrix V
	 */
	inline std::vector<std::vector<double> > getMatrixV() const { return polytope_V; }

	/**
	 * Get vector h in Wq [sense] h
	 * @return right hand side vector h
	 */
	inline std::vector<double> getMatrixH() const { return polytope_h; }

	/**
	 * Get vector h in Wq [sense] h
	 * @return right hand side vector h
	 */
	inline std::vector<char> getMatrixSense() const { return polytope_sense; }
};

typedef const UncertaintySet* UNCSetCPtr;
typedef UncertaintySet* UNCSetPtr;


#endif


