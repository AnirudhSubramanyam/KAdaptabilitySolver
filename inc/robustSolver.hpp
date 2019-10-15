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


#ifndef KADAPTABLESOLVER_HPP
#define KADAPTABLESOLVER_HPP

#include "problemInfo.hpp"
#include <ilcplex/cplexx.h>
#include <vector>



/**
 * Class to solve all CPLEX-based computational examples
 * in K-adaptability paper.
 * 
 * Note that the problem is stored in pInfo (temporary variable),
 * which must be in the same scope as the solver object.
 */
class KAdaptableSolver {
protected:

	/** Pointer to const information */
	KAdaptableInfo *pInfo = NULL;

	/** Best known solution for K-adaptable problem (in MILP representation) */
	std::vector<double> xsol;






public:
	/**
	 * Delete default constructor
	 */
	KAdaptableSolver() = delete;


	/**
	 * Constructor takes a specific instance as input
	 */
	KAdaptableSolver(const KAdaptableInfo& pInfoData);


	/**
	 * Define pre-C++0x methods
	 */
	~KAdaptableSolver();
	KAdaptableSolver (const KAdaptableSolver&);
	KAdaptableSolver& operator=(const KAdaptableSolver&);


	/**
	 * Delete C++11 methods
	 */
	KAdaptableSolver (KAdaptableSolver&&) = delete;
	KAdaptableSolver& operator=(KAdaptableSolver&&) = delete;


	/**
	 * Reset solver with information about a new problem
	 * @param pInfoData reference to a new problem instance
	 */
	void setInfo(const KAdaptableInfo& pInfoData);
	

	/**
	 * Reset the solver with a new problem
	 * @param pInfoData reference to a new problem instance
	 * @param K         upper bound on # of 2nd-stage policies to be supported
	 */
	void reset(const KAdaptableInfo& pInfoData, const unsigned int K = 1);






public:

	/**
	 * Set the best known solution for K-adaptable problem
	 * @param x candidate best known solution for K-adaptable problem
	 * @param K # of second-stage policies
	 */
	void setX(const std::vector<double>& x, const unsigned int K);

	/**
	 * Return # of 2nd-stage policies in candidate solution
	 * @param  x candidate solution
	 * @return   # of 2nd-stage policies in x
	 */
	unsigned int getNumPolicies(const std::vector<double>& x) const;

	/**
	 * Resize the K'-adaptable solution x to the size of a K-adaptable solution
	 * @param x candidate solution
	 * @param K # of 2nd-stage policies
	 */
	void resizeX(std::vector<double>& x, const unsigned int K) const;

	/**
	 * Remove policy number k of K-adaptable solution x (remove only 2nd-stage variables)
	 * @param  x candidate solution
	 * @param  K # of 2nd-stage policies
	 * @param  k policy number of interest
	 */
	void removeXPolicy(std::vector<double>& x, const unsigned int K, const unsigned int k) const;

	/**
	 * Return policy number k of K-adaptable solution x (1st-stage and 2nd-stage variables)
	 * @param  x candidate solution
	 * @param  K # of 2nd-stage policies
	 * @param  k policy number of interest
	 * @return   policy k of candidate solution
	 */
	const std::vector<double> getXPolicy(const std::vector<double>& x, const unsigned int K, const unsigned int k) const;






public:
	/** 
	 * Add 1st-stage variables and constraints (not involving uncertain parameters)
	 * [Does not check if model is of the correct size]
	 * 
	 * @param env solver environment object
	 * @param lp  solver model object
	 */
	void updateX(CPXCENVptr env, CPXLPptr lp) const;

	/**
	 * Add 1st-stage constraints involving uncertain parameters
	 * [Does not check if model is of the correct size]
	 * 
	 * @param env solver environment object
	 * @param lp   solver model object
	 * @param q    candidate scenario (if empty, then attempt to reformulate via dualization)
	 * @param lazy add constraints in a lazy manner
	 */
	void updateXQ(CPXCENVptr env, CPXLPptr lp, const std::vector<double>& q = {}, const bool lazy = false) const;
	
	/**
	 * Add 2nd-stage variables and constraints (not involving uncertain parameters) for policy k
	 * [Does not check if model is of the correct size]
	 * 
	 * @param env solver environment object
	 * @param lp  solver model object
	 * @param k   policy to which to add 2nd-stage variables and constraints
	 */
	void updateY(CPXCENVptr env, CPXLPptr lp, const unsigned int k) const;

	/**
	 * Add 2nd-stage constraints involving uncertain parameters in policy k
	 * [Does not check if model is of the correct size]
	 * 
	 * @param env solver environment object
	 * @param lp   solver model object
	 * @param k    policy to which to add constraints
	 * @param q    candidate scenario (if empty, then attempt to reformulate via dualization)
	 * @param lazy add constraints in a lazy manner
	 */
	void updateYQ(CPXCENVptr env, CPXLPptr lp, const unsigned int k, const std::vector<double>& q = {}, const bool lazy = false) const;

	/**
	 * Get all constraints involving uncertain parameters and first-stage variables using realization q
	 * @param q       candidate scenario
	 * @param rcnt    # of constraints will be returned here (size of rhs, sense, rmatbeg)
	 * @param nzcnt   # of non-zeros across all constraints (size of rmatind, rmatval)
	 * @param rhs     rhs of each constraint
	 * @param sense   sense of each constraint
	 * @param rmatbeg rmatbeg[i] = starting position of constraint i in rmatind/rmatval, rmatbeg[0] = 0
	 * @param rmatind linear indices of variables appearing in constraints
	 * @param rmatval coefficients of variables appearing in constraints
	 */
	void getXQ_fixedQ(const std::vector<double>& q, CPXDIM& rcnt, CPXNNZ& nzcnt, std::vector<double>& rhs, std::vector<char>& sense, std::vector<CPXNNZ>& rmatbeg, std::vector<CPXDIM>& rmatind, std::vector<double>& rmatval) const;

	/**
	 * Get all constraints involving uncertain parameters and first-stage variables using decisions x
	 * @param x       candidate solution
	 * @param rcnt    # of constraints will be returned here (size of rhs, sense, rmatbeg)
	 * @param nzcnt   # of non-zeros across all constraints (size of rmatind, rmatval)
	 * @param rhs     rhs of each constraint
	 * @param sense   sense of each constraint
	 * @param rmatbeg rmatbeg[i] = starting position of constraint i in rmatind/rmatval, rmatbeg[0] = 0
	 * @param rmatind linear indices of parameters appearing in constraints
	 * @param rmatval coefficients of parameters appearing in constraints
	 */
	void getXQ_fixedX(const std::vector<double>& x, CPXDIM& rcnt, CPXNNZ& nzcnt, std::vector<double>& rhs, std::vector<char>& sense, std::vector<CPXNNZ>& rmatbeg, std::vector<CPXDIM>& rmatind, std::vector<double>& rmatval) const;

	/**
	 * Get all constraints involving uncertain parameters in policy k using realization q
	 * @param k       policy for which to get constraints
	 * @param q       candidate scenario
	 * @param rcnt    # of constraints will be returned here (size of rhs, sense, rmatbeg)
	 * @param nzcnt   # of non-zeros across all constraints (size of rmatind, rmatval)
	 * @param rhs     rhs of each constraint
	 * @param sense   sense of each constraint
	 * @param rmatbeg rmatbeg[i] = starting position of constraint i in rmatind/rmatval, rmatbeg[0] = 0
	 * @param rmatind linear indices of variables appearing in constraints
	 * @param rmatval coefficients of variables appearing in constraints
	 */
	void getYQ_fixedQ(const unsigned int k, const std::vector<double>& q, CPXDIM& rcnt, CPXNNZ& nzcnt, std::vector<double>& rhs, std::vector<char>& sense, std::vector<CPXNNZ>& rmatbeg, std::vector<CPXDIM>& rmatind, std::vector<double>& rmatval) const;

	/**
	 * Get all constraints involving uncertain parameters in policy k using decisions x
	 * @param k       policy for which to get constraints
	 * @param x       candidate solution
	 * @param rcnt    # of constraints will be returned here (size of rhs, sense, rmatbeg)
	 * @param nzcnt   # of non-zeros across all constraints (size of rmatind, rmatval)
	 * @param rhs     rhs of each constraint
	 * @param sense   sense of each constraint
	 * @param rmatbeg rmatbeg[i] = starting position of constraint i in rmatind/rmatval, rmatbeg[0] = 0
	 * @param rmatind linear indices of parameters appearing in constraints
	 * @param rmatval coefficients of parameters appearing in constraints
	 */
	void getYQ_fixedX(const unsigned int k, const std::vector<double>& x, CPXDIM& rcnt, CPXNNZ& nzcnt, std::vector<double>& rhs, std::vector<char>& sense, std::vector<CPXNNZ>& rmatbeg, std::vector<CPXDIM>& rmatind, std::vector<double>& rmatval) const;






public:
	
	/**
	 * Get a lower bound on the fully adaptive solution value
	 * @return lower bound on fully adaptive solution value
	 */
	double getLowerBound() const;

	/**
	 * Get the worst-case deterministic objective value
	 * @param  q scenario which generates the worst-case deterministic value is returned here
	 * @return   worst-case deterministic objective value
	 */
	double getLowerBound_2(std::vector<double>& q) const;

	/**
	 * Get the worst-case objective value of the given solution
	 * @param  x candidate solution
	 * @param  K # of 2nd-stage policies
	 * @param  q scenario which generates the worst-case objective value is returned here
	 * @return   worst-case objective value
	 */
	double getWorstCase(const std::vector<double>& x, const unsigned int K, std::vector<double>& q) const;






public:

	/**
	 * Check if the given solution is 'structurally' feasible for the K-adaptable problem
	 * Structural feasibility = satisfaction of all uncertainty-independent constraints.
	 * This includes bounds and constraints linking first- and second-stage variables.
	 * 
	 * @param  x candidate k-adaptable solution with k <= K
	 * @param  K  # of 2nd-stage policies
	 * @param  xx structurally feasible solution with (possibly) reduced number of policies will be returned here (if feasible)
	 * @return    true if feasible
	 */
	bool feasible_DET_K(const std::vector<double>& x, const unsigned int K, std::vector<double>& xx) const;

	/**
	 * Check if the given solution satisfies 1st-stage constraints involving uncertain parameters
	 * [to be used by internal routines only]
	 * 
	 * @param  x candidate solution
	 * @param  q scenario from the uncertainty set will be returned here
	 * @return   true if solution is feasible
	 */
	bool feasible_XQ(const std::vector<double>& x, std::vector<double>& q) const;

	/**
	 * Check if the given solution satisfies 1st-stage constraints involving uncertain parameter vectors contained in 'samples'
	 * [to be used by internal routines only]
	 * 
	 * @param  x       candidate solution
	 * @param  samples collection of scenarios to check feasibility against
	 * @param  label   sequence number of violating scenario in samples will be returned here (if not feasible)
	 * @return         true if solution is feasible
	 */
	bool feasible_XQ(const std::vector<double>& x, const std::vector<std::vector<double> >& samples, int & label) const;

	/**
	 * Check if the given K-adaptable solution satisfies 2nd-stage constraints involving uncertain parameters
	 * [to be used by internal routines only]
	 * 
	 * @param  x candidate solution
	 * @param  K # of 2nd-stage policies
	 * @param  q scenario from the uncertainty set will be returned here
	 * @return   true if solution is feasible
	 */
	bool feasible_YQ(const std::vector<double>& x, const unsigned int K, std::vector<double>& q) const;

	/**
	 * Check if the given K-adaptable solution satisfies 2nd-stage constraints involving uncertain parameter vectors contained in 'samples'
	 * [to be used by internal routines only]
	 * 
	 * @param  x       candidate solution
	 * @param  K       # of 2nd-stage policies
	 * @param  samples collection of scenarios to check feasibility against
	 * @param  label   sequence number of violating scenario in samples will be returned here (if not feasible)
	 * @return         true if feasible
	 */
	bool feasible_YQ(const std::vector<double>& x, const unsigned int K, const std::vector<std::vector<double> >& samples, int & label, bool heur = 0) const;






public:

	/**
	 * Check if the given solution is deterministically feasible for the given scenario
	 * @param  x candidate solution
	 * @param  q candidate scenario
	 * @return   true if feasible
	 */
	bool feasible_DET(const std::vector<double>& x, const std::vector<double>& q) const;

	/**
	 * Check if the given solution is static robust feasible
	 * @param  x candidate solution
	 * @param  q scenario from the uncertainty set will be returned here (if not feasible)
	 * @return   true if feasible
	 */
	bool feasible_SRO(const std::vector<double>& x, std::vector<double>& q) const;

	/**
	 * Check if the given solution is feasible for the collection of scenarios contained in 'samples'
	 * @param  x       candidate solution
	 * @param  samples collection of scenarios to check feasibility against
	 * @param  label   sequence number of violating scenario in samples will be returned here (if not feasible)
	 * @return         true if feasible
	 */
	bool feasible_SRO(const std::vector<double>& x, const std::vector<std::vector<double> >& samples, int & label) const;

	/**
	 * Check if the given solution is feasible for the K-adaptable problem
	 * @param  x candidate k-adaptable solution with k <= K
	 * @param  K # of 2nd-stage policies
	 * @param  q scenario from the uncertainty set will be returned here (if not feasible)
	 * @return   true if feasible
	 */
	bool feasible_KAdaptability(const std::vector<double>& x, const unsigned int K, std::vector<double>& q) const;

	/**
	 * Check if the given solution is K-adaptable feasible for the collection of scenarios contained in 'samples'
	 * @param  x       candidate k-adaptable solution with k <= K
	 * @param  K       # of 2nd-stage policies
	 * @param  samples collection of scenarios to check feasibility against
	 * @param  label   sequence number of violating scenario in samples will be returned here (if not feasible)
	 * @return         true if feasible
	 */
	bool feasible_KAdaptability(const std::vector<double>& x, const unsigned int K, const std::vector<std::vector<double> >& samples, int & label) const;






public:
	/**
	 * Solve the deterministic problem for the given scenario
	 * @param  q scenario from the uncertainty set
	 * @param  x optimal solution (if any) will be returned here
	 * @return   solve status (non-zero value indicates unsuccessful termination)
	 */
	int solve_DET(const std::vector<double>& q, std::vector<double>& x) const;

	/**
	 * Solve the static robust problem via duality-based reformulation
	 * @param  x optimal solution (if any) will be returned here
	 * @return   solve status (non-zero value indicates unsuccessful termination)
	 */
	int solve_SRO_duality(std::vector<double>& x) const;

	/**
	 * Solve the static robust problem via cutting plane method implemented using solver callbacks
	 * @param  x optimal solution (if any) will be returned here
	 * @return   solve status (non-zero value indicates unsuccessful termination)
	 */
	int solve_SRO_cuttingPlane(std::vector<double>& x);

	/**
	 * Solve the scenario-based static robust problem (uncertainty set is a finite set of points)
	 * @param  samples collection of scenarios to be robust against
	 * @param  x       optimal solution (if any) will be returned here
	 * @return         solve status (non-zero value indicates unsuccessful termination)
	 */
	int solve_ScSRO(const std::vector<std::vector<double> >& samples, std::vector<double>& x) const;






public:
	/** Library of samples (temporary var -- to be used by solve_KAdaptability() only) */
	std::vector<std::vector<double> > bb_samples;

	/** Feasible (ideally optimal) static robust solution -- to be used by solve_KAdaptability() only) */
	std::vector<double> xstatic;

	/** (K - 1)-adaptable solution -- to be used by solve_KAdaptability() in heuristic_mode only */
	std::vector<double> xfix;

	/** Run in heuristic mode -- to be used by solve_KAdaptability() only */
	bool heuristic_mode;

	/** Current value of K in K-Adaptability -- to be used by solve_KAdaptability() only) */
	unsigned int NK;

	/**
	 * Get best solution found
	 * @return best solution found
	 */
	inline const std::vector<double>& getX() const {
		return xsol;
	}

	/**
	 * Get nominal scenario
	 * @return nominal scenario from uncertainty set
	 */
	inline std::vector<double> getNominal() const {
		return pInfo->getNominal();
	}

	/**
	 * Is the underlying problem a pure LP?
	 * @return true if the problem is continuous
	 */
	inline bool isContinuous() const {
		return pInfo->isContinuous();
	}

	/**
	 * Indicates if problem has 2nd-stage decisions only
	 * @return true if problem has 2nd-stage decisions only
	 */
	inline bool isSecondStageOnly() const {
		return pInfo->isSecondStageOnly();
	}

	/**
	 * Indicates if problem has objective uncertainty only
	 * @return true if problem has objective uncertainty only
	 */
	inline bool hasObjectiveUncOnly() const {
		return pInfo->hasObjectiveUncOnly();
	}

	/**
	 * Solve the K-adaptability problem using solver callbacks
	 * @param  K # of 2nd-stage policies
	 * @param  h run in heuristic mode
	 * @param  x optimal solution (if any) will be returned here
	 * @return   solve status (non-zero value indicates unsuccessful termination)
	 */
	int solve_KAdaptability(const unsigned int K, const bool h, std::vector<double>& x);

	/**
	 * Solve the separation problem arising in solve_KAdaptability()
	 * @param  x     candidate solution
	 * @param  K     # of 2nd-stage policies
	 * @param  label sequence number of violating scenario in bb_samples will be returned here
	 * @param  heur  separation is not guaranteed to be exact if this parameter is true
	 * @return       true if separating scenario was found, false otherwise
	 */
	bool solve_separationProblem(const std::vector<double>& x, const unsigned int K, int& label, bool heur = 0);

	/**
	 * Solve the min-max-min robust combinatorial optimization problem (as in Buchheim & Kurtz)
	 * @param  K # of 2nd-stage policies
	 * @param  x optimal solution (if any) will be returned here
	 * @return   solve status (non-zero value indicates unsuccessful termination)
	 */
	int solve_KAdaptability_minMaxMin(const unsigned int K, std::vector<double>& x);

};

#endif