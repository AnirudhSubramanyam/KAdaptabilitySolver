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


#ifndef KNP_INSTANCE_HPP
#define KNP_INSTANCE_HPP


#include <vector>
#include <cmath>
#include <cassert>
#include <string>
#include <algorithm> //std::sort, std::random_shuffle
#include <random>
#include <array>

#define KNP_NFACTORS 4

/**
 * Class represesnting an instance of the Knapsack Problem
 * (as defined in the K-Adaptability paper)
 * 
 * Note: All arrays and matrices are 1-indexed
 */
struct KNP {
	/** # of items */
	int N;

	/** Cost of items */
	std::vector<double> cost;

	/** Budget (size of knapsack) */
	double B;

	/** Profit of items */
	std::vector<double> profit;

	/** Late-stage profit penalty */
	double theta;

	/** Option of taking loans */
	bool loan;

	/** Cost of borrowing in the first-stage */
	double ell;

	/** Cost of borrowing in the second-stage */
	double lambda;

	/** Risk factor coefficients for costs */
	std::vector<std::array<double, 1 + KNP_NFACTORS> > phi;

	/** Risk factor coefficients for profits */
	std::vector<std::array<double, 1 + KNP_NFACTORS> > ksi;

	/** Solution file name */
	std::string solfilename;

};




/**
 * Construct an instance of the Knapsack Problem
 * @param data          (reference to) instance of Knapsack Problem
 * @param n             number of items
 * @param seed          seed that is unique to this instance
 * @param mixed_integer true if mixed-integer knapsack problem, otherwise pure integer
 */
static inline void gen_KNP(KNP& data, unsigned int n, int seed = 1, bool mixed_integer = false) {

	if (n < 5) {
		fprintf(stderr, "warning: N < 5 in Knapsack Problem. Setting N = 5.\n");
		n = 5;
	}

	// seed for re-producibility
	std::default_random_engine gen (1111 + seed);
	std::uniform_real_distribution<double> interval_1 (0.0, 1.0);
	std::uniform_real_distribution<double> interval_10 (0.0, 10.0);

	// # of items
	data.N = n;
	data.B = 0.0;
	data.theta = 0.8;
	data.loan = mixed_integer;
	data.ell = 0.12;
	data.lambda = 1.2;

	// Cost of items
	data.cost.assign(1 + data.N, 0.0);
	for (n = 1; (int)n <= data.N; n++) {
		data.cost[n] = interval_10(gen);
		data.B += data.cost[n];
	}

	// Budget (size of knapsack)
	data.B /= 2.0;


	// Profit of items
	data.profit.assign(1 + data.N, 0.0);
	for (n = 1; (int)n <= data.N; n++) {
		data.profit[n] = data.cost[n] / 5.0;
	}


	// Risk factor coefficients
	data.phi.resize(1 + data.N);
	data.ksi.resize(1 + data.N);
	for (n = 1; (int)n <= data.N; n++) {
		// generate F-1 numbers between 0 and 1
		std::array<double, KNP_NFACTORS - 1> xn, yn;
		for (size_t f = 0; f < xn.size(); f++) {
			xn[f] = interval_1(gen);
			yn[f] = interval_1(gen);
		}

		// sort these numbers
		std::sort(xn.begin(), xn.end());
		std::sort(yn.begin(), yn.end());

		// generate F numbers such that their sum = 1
		data.phi[n][0] = 0.0;
		data.phi[n][1] = xn[0];
		data.ksi[n][0] = 0.0;
		data.ksi[n][1] = yn[0];
		for (size_t f = 2; f < KNP_NFACTORS; f++) {
			data.phi[n][f] = xn[f-1] - xn[f-2];
			data.ksi[n][f] = yn[f-1] - yn[f-2];
		}
		data.phi[n][KNP_NFACTORS] = 1 - xn[KNP_NFACTORS-2];
		data.ksi[n][KNP_NFACTORS] = 1 - yn[KNP_NFACTORS-2];

		// randomize in order to eliminate bias
		std::random_shuffle(data.phi[n].begin() + 1, data.phi[n].end());
		std::random_shuffle(data.ksi[n].begin() + 1, data.ksi[n].end());

		// check (just to be safe)
		#ifndef NDEBUG
		double check_phi = 0.0;
		double check_ksi = 0.0;
		for (size_t f = 1; f <= KNP_NFACTORS; f++) {
			check_phi += data.phi[n][f];
			check_ksi += data.ksi[n][f];
		}
		assert(std::abs(1 - check_phi) <= 1.E-12);
		assert(std::abs(1 - check_ksi) <= 1.E-12);
		#endif
	}

	// set solution file name
	data.solfilename = "KNP-n" + std::to_string(data.N) + "-s" + std::to_string(seed) + "-t" + std::to_string(data.loan) + ".opt";
}
#undef KNP_NFACTORS
#endif