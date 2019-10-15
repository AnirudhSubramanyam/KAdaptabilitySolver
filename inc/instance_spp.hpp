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


#ifndef SPP_INSTANCE_HPP
#define SPP_INSTANCE_HPP


#include <vector>
#include <cmath>
#include <string>
#include <cassert>
#include <random>

/**
 * Class represesnting an instance of the Shortest Path Problem
 * (as defined in the K-Adaptability paper)
 * 
 * Note: All arrays and matrices are 1-indexed
 */
struct SPP {
	/** # of nodes in graph */
	int N;

	/** # of arcs in graph */
	int A;

	/** Arc-weights (an entry of -1 indicates the arc does not exist) */
	std::vector<std::vector<double> > costs;

	/** start node */
	int src;

	/** end node */
	int tgt;

	/** Solution file name */
	std::string solfilename;

};




/**
 * Construct an instance of the Shortest Path Problem
 * @param data (reference to) instance of Shortest Path Problem
 * @param n    number of nodes
 * @param seed seed that is unique to this instance
 */
static inline void gen_SPP(SPP& data, unsigned int n, int seed = 1) {
	
	if (n < 20) {
		fprintf(stderr, "warning: N < 20 in Shortest Path Problem. Setting N = 20.\n");
		n = 20;
	}

	// seed for re-producibility
	std::default_random_engine gen (1111 + seed);

	// # of nodes, # of arcs
	data.N = n;
	data.A = n * (n - 1);

	// costs(i,j) = weight of arc from i to j for i, j = 1, ..., N
	data.costs.assign(1 + data.N, std::vector<double>(1 + data.N, -1.0));

	// coordinates(i) = co-ordinates of the i'th node
	std::vector<std::vector<double> > coordinates (1);
	std::uniform_real_distribution<double> interval (0.0, 10.0);
	for (n = 1; (int)n <= data.N; ++n) {
		double xn = interval(gen);
		double yn = interval(gen);
		coordinates.emplace_back(std::vector<double>{xn, yn});
	}

	// compute arc-weights
	for (int i = 1; i <= data.N; ++i) for (int j = i + 1; j <= data.N; ++j) {
		double xd   = std::abs(coordinates[i][0] - coordinates[j][0]);
		double yd   = std::abs(coordinates[i][1] - coordinates[j][1]);
		data.costs[i][j] = std::sqrt((xd * xd) + (yd * yd));
		data.costs[j][i] = data.costs[i][j];
	}

	// delete 70\% of the arcs
	int D = (7 * data.A) / 10;
	for (int count = 1; count <= D; ++count) {
		
		// Compute max element in cost matrix
		double max = -1.0;
		int s = 0;
		int t = 0;
		for (int i = 1; i <= data.N; ++i) {
			data.costs[i][i] = -1.0;
			for (int j = 1; j <= data.N; ++j) if (i != j) {
				if (data.costs[i][j] > max) {
					max = data.costs[i][j];
					s = i;
					t = j;
				}
			}
		}
		assert(max >= 0.0);
		assert(s > 0 && t > 0);
		assert(s != t);

		// Kill this arc by setting cost = -1
		data.costs[s][t] = -1.0;
		data.A--;

		// Source and terminal nodes
		if (count == 1) {
			data.src = s;
			data.tgt = t;
		}
	}
	assert(data.src >= 1 && data.src <= data.N);
	assert(data.tgt >= 1 && data.tgt <= data.N);
	assert(data.src != data.tgt);

	// set solution file name
	data.solfilename = "SPP-n" + std::to_string(data.N) + "-s" + std::to_string(seed) + ".opt";
}

#endif