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


#ifndef PSP_INSTANCE_HPP
#define PSP_INSTANCE_HPP


#include <vector>
#include <cmath>
#include <string>
#include <random>
#include <iterator>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iostream>
#include <algorithm> //std::sort, std::random_shuffle
#include <array>

#define PSP_NFACTORS 4

/**
 * Class represesnting an instance of the Project Scheduling Problem
 * (as defined in the K-Adaptability paper)
 * 
 * Note: All arrays and matrices are 1-indexed
 */
struct PSP {
	/** # of nodes in project graph = # of activities including dummy start and end nodes */
	int N;

	/** Budget of resources */
	double B;

	/** Use affine decision rules */
	bool affine;

	/** Special case of project graphs */
	bool special;

	/** Activity durations */
	std::vector<double> duration;

	/** All nodes in the graph which succeed this node */
	std::vector<std::vector<int> > successor;

	/** Risk factor coefficients */
	std::vector<std::array<double, 1 + PSP_NFACTORS> > phi;

	/** Solution file name */
	std::string solfilename;

};



/**
 * Construct an instance of the Project Scheduling Problem
 * See Example 2.2 in "Robust Resource Allocations in Temporal Networks", Math. Prog. (2012)
 * @param data (reference to) instance of Project Scheduling Problem
 * @param k    parameter to generate instances of Example 2.2
 * @param seed seed that is unique to this instance
 */
static inline void gen_PSP(PSP& data, unsigned int k, int seed = 1) {
	
	if (k == 0) {
		fprintf(stderr, "warning: k = 0 in Resource Allocation Problem (special case). Setting k = 1.\n");
		k = 1;
	}

	// initialize
	data.N = (3 * k) + 1;
	data.B = 0;
	data.affine = (seed != 0);
	data.special = true;
	data.duration.resize(1 + data.N, 0.0);
	data.successor.resize(1 + data.N);
	data.phi.resize(data.N);
	data.solfilename = "PSP-n" + std::to_string(data.N) + "-s" + std::to_string(seed) + "-d0.opt";

	// generate project graph
	for (unsigned int l = 0; l < k; l++) for (unsigned int p = 2; p <= 3; p++) {
		data.successor[(3 * l) + 1].emplace_back((3 * l) + p);
		data.successor[(3 * l) + p].emplace_back((3 * l) + 4);
	}

	// set durations
	for (unsigned int l = 0; l < k; l++) {
		data.duration[(3 * l) + 1] = 0;
		data.duration[(3 * l) + 2] = 1;
		data.duration[(3 * l) + 3] = -1;
	}
	data.duration[data.N] = 0;

	// dummy end node has no successors
	assert(data.successor[data.N].empty());
}
#undef PSP_NFACTORS
#endif