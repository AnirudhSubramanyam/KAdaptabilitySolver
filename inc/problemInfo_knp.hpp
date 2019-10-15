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


#ifndef KADAPTABLEINFO_KNP_HPP
#define KADAPTABLEINFO_KNP_HPP

#include "instance_knp.hpp"
#include "problemInfo.hpp"

class KAdaptableInfo_KNP : public KAdaptableInfo {
protected:
	/** Specific instance data */
	KNP data;

	/**
	 * Define 1st-stage and 2nd-stage variables
	 */
	void makeVars() override;

	/**
	 * Define uncertainty set
	 */
	void makeUncSet() override;

	/**
	 * Define constraints involving 1st-stage variables only
	 */
	void makeConsX() override;
	
	/**
	 * Define constraints involving 2nd-stage variables
	 * @param k policy number to make constraints for
	 */
	void makeConsY(unsigned int k = 0) override;

public:

	/**
	 * Set data of the instance that this object refers to
	 * @param data problem instance
	 */
	void setInstance(const KNP& data);

	/**
	 * Virtual (default) constructor
	 */
	inline KAdaptableInfo_KNP* create() const override {
		return new KAdaptableInfo_KNP();
	}

	/**
	 * Virtual (copy) constructor
	 */
	inline KAdaptableInfo_KNP* clone() const override {
		return new KAdaptableInfo_KNP (*this);
	}
	


};

#endif