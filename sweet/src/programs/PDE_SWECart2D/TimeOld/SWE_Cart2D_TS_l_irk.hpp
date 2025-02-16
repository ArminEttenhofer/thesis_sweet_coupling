/*
 * SWE_Cart2D_TS_l_irk.hpp
 *
 *  Created on: 29 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWECART2D_TIME_SWE_CART2D_TS_L_IRK_HPP
#define PROGRAMS_PDE_SWECART2D_TIME_SWE_CART2D_TS_L_IRK_HPP

#include <programs/PDE_SWECart2D/TimeOld/PDESWECart2DTS_BaseInterface.hpp>
#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKCart2DData.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>

#if !SWEET_USE_CART2D_SPECTRAL_SPACE
	#include <sweet/Data/Cart2D/Cart2DOperatorsComplex.hpp>
#endif


class SWE_Cart2D_TS_l_irk	: public PDESWECart2DTS_BaseInterface
{
#if !SWEET_USE_CART2D_SPECTRAL_SPACE
	sweet::Data::Cart2D::OperatorsComplex opComplex;
#endif

	int timestepping_order;

public:
	bool setup(
			sweet::Data::Cart2D::Operators *io_ops,
			int i_order
	);

	bool setup(
			sweet::Data::Cart2D::Operators *io_ops
	);

	void runTimestep(
			sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

	virtual ~SWE_Cart2D_TS_l_irk() {}
};

#endif
