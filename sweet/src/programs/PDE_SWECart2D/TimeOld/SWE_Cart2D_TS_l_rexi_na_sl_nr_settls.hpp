/*
 * SWE_Cart2D_TS_l_rexi_na_sl_nd_settls.hpp
 *
 *  Created on: 29 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_cart2d.cpp
 *					which was also written by Pedro Peixoto
 */

#ifndef PROGRAMS_PDE_SWECART2D_TIME_SWE_CART2D_TS_L_REXI_NA_SL_NR_SETTLS_HPP
#define PROGRAMS_PDE_SWECART2D_TIME_SWE_CART2D_TS_L_REXI_NA_SL_NR_SETTLS_HPP

#include <programs/PDE_SWECart2D/TimeOld/PDESWECart2DTS_BaseInterface.hpp>
#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_direct.hpp>
#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_rexi.hpp>
#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKCart2DData.hpp>
#include <sweet/Data/Cart2D/DataSampler.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2D/Operators.hpp>
#include <sweet/SemiLagrangian/Cart2D.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>



class SWE_Cart2D_TS_l_rexi_na_sl_nr_settls	: public PDESWECart2DTS_BaseInterface
{



	SWE_Cart2D_TS_l_rexi ts_l_rexi;

	//int with_linear_div_only;

	sweet::SemiLagrangian::Cart2D semiLagrangian;
	sweet::Data::Cart2D::DataSampler sampler2D;

	//Previous values (t_n-1)
	sweet::Data::Cart2D::DataSpectral h_prev, u_prev, v_prev;

	// Arrival points for semi-lag
	sweet::Data::Vector::Vector<double> posx_a, posy_a;

	// Departure points for semi-lag
	sweet::Data::Vector::Vector<double> posx_d, posy_d;

	bool use_only_linear_divergence;

public:
	bool setup(
			sweet::Data::Cart2D::Operators *io_ops
	);

	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	);

	void runTimestep(
			sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

#if ( SWEET_PARAREAL && SWEET_PARAREAL_CART2D ) || ( SWEET_XBRAID && SWEET_XBRAID_CART2D )
	void set_previous_solution(
				sweet::Data::Cart2D::DataSpectral &i_h_prev,
				sweet::Data::Cart2D::DataSpectral &i_u_prev,
				sweet::Data::Cart2D::DataSpectral &i_v_prev
	) override
	{
		h_prev = i_h_prev;
		u_prev = i_u_prev;
		v_prev = i_v_prev;
	}
#endif

	virtual ~SWE_Cart2D_TS_l_rexi_na_sl_nr_settls() {}
};

#endif
