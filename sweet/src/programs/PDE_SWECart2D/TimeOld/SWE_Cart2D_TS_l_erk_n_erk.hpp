/*
 * SWE_Cart2D_TS_ln_erk.hpp
 *
 *  Created on: 29 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWECART2D_TIME_SWE_CART2D_TS_L_ERK_N_ERK_HPP
#define PROGRAMS_PDE_SWECART2D_TIME_SWE_CART2D_TS_L_ERK_N_ERK_HPP

#include <programs/PDE_SWECart2D/TimeOld/PDESWECart2DTS_BaseInterface.hpp>
#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKCart2DData.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>

class SWE_Cart2D_TS_l_erk_n_erk	: public PDESWECart2DTS_BaseInterface
{
	int timestepping_order;
	int timestepping_order2;
	bool use_only_linear_divergence;

	sweet::DEPRECATED_TimesteppingExplicitRKCart2DData timestepping_rk_linear;
	sweet::DEPRECATED_TimesteppingExplicitRKCart2DData timestepping_rk_nonlinear;

private:
	void euler_timestep_update_linear(
			const sweet::Data::Cart2D::DataSpectral &i_h,	//!< prognostic variables
			const sweet::Data::Cart2D::DataSpectral &i_u,	//!< prognostic variables
			const sweet::Data::Cart2D::DataSpectral &i_v,	//!< prognostic variables

			sweet::Data::Cart2D::DataSpectral &o_h_t,	//!< time updates
			sweet::Data::Cart2D::DataSpectral &o_u_t,	//!< time updates
			sweet::Data::Cart2D::DataSpectral &o_v_t,	//!< time updates

			double i_simulation_timestamp
	);


private:
	void euler_timestep_update_nonlinear(
			const sweet::Data::Cart2D::DataSpectral &i_h,	//!< prognostic variables
			const sweet::Data::Cart2D::DataSpectral &i_u,	//!< prognostic variables
			const sweet::Data::Cart2D::DataSpectral &i_v,	//!< prognostic variables

			sweet::Data::Cart2D::DataSpectral &o_h_t,	//!< time updates
			sweet::Data::Cart2D::DataSpectral &o_u_t,	//!< time updates
			sweet::Data::Cart2D::DataSpectral &o_v_t,	//!< time updates

			double i_simulation_timestamp
	);


public:
	bool setup(
			sweet::Data::Cart2D::Operators *io_ops
	);

	void runTimestep(
			sweet::Data::Cart2D::DataSpectral &io_h_pert,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,		//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,		//!< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

	virtual ~SWE_Cart2D_TS_l_erk_n_erk() {}
};

#endif
