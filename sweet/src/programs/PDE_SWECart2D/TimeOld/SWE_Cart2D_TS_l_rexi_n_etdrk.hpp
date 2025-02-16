/*
 * SWE_Cart2D_TS_l_rexi_n_etdrk.hpp
 *
 *  Created on: 29 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWECART2D_TIME_SWE_CART2D_TS_L_REXI_N_ETDRK_HPP
#define PROGRAMS_PDE_SWECART2D_TIME_SWE_CART2D_TS_L_REXI_N_ETDRK_HPP

#include <programs/PDE_SWECart2D/TimeOld/PDESWECart2DTS_BaseInterface.hpp>
#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_rexi.hpp>
#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKCart2DData.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>


class SWE_Cart2D_TS_l_rexi_n_etdrk	:
		public PDESWECart2DTS_BaseInterface
{
	SWE_Cart2D_TS_l_rexi ts_phi0_rexi;
	SWE_Cart2D_TS_l_rexi ts_phi1_rexi;
	SWE_Cart2D_TS_l_rexi ts_phi2_rexi;

	SWE_Cart2D_TS_l_rexi ts_ups0_rexi;
	SWE_Cart2D_TS_l_rexi ts_ups1_rexi;
	SWE_Cart2D_TS_l_rexi ts_ups2_rexi;
	SWE_Cart2D_TS_l_rexi ts_ups3_rexi;

	int timestepping_order;
	bool use_only_linear_divergence;


public:
	bool setup(
			sweet::Data::Cart2D::Operators *io_ops
	);

	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	);

	void euler_timestep_update_nonlinear(
			const sweet::Data::Cart2D::DataSpectral &i_h,	//!< prognostic variables
			const sweet::Data::Cart2D::DataSpectral &i_u,	//!< prognostic variables
			const sweet::Data::Cart2D::DataSpectral &i_v,	//!< prognostic variables

			sweet::Data::Cart2D::DataSpectral &o_h_t,	//!< time updates
			sweet::Data::Cart2D::DataSpectral &o_u_t,	//!< time updates
			sweet::Data::Cart2D::DataSpectral &o_v_t,	//!< time updates

			double i_timestamp
	);

	void runTimestep(
			sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Cart2D_TS_l_rexi_n_etdrk();
};

#endif
