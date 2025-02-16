/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_cart2d.cpp
 *					which was also written by Pedro Peixoto
 */

#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_rexi_n_erk.hpp>
#include <sweet/Data/Cart2DComplex/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>



bool SWE_Cart2D_TS_l_rexi_n_erk::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{

	ts_l_rexi.shackRegistration(io_shackDict);
	PDESWECart2DTS_BaseInterface::shackRegistration(io_shackDict);

	return true;
}

/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Cart2D_TS_l_rexi_n_erk::euler_timestep_update_nonlinear(
		const sweet::Data::Cart2D::DataSpectral &i_h,	//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_u,	//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_v,	//!< prognostic variables

		sweet::Data::Cart2D::DataSpectral &o_h_t,	//!< time updates
		sweet::Data::Cart2D::DataSpectral &o_u_t,	//!< time updates
		sweet::Data::Cart2D::DataSpectral &o_v_t,	//!< time updates

		double i_timestamp
)
{
	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */
	//o_h_t = -ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);
	o_u_t = -i_u*ops->diff_c_x(i_u) - i_v*ops->diff_c_y(i_u);
	o_v_t = -i_u*ops->diff_c_x(i_v) - i_v*ops->diff_c_y(i_v);
	if (use_only_linear_divergence) //only nonlinear advection left to solve
		o_h_t = - (i_u*ops->diff_c_x(i_h) + i_v*ops->diff_c_y(i_h));
	else //full nonlinear equation on h
		o_h_t = -ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);

}


void SWE_Cart2D_TS_l_rexi_n_erk::runTimestep(
		sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETErrorFatal("SWE_Cart2D_TS_l_rexi_n_erk: Only constant time step size allowed");

	if (timestepping_order_nonlinear == 1)
	{
		ts_l_rexi.runTimestep(
				io_h, io_u, io_v,
				i_dt,
				i_simulation_timestamp
			);

		// standard time stepping
		timestepping_rk.runTimestep(
				this,
				&SWE_Cart2D_TS_l_rexi_n_erk::euler_timestep_update_nonlinear,	//!< pointer to function to compute euler time step updates
				io_h, io_u, io_v,
				i_dt,
				timestepping_order_nonlinear,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order_nonlinear == 2)
	{
		ts_l_rexi.runTimestep(
				io_h, io_u, io_v,
				i_dt*0.5,
				i_simulation_timestamp
			);

		// standard time stepping
		timestepping_rk.runTimestep(
				this,
				&SWE_Cart2D_TS_l_rexi_n_erk::euler_timestep_update_nonlinear,	//!< pointer to function to compute euler time step updates
				io_h, io_u, io_v,
				i_dt,
				timestepping_order_nonlinear,
				i_simulation_timestamp
			);

		ts_l_rexi.runTimestep(
				io_h, io_u, io_v,
				i_dt*0.5,
				i_simulation_timestamp
			);
	}
	else
	{
		SWEETErrorFatal("SWE_Cart2D_TS_l_rexi_n_erk: Explicit erk order not implemented for this scheme, please set --timestepping-order2 to 1 or 2.");
	}
}



/*
 * Setup
 */
bool SWE_Cart2D_TS_l_rexi_n_erk::setup(
	sweet::Data::Cart2D::Operators *io_ops
)
{
	PDESWECart2DTS_BaseInterface::setup(io_ops);

	use_only_linear_divergence = shackPDESWECart2D->use_only_linear_divergence;

	ts_l_rexi.setup(io_ops, "phi0");

	timestepping_order_nonlinear = shackPDESWETimeDisc->timestepping_order;
	timestepping_rk.setupBuffers(ops->cart2DDataConfig, timestepping_order_nonlinear);

	if (shackCart2DDataOps->space_grid_use_c_staggering)
		SWEETErrorFatal("Staggering not supported for l_rexi_n_erk");

	return true;
}

