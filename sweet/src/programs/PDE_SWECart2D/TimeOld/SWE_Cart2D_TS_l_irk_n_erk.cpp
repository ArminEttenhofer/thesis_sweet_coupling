/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_cart2d.cpp
 *					which was also written by Pedro Peixoto
 */


#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_irk_n_erk.hpp>



bool SWE_Cart2D_TS_l_irk_n_erk::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	ts_l_irk.shackRegistration(io_shackDict);
	PDESWECart2DTS_BaseInterface::shackRegistration(io_shackDict);

	return true;
}

/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Cart2D_TS_l_irk_n_erk::euler_timestep_update_nonlinear(
		const sweet::Data::Cart2D::DataSpectral &i_h,	//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_u,	//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_v,	//!< prognostic variables

		sweet::Data::Cart2D::DataSpectral &o_h_t,	//!< time updates
		sweet::Data::Cart2D::DataSpectral &o_u_t,	//!< time updates
		sweet::Data::Cart2D::DataSpectral &o_v_t	//!< time updates
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


void SWE_Cart2D_TS_l_irk_n_erk::runTimestep(
		sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETErrorFatal("SWE_Cart2D_TS_l_irk_n_erk: Only constant time step size allowed");

	sweet::Data::Cart2D::DataSpectral h_linear_t1 = io_h;
	sweet::Data::Cart2D::DataSpectral u_linear_t1 = io_u;
	sweet::Data::Cart2D::DataSpectral v_linear_t1 = io_v;

	ts_l_irk.runTimestep(
			h_linear_t1, u_linear_t1, v_linear_t1,
			i_dt,
			i_simulation_timestamp
		);

	// compute non-linear tendencies at half time step
	sweet::Data::Cart2D::DataSpectral h_dt_nonlinear(ops->cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral u_dt_nonlinear(ops->cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral v_dt_nonlinear(ops->cart2DDataConfig);

	// standard time stepping
	euler_timestep_update_nonlinear(
			io_h, io_u, io_v,
			h_dt_nonlinear, u_dt_nonlinear, v_dt_nonlinear
		);

	io_h = h_linear_t1 + h_dt_nonlinear*i_dt;
	io_u = u_linear_t1 + u_dt_nonlinear*i_dt;
	io_v = v_linear_t1 + v_dt_nonlinear*i_dt;
}



/*
 * Setup
 */
bool SWE_Cart2D_TS_l_irk_n_erk::setup(
	sweet::Data::Cart2D::Operators *io_ops
)
{
	PDESWECart2DTS_BaseInterface::setup(io_ops);

	SWEET_ASSERT(io_ops != nullptr);
	SWEET_ASSERT(shackPDESWETimeDisc != nullptr);
	SWEET_ASSERT(shackPDESWECart2D != nullptr);

	timestepping_order_linear = shackPDESWETimeDisc->timestepping_order;
	use_only_linear_divergence = shackPDESWECart2D->use_only_linear_divergence;

	ts_l_irk.setup(io_ops, timestepping_order_linear);

	if (shackCart2DDataOps->space_grid_use_c_staggering)
		SWEETErrorFatal("Staggering not supported for l_irk_n_erk");


	if (timestepping_order_linear != 1)
		SWEETErrorFatal("SWE_Cart2D_TS_l_irk_n_erk: Only 1st order TS supported with this implementation. Please set --timestepping-order 1.");

	timestepping_order_nonlinear = timestepping_order_linear;
	timestepping_rk.setupBuffers(ops->cart2DDataConfig, timestepping_order_nonlinear);
	return true;
}

