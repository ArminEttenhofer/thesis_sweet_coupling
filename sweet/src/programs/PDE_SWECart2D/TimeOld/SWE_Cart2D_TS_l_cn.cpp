/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_cart2d.cpp
 *					which was also written by Pedro Peixoto
 *
 */

#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_cn.hpp>



bool SWE_Cart2D_TS_l_cn::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	ts_l_erk.shackRegistration(io_shackDict);
	ts_l_irk.shackRegistration(io_shackDict);
	PDESWECart2DTS_BaseInterface::shackRegistration(io_shackDict);

	return true;
}



/*
 * This is a Crank-Nicolson scheme for linear equation
 *
 * 1) Takes an explicit 1/2 step (or whatever controlled by crank_nicolson_damping_factor) of explicit euler
 * 2) Then, takes an implicit 1/2 step (or whatever controlled by crank_nicolson_damping_factor) with and implicit euler scheme
 *
 * With explicit Euler, may be viewed as classic CN:
 *
 * U_t = L U(n)
 *
 *
 * (U(n+1) - U(n)) / dt = 0.5*(L U(n+1) + L U(n))
 *
 * <=> U(n+1) - U(n) =  dt/2 *(L U(n+1) + L U(n))
 *
 * <=> (1-dt/2L) U(n+1)= (1 + dt/2 L) U(n)
 *     ---------------   -----------------
 *           |                  |
 *  (1/2 implicit euler)  (1/2 explicit euler)
 *
 * Comment added by P. Peixoto on 4 Sept 2017
 */

void SWE_Cart2D_TS_l_cn::runTimestep(
		sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETErrorFatal("SWE_Cart2D_TS_l_cn: Only constant time step size allowed (Please set --dt)");


	sweet::Data::Cart2D::DataSpectral h_linear_t1 = io_h;
	sweet::Data::Cart2D::DataSpectral u_linear_t1 = io_u;
	sweet::Data::Cart2D::DataSpectral v_linear_t1 = io_v;

	ts_l_erk.runTimestep(
			io_h,
			io_u,
			io_v,
			i_dt*(1.0-crank_nicolson_damping_factor),
			i_simulation_timestamp
		);

	ts_l_irk.runTimestep(
			io_h, io_u, io_v,
			i_dt*crank_nicolson_damping_factor,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
bool SWE_Cart2D_TS_l_cn::setup(
		sweet::Data::Cart2D::Operators *io_ops
)
{
	PDESWECart2DTS_BaseInterface::setup(io_ops);

	if (shackCart2DDataOps->space_grid_use_c_staggering)
		SWEETErrorFatal("Staggering not supported for l_cn");

	crank_nicolson_damping_factor = 0.5;

	ts_l_irk.setup(io_ops, 1);
	ts_l_erk.setup(io_ops, 1);

	return true;
}
