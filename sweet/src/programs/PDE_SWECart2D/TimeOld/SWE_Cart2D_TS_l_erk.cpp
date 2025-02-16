/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_cart2d.cpp
 *					which was also written by Pedro Peixoto
 */

#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_erk.hpp>


/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Cart2D_TS_l_erk::euler_timestep_update(
		const sweet::Data::Cart2D::DataSpectral &i_h,	//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_u,	//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_v,	//!< prognostic variables

		sweet::Data::Cart2D::DataSpectral &o_h_t,	//!< time updates
		sweet::Data::Cart2D::DataSpectral &o_u_t,	//!< time updates
		sweet::Data::Cart2D::DataSpectral &o_v_t,	//!< time updates

		double i_simulation_timestamp
)
{
	// A- grid method
	if (!shackCart2DDataOps->space_grid_use_c_staggering)
	{
		/*
		 * linearized non-conservative (advective) formulation:
		 *
		 * h_t = -h0*u_x - h0*v_ym
		 * u_t = -g * h_x + f*v
		 * v_t = -g * h_y - f*u
		 */

#if 1
		o_u_t = -shackPDESWECart2D->gravitation*ops->diff_c_x(i_h) + shackPDESWECart2D->cart2d_rotating_f0*i_v;
		o_v_t = -shackPDESWECart2D->gravitation*ops->diff_c_y(i_h) - shackPDESWECart2D->cart2d_rotating_f0*i_u;

		// standard update
		o_h_t = -(ops->diff_c_x(i_u) + ops->diff_c_y(i_v))*shackPDESWECart2D->h0;
#else

	#if 0
		// U-only
		o_u_t = -shackPDESWECart2D->gravitation*ops->diff_c_x(i_h) + shackPDESWECart2D->cart2d_rotating_f0*i_v;
		//o_v_t.grid_setZero();
		o_v_t = - shackPDESWECart2D->cart2d_rotating_f0*i_u;

		// standard update
		o_h_t = -(ops->diff_c_x(i_u))*shackPDESWECart2D->h0;

	#else
		// V-only
		//o_u_t.spectral_setZero();
		o_u_t = +shackPDESWECart2D->cart2d_rotating_f0*i_v;
		o_v_t = -shackPDESWECart2D->gravitation*ops->diff_c_y(i_h) - shackPDESWECart2D->cart2d_rotating_f0*i_u;// - shackDict.sim.f0*i_u;

		// standard update
		o_h_t = -(ops->diff_c_y(i_v))*shackPDESWECart2D->h0;
	#endif

#endif
	}
	else // shackDict.disc.use_staggering = true
	{
		// STAGGERED GRID

		/*
		 * Sadourny energy conserving scheme
		 *
		 * Note, that this grid does not follow the formulation
		 * in the paper of Robert Sadourny, but looks as follows:
		 *
		 *              ^
		 *              |
		 *       ______v0,1_____
		 *       |             |
		 *       |			   |
		 *       |             |
		 *  u0,0 |->  H/P0,0   |u1,0 ->
		 *(0,0.5)|			   |
		 *       |      ^      |
		 *   q0,0|______|______|
		 * (0,0)      v0,0
		 *           (0.5,0)
		 *
		 * V_t + q N x (P V) + grad( g P + 1/2 V*V) = 0
		 * P_t + div(P V) = 0
		 */

		sweet::Data::Cart2D::DataSpectral H = shackPDESWECart2D->gravitation*i_h;// + 0.5*(ops->avg_f_x(i_u*i_u) + ops->avg_f_y(i_v*i_v));

		sweet::Data::Cart2D::DataGrid o_u_t_phys(o_u_t.cart2DDataConfig);
		sweet::Data::Cart2D::DataGrid o_v_t_phys(o_v_t.cart2DDataConfig);
		o_u_t_phys = ops->avg_f_y(shackPDESWECart2D->cart2d_rotating_f0*ops->avg_b_x(i_v.toGrid())) - ops->diff_b_x(H).toGrid();
		o_v_t_phys = -ops->avg_f_x(shackPDESWECart2D->cart2d_rotating_f0*ops->avg_b_y(i_u.toGrid())) - ops->diff_b_y(H).toGrid();
		o_u_t.loadCart2DDataGrid(o_u_t_phys);
		o_v_t.loadCart2DDataGrid(o_v_t_phys);

		/*
		 * P UPDATE
		 */
		o_h_t = -ops->diff_f_x(shackPDESWECart2D->h0*i_u) - ops->diff_f_y(shackPDESWECart2D->h0*i_v);
	}
}



void SWE_Cart2D_TS_l_erk::runTimestep(
		sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETErrorFatal("SWE_Cart2D_TS_l_erk: Only constant time step size allowed (please set --dt)");

	// standard time stepping
	timestepping_rk.runTimestep(
			this,
			&SWE_Cart2D_TS_l_erk::euler_timestep_update,	//!< pointer to function to compute euler time step updates
			io_h, io_u, io_v,
			i_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
bool SWE_Cart2D_TS_l_erk::setup(
		sweet::Data::Cart2D::Operators *io_ops,
		int i_order
)
{
	PDESWECart2DTS_BaseInterface::setup(io_ops);

	SWEET_ASSERT(i_order > 0);
	timestepping_order = i_order;

	timestepping_rk.setupBuffers(ops->cart2DDataConfig, timestepping_order);

	//if (shackDict.disc.use_staggering)
	//	SWEETErrorFatal("Staggering not supported for l_erk");
	return true;
}

bool SWE_Cart2D_TS_l_erk::setup(
		sweet::Data::Cart2D::Operators *io_ops
)
{
	return setup(io_ops, shackPDESWETimeDisc->timestepping_order);
}
