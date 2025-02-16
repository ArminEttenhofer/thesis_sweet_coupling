/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *
 *  	2017-07-13: Updated and validated by P. Peixoto.
 *  	2017-05-29: Based on source swe_cart2d.cpp
 *					which was also written by Pedro Peixoto
 */

#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_ln_erk.hpp>


/*
 * Main routine for method to be used in case of finite differences
 *
 * - A-Grid with spectral or FD spatial discretizations
 * - C-Grid with energy conserving FD scheme for spatial discretizations
 *
 */
void SWE_Cart2D_TS_ln_erk::euler_timestep_update(
		const sweet::Data::Cart2D::DataSpectral &i_h,	//!< prognostic variables (perturbed part of height)
		const sweet::Data::Cart2D::DataSpectral &i_u,	//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_v,	//!< prognostic variables

		sweet::Data::Cart2D::DataSpectral &o_h_t,	//!< time updates
		sweet::Data::Cart2D::DataSpectral &o_u_t,	//!< time updates
		sweet::Data::Cart2D::DataSpectral &o_v_t,	//!< time updates

		double i_simulation_timestamp
)
{
	// A-grid method
	if (!shackCart2DDataOps->space_grid_use_c_staggering)
	{
		/*
		 * non-conservative (advective) formulation:
		 *
		 *	h_t = -(u*h)_x - (v*h)_y
		 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
		 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
		 */

		sweet::Data::Cart2D::DataSpectral total_h = i_h + shackPDESWECart2D->h0;

		o_u_t = -shackPDESWECart2D->gravitation*ops->diff_c_x(total_h) - i_u*ops->diff_c_x(i_u) - i_v*ops->diff_c_y(i_u);
		o_v_t = -shackPDESWECart2D->gravitation*ops->diff_c_y(total_h) - i_u*ops->diff_c_x(i_v) - i_v*ops->diff_c_y(i_v);

		o_u_t += shackPDESWECart2D->cart2d_rotating_f0*i_v;
		o_v_t -= shackPDESWECart2D->cart2d_rotating_f0*i_u;

		// standard update
		/*
		 * P UPDATE
		 */
		if (!use_only_linear_divergence){ //full nonlinear divergence
			// standard update
			//o_h_t = -ops->diff_f_x(U) - ops->diff_f_y(V);
			o_h_t = -ops->diff_c_x(i_u*total_h) - ops->diff_c_y(i_v*total_h);
		}
		else // use linear divergence
		{
			//o_h_t = -ops->diff_f_x(shackPDESWECart2D->h0*i_u) - ops->diff_f_y(shackPDESWECart2D->h0*i_v);
			o_h_t = -i_u*ops->diff_c_x(total_h) - i_v*ops->diff_c_y(total_h) + //nonlinear adv
					-ops->diff_c_x(i_u*shackPDESWECart2D->h0) - ops->diff_c_y(i_v*shackPDESWECart2D->h0); //linear div
		}

	}
	else // shackDict.disc.use_staggering = true
	{
		// STAGGERED GRID

		sweet::Data::Cart2D::DataSpectral U(i_h.cart2DDataConfig); // U flux
		sweet::Data::Cart2D::DataSpectral V(i_h.cart2DDataConfig); // V flux
		sweet::Data::Cart2D::DataSpectral H(i_h.cart2DDataConfig); //Bernoulli potential

		sweet::Data::Cart2D::DataGrid U_phys(i_h.cart2DDataConfig); // U flux
		sweet::Data::Cart2D::DataGrid V_phys(i_h.cart2DDataConfig); // V flux
		sweet::Data::Cart2D::DataGrid H_phys(i_h.cart2DDataConfig); //Bernoulli potential

		sweet::Data::Cart2D::DataGrid i_u_phys = i_u.toGrid();
		sweet::Data::Cart2D::DataGrid i_v_phys = i_v.toGrid();

		sweet::Data::Cart2D::DataGrid total_h_phys = i_h.toGrid() + shackPDESWECart2D->h0;
		sweet::Data::Cart2D::DataSpectral total_h(i_h.cart2DDataConfig);
		total_h.loadCart2DDataGrid(total_h_phys);


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
		/*
		 * U and V updates
		 */


		U_phys = ops->avg_b_x(total_h_phys)*i_u_phys;
		V_phys = ops->avg_b_y(total_h_phys)*i_v_phys;
		H_phys = shackPDESWECart2D->gravitation*total_h_phys + 0.5*(ops->avg_f_x(i_u_phys*i_u_phys) + ops->avg_f_y(i_v_phys*i_v_phys));

		U.loadCart2DDataGrid(U_phys);
		V.loadCart2DDataGrid(V_phys);
		H.loadCart2DDataGrid(H_phys);


		// Potential vorticity
		sweet::Data::Cart2D::DataGrid total_h_pv_phys = total_h_phys;
		sweet::Data::Cart2D::DataSpectral total_h_pv = total_h_phys(i_h.cart2DDataConfig);
		total_h_pv_phys = ops->avg_b_x(ops->avg_b_y(total_h_phys));
		total_h_pv.loadCart2DDataGrid(total_h_pv_phys);

#if 0
		if (total_h_pv.reduce_min() < 0.00000001)
		{
			std::cerr << "Test case not adequate for vector invariant formulation. Null or negative water height" << std::endl;
			std::cerr << "Min h_pv   : " << total_h_pv.reduce_min() << std::endl;
			std::cerr << "Min h_total: " << total_h.reduce_min() << std::endl;
			std::cerr << "Min h_pert : " << i_h.reduce_min() << std::endl;
			SWEETErrorFatal("SWE_Cart2D_TS_ln_erk: Methods unstable or inadequate for vector invariant swe");;
		}
#endif

		sweet::Data::Cart2D::DataSpectral q = (ops->diff_b_x(i_v) - ops->diff_b_y(i_u) + shackPDESWECart2D->cart2d_rotating_f0) / total_h_pv;
		sweet::Data::Cart2D::DataGrid q_phys = q.toGrid();

		// u, v tendencies
		// Energy conserving scheme
		sweet::Data::Cart2D::DataGrid o_u_t_phys = ops->avg_f_y(q_phys*ops->avg_b_x(V_phys)) - ops->diff_b_x(H).toGrid();
		sweet::Data::Cart2D::DataGrid o_v_t_phys = -ops->avg_f_x(q_phys*ops->avg_b_y(U_phys)) - ops->diff_b_y(H).toGrid();
		o_u_t.loadCart2DDataGrid(o_u_t_phys);
		o_v_t.loadCart2DDataGrid(o_v_t_phys);

		/*
		 * P UPDATE
		 */
		if (!use_only_linear_divergence){ //full nonlinear divergence
			// standard update
			o_h_t = -ops->diff_f_x(U) - ops->diff_f_y(V);
		}
		else // use linear divergence
		{
			o_h_t = -i_u*ops->diff_f_x(total_h) - i_v*ops->diff_f_y(total_h) + //nonlinear adv
					-ops->diff_f_x(i_u*shackPDESWECart2D->h0) - ops->diff_f_y(i_v*shackPDESWECart2D->h0); //linear div
			//o_h_t = -ops->diff_f_x(shackPDESWECart2D->h0*i_u) - ops->diff_f_y(shackPDESWECart2D->h0*i_v);
		}
	}
}



void SWE_Cart2D_TS_ln_erk::runTimestep(
		sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETErrorFatal("SWE_Cart2D_TS_ln_erk: Only constant time step size allowed");

	// standard time stepping
	timestepping_rk.runTimestep(
			this,
			&SWE_Cart2D_TS_ln_erk::euler_timestep_update,	//!< pointer to function to compute euler time step updates
			io_h, io_u, io_v,
			i_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
bool SWE_Cart2D_TS_ln_erk::setup(
		sweet::Data::Cart2D::Operators *io_ops
)
{
	PDESWECart2DTS_BaseInterface::setup(io_ops);

	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	use_only_linear_divergence = shackPDESWECart2D->use_only_linear_divergence;

	return true;
}

