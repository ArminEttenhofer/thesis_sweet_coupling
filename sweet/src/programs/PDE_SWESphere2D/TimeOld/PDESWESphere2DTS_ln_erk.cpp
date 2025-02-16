/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere2DTS_ln_erk.hpp"


bool PDESWESphere2DTS_ln_erk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup(io_ops, shackPDESWETimeDisc->timestepping_order);
}


bool PDESWESphere2DTS_ln_erk::setup(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_order	//!< order of RK time stepping method
)
{
	ops = io_ops,
	timestepping_order = i_order;

	setupFG();

	return true;
}


/*
 * Main routine for method to be used in case of finite differences
 */
void PDESWESphere2DTS_ln_erk::euler_timestep_update_pert(
		const sweet::Data::Sphere2D::DataSpectral &i_phi_pert,
		const sweet::Data::Sphere2D::DataSpectral &i_vrt,
		const sweet::Data::Sphere2D::DataSpectral &i_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_pert_t,	//!< time updates
		sweet::Data::Sphere2D::DataSpectral &o_vrt_t,	//!< time updates
		sweet::Data::Sphere2D::DataSpectral &o_div_t,	//!< time updates

		double i_simulation_timestamp
)
{
	double gh0 = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0;

	/*
	 * NON-LINEAR SWE
	 *
	 * See
	 * 	Williamson, David L., Drake, John B., Hack, James J., Jakob, Rudiger, & Swarztrauber, Paul N. (1992).
	 * 	A standard test set for numerical approximations to the shallow water equations in spherical geometry.
	 * 	Journal of Computational Physics, 102(1), 211–224. https://doi.org/10.1016/S0021-9991(05)80016-6
	 *
	 * "2.3 Vorticity/Divergence Form"
	 */

	/*
	 * See documentation in [sweet]/doc/swe/swe_sphere2d_formulation/
	 */
	sweet::Data::Sphere2D::DataGrid phi_pert_phys = i_phi_pert.toGrid();

	/*
	 * Step 1a
	 */
	sweet::Data::Sphere2D::DataGrid ug, vg;
	ops->vrtdiv_2_uv(i_vrt, i_div, ug, vg);

	/*
	 * Step 1b
	 */
	sweet::Data::Sphere2D::DataGrid vrtg = i_vrt.toGrid();

	/*
	 * Step 1c
	 */

	using namespace sweet;

	// left part of eq. (19)
	sweet::Data::Sphere2D::DataGrid u_nl = ug*(vrtg+fg);

	// left part of eq. (20)
	sweet::Data::Sphere2D::DataGrid v_nl = vg*(vrtg+fg);

	/*
	 * Step 1d
	 */
	// Eq. (21) & left part of Eq. (22)
	ops->uv_2_vrtdiv(u_nl, v_nl, o_div_t, o_vrt_t);


	/*
	 * Step 1e
	 */
	o_vrt_t *= -1.0;

	/*
	 * Step 1f
	 */
	// Right part of Eq. (22)
	sweet::Data::Sphere2D::DataGrid tmpg = 0.5*(ug*ug+vg*vg);

	sweet::Data::Sphere2D::DataSpectral e = phi_pert_phys+tmpg;

	/*
	 * Step 1g
	 */
	o_div_t -= ops->laplace(e);

	/*
	 * Compute Phi geopotential tendencies
	 */

	/*
	 * Step 2a
	 */
	u_nl = ug*(phi_pert_phys + gh0);
	v_nl = vg*(phi_pert_phys + gh0);

	ops->uv_2_vrtdiv(u_nl,v_nl, e, o_phi_pert_t);

	o_phi_pert_t *= -1.0;


}



void PDESWESphere2DTS_ln_erk::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi,
		sweet::Data::Sphere2D::DataSpectral &io_vort,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		SWEETErrorFatal("Only constant time step size allowed");

	// standard time stepping
	timestepping_rk.runTimestep(
			this,
			&PDESWESphere2DTS_ln_erk::euler_timestep_update_pert,	//!< pointer to function to compute euler time step updates
			io_phi, io_vort, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}


PDESWESphere2DTS_ln_erk::PDESWESphere2DTS_ln_erk()
{
}



PDESWESphere2DTS_ln_erk::~PDESWESphere2DTS_ln_erk()
{
}

