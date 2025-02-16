/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere2DTS_l_erk.hpp"

bool PDESWESphere2DTS_l_erk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order);
}


bool PDESWESphere2DTS_l_erk::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_order	//!< order of RK time stepping method
)
{
	ops = io_ops;

	setupFG();

	timestepping_order = i_order;
	return true;
}



void PDESWESphere2DTS_l_erk::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	// standard time stepping
	timestepping_rk.runTimestep(
			this,
			&PDESWESphere2DTS_l_erk::euler_timestep_update,	//!< pointer to function to compute euler time step updates
			io_phi_pert, io_vrt, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}


/*
 * Main routine for method to be used in case of finite differences
 */
void PDESWESphere2DTS_l_erk::euler_timestep_update(
		const sweet::Data::Sphere2D::DataSpectral &i_phi_pert,
		const sweet::Data::Sphere2D::DataSpectral &i_vort,
		const sweet::Data::Sphere2D::DataSpectral &i_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_pert_t,	//!< time updates
		sweet::Data::Sphere2D::DataSpectral &o_vort_t,	//!< time updates
		sweet::Data::Sphere2D::DataSpectral &o_div_t,	//!< time updates

		double i_simulation_timestamp
)
{
	if (!shackPDESWESphere2D->sphere2d_use_fsphere2D)
	{

		double gh0 = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0;

		/*
		 * See documentation in [sweet]/doc/swe/swe_sphere2d_formulation/
		 */
		/*
		 * Step 1a
		 */
		sweet::Data::Sphere2D::DataGrid ug(i_phi_pert.sphere2DDataConfig);
		sweet::Data::Sphere2D::DataGrid vg(i_phi_pert.sphere2DDataConfig);
		ops->vrtdiv_2_uv(i_vort, i_div, ug, vg);

		/*
		 * Step 1b
		 */
		sweet::Data::Sphere2D::DataGrid tmpg1 = ug*fg;
		sweet::Data::Sphere2D::DataGrid tmpg2 = vg*fg;

		/*
		 * Step 1c
		 */
		ops->uv_2_vrtdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

		/*
		 * Step 1d
		 */
		o_vort_t *= -1.0;

		/*
		 * Step 1e
		 */
		o_div_t += -ops->laplace(i_phi_pert);

		/*
		 * DIV on velocity field
		 */
		o_phi_pert_t = (-gh0)*i_div;
	}
	else
	{
		double gh = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0;

		o_div_t = -ops->laplace(i_phi_pert);

		o_vort_t = -shackPDESWESphere2D->sphere2d_fsphere2d_f0*i_div;
		o_div_t += shackPDESWESphere2D->sphere2d_fsphere2d_f0*i_vort;

		o_phi_pert_t = -gh*i_div;
	}
}




PDESWESphere2DTS_l_erk::PDESWESphere2DTS_l_erk()
{
}



PDESWESphere2DTS_l_erk::~PDESWESphere2DTS_l_erk()
{
}

