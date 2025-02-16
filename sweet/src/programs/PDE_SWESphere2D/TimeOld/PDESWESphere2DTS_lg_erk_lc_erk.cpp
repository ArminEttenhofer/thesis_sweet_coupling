/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere2DTS_lg_erk_lc_erk.hpp"



bool PDESWESphere2DTS_lg_erk_lc_erk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order);
}


bool PDESWESphere2DTS_lg_erk_lc_erk::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_order	//!< order of RK time stepping method
)
{
	ops = io_ops;

	timestepping_order = i_order;
	setupFG();

	return true;
}


void PDESWESphere2DTS_lg_erk_lc_erk::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		timestepping_rk_linear.runTimestep(
				this,
				&PDESWESphere2DTS_lg_erk_lc_erk::euler_timestep_update_lg,	//!< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);

		timestepping_rk_nonlinear.runTimestep(
				this,
				&PDESWESphere2DTS_lg_erk_lc_erk::euler_timestep_update_lc,	//!< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{
		// HALF time step for linear part
		timestepping_rk_linear.runTimestep(
				this,
				&PDESWESphere2DTS_lg_erk_lc_erk::euler_timestep_update_lg,	//!< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				timestepping_order,		//! This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// FULL time step for non-linear part
		timestepping_rk_nonlinear.runTimestep(
				this,
				&PDESWESphere2DTS_lg_erk_lc_erk::euler_timestep_update_lc,	//!< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,		//! This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// HALF time step for linear part
		timestepping_rk_linear.runTimestep(
				this,
				&PDESWESphere2DTS_lg_erk_lc_erk::euler_timestep_update_lg,	//!< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				timestepping_order,		//! This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);
	}
	else
	{
		SWEETErrorFatal("This time stepping order is not yet supported!");
	}
}





/*
 * Main routine for method to be used in case of finite differences
 */
void PDESWESphere2DTS_lg_erk_lc_erk::euler_timestep_update(
		const sweet::Data::Sphere2D::DataSpectral &i_phi,
		const sweet::Data::Sphere2D::DataSpectral &i_vort,
		const sweet::Data::Sphere2D::DataSpectral &i_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_t,	//!< time updates
		sweet::Data::Sphere2D::DataSpectral &o_vort_t,	//!< time updates
		sweet::Data::Sphere2D::DataSpectral &o_div_t,	//!< time updates

		double i_simulation_timestamp
)
{
	if (shackPDESWESphere2D->sphere2d_use_fsphere2D)
	{
		double gh0 = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0;

		o_phi_t = -gh0*i_div;
		o_div_t = -ops->laplace(i_phi);

		o_vort_t = -shackPDESWESphere2D->sphere2d_fsphere2d_f0*i_div;
		o_div_t += shackPDESWESphere2D->sphere2d_fsphere2d_f0*i_vort;
	}
	else
	{
#if 0
		double gh = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0;

		o_phi_t = -gh*i_div;
		o_div_t = -ops->laplace(i_phi);

		/*
		 * This doesn't converge to the reference solution
		 */
		sweet::Data::Sphere2D::DataSpectral f(fg);
		o_vort_t = -f*i_div;
		o_div_t += f*i_vort;

#else

		double gh = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0;

		/*
		 * Apply Coriolis Effect in physical VELOCITY space
		 */
		sweet::Data::Sphere2D::DataGrid ug(i_phi.sphere2DDataConfig);
		sweet::Data::Sphere2D::DataGrid vg(i_phi.sphere2DDataConfig);
		ops->vrtdiv_2_uv(i_vort, i_div, ug, vg);

		sweet::Data::Sphere2D::DataGrid tmpg1 = ug*fg;
		sweet::Data::Sphere2D::DataGrid tmpg2 = vg*fg;

		ops->uv_2_vrtdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

		o_vort_t *= -1.0;
		o_div_t += -ops->laplace(i_phi);

		/*
		 * DIV on velocity field
		 */
		o_phi_t = (-gh)*i_div;
#endif
	}
}




void PDESWESphere2DTS_lg_erk_lc_erk::euler_timestep_update_lg(
		const sweet::Data::Sphere2D::DataSpectral &i_phi,
		const sweet::Data::Sphere2D::DataSpectral &i_vort,
		const sweet::Data::Sphere2D::DataSpectral &i_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_t,	//!< time updates
		sweet::Data::Sphere2D::DataSpectral &o_vort_t,	//!< time updates
		sweet::Data::Sphere2D::DataSpectral &o_div_t,	//!< time updates

		double i_simulation_timestamp
)
{
	double gh0 = shackPDESWESphere2D->gravitation*shackPDESWESphere2D->h0;


	/*
	 * See documentation in [sweet]/doc/swe/swe_sphere2d_formulation/
	 * Section "lg_erk"
	 */


	/*
	 * step 1a
	 */
	o_vort_t.spectral_setZero();

	/*
	 * step 1b
	 */
	o_div_t = -ops->laplace(i_phi);

	/*
	 * step 2a
	 */
	o_phi_t = -gh0*i_div;
}



void PDESWESphere2DTS_lg_erk_lc_erk::euler_timestep_update_lc(
		const sweet::Data::Sphere2D::DataSpectral &i_phi,
		const sweet::Data::Sphere2D::DataSpectral &i_vort,
		const sweet::Data::Sphere2D::DataSpectral &i_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_t,	//!< time updates
		sweet::Data::Sphere2D::DataSpectral &o_vort_t,	//!< time updates
		sweet::Data::Sphere2D::DataSpectral &o_div_t,	//!< time updates

		double i_simulation_timestamp
)
{
	/*
	 * See documentation in [sweet]/doc/swe/swe_sphere2d_formulation/
	 * Section "lc_erk"
	 */

	/*
	 * step 1a
	 */
	sweet::Data::Sphere2D::DataGrid ug(i_phi.sphere2DDataConfig);
	sweet::Data::Sphere2D::DataGrid vg(i_phi.sphere2DDataConfig);
	ops->vrtdiv_2_uv(i_vort, i_div, ug, vg);

	/*
	 * step 1b
	 */
	sweet::Data::Sphere2D::DataGrid tmp_u = ug*fg;
	sweet::Data::Sphere2D::DataGrid tmp_v = vg*fg;

	/*
	 * step 1c
	 */
	ops->uv_2_vrtdiv(tmp_u, tmp_v, o_div_t, o_vort_t);

	/*
	 * step 1d
	 */
	o_vort_t *= -1.0;


	/*
	 * step 1e
	 * Nothing to do
	 */

	/*
	 * step 2a
	 * Zero tendencies
	 */
	o_phi_t.spectral_setZero();

}



/**
 * This routine is used by other time step implementations
 */
void PDESWESphere2DTS_lg_erk_lc_erk::euler_timestep_update_lc(
		sweet::Data::Sphere2D::DataSpectral &io_phi,
		sweet::Data::Sphere2D::DataSpectral &io_vort,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_dt,
		double i_simulation_timestamp
)
{
	sweet::Data::Sphere2D::DataSpectral tmp_phi(io_phi.sphere2DDataConfig);
	sweet::Data::Sphere2D::DataSpectral tmp_vort(io_vort.sphere2DDataConfig);
	sweet::Data::Sphere2D::DataSpectral tmp_div(io_div.sphere2DDataConfig);

	euler_timestep_update_lc(
			io_phi,
			io_vort,
			io_div,

			tmp_phi,
			tmp_vort,
			tmp_div,

			i_simulation_timestamp
		);

	io_phi += i_dt*tmp_phi;
	io_vort += i_dt*tmp_vort;
	io_div += i_dt*tmp_div;
}



PDESWESphere2DTS_lg_erk_lc_erk::PDESWESphere2DTS_lg_erk_lc_erk()
{
}



PDESWESphere2DTS_lg_erk_lc_erk::~PDESWESphere2DTS_lg_erk_lc_erk()
{
}

