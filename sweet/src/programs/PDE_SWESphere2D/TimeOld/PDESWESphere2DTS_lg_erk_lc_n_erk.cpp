/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <sweet/Tools/StopwatchBox.hpp>
#include "PDESWESphere2DTS_lg_erk_lc_n_erk.hpp"




bool PDESWESphere2DTS_lg_erk_lc_n_erk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	int version_id = 0;
	if (timestepping_method == "lg_exp_lc_n_erk_ver1")
		version_id = 1;

	return setup(
			io_ops,
			timestepping_order,
			version_id
		);
}

bool PDESWESphere2DTS_lg_erk_lc_n_erk::setup(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_order,	//!< order of RK time stepping method
		int i_version_id
)
{
	ops = io_ops;

	timestepping_order = i_order;

	version_id = i_version_id;

	setupFG();
	return true;
}



void PDESWESphere2DTS_lg_erk_lc_n_erk::runTimestep(
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
				&PDESWESphere2DTS_lg_erk_lc_n_erk::euler_timestep_update_linear,	//!< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);

		timestepping_rk_nonlinear.runTimestep(
				this,
				&PDESWESphere2DTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	//!< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{

		if (version_id == 0)
		{
			// HALF time step for linear part
			timestepping_rk_linear.runTimestep(
					this,
					&PDESWESphere2DTS_lg_erk_lc_n_erk::euler_timestep_update_linear,	//!< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					this,
					&PDESWESphere2DTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	//!< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					timestepping_order,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_rk_linear.runTimestep(
					this,
					&PDESWESphere2DTS_lg_erk_lc_n_erk::euler_timestep_update_linear,	//!< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					this,
					&PDESWESphere2DTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	//!< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);


			// FULL time step for linear part
			timestepping_rk_linear.runTimestep(
					this,
					&PDESWESphere2DTS_lg_erk_lc_n_erk::euler_timestep_update_linear,	//!< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					timestepping_order,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					this,
					&PDESWESphere2DTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	//!< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);
		}
		else
		{
			SWEETErrorFatal("Invalid verison id");
		}
	}
	else
	{
		SWEETErrorFatal("Not yet supported!");
	}
}



void PDESWESphere2DTS_lg_erk_lc_n_erk::euler_timestep_update_linear(
		const sweet::Data::Sphere2D::DataSpectral &i_phi,
		const sweet::Data::Sphere2D::DataSpectral &i_vrt,
		const sweet::Data::Sphere2D::DataSpectral &i_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_t,
		sweet::Data::Sphere2D::DataSpectral &o_vrt_t,
		sweet::Data::Sphere2D::DataSpectral &o_div_t,

		double i_simulation_timestamp
)
{
	/*
	 * LINEAR
	 */
	double gh = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0;

	o_phi_t = -gh*i_div;
	o_div_t = -ops->laplace(i_phi);
	o_vrt_t.spectral_setZero();
}


void PDESWESphere2DTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n(
		const sweet::Data::Sphere2D::DataSpectral &i_phi,
		const sweet::Data::Sphere2D::DataSpectral &i_vort,
		const sweet::Data::Sphere2D::DataSpectral &i_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_t,
		sweet::Data::Sphere2D::DataSpectral &o_vort_t,
		sweet::Data::Sphere2D::DataSpectral &o_div_t,

		double i_simulation_timestamp
)
{
#if SWEET_BENCHMARK_TIMINGS
	sweet::Tools::StopwatchBox::getInstance().main_timestepping_nonlinearities.start();
#endif

	/*
	 * NON-LINEAR
	 *
	 * Follows Hack & Jakob formulation
	 */

	sweet::Data::Sphere2D::DataGrid ug(i_phi.sphere2DDataConfig);
	sweet::Data::Sphere2D::DataGrid vg(i_phi.sphere2DDataConfig);

	sweet::Data::Sphere2D::DataGrid vrtg = i_vort.toGrid();
	sweet::Data::Sphere2D::DataGrid divg = i_div.toGrid();
	ops->vrtdiv_2_uv(i_vort, i_div, ug, vg);

	sweet::Data::Sphere2D::DataGrid phig = i_phi.toGrid();

	using namespace sweet;
	sweet::Data::Sphere2D::DataGrid tmpg1 = ug*(vrtg+fg);
	sweet::Data::Sphere2D::DataGrid tmpg2 = vg*(vrtg+fg);

	ops->uv_2_vrtdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

	o_vort_t *= -1.0;

	tmpg1 = ug*phig;
	tmpg2 = vg*phig;

	sweet::Data::Sphere2D::DataSpectral tmpspec(i_phi.sphere2DDataConfig);
	ops->uv_2_vrtdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

	o_phi_t *= -1.0;

	sweet::Data::Sphere2D::DataGrid tmpg(i_phi.sphere2DDataConfig);
	tmpg = 0.5*(ug*ug+vg*vg);

	tmpspec = tmpg;
	o_div_t += -ops->laplace(tmpspec);


#if SWEET_BENCHMARK_TIMINGS
	sweet::Tools::StopwatchBox::getInstance().main_timestepping_nonlinearities.stop();
#endif
}



/**
 * This routine is used by other time step implementations
 */
void PDESWESphere2DTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n(
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

	euler_timestep_update_lc_n(
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





PDESWESphere2DTS_lg_erk_lc_n_erk::PDESWESphere2DTS_lg_erk_lc_n_erk()
{
}



PDESWESphere2DTS_lg_erk_lc_n_erk::~PDESWESphere2DTS_lg_erk_lc_n_erk()
{
}

