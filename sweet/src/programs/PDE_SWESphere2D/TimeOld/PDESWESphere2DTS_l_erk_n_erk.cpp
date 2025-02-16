/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <sweet/Tools/StopwatchBox.hpp>
#include "PDESWESphere2DTS_l_erk_n_erk.hpp"


bool PDESWESphere2DTS_l_erk_n_erk::setup_auto(
	const std::string &i_timestepping_method,
	sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	setup_main(	io_ops,
			shackPDESWETimeDisc->timestepping_order,
			shackPDESWETimeDisc->timestepping_order2
		);

	return true;
}


bool PDESWESphere2DTS_l_erk_n_erk::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_order,	//!< order of RK time stepping method for non-linear parts
		int i_order2	//!< order of RK time stepping method for non-linear parts
)
{
	ops = io_ops;

	timestepping_order = i_order;
	timestepping_order2 = i_order2;

	setupFG();

	return true;
}



bool PDESWESphere2DTS_l_erk_n_erk::implementsTimesteppingMethod(
		const std::string &i_timestepping_method
)
{
	timestepping_method = i_timestepping_method;
	return i_timestepping_method == "l_erk_n_erk";
}

std::string PDESWESphere2DTS_l_erk_n_erk::getIDString()
{
	return "l_erk_n_erk";
}


void PDESWESphere2DTS_l_erk_n_erk::runTimestep(
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
				&PDESWESphere2DTS_l_erk_n_erk::euler_timestep_update_linear,	//!< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);

		timestepping_rk_nonlinear.runTimestep(
				this,
				&PDESWESphere2DTS_l_erk_n_erk::euler_timestep_update_nonlinear,	//!< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order2,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{
		// HALF time step for linear part
		timestepping_rk_linear.runTimestep(
				this,
				&PDESWESphere2DTS_l_erk_n_erk::euler_timestep_update_linear,	//!< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				timestepping_order,		//! This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// FULL time step for non-linear part
		timestepping_rk_nonlinear.runTimestep(
				this,
				&PDESWESphere2DTS_l_erk_n_erk::euler_timestep_update_nonlinear,	//!< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order2,		//! This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// HALF time step for linear part
		timestepping_rk_linear.runTimestep(
				this,
				&PDESWESphere2DTS_l_erk_n_erk::euler_timestep_update_linear,	//!< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				timestepping_order,		//! This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 0)
	{
		SWEETErrorFatal("Please specify the timestepping order via --timestepping-order=[int]");
	}
	else
	{
		SWEETErrorFatal("programs/swe_sphere2D/PDESWESphere2DTS_l_erk_n_erk.cpp: This order is not yet supported!");
	}
}



void PDESWESphere2DTS_l_erk_n_erk::euler_timestep_update_linear(
	const sweet::Data::Sphere2D::DataSpectral &i_phi,
	const sweet::Data::Sphere2D::DataSpectral &i_vrt,
	const sweet::Data::Sphere2D::DataSpectral &i_div,

	sweet::Data::Sphere2D::DataSpectral &o_phi_t,	//!< time updates
	sweet::Data::Sphere2D::DataSpectral &o_vrt_t,	//!< time updates
	sweet::Data::Sphere2D::DataSpectral &o_div_t,	//!< time updates

	double i_simulation_timestamp
)
{
	double gh0 = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0;

	if (!shackPDESWESphere2D->sphere2d_use_fsphere2D)
	{

		/*
		 * See documentation in [sweet]/doc/swe/swe_sphere2d_formulation/
		 */
		/*
		 * Step 1a
		 */
		sweet::Data::Sphere2D::DataGrid ug(i_phi.sphere2DDataConfig);
		sweet::Data::Sphere2D::DataGrid vg(i_phi.sphere2DDataConfig);
		ops->vrtdiv_2_uv(i_vrt, i_div, ug, vg);

		/*
		 * Step 1b
		 */
		using namespace sweet;
		sweet::Data::Sphere2D::DataGrid tmpg1 = ug*fg;
		sweet::Data::Sphere2D::DataGrid tmpg2 = vg*fg;

		/*
		 * Step 1c
		 */
		ops->uv_2_vrtdiv(tmpg1, tmpg2, o_div_t, o_vrt_t);

		/*
		 * Step 1d
		 */
		o_vrt_t *= -1.0;

		/*
		 * Step 1e
		 */
		o_div_t += -ops->laplace(i_phi);

		/*
		 * DIV on velocity field
		 */
		o_phi_t = (-gh0)*i_div;
	}
	else
	{
		o_div_t = -ops->laplace(i_phi);

		o_vrt_t = -shackPDESWESphere2D->sphere2d_fsphere2d_f0*i_div;
		o_div_t += shackPDESWESphere2D->sphere2d_fsphere2d_f0*i_vrt;

		o_phi_t = (-gh0)*i_div;
	}
}



void PDESWESphere2DTS_l_erk_n_erk::euler_timestep_update_nonlinear(
		const sweet::Data::Sphere2D::DataSpectral &i_phi,
		const sweet::Data::Sphere2D::DataSpectral &i_vrt,
		const sweet::Data::Sphere2D::DataSpectral &i_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_dt,	//!< time updates
		sweet::Data::Sphere2D::DataSpectral &o_vrt_dt,	//!< time updates
		sweet::Data::Sphere2D::DataSpectral &o_div_dt,	//!< time updates

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

	sweet::Data::Sphere2D::DataGrid vrtg = i_vrt.toGrid();
	sweet::Data::Sphere2D::DataGrid divg = i_div.toGrid();
	ops->vrtdiv_2_uv(i_vrt, i_div, ug, vg);

	sweet::Data::Sphere2D::DataGrid phig = i_phi.toGrid();

	sweet::Data::Sphere2D::DataGrid tmpg1 = ug*(vrtg/*+fg*/);
	sweet::Data::Sphere2D::DataGrid tmpg2 = vg*(vrtg/*+fg*/);

	ops->uv_2_vrtdiv(tmpg1, tmpg2, o_div_dt, o_vrt_dt);

	o_vrt_dt *= -1.0;

	tmpg1 = ug*phig;
	tmpg2 = vg*phig;

	sweet::Data::Sphere2D::DataSpectral tmpspec(i_phi.sphere2DDataConfig);
	ops->uv_2_vrtdiv(tmpg1,tmpg2, tmpspec, o_phi_dt);

	o_phi_dt *= -1.0;

	sweet::Data::Sphere2D::DataGrid tmpg = 0.5*(ug*ug+vg*vg);

	tmpspec = /*phig+*/tmpg;

	o_div_dt += -ops->laplace(tmpspec);


#if SWEET_BENCHMARK_TIMINGS
	sweet::Tools::StopwatchBox::getInstance().main_timestepping_nonlinearities.stop();
#endif
}



/**
 * This routine is used by other time step implementations
 */
void PDESWESphere2DTS_l_erk_n_erk::euler_timestep_update_nonlinear(
		sweet::Data::Sphere2D::DataSpectral &io_phi,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_dt,
		double i_simulation_timestamp
)
{
#if SWEET_BENCHMARK_TIMINGS
	sweet::Tools::StopwatchBox::getInstance().main_timestepping_nonlinearities.start();
#endif

	sweet::Data::Sphere2D::DataSpectral tmp_phi(io_phi.sphere2DDataConfig);
	sweet::Data::Sphere2D::DataSpectral tmp_vrt(io_vrt.sphere2DDataConfig);
	sweet::Data::Sphere2D::DataSpectral tmp_div(io_div.sphere2DDataConfig);

	euler_timestep_update_nonlinear(
			io_phi,
			io_vrt,
			io_div,

			tmp_phi,
			tmp_vrt,
			tmp_div,

			i_simulation_timestamp
		);

	io_phi += i_dt*tmp_phi;
	io_vrt += i_dt*tmp_vrt;
	io_div += i_dt*tmp_div;

#if SWEET_BENCHMARK_TIMINGS
	sweet::Tools::StopwatchBox::getInstance().main_timestepping_nonlinearities.stop();
#endif
}



PDESWESphere2DTS_l_erk_n_erk::PDESWESphere2DTS_l_erk_n_erk()
{
}



PDESWESphere2DTS_l_erk_n_erk::~PDESWESphere2DTS_l_erk_n_erk()
{
}

