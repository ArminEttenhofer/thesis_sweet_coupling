/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere2DTS_lg_exp_lc_n_erk.hpp"



void PDESWESphere2DTS_lg_exp_lc_n_erk::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi,
		sweet::Data::Sphere2D::DataSpectral &io_vort,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		sweet::Data::Sphere2D::DataSpectral tmp_phi = io_phi;
		sweet::Data::Sphere2D::DataSpectral tmp_vort = io_vort;
		sweet::Data::Sphere2D::DataSpectral tmp_div = io_div;

		// first order IRK for linear
		timestepping_lg_rexi.runTimestep(
				io_phi, io_vort, io_div,
				i_dt,
				i_simulation_timestamp
			);

		sweet::Data::Sphere2D::DataSpectral phi_dt(io_phi.sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral vort_dt(io_vort.sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral div_dt(io_div.sphere2DDataConfig);

		// first order explicit for non-linear
		timestepping_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				tmp_phi, tmp_vort, tmp_div,
				phi_dt, vort_dt, div_dt,
				i_simulation_timestamp
			);


		io_phi += i_dt*phi_dt;
		io_vort += i_dt*vort_dt;
		io_div += i_dt*div_dt;
	}
	else if (timestepping_order == 2 || timestepping_order == 4)
	{
		if (version_id == 0)
		{
			// HALF time step for linear part
			timestepping_lg_rexi.runTimestep(
					io_phi, io_vort, io_div,
					i_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_lg_erk_lc_n_erk,
					&PDESWESphere2DTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	//!< pointer to function to compute Euler time step updates
					io_phi, io_vort, io_div,
					i_dt,
					timestepping_order2,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_lg_rexi.runTimestep(
					io_phi, io_vort, io_div,
					i_dt*0.5,
					i_simulation_timestamp+i_dt*0.5	/* TODO: CHECK THIS, THIS MIGHT BE WRONG!!! */
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_lg_erk_lc_n_erk,
					&PDESWESphere2DTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	//!< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt*0.5,
					timestepping_order2,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_lg_rexi.runTimestep(
					io_phi, io_vort, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_lg_erk_lc_n_erk,
					&PDESWESphere2DTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	//!< pointer to function to compute euler time step updates
					io_phi, io_vort, io_div,
					i_dt*0.5,
					timestepping_order2,		//! This must be 2nd order accurate to get overall 2nd order accurate method
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


bool PDESWESphere2DTS_lg_exp_lc_n_erk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	int version = 0;
	if (timestepping_method == "lg_exp_lc_n_erk_ver1")
		version = 1;

	return setup_main(
			io_ops,
			shackExpIntegration,
			shackPDESWETimeDisc->timestepping_order,
			shackPDESWETimeDisc->timestepping_order2,
			shackTimestepControl->currentTimestepSize,
			version
		);
}


/*
 * Setup
 */
bool PDESWESphere2DTS_lg_exp_lc_n_erk::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		sweet::ExpIntegration::Shack *i_shackExpIntegration,
		int i_timestepping_order,
		int i_timestepping_order2,
		double i_timestepSize,
		int i_version_id
)
{
	ops = io_ops;
	version_id = i_version_id;

	timestepping_order = i_timestepping_order;
	timestepping_order2 = i_timestepping_order2;

	timestep_size = shackTimestepControl->currentTimestepSize;

	if (timestepping_order2 == 1)
	{
		timestepping_lg_rexi.setup_variant_10(
				ops,
				i_shackExpIntegration,
				"phi0",
				i_timestepSize,
				false,
				true
			);
	}
	else if (timestepping_order2 == 2 || timestepping_order2 == 4)
	{
		if (version_id == 0)
		{
			timestepping_lg_rexi.setup_variant_10(
					ops,
					i_shackExpIntegration,
					"phi0",
					i_timestepSize*0.5,
					false,
					true
			);
		}
		else if (version_id == 1)
		{
			timestepping_lg_rexi.setup_variant_10(
					ops,
					i_shackExpIntegration,
					"phi0",
					i_timestepSize,
					false,
					true
			);
		}
		else
		{
			SWEETErrorFatal("Invalid version id");
		}
	}
	else
	{
		SWEETErrorFatal("Invalid timestepping order");
	}


	//
	// Only request 1st order time stepping methods for irk and erk
	// These 1st order methods will be combined to higher-order methods in this class
	//
	timestepping_lg_erk_lc_n_erk.setup(ops, 1, -1);

	return true;
}



PDESWESphere2DTS_lg_exp_lc_n_erk::PDESWESphere2DTS_lg_exp_lc_n_erk()
{
}



PDESWESphere2DTS_lg_exp_lc_n_erk::~PDESWESphere2DTS_lg_exp_lc_n_erk()
{
}

