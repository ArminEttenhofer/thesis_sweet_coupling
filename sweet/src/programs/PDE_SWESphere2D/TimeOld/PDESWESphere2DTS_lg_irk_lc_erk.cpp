/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere2DTS_lg_irk_lc_erk.hpp"



bool PDESWESphere2DTS_lg_irk_lc_erk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	if (timestepping_method == "lg_irk_lc_erk_ver1")
		return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, 1);
	else
		return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, 0);
}



bool PDESWESphere2DTS_lg_irk_lc_erk::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_timestepping_order,	//!< order of RK time stepping method
		int i_version_id
)
{
	ops = io_ops;

	version_id = i_version_id;

	timestepping_order = i_timestepping_order;
	timestep_size = shackTimestepControl->currentTimestepSize;

	if (timestepping_order == 1)
	{
		timestepping_lg_irk.setup(
				io_ops,
				1,
				timestep_size*0.5
		);
	}
	else if (timestepping_order == 2)
	{
		if (version_id == 0)
		{
			timestepping_lg_irk.setup(
					io_ops,
					2,
					timestep_size*0.5
			);
		}
		else if (version_id == 1)
		{
			timestepping_lg_irk.setup(
					ops,
					2,
					timestep_size
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
	timestepping_lg_erk_lc_erk.setup_main(ops, 1);

	return true;
}




void PDESWESphere2DTS_lg_irk_lc_erk::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		if (version_id == 0)
		{
			// first order IRK for linear
			timestepping_lg_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);


			// first order explicit for non-linear
			timestepping_lg_erk_lc_erk.euler_timestep_update_lc(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);
		}
		else if (version_id == 0)
		{
			// first order explicit for non-linear
			timestepping_lg_erk_lc_erk.euler_timestep_update_lc(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// first order IRK for linear
			timestepping_lg_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);
		}
	}
	else if (timestepping_order == 2)
	{
		if (version_id == 0)
		{
			// HALF time step for linear part
			timestepping_lg_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_lg_erk_lc_erk,
					&PDESWESphere2DTS_lg_erk_lc_erk::euler_timestep_update_lc,	//!< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					timestepping_order,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_lg_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp+i_fixed_dt*0.5	/* TODO: CHECK THIS, THIS MIGHT BE WRONG!!! */
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_lg_erk_lc_erk,
					&PDESWESphere2DTS_lg_erk_lc_erk::euler_timestep_update_lc,	//!< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_lg_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_lg_erk_lc_erk,
					&PDESWESphere2DTS_lg_erk_lc_erk::euler_timestep_update_lc,	//!< pointer to function to compute euler time step updates
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




PDESWESphere2DTS_lg_irk_lc_erk::PDESWESphere2DTS_lg_irk_lc_erk()
{
}



PDESWESphere2DTS_lg_irk_lc_erk::~PDESWESphere2DTS_lg_irk_lc_erk()
{
}

