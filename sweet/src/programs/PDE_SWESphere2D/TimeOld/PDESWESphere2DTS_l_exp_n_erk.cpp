/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere2DTS_l_exp_n_erk.hpp"



bool PDESWESphere2DTS_l_exp_n_erk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	int _version_id = 0;

	if (timestepping_method == "l_exp_n_erk_ver1")
		_version_id = 1;

	return setup_main(
			io_ops,
			shackExpIntegration,
			shackExpIntegration->exp_method,
			shackPDESWETimeDisc->timestepping_order,
			shackPDESWETimeDisc->timestepping_order2,
			shackTimestepControl->currentTimestepSize,
			shackPDESWESphere2D->sphere2d_use_fsphere2D,
			_version_id,
			shackExpIntegration->sphere2d_solver_preallocation
		);
}


bool PDESWESphere2DTS_l_exp_n_erk::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		sweet::ExpIntegration::Shack *i_shackExpIntegration,
		const std::string &i_exp_method,
		int i_order,	//!< order of RK time stepping method
		int i_order2,	//!< order of RK time stepping method of non-linear parts
		double i_timestepSize,
		bool i_use_f_sphere2D,
		int i_version_id,
		bool i_use_rexi_sphere2d_solver_preallocation
)
{
	ops = io_ops;

	timestepping_order = i_order;
	timestepping_order2 = i_order2;

	timestep_size = i_timestepSize;

	version_id = i_version_id;

	if (timestepping_order != timestepping_order2)
		SWEETErrorFatal("Mismatch of timestepping method orders, they should be equal");


	if (timestepping_order2 == 1)
	{
		timestepping_l_rexi.setup_variant_100(
				io_ops,
				i_shackExpIntegration,
				"phi0",
				i_exp_method,
				i_timestepSize,
				i_use_f_sphere2D,
				false,
				timestepping_order,
				i_use_rexi_sphere2d_solver_preallocation
		);
	}
	else if (timestepping_order2 == 2 || timestepping_order2 == 4)
	{
		if (version_id == 0)
		{
			timestepping_l_rexi.setup_variant_100(
					io_ops,
					i_shackExpIntegration,
					"phi0",
					i_exp_method,
					i_timestepSize*0.5,
					i_use_f_sphere2D,
					false,
					timestepping_order,
					i_use_rexi_sphere2d_solver_preallocation
			);
		}
		else if (version_id == 1)
		{
			timestepping_l_rexi.setup_variant_100(
					io_ops,
					i_shackExpIntegration,
					"phi0",
					i_exp_method,
					i_timestepSize,
					i_use_f_sphere2D,
					false,
					timestepping_order,
					i_use_rexi_sphere2d_solver_preallocation
			);
		}
		else
		{
			SWEETErrorFatal("Invalid version");
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
	timestepping_l_erk_n_erk.setup_main(io_ops, 1, 1);
	return true;
}




void PDESWESphere2DTS_l_exp_n_erk::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order2 == 1)
	{
		if (version_id == 0)
		{
			// first order REXI for linear part
			timestepping_l_rexi.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// first order explicit for non-linear
			timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);
		}
		else if (version_id == 1)
		{
			// first order explicit for non-linear
			timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// first order REXI for linear part
			timestepping_l_rexi.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);
		}
		else
		{
			SWEETErrorFatal("Invalid version id");
		}
	}
	else if (timestepping_order2 == 2 || timestepping_order2)
	{
		if (version_id == 0)
		{
			// HALF time step for linear part
			timestepping_l_rexi.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_l_erk_n_erk,
					&PDESWESphere2DTS_l_erk_n_erk::euler_timestep_update_nonlinear,	//!< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					timestepping_order2,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_l_rexi.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp+i_fixed_dt*0.5	/* TODO: CHECK THIS, THIS MIGHT BE WRONG!!! */
				);

		}
		else if (version_id == 1)
		{

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_l_erk_n_erk,
					&PDESWESphere2DTS_l_erk_n_erk::euler_timestep_update_nonlinear,	//!< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order2,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_l_rexi.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_l_erk_n_erk,
					&PDESWESphere2DTS_l_erk_n_erk::euler_timestep_update_nonlinear,	//!< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order2,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);
		}
		else
		{
			SWEETErrorFatal("Invalid version id");
		}
	}
	else
	{
		SWEETErrorFatal("Not yet supported!");
	}
}



PDESWESphere2DTS_l_exp_n_erk::PDESWESphere2DTS_l_exp_n_erk()	:
		version_id(0)
{
}



PDESWESphere2DTS_l_exp_n_erk::~PDESWESphere2DTS_l_exp_n_erk()
{
}

