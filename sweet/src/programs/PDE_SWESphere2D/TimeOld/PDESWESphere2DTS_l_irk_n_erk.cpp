/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere2DTS_l_irk_n_erk.hpp"




bool PDESWESphere2DTS_l_irk_n_erk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	if (
		timestepping_method == "l_irk_n_erk" ||
		timestepping_method == "l_irk_n_erk_ver0" ||
		timestepping_method == "l_cn_n_erk" ||
		timestepping_method == "l_cn_n_erk_ver0"
	)
		return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, shackPDESWETimeDisc->timestepping_order2, 0);
	else
		return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, shackPDESWETimeDisc->timestepping_order2, 1);
}



/*
 * Setup
 */
bool PDESWESphere2DTS_l_irk_n_erk::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_order,	//!< order of RK time stepping method
		int i_order2,	//!< order of RK time stepping method for non-linear parts
		int i_version_id
)
{
	ops = io_ops;

	if (i_order2 < 0)
		i_order2 = i_order;

	if (i_order != i_order2)
		SWEETErrorFatal("Orders of 1st and 2nd one must match");

	version_id = i_version_id;

	timestepping_order = i_order;
	timestepping_order2 = i_order2;
	timestep_size = shackTimestepControl->currentTimestepSize;

	if (timestepping_order == 1)
	{
		timestepping_l_irk.setup(
			ops,
			1,
			timestep_size
		);
	}
	else if (timestepping_order == 2)
	{
		if (version_id == 0)
		{
			timestepping_l_irk.setup_main(
					ops,
					2,
					timestep_size*0.5,
					shackPDESWETimeDisc->timestepping_crank_nicolson_filter,
					false
			);
		}
		else if (version_id == 1)
		{
			timestepping_l_irk.setup_main(
					ops,
					2,
					timestep_size,
					shackPDESWETimeDisc->timestepping_crank_nicolson_filter,
					false
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
	timestepping_l_erk_n_erk.setup_main(ops, 1, 1);

	return true;
}


bool PDESWESphere2DTS_l_irk_n_erk::implementsTimesteppingMethod(
		const std::string &i_timestepping_method
)
{
	timestepping_method = i_timestepping_method;

	if (
		i_timestepping_method == "l_irk_n_erk"	|| i_timestepping_method == "l_irk_n_erk_ver0"	||
		i_timestepping_method == "l_cn_n_erk"	|| i_timestepping_method == "l_cn_n_erk_ver0"		||
		i_timestepping_method == "l_irk_n_erk_ver1"
	)
		return true;

	return false;
}


std::string PDESWESphere2DTS_l_irk_n_erk::getIDString()
{
	std::string s = "l_irk_n_erk_ver";

	if (version_id == 0)
		s += "0";
	else if (version_id == 1)
		s += "1";
	else
		SWEETErrorFatal("Version ID");

	return s;
}


void PDESWESphere2DTS_l_irk_n_erk::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi,
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
			timestepping_l_irk.runTimestep(
					io_phi, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// first order explicit for non-linear
			timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
					io_phi, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);
		}
		else
		{
			// first order explicit for non-linear
			timestepping_l_erk_n_erk.euler_timestep_update_nonlinear(
					io_phi, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// first order IRK for linear
			timestepping_l_irk.runTimestep(
					io_phi, io_vrt, io_div,
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
			timestepping_l_irk.runTimestep(
					io_phi, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_l_erk_n_erk,
					&PDESWESphere2DTS_l_erk_n_erk::euler_timestep_update_nonlinear,	//!< pointer to function to compute euler time step updates
					io_phi, io_vrt, io_div,
					i_fixed_dt,
					timestepping_order2,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_l_irk.runTimestep(
					io_phi, io_vrt, io_div,
					i_fixed_dt*0.5,
					i_simulation_timestamp+i_fixed_dt*0.5	/* TODO: CHECK THIS, THIS MIGHT BE WRONG!!! However, we only have autonomous simulations so far */
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_l_erk_n_erk,
					&PDESWESphere2DTS_l_erk_n_erk::euler_timestep_update_nonlinear,	//!< pointer to function to compute euler time step updates
					io_phi, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order2,		//! This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_l_irk.runTimestep(
					io_phi, io_vrt, io_div,
					i_fixed_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.runTimestep(
					&timestepping_l_erk_n_erk,
					&PDESWESphere2DTS_l_erk_n_erk::euler_timestep_update_nonlinear,	//!< pointer to function to compute euler time step updates
					io_phi, io_vrt, io_div,
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





PDESWESphere2DTS_l_irk_n_erk::PDESWESphere2DTS_l_irk_n_erk()	:
		version_id(0),
		timestepping_order(-1)
{

}



PDESWESphere2DTS_l_irk_n_erk::~PDESWESphere2DTS_l_irk_n_erk()
{
}

