/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere2DTS_lg_irk_lc_na_erk_vd.hpp"



bool PDESWESphere2DTS_lg_irk_lc_na_erk_vd::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	if (
		i_timestepping_method == "lg_irk_lc_na_erk_vd" || i_timestepping_method == "lg_irk_lc_na_erk_vd_ver0" ||
		i_timestepping_method == "lg_irk_lc_na_erk_vd_ver1"
	)
		return true;

	return false;
}



bool PDESWESphere2DTS_lg_irk_lc_na_erk_vd::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	if (
		timestepping_method == "lg_irk_lc_na_erk_vd" ||
		timestepping_method == "lg_irk_lc_na_erk_vd_ver0"
	)
	{
		return setup(io_ops, timestepping_order, timestepping_order2, 0);
	}
	else if (
			timestepping_method == "lg_irk_lc_na_erk_vd_ver1"
		)
	{
		return setup(io_ops, timestepping_order,timestepping_order2, 1);
	}
	else
	{
		SWEETErrorFatal("Not implemented");
	}

	return false;
}


bool PDESWESphere2DTS_lg_irk_lc_na_erk_vd::setup(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_timestepping_order,	//!< order of RK time stepping method
		int i_timestepping_order2,	//!< order of RK time stepping method for non-linear parts
		int i_version_id
)
{
	ops = io_ops;

	version_id = i_version_id;

	timestepping_order = i_timestepping_order;
	timestepping_order2 = i_timestepping_order2;

	if (i_timestepping_order2 < 0)
		i_timestepping_order2 = i_timestepping_order;

	if (i_timestepping_order != i_timestepping_order2)
		SWEETErrorFatal("Orders of 1st and 2nd one must match");

	timestep_size = shackTimestepControl->currentTimestepSize;

	if (timestepping_order == 1)
	{
		timestepping_lg_irk.setup(
			ops,
			1,
			timestep_size
		);
	}
	else if (timestepping_order == 2)
	{
		if (version_id == 0)
		{
			timestepping_lg_irk.setup(
				ops,
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
			SWEETErrorFatal("Invalid version");
		}
	}
	else
	{
		SWEETErrorFatal("Invalid timestepping order");
	}


	// Only NA part!
	timestepping_ln_erk_split_vd.setup_main(
			ops,
			i_timestepping_order2,
			false, true, true, false, false
		);

	return true;
}




std::string PDESWESphere2DTS_lg_irk_lc_na_erk_vd::getIDString()
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



void PDESWESphere2DTS_lg_irk_lc_na_erk_vd::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_dt,
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
					i_dt,
					i_simulation_timestamp
				);

			// first order explicit for non-linear
			timestepping_ln_erk_split_vd.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);
		}
		else
		{
			// first order explicit for non-linear
			timestepping_ln_erk_split_vd.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// first order IRK for linear
			timestepping_lg_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
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
					i_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_ln_erk_split_vd.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_lg_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt*0.5,
					i_simulation_timestamp+i_dt*0.5
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_ln_erk_split_vd.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for linear part
			timestepping_lg_irk.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt,
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_ln_erk_split_vd.runTimestep(
					io_phi_pert, io_vrt, io_div,
					i_dt*0.5,
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




PDESWESphere2DTS_lg_irk_lc_na_erk_vd::PDESWESphere2DTS_lg_irk_lc_na_erk_vd()
{

}



PDESWESphere2DTS_lg_irk_lc_na_erk_vd::~PDESWESphere2DTS_lg_irk_lc_na_erk_vd()
{
}

