/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere2DTS_l_erk_na_erk_vd.hpp"


bool PDESWESphere2DTS_l_erk_na_erk_vd::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, shackPDESWETimeDisc->timestepping_order2);
}

bool PDESWESphere2DTS_l_erk_na_erk_vd::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_order,	//!< order of RK time stepping method for non-linear parts
		int i_order2	//!< order of RK time stepping method for non-linear parts
)
{
	ops = io_ops;

	timestepping_order = i_order;
	timestepping_order2 = i_order2;

	if (timestepping_order != timestepping_order2)
		SWEETErrorFatal("Timestepping orders must match");

	l_erk_split_vd.setup_main(io_ops, timestepping_order, true, true, false, false, false);
	na_erk_split_vd.setup_main(io_ops, timestepping_order, false, false, true, false, false);

	return true;
}




bool PDESWESphere2DTS_l_erk_na_erk_vd::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	timestepping_method = i_timestepping_method;
	return i_timestepping_method == "l_erk_na_erk_vd";
}


std::string PDESWESphere2DTS_l_erk_na_erk_vd::getIDString()
{
	return "l_erk_na_erk_vd";
}


void PDESWESphere2DTS_l_erk_na_erk_vd::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		l_erk_split_vd.runTimestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		na_erk_split_vd.runTimestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{

		na_erk_split_vd.runTimestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				i_simulation_timestamp
			);

		l_erk_split_vd.runTimestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		na_erk_split_vd.runTimestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				i_simulation_timestamp
			);
	}
	else
	{
		SWEETErrorFatal("Not yet supported!");
	}
}



PDESWESphere2DTS_l_erk_na_erk_vd::PDESWESphere2DTS_l_erk_na_erk_vd()
{
}



PDESWESphere2DTS_l_erk_na_erk_vd::~PDESWESphere2DTS_l_erk_na_erk_vd()
{
}

