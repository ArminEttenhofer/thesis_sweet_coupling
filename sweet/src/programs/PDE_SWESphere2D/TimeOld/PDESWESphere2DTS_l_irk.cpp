/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <programs/PDE_SWESphere2D/TimeOld/DEPRECATEDTimeStepSizeChanged.hpp>
#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include "PDESWESphere2DTS_l_irk.hpp"

#include <complex>
#include "../TimeHelpers/SphBandedMatrix_GridReal.hpp"


bool PDESWESphere2DTS_l_irk::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	if (
			i_timestepping_method == "l_irk"	||
			i_timestepping_method == "lg_irk"	||
			false
	)
	{
		timestepping_method = i_timestepping_method;
		return true;
	}

	return false;
}


std::string PDESWESphere2DTS_l_irk::getIDString()
{
	return timestepping_method;
}


bool PDESWESphere2DTS_l_irk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	if (shackPDESWESphere2D->sphere2d_use_fsphere2D)
		SWEETErrorFatal("TODO: Not yet supported");

	if (timestepping_method == "l_irk")
	{
		return setup_main(
				io_ops,
				shackPDESWETimeDisc->timestepping_order,
				shackTimestepControl->currentTimestepSize,
				0.5,
				false
			);
	}
	else if (timestepping_method == "lg_irk")
	{
		return setup_main(
				io_ops,
				shackPDESWETimeDisc->timestepping_order,
				shackTimestepControl->currentTimestepSize,
				0.5,
				true
			);
	}
	else
	{
		SWEETErrorFatal("ERROR");
	}
	return false;
}

bool PDESWESphere2DTS_l_irk::setup(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_timestep_order,
		double i_timestepSize
)
{
	ops = io_ops;

	return setup_main(
			io_ops,
			i_timestep_order,
			i_timestepSize,
			shackPDESWETimeDisc->timestepping_crank_nicolson_filter,
			false
		);
}


bool PDESWESphere2DTS_l_irk::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_timestepping_order,
		double i_timestepSize,
		double i_crank_nicolson_damping_factor,
		bool i_no_coriolis
)
{
	ops = io_ops;

	timestepping_order = i_timestepping_order;
	timestep_size = i_timestepSize;
	crank_nicolson_damping_factor = i_crank_nicolson_damping_factor;

	use_f_sphere2D = shackPDESWESphere2D->sphere2d_use_fsphere2D;
	no_coriolis = i_no_coriolis;

	if (timestepping_order != 1 && timestepping_order != 2)
		SWEETErrorFatal("Only 1st and 2nd order supported!");

	clear();

	if (no_coriolis)
	{
		swe_sphere2d_ts_lg_erk.shackRegistration(this);
		swe_sphere2d_ts_lg_erk.setup_main(ops, 1);
	}
	else
	{
		if (use_f_sphere2D)
		{
			f0 = shackPDESWESphere2D->sphere2d_fsphere2d_f0;
			two_coriolis = 0.0;
		}
		else
		{
			f0 = 0.0;
			two_coriolis = 2.0*shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;
		}

		swe_sphere2d_ts_l_erk.shackRegistration(this);
		swe_sphere2d_ts_l_erk.setup_main(ops, 1);
	}


	sphere2d_radius = shackSphere2DDataOps->sphere2d_radius;

	update_coefficients(i_timestepSize);
	return true;
}


void PDESWESphere2DTS_l_irk::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (DEPRECATED_TimeStepSizeChanged::is_changed(timestep_size, i_fixed_dt, true))
		update_coefficients(i_fixed_dt);


	/*
	 * Crank-Nicolson method:
	 *
	 * (U(t+1) - q dt F(U(t+1))) = (U(t) + q dt F(U(t)))
	 *
	 * with q the CN damping facor with no damping for q=0.5
	 */
	if (timestepping_order == 2)
	{
		/*
		 * Explicit Euler
		 */
		if (no_coriolis)
			swe_sphere2d_ts_lg_erk.runTimestep(io_phi, io_vrt, io_div, dt_explicit, i_simulation_timestamp);
		else
			swe_sphere2d_ts_l_erk.runTimestep(io_phi, io_vrt, io_div, dt_explicit, i_simulation_timestamp);
	}


	double gh0 = shackPDESWESphere2D->gravitation*shackPDESWESphere2D->h0;

	if (no_coriolis)
	{
		sweet::Data::Sphere2D::DataSpectral rhs = io_div + ops->implicit_L(io_phi, dt_implicit);
		sweet::Data::Sphere2D::DataSpectral div1 = ops->implicit_helmholtz(rhs, -gh0*dt_implicit*dt_implicit, sphere2d_radius);
		sweet::Data::Sphere2D::DataSpectral phi1 = io_phi - dt_implicit*gh0*div1;
		sweet::Data::Sphere2D::DataSpectral vrt1 = io_vrt;

		io_phi = phi1;
		io_vrt = vrt1;
		io_div = div1;
	}
	else
	{
		double dt_two_omega = dt_implicit*2.0*shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;

		sweet::Data::Sphere2D::DataSpectral rhs = io_div + ops->implicit_FJinv(io_vrt, dt_two_omega) + ops->implicit_L(io_phi, dt_implicit);
		sweet::Data::Sphere2D::DataSpectral div1 = sphSolverDiv.solve(rhs);

		sweet::Data::Sphere2D::DataSpectral phi1 = io_phi - dt_implicit*gh0*div1;
		sweet::Data::Sphere2D::DataSpectral vrt1 = ops->implicit_Jinv(io_vrt - ops->implicit_F(div1, dt_two_omega), dt_two_omega);


		io_phi = phi1;
		io_vrt = vrt1;
		io_div = div1;
	}
}

void PDESWESphere2DTS_l_irk::solveImplicit(
		sweet::Data::Sphere2D::DataSpectral &io_phi,	//!< rhs variables
		sweet::Data::Sphere2D::DataSpectral &io_vrt,	//!< rhs variables
		sweet::Data::Sphere2D::DataSpectral &io_div,	//!< rhs variables

		double dt
)
{
	double gh0 = shackPDESWESphere2D->gravitation*shackPDESWESphere2D->h0;

	if (no_coriolis)
	{
		SWEETErrorFatal("Not supported yet");
	}
	else
	{
		double dt_two_omega = dt*2.0*shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;

		// Implicit update using explicit evaluation for implicit_FJinv and implicit_L (with or without NL update)
		sweet::Data::Sphere2D::DataSpectral rhs_div = io_div + ops->implicit_FJinv(io_vrt, dt_two_omega) + ops->implicit_L(io_phi, dt);
		sweet::Data::Sphere2D::DataSpectral div1 = sphSolverDiv.solve(rhs_div);

		// Update for phi using implicit update for div
		sweet::Data::Sphere2D::DataSpectral phi1 = io_phi - dt*gh0*div1;

		// Decoupled implicit update, using the implicit update for div 
		sweet::Data::Sphere2D::DataSpectral vrt1 = ops->implicit_Jinv(io_vrt - ops->implicit_F(div1, dt_two_omega), dt_two_omega);

		io_phi = phi1;
		io_vrt = vrt1;
		io_div = div1;
	}
}



void PDESWESphere2DTS_l_irk::update_coefficients(double i_timestepSize)
{
	timestep_size = i_timestepSize;

	double 	gh0 = shackPDESWESphere2D->gravitation*shackPDESWESphere2D->h0;

	if (timestepping_order == 1)
	{
		dt_explicit = -666.0;  // Yeah !!!
		dt_implicit = timestep_size;
	}
	else
	{
		dt_explicit = timestep_size*(1.0-crank_nicolson_damping_factor);
		dt_implicit = timestep_size*crank_nicolson_damping_factor;
	}

	if (!no_coriolis)
	{
		if (!use_f_sphere2D)
		{
			double dt_two_omega = dt_implicit*two_coriolis;

			sphSolverDiv.setup(ops->sphere2DDataConfig, 4);
			sphSolverDiv.solver_addComponent_implicit_J(dt_two_omega);
			sphSolverDiv.solver_addComponent_implicit_FJinvF(dt_two_omega);
			sphSolverDiv.solver_addComponent_implicit_L(gh0*dt_implicit, dt_implicit, sphere2d_radius);
		}
	}
}



PDESWESphere2DTS_l_irk::PDESWESphere2DTS_l_irk()
{
}


void PDESWESphere2DTS_l_irk::clear()
{
}




PDESWESphere2DTS_l_irk::~PDESWESphere2DTS_l_irk()
{
	clear();
}
