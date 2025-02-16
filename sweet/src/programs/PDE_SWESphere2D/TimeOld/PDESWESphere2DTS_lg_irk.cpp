/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */


#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include "PDESWESphere2DTS_lg_irk.hpp"

#include <complex>
#include <sweet/Error/Fatal.hpp>
#include "../TimeHelpers/SphBandedMatrix_GridReal.hpp"



bool PDESWESphere2DTS_lg_irk::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	/*
	 * Supported directly in l_irk, not in this class anymore
	 */
	timestepping_method = i_timestepping_method;
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	if (i_timestepping_method == "lg_irk_DEPRECATED")
		return true;

	return false;
}



std::string PDESWESphere2DTS_lg_irk::getIDString()
{
	return "lg_irk_DEPRECATED";
}



bool PDESWESphere2DTS_lg_irk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	if (shackPDESWESphere2D->sphere2d_use_fsphere2D)
		SWEETErrorFatal("TODO: Not yet supported");

	return setup_main(
			io_ops,
			shackPDESWETimeDisc->timestepping_order,
			shackTimestepControl->currentTimestepSize,
			shackPDESWETimeDisc->timestepping_crank_nicolson_filter
		);
}


bool PDESWESphere2DTS_lg_irk::setup(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_timestep_order,
		double i_timestepSize
)
{
	return setup_main(
			io_ops,
			i_timestep_order,
			i_timestepSize,
			shackPDESWETimeDisc->timestepping_crank_nicolson_filter
		);
}



bool PDESWESphere2DTS_lg_irk::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_timestep_order,
		double i_timestepSize,
		double i_crank_nicolson_damping_factor
)
{
	ops = io_ops;

	timestepping_order = i_timestep_order;
	timestep_size = i_timestepSize;

	if (i_timestep_order == 1)
	{
		// set this to 1 to ignore it
		crank_nicolson_damping_factor = 1.0;
	}
	else if (i_timestep_order == 2)
	{
		crank_nicolson_damping_factor = i_crank_nicolson_damping_factor;
		lg_erk.setup_main(ops, 1);
	}
	else
	{
		SWEETErrorFatal("Only 1st and 2nd order IRK supported so far with this implementation! Use l_cn if you want to have 2nd order Crank-Nicolson!");
	}


	update_coefficients();

	r = shackSphere2DDataOps->sphere2d_radius;
	inv_r = 1.0/r;

	gh = shackPDESWESphere2D->gravitation*shackPDESWESphere2D->h0;

	return true;
}




void PDESWESphere2DTS_lg_irk::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,		
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		SWEETErrorFatal("Only constant time step size allowed");

	if (std::abs(timestep_size - i_fixed_dt)/std::max(timestep_size, i_fixed_dt) > 1e-10)
	{
		//std::cout << "Warning: Reducing time step size from " << i_fixed_dt << " to " << timestep_size << std::endl;
		timestep_size = i_fixed_dt;

		update_coefficients();
	}

	if (timestepping_order == 2)
	{
		/*
		 * Execute a half ERK time step first for 2nd order
		 */
		lg_erk.runTimestep(io_phi_pert, io_vrt, io_div, i_fixed_dt*(1.0-crank_nicolson_damping_factor), i_simulation_timestamp);
	}

	sweet::Data::Sphere2D::DataSpectral phi0 = io_phi_pert;
	sweet::Data::Sphere2D::DataSpectral vrt0 = io_vrt;
	sweet::Data::Sphere2D::DataSpectral div0 = io_div;

	sweet::Data::Sphere2D::DataSpectral phi(ops->sphere2DDataConfig);
	sweet::Data::Sphere2D::DataSpectral vort(ops->sphere2DDataConfig);
	sweet::Data::Sphere2D::DataSpectral div(ops->sphere2DDataConfig);

	{
		sweet::Data::Sphere2D::DataSpectral rhs = gh*div0 + alpha*phi0;
		phi = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);
		io_phi_pert = phi*beta;

		rhs = alpha*div0 + ops->laplace(phi0);
		div = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);
		io_div = div*beta;

		io_vrt = vrt0;
	}
}



PDESWESphere2DTS_lg_irk::PDESWESphere2DTS_lg_irk()
{
}


void PDESWESphere2DTS_lg_irk::update_coefficients()
{
	alpha = -1.0/timestep_size;
	beta = -1.0/timestep_size;

	{
		/*
		 * Crank-Nicolson method:
		 *
		 * (U(t+1) - q dt F(U(t+1))) = (U(t) + q dt F(U(t)))
		 *
		 * with q the CN damping facor with no damping for q=0.5
		 */

		// scale dt by the damping factor to reuse solver structure

		alpha /= crank_nicolson_damping_factor;
		beta /= crank_nicolson_damping_factor;
	}
}



PDESWESphere2DTS_lg_irk::~PDESWESphere2DTS_lg_irk()
{
}
