/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */



#include "PDESWESphere2DTS_lg_exp_lc_taylor.hpp"



bool PDESWESphere2DTS_lg_exp_lc_taylor::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order);
}




bool PDESWESphere2DTS_lg_exp_lc_taylor::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_order
)
{
	ops = io_ops;

	if (shackExpIntegration->taylor_num_expansions <= 0)
		SWEETErrorFatal("This time stepping method requires setting the number of Taylor expansions");

	timestepping_order = i_order;

	if (timestepping_order != 1 && timestepping_order != 2)
		SWEETErrorFatal("Only 1st and 2nd order time stepping order supported");

	timestepping_lg_exp.setup_main(io_ops, "phi0");

	if (shackExpIntegration->exp_method != "ss_taylor")
		SWEETErrorFatal("Use --exp-method=ss_taylor to use this time stepper [TODO: This might not make sense after some internal changes]");

	return true;
}



void PDESWESphere2DTS_lg_exp_lc_taylor::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		run_timestep_lc(io_phi_pert, io_vrt, io_div, i_fixed_dt*0.5, i_simulation_timestamp);
		run_timestep_lg(io_phi_pert, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
		run_timestep_lc(io_phi_pert, io_vrt, io_div, i_fixed_dt*0.5, i_simulation_timestamp);
	}
	else if (timestepping_order == 2)
	{
		run_timestep_lc(io_phi_pert, io_vrt, io_div, i_fixed_dt*0.5, i_simulation_timestamp);
		run_timestep_lg(io_phi_pert, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
		run_timestep_lc(io_phi_pert, io_vrt, io_div, i_fixed_dt*0.5, i_simulation_timestamp);
	}
	else
	{
		SWEETErrorFatal("This time stepping order is not yet supported!");
	}
}



void PDESWESphere2DTS_lg_exp_lc_taylor::run_timestep_lg(
		sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	timestepping_lg_exp.runTimestep(
			io_phi_pert,
			io_vrt,
			io_div,
			i_fixed_dt,
			i_simulation_timestamp
		);
}



void PDESWESphere2DTS_lg_exp_lc_taylor::run_timestep_lc(
		sweet::Data::Sphere2D::DataSpectral &io_phi,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_dt,
		double i_simulation_timestamp
)
{
	/*
	 * Strang split version
	 *
	 * Treat $L_g$ terms exponentially
	 * Treat $L_c$ terms with Taylor expansion
	 */
	sweet::Data::Sphere2D::DataSpectral phi_Lc_pow_i_prev = io_phi;
	sweet::Data::Sphere2D::DataSpectral vrt_Lc_pow_i_prev = io_vrt;
	sweet::Data::Sphere2D::DataSpectral div_Lc_pow_i_prev = io_div;

	double taylor_coeff = 1.0;

	for (int i = 1; i <= shackExpIntegration->taylor_num_expansions; i++)
	{
		// Dummy output variables
		sweet::Data::Sphere2D::DataSpectral phi_Lc_pow_i(io_phi.sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral vrt_Lc_pow_i(io_phi.sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral div_Lc_pow_i(io_phi.sphere2DDataConfig, 0);

		// Time tendencies
		timestepping_lc_erk.euler_timestep_update_lc_spectral_only(
				phi_Lc_pow_i_prev, vrt_Lc_pow_i_prev, div_Lc_pow_i_prev,
				phi_Lc_pow_i, vrt_Lc_pow_i, div_Lc_pow_i,
				i_simulation_timestamp
		);

		// Taylor coefficients
		taylor_coeff *= i_dt/(double)i;

		io_phi += taylor_coeff*phi_Lc_pow_i;
		io_vrt += taylor_coeff*vrt_Lc_pow_i;
		io_div += taylor_coeff*div_Lc_pow_i;

		phi_Lc_pow_i_prev.swap(phi_Lc_pow_i);
		vrt_Lc_pow_i_prev.swap(vrt_Lc_pow_i);
		div_Lc_pow_i_prev.swap(div_Lc_pow_i);
	}
}



PDESWESphere2DTS_lg_exp_lc_taylor::PDESWESphere2DTS_lg_exp_lc_taylor()
{
}



PDESWESphere2DTS_lg_exp_lc_taylor::~PDESWESphere2DTS_lg_exp_lc_taylor()
{
}

