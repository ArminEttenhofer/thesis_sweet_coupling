/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere2DTS_lg_exp_lc_n_etd_uv.hpp"


bool PDESWESphere2DTS_lg_exp_lc_n_etd_uv::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	timestepping_method = i_timestepping_method;

	if (	i_timestepping_method == "lg_exp_lc_n_etd_uv"	||
			i_timestepping_method == "lg_exp_lc_na_nr_etd_uv"	||
			i_timestepping_method == "lg_exp_lc_na_etd_uv"	||
			i_timestepping_method == "lg_exp_lc_nr_etd_uv"	||
			i_timestepping_method == "lg_exp_lc_etd_uv"	||
			false
	)
		return true;

	return false;
}


std::string PDESWESphere2DTS_lg_exp_lc_n_etd_uv::getIDString()
{
	return "lg_exp_lc_n_etd";
}


bool PDESWESphere2DTS_lg_exp_lc_n_etd_uv::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	if (shackPDESWESphere2D->sphere2d_use_fsphere2D)
		SWEETErrorFatal("TODO: Not yet supported");

	bool _with_na = false;
	bool _with_nr = false;

	if (	timestepping_method == "lg_exp_lc_n_etd_uv"	||
			timestepping_method == "lg_exp_lc_na_nr_etd_uv")
	{
		_with_na = true;
		_with_nr = true;
	}
	else if (timestepping_method == "lg_exp_lc_na_etd_uv")
	{
		_with_na = true;
		_with_nr = false;
	}
	else if (timestepping_method == "lg_exp_lc_nr_etd_uv")
	{
		_with_na = false;
		_with_nr = true;
	}
	else if (timestepping_method == "lg_exp_lc_etd_uv")
	{
		_with_na = false;
		_with_nr = false;
	}
	else
	{
		SWEETErrorFatal("Unknown TM");
	}

	return setup_main(
			ops,
			shackExpIntegration,
			shackPDESWETimeDisc->timestepping_order,
			shackTimestepControl->currentTimestepSize,
			_with_na,
			_with_nr
		);
}



bool PDESWESphere2DTS_lg_exp_lc_n_etd_uv::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		sweet::ExpIntegration::Shack *i_shackExpIntegration,
		int i_timestepping_order,
		double i_timestepSize,
		bool i_with_na,
		bool i_with_nr
)
{
	ops = io_ops;
	timestepping_order = i_timestepping_order;

	with_na = i_with_na;
	with_nr = i_with_nr;

	ts_ln_erk_split_uv.setup_main(ops, i_timestepping_order, true, true, true, true, false);

	if (timestepping_order == 0 || timestepping_order == 1)
	{
		ts_phi0_exp.setup_variant_50(ops, i_shackExpIntegration, "phi0", i_timestepSize, false, true, timestepping_order);	/* NO Coriolis */
		ts_phi1_exp.setup_variant_50(ops, i_shackExpIntegration, "phi1", i_timestepSize, false, true, timestepping_order);
	}
	else if (timestepping_order == 2)
	{
		ts_phi0_exp.setup_variant_50(ops, i_shackExpIntegration, "phi0", i_timestepSize, false, true, timestepping_order);	/* NO Coriolis */
		ts_phi1_exp.setup_variant_50(ops, i_shackExpIntegration, "phi1", i_timestepSize, false, true, timestepping_order);
		ts_phi2_exp.setup_variant_50(ops, i_shackExpIntegration, "phi2", i_timestepSize, false, true, timestepping_order);

		NU_phi_prev.setup(ops->sphere2DDataConfig);
		NU_vrt_prev.setup(ops->sphere2DDataConfig);
		NU_div_prev.setup(ops->sphere2DDataConfig);
	}
	else if (timestepping_order == 3)
	{
		ts_phi0_exp.setup_variant_50(ops, i_shackExpIntegration, "phi0", i_timestepSize, false, true, timestepping_order);	/* NO Coriolis */
		ts_phi1_exp.setup_variant_50(ops, i_shackExpIntegration, "phi1", i_timestepSize, false, true, timestepping_order);
		ts_phi2_exp.setup_variant_50(ops, i_shackExpIntegration, "phi2", i_timestepSize, false, true, timestepping_order);
		ts_phi3_exp.setup_variant_50(ops, i_shackExpIntegration, "phi3", i_timestepSize, false, true, timestepping_order);

		NU_phi_prev.setup(ops->sphere2DDataConfig);
		NU_vrt_prev.setup(ops->sphere2DDataConfig);
		NU_div_prev.setup(ops->sphere2DDataConfig);

		NU_phi_prev_2.setup(ops->sphere2DDataConfig);
		NU_vrt_prev_2.setup(ops->sphere2DDataConfig);
		NU_div_prev_2.setup(ops->sphere2DDataConfig);
	}
	else
	{
		SWEETErrorFatal("TODO: This order is not implemented, yet!");
	}

	return true;
}



void PDESWESphere2DTS_lg_exp_lc_n_etd_uv::printHelp()
{
	std::cout << "	Exponential ETD:" << std::endl;
	std::cout << "		+ lg_exp_lc_n_etd_uv" << std::endl;
	std::cout << "		+ lg_exp_lc_na_nr_etd_uv" << std::endl;
	std::cout << "		+ lg_exp_lc_na_etd_uv" << std::endl;
	std::cout << "		+ lg_exp_lc_nr_etd_uv" << std::endl;
	std::cout << "		+ lg_exp_lc_etd_uv" << std::endl;
}




void PDESWESphere2DTS_lg_exp_lc_n_etd_uv::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_U_phi,
		sweet::Data::Sphere2D::DataSpectral &io_U_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_U_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	const sweet::Data::Sphere2D::Config *sphere2DDataConfig = io_U_phi.sphere2DDataConfig;

	if (timestepping_order == 1)
	{
		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}
		 * 			+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 */

		sweet::Data::Sphere2D::DataSpectral phi0_Un_phi(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi0_Un_vrt(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi0_Un_div(sphere2DDataConfig);

		ts_phi0_exp.runTimestep(
				io_U_phi, io_U_vrt, io_U_div,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		sweet::Data::Sphere2D::DataSpectral F_phi(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral F_vrt(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral F_div(sphere2DDataConfig, 0);

		ts_ln_erk_split_uv.euler_timestep_update_lc(
				io_U_phi, io_U_vrt, io_U_div,
				F_phi, F_vrt, F_div,
				i_simulation_timestamp
		);


		if (with_na)
		{
			ts_ln_erk_split_uv.euler_timestep_update_na(
					io_U_phi, io_U_vrt, io_U_div,
					F_phi, F_vrt, F_div,
					i_simulation_timestamp
			);
		}

		if (with_nr)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					io_U_phi, io_U_vrt, io_U_div,
					F_phi, F_vrt, F_div,
					i_simulation_timestamp
			);
		}


		sweet::Data::Sphere2D::DataSpectral phi1_FUn_phi(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi1_FUn_vrt(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi1_FUn_div(sphere2DDataConfig);

		ts_phi1_exp.runTimestep(
				F_phi, F_vrt, F_div,
				phi1_FUn_phi, phi1_FUn_vrt, phi1_FUn_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		io_U_phi = phi0_Un_phi + i_fixed_dt*phi1_FUn_phi;
		io_U_vrt = phi0_Un_vrt + i_fixed_dt*phi1_FUn_vrt;
		io_U_div = phi0_Un_div + i_fixed_dt*phi1_FUn_div;
	}
	else if (timestepping_order == 2)
	{
		const sweet::Data::Sphere2D::DataSpectral &U_phi = io_U_phi;
		const sweet::Data::Sphere2D::DataSpectral &U_vrt = io_U_vrt;
		const sweet::Data::Sphere2D::DataSpectral &U_div = io_U_div;

		if (i_simulation_timestamp == 0)
		{
			/*
			 * First time step:
			 * Simply backup existing fields for multi-step parts of this algorithm.
			 */
			NU_phi_prev.spectral_setZero();
			NU_vrt_prev.spectral_setZero();
			NU_div_prev.spectral_setZero();

			ts_ln_erk_split_uv.euler_timestep_update_lc(
					U_phi, U_vrt, U_div,
					NU_phi_prev, NU_vrt_prev, NU_div_prev,
					i_simulation_timestamp
			);

			if (with_na)
			{
				ts_ln_erk_split_uv.euler_timestep_update_na(
						U_phi, U_vrt, U_div,
						NU_phi_prev, NU_vrt_prev, NU_div_prev,
						i_simulation_timestamp
				);
			}

			if (with_nr)
			{
				ts_ln_erk_split_uv.euler_timestep_update_nr(
						U_phi, U_vrt, U_div,
						NU_phi_prev, NU_vrt_prev, NU_div_prev,
						i_simulation_timestamp
				);
			}

		}


		/*
		 * u_{n+1} = \varphi_{0} (\Delta t L )u_{n} + \Delta t K
		 *
		 * K = \left[
		 * 		\varphi_{1}(\Delta t L)F(u_{n}) +
		 * 		\varphi_{2}(\Delta t L) (F(u_{n})-F(u_{n-1}))
		 * 	\right]
		 */

		/*
		 * phi0
		 */
		sweet::Data::Sphere2D::DataSpectral phi0_phi(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi0_vrt(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi0_div(sphere2DDataConfig);

		ts_phi0_exp.runTimestep(
				U_phi, U_vrt, U_div,
				phi0_phi, phi0_vrt, phi0_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		/*
		 * phi1
		 */
		sweet::Data::Sphere2D::DataSpectral F_phi(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral F_vrt(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral F_div(sphere2DDataConfig, 0);

		ts_ln_erk_split_uv.euler_timestep_update_lc(
				io_U_phi, io_U_vrt, io_U_div,
				F_phi, F_vrt, F_div,
				i_simulation_timestamp
		);

		if (with_na)
		{
			ts_ln_erk_split_uv.euler_timestep_update_na(
					io_U_phi, io_U_vrt, io_U_div,
					F_phi, F_vrt, F_div,
					i_simulation_timestamp
			);
		}

		if (with_nr)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					io_U_phi, io_U_vrt, io_U_div,
					F_phi, F_vrt, F_div,
					i_simulation_timestamp
			);
		}

		sweet::Data::Sphere2D::DataSpectral phi1_phi(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi1_vrt(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi1_div(sphere2DDataConfig);

		ts_phi1_exp.runTimestep(
				F_phi, F_vrt, F_div,
				phi1_phi, phi1_vrt, phi1_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		/*
		 * phi2
		 */

		sweet::Data::Sphere2D::DataSpectral phi2_phi(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi2_vrt(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi2_div(sphere2DDataConfig);

		ts_phi2_exp.runTimestep(
				F_phi - NU_phi_prev,
				F_vrt - NU_vrt_prev,
				F_div - NU_div_prev,

				phi2_phi,
				phi2_vrt,
				phi2_div,

				i_fixed_dt,
				i_simulation_timestamp
			);

		io_U_phi = phi0_phi + i_fixed_dt*(phi1_phi + phi2_phi);
		io_U_vrt = phi0_vrt + i_fixed_dt*(phi1_vrt + phi2_vrt);
		io_U_div = phi0_div + i_fixed_dt*(phi1_div + phi2_div);

		/*
		 * Backup nonlinear evaluations
		 */
		{
			NU_phi_prev = F_phi;
			NU_vrt_prev = F_vrt;
			NU_div_prev = F_div;
		}
	}
	else if (timestepping_order == 3)
	{
		const sweet::Data::Sphere2D::DataSpectral &U_phi = io_U_phi;
		const sweet::Data::Sphere2D::DataSpectral &U_vrt = io_U_vrt;
		const sweet::Data::Sphere2D::DataSpectral &U_div = io_U_div;

		if (i_simulation_timestamp == 0)
		{
			/*
			 * First time step:
			 * Simply backup existing fields for multi-step parts of this algorithm.
			 */
			NU_phi_prev.spectral_setZero();
			NU_vrt_prev.spectral_setZero();
			NU_div_prev.spectral_setZero();

			ts_ln_erk_split_uv.euler_timestep_update_lc(
					U_phi, U_vrt, U_div,
					NU_phi_prev, NU_vrt_prev, NU_div_prev,
					i_simulation_timestamp
			);

			if (with_na)
			{
				ts_ln_erk_split_uv.euler_timestep_update_na(
						U_phi, U_vrt, U_div,
						NU_phi_prev, NU_vrt_prev, NU_div_prev,
						i_simulation_timestamp
				);
			}

			if (with_nr)
			{
				ts_ln_erk_split_uv.euler_timestep_update_nr(
						U_phi, U_vrt, U_div,
						NU_phi_prev, NU_vrt_prev, NU_div_prev,
						i_simulation_timestamp
				);
			}

			NU_phi_prev_2 = NU_phi_prev;
			NU_vrt_prev_2 = NU_vrt_prev;
			NU_div_prev_2 = NU_div_prev;
		}


		/*
		 * Compute nonlinear tendencies
		 */
		sweet::Data::Sphere2D::DataSpectral NU_phi(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral NU_vrt(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral NU_div(sphere2DDataConfig, 0);

		ts_ln_erk_split_uv.euler_timestep_update_lc(
				io_U_phi, io_U_vrt, io_U_div,
				NU_phi, NU_vrt, NU_div,
				i_simulation_timestamp
		);

		if (with_na)
		{
			ts_ln_erk_split_uv.euler_timestep_update_na(
					io_U_phi, io_U_vrt, io_U_div,
					NU_phi, NU_vrt, NU_div,
					i_simulation_timestamp
			);
		}

		if (with_nr)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					io_U_phi, io_U_vrt, io_U_div,
					NU_phi, NU_vrt, NU_div,
					i_simulation_timestamp
			);
		}


		/*
		 * u_{n+1} = \varphi_{0} (\Delta t L )u_{n} + \Delta t K
		 *
		 * K = \left[
		 * 		\varphi_{1}(\Delta t L)F(u_{n}) +
		 * 		\varphi_{2}(\Delta t L) (F(u_{n})-F(u_{n-1}))
		 * 	\right]
		 */

		/*
		 * phi0
		 */
		sweet::Data::Sphere2D::DataSpectral phi0_phi(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi0_vrt(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi0_div(sphere2DDataConfig);

		ts_phi0_exp.runTimestep(
				U_phi, U_vrt, U_div,
				phi0_phi, phi0_vrt, phi0_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		/*
		 * phi1
		 */
		sweet::Data::Sphere2D::DataSpectral phi1_phi(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi1_vrt(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi1_div(sphere2DDataConfig);

		ts_phi1_exp.runTimestep(
				NU_phi, NU_vrt, NU_div,
				phi1_phi, phi1_vrt, phi1_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		/*
		 * phi2
		 */
		sweet::Data::Sphere2D::DataSpectral phi2_phi(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi2_vrt(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi2_div(sphere2DDataConfig);

		ts_phi2_exp.runTimestep(
				3.0/2.0*NU_phi - 2.0*NU_phi_prev + 1.0/2.0*NU_phi_prev_2,
				3.0/2.0*NU_vrt - 2.0*NU_vrt_prev + 1.0/2.0*NU_vrt_prev_2,
				3.0/2.0*NU_div - 2.0*NU_div_prev + 1.0/2.0*NU_div_prev_2,

				phi2_phi,
				phi2_vrt,
				phi2_div,

				i_fixed_dt,
				i_simulation_timestamp
			);


		/*
		 * phi3
		 */
		sweet::Data::Sphere2D::DataSpectral phi3_phi(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi3_vrt(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral phi3_div(sphere2DDataConfig);

		ts_phi3_exp.runTimestep(
				1.0/2.0*NU_phi - NU_phi_prev + 1.0/2.0*NU_phi_prev_2,
				1.0/2.0*NU_vrt - NU_vrt_prev + 1.0/2.0*NU_vrt_prev_2,
				1.0/2.0*NU_div - NU_div_prev + 1.0/2.0*NU_div_prev_2,

				phi3_phi,
				phi3_vrt,
				phi3_div,

				i_fixed_dt,
				i_simulation_timestamp
			);

		/*
		 * Compute full time step
		 */
		io_U_phi = phi0_phi + i_fixed_dt*(phi1_phi + phi2_phi + phi3_phi);
		io_U_vrt = phi0_vrt + i_fixed_dt*(phi1_vrt + phi2_vrt + phi3_vrt);
		io_U_div = phi0_div + i_fixed_dt*(phi1_div + phi2_div + phi3_div);

		/*
		 * Backup nonlinear evaluations
		 */
		{
			NU_phi_prev_2 = NU_phi_prev;
			NU_vrt_prev_2 = NU_vrt_prev;
			NU_div_prev_2 = NU_div_prev;

			NU_phi_prev = NU_phi;
			NU_vrt_prev = NU_vrt;
			NU_div_prev = NU_div;
		}
	}
	else
	{
		SWEETErrorFatal("TODO: This order is not implemented, yet!");
	}
}



PDESWESphere2DTS_lg_exp_lc_n_etd_uv::PDESWESphere2DTS_lg_exp_lc_n_etd_uv()
{
}

PDESWESphere2DTS_lg_exp_lc_n_etd_uv::~PDESWESphere2DTS_lg_exp_lc_n_etd_uv()
{
}

