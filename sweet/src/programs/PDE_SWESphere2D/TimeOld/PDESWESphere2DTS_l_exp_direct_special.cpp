/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */


#include "PDESWESphere2DTS_l_exp_direct_special.hpp"



bool PDESWESphere2DTS_l_exp_direct_special::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	if (timestepping_method == "lg_exp_special")
		return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, false, "phi0");
	else if (timestepping_method == "l_exp_special")
		return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, true, "phi0");

	return false;
}




bool PDESWESphere2DTS_l_exp_direct_special::setup_main(
	const sweet::Data::Sphere2D::Operators *io_ops,
	int i_order,
	bool i_use_coriolis,	//!< Include Coriolis term
	const std::string &i_function_name
)
{
	ops = io_ops;

	timestepping_order = i_order;
	use_coriolis = i_use_coriolis;

	if (use_coriolis)
		timestepping_method = "l_exp_special";
	else
		timestepping_method = "lg_exp_special";


	if (timestepping_order != 1 && timestepping_order != 2 && timestepping_order != 4)
	{
		std::ostringstream ss;
		ss << "Time stepping order " << i_order << " requested";
		SWEETErrorFatal(ss.str());
	}

	if (use_coriolis)
	{
		if (i_function_name != "phi0")
		{
			SWEETErrorFatal("Only phi0 functions supported for the ETDnRK lg/lc scheme");
		}

		timestepping_lg_exp_phi0.setup_main(io_ops, "phi0");
		timestepping_lg_exp_phi1.setup_main(io_ops, "phi1");

		if (timestepping_order >= 2)
		{
			timestepping_lg_exp_phi2.setup_main(io_ops, "phi2");

			if (timestepping_order >= 4)
			{
				timestepping_lg_exp_ups1.setup_main(io_ops, "ups1");
				timestepping_lg_exp_ups2.setup_main(io_ops, "ups2");
				timestepping_lg_exp_ups3.setup_main(io_ops, "ups3");
			}
		}
	}
	else
	{
		timestepping_lg_exp_phi0.setup_main(io_ops, i_function_name);
	}

	return true;
}


bool PDESWESphere2DTS_l_exp_direct_special::implementsTimesteppingMethod(
	const std::string &i_timestepping_method
)
{
	timestepping_method = i_timestepping_method;

	if (timestepping_method == "l_exp_special")
		return true;

	if (timestepping_method == "lg_exp_special")
		return true;

	return false;
}


bool PDESWESphere2DTS_l_exp_direct_special::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	PDESWESphere2DTS_BaseInterface::shackRegistration(io_shackDict);

	timestepping_lg_exp_phi0.shackRegistration(io_shackDict);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timestepping_lg_exp_phi0);

	timestepping_lg_exp_phi1.shackRegistration(io_shackDict);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timestepping_lg_exp_phi1);

	timestepping_lg_exp_phi2.shackRegistration(io_shackDict);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timestepping_lg_exp_phi2);

	timestepping_lg_exp_ups1.shackRegistration(io_shackDict);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timestepping_lg_exp_ups1);

	timestepping_lg_exp_ups2.shackRegistration(io_shackDict);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timestepping_lg_exp_ups2);

	timestepping_lg_exp_ups3.shackRegistration(io_shackDict);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timestepping_lg_exp_ups3);

	timestepping_lc_erk.shackRegistration(io_shackDict);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timestepping_lc_erk);
	return true;
}


PDESWESphere2DTS_l_exp_direct_special::PDESWESphere2DTS_l_exp_direct_special()	:
		use_coriolis(false)
{
}



void PDESWESphere2DTS_l_exp_direct_special::euler_timestep_store_update_lc(
		const sweet::Data::Sphere2D::DataSpectral &i_phi_pert,
		const sweet::Data::Sphere2D::DataSpectral &i_vrt,
		const sweet::Data::Sphere2D::DataSpectral &i_div,
		sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &o_vrt,
		sweet::Data::Sphere2D::DataSpectral &o_div,
		double i_simulation_timestamp
)
{
	/*
	 * We have this special function for the lc updates since we need
	 * to initialize the output variables with 0 values and this makes
	 * it more comfortable.
	 */
	o_phi_pert.spectral_setZero();
	o_vrt.spectral_setZero();
	o_div.spectral_setZero();

	//timestepping_lc_erk.euler_timestep_update_lc_spectral_only(
	timestepping_lc_erk.euler_timestep_update_lc(
			i_phi_pert, i_vrt, i_div,
			o_phi_pert, o_vrt, o_div,
			i_simulation_timestamp
		);
}



void PDESWESphere2DTS_l_exp_direct_special::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	const sweet::Data::Sphere2D::Config *sphere2DDataConfig = io_phi_pert.sphere2DDataConfig;

	if (!use_coriolis)
	{
		timestepping_lg_exp_phi0.runTimestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		return;
	}

	if (timestepping_order == 1)
	{
		/*
		 * U_{1} =	\psi_{0}( \Delta t L ) U_{0}
		 * 			+ \Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 */

		sweet::Data::Sphere2D::DataSpectral phi0_Un_phi_pert(sphere2DDataConfig), phi0_Un_vrt(sphere2DDataConfig), phi0_Un_div(sphere2DDataConfig);
		timestepping_lg_exp_phi0.runTimestep(
				io_phi_pert, io_vrt, io_div,
				phi0_Un_phi_pert, phi0_Un_vrt, phi0_Un_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		sweet::Data::Sphere2D::DataSpectral FUn_phi_pert(sphere2DDataConfig), FUn_vrt(sphere2DDataConfig), FUn_div(sphere2DDataConfig);
		euler_timestep_store_update_lc(
				io_phi_pert, io_vrt, io_div,
				FUn_phi_pert, FUn_vrt, FUn_div,
				i_simulation_timestamp
			);

		sweet::Data::Sphere2D::DataSpectral phi1_FUn_phi_pert(sphere2DDataConfig), phi1_FUn_vrt(sphere2DDataConfig), phi1_FUn_div(sphere2DDataConfig);
		timestepping_lg_exp_phi1.runTimestep(
				FUn_phi_pert, FUn_vrt, FUn_div,
				phi1_FUn_phi_pert, phi1_FUn_vrt, phi1_FUn_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		io_phi_pert = phi0_Un_phi_pert + i_fixed_dt*phi1_FUn_phi_pert;
		io_vrt = phi0_Un_vrt + i_fixed_dt*phi1_FUn_vrt;
		io_div = phi0_Un_div + i_fixed_dt*phi1_FUn_div;
		return;
	}

	if (timestepping_order == 2)
	{
		/*
		 * A_{n}=\psi_{0}(\Delta tL)U_{n}+\Delta t\psi_{1}(\Delta tL)F(U_{n})
		 */

		sweet::Data::Sphere2D::DataSpectral phi0_Un_h(sphere2DDataConfig), phi0_Un_u(sphere2DDataConfig), phi0_Un_v(sphere2DDataConfig);

		timestepping_lg_exp_phi0.runTimestep(
				io_phi_pert, io_vrt, io_div,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_fixed_dt,
				i_simulation_timestamp
			);

		sweet::Data::Sphere2D::DataSpectral FUn_h(sphere2DDataConfig), FUn_u(sphere2DDataConfig), FUn_v(sphere2DDataConfig);

		euler_timestep_store_update_lc(
				io_phi_pert, io_vrt, io_div,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);

		sweet::Data::Sphere2D::DataSpectral phi1_FUn_h(sphere2DDataConfig), phi1_FUn_u(sphere2DDataConfig), phi1_FUn_v(sphere2DDataConfig);

		timestepping_lg_exp_phi1.runTimestep(
				FUn_h, FUn_u, FUn_v,
				phi1_FUn_h, phi1_FUn_u, phi1_FUn_v,
				i_fixed_dt,
				i_simulation_timestamp
			);

		sweet::Data::Sphere2D::DataSpectral A_h = phi0_Un_h + i_fixed_dt*phi1_FUn_h;
		sweet::Data::Sphere2D::DataSpectral A_u = phi0_Un_u + i_fixed_dt*phi1_FUn_u;
		sweet::Data::Sphere2D::DataSpectral A_v = phi0_Un_v + i_fixed_dt*phi1_FUn_v;

		/*
		 * U_{n+1} = A_{n}+ \Delta t \psi_{2}(\Delta tL)
		 * 				\left(F(A_{n},t_{n}+\Delta t)-F(U_{n})\right)
		 */

		sweet::Data::Sphere2D::DataSpectral FAn_h(sphere2DDataConfig), FAn_u(sphere2DDataConfig), FAn_v(sphere2DDataConfig);

		euler_timestep_store_update_lc(
				A_h, A_u, A_v,
				FAn_h, FAn_u, FAn_v,
				i_simulation_timestamp
		);


		sweet::Data::Sphere2D::DataSpectral phi2_X_h(sphere2DDataConfig), phi2_X_u(sphere2DDataConfig), phi2_X_v(sphere2DDataConfig);

		timestepping_lg_exp_phi2.runTimestep(
				FAn_h - FUn_h,
				FAn_u - FUn_u,
				FAn_v - FUn_v,

				phi2_X_h,
				phi2_X_u,
				phi2_X_v,

				i_fixed_dt,
				i_simulation_timestamp
			);

		io_phi_pert = A_h + i_fixed_dt*phi2_X_h;
		io_vrt = A_u + i_fixed_dt*phi2_X_u;
		io_div = A_v + i_fixed_dt*phi2_X_v;
		return;
	}

	if (timestepping_order == 4)
	{
		double dt = i_fixed_dt;
		double dt_half = dt*0.5;

		/*
		 * Precompute commonly used terms
		 */
		sweet::Data::Sphere2D::DataSpectral phi0_Un_half_phi_pert(sphere2DDataConfig), phi0_Un_half_vrt(sphere2DDataConfig), phi0_Un_half_div(sphere2DDataConfig);

		timestepping_lg_exp_phi0.runTimestep(
				io_phi_pert, io_vrt, io_div,
				phi0_Un_half_phi_pert, phi0_Un_half_vrt, phi0_Un_half_div,
				dt_half,
				i_simulation_timestamp
			);

		sweet::Data::Sphere2D::DataSpectral FUn_phi_pert(sphere2DDataConfig), FUn_vrt(sphere2DDataConfig), FUn_div(sphere2DDataConfig);
		euler_timestep_store_update_lc(
				io_phi_pert, io_vrt, io_div,
				FUn_phi_pert, FUn_vrt, FUn_div,
				i_simulation_timestamp
		);

		/*
		 * A_{n} = \phi_{0}\left(0.5\Delta tL\right)U_{n}+0.5\Delta t\phi_{1}\left(0.5\Delta tL\right)F(U_{n},t_{n})
		 */
		sweet::Data::Sphere2D::DataSpectral phi1_phi_pert(sphere2DDataConfig), phi1_vrt(sphere2DDataConfig), phi1_div(sphere2DDataConfig);
		timestepping_lg_exp_phi1.runTimestep(
				FUn_phi_pert, FUn_vrt, FUn_div,
				phi1_phi_pert, phi1_vrt, phi1_div,
				dt_half,
				i_simulation_timestamp
			);

		sweet::Data::Sphere2D::DataSpectral A_phi_pert = phi0_Un_half_phi_pert + dt_half*phi1_phi_pert;
		sweet::Data::Sphere2D::DataSpectral A_vrt = phi0_Un_half_vrt + dt_half*phi1_vrt;
		sweet::Data::Sphere2D::DataSpectral A_div = phi0_Un_half_div + dt_half*phi1_div;


		/*
		 * B_{n} = \phi_{0}\left(0.5\Delta tL\right)U_{n}+0.5\Delta t\phi_{1}\left(0.5\Delta tL\right)F(A_{n},t_{n}+0.5\Delta t)
		 */
		sweet::Data::Sphere2D::DataSpectral FAn_phi_pert(sphere2DDataConfig), FAn_vrt(sphere2DDataConfig), FAn_div(sphere2DDataConfig);

		euler_timestep_store_update_lc(
				A_phi_pert, A_vrt, A_div,
				FAn_phi_pert, FAn_vrt, FAn_div,
				i_simulation_timestamp + dt_half
		);

		timestepping_lg_exp_phi1.runTimestep(
				FAn_phi_pert, FAn_vrt, FAn_div,
				phi1_phi_pert, phi1_vrt, phi1_div,
				dt_half,
				i_simulation_timestamp
			);

		sweet::Data::Sphere2D::DataSpectral B_phi_pert = phi0_Un_half_phi_pert + dt_half*phi1_phi_pert;
		sweet::Data::Sphere2D::DataSpectral B_vrt = phi0_Un_half_vrt + dt_half*phi1_vrt;
		sweet::Data::Sphere2D::DataSpectral B_div = phi0_Un_half_div + dt_half*phi1_div;


		/*
		 * C_{n} = \phi_{0}\left(0.5\Delta tL\right)A_{n}+0.5\Delta t\phi_{1}\left(0.5\Delta tL\right)\left(2F(B_{n},t_{n}+0.5\Delta t)-F(U_{n},t_{n})\right)
		 */
		sweet::Data::Sphere2D::DataSpectral FBn_phi_pert(sphere2DDataConfig), FBn_vrt(sphere2DDataConfig), FBn_div(sphere2DDataConfig);

		euler_timestep_store_update_lc(
				B_phi_pert, B_vrt, B_div,
				FBn_phi_pert, FBn_vrt, FBn_div,
				i_simulation_timestamp + dt_half
		);

		timestepping_lg_exp_phi1.runTimestep(
				2.0*FBn_phi_pert - FUn_phi_pert,
				2.0*FBn_vrt - FUn_vrt,
				2.0*FBn_div - FUn_div,
				phi1_phi_pert,	phi1_vrt,	phi1_div,
				dt_half,
				i_simulation_timestamp
			);

		sweet::Data::Sphere2D::DataSpectral phi0_An_phi_pert(sphere2DDataConfig), phi0_An_vrt(sphere2DDataConfig), phi0_An_div(sphere2DDataConfig);
		timestepping_lg_exp_phi0.runTimestep(
				A_phi_pert, A_vrt, A_div,
				phi0_An_phi_pert, phi0_An_vrt, phi0_An_div,
				dt_half,
				i_simulation_timestamp
			);

		sweet::Data::Sphere2D::DataSpectral C_phi_pert = phi0_An_phi_pert + dt_half*phi1_phi_pert;
		sweet::Data::Sphere2D::DataSpectral C_vrt = phi0_An_vrt + dt_half*phi1_vrt;
		sweet::Data::Sphere2D::DataSpectral C_div = phi0_An_div + dt_half*phi1_div;



		/*
		 * R0 - R3
		 */
		sweet::Data::Sphere2D::DataSpectral FCn_phi_pert(sphere2DDataConfig), FCn_vrt(sphere2DDataConfig), FCn_div(sphere2DDataConfig);

		euler_timestep_store_update_lc(
				C_phi_pert, C_vrt, C_div,
				FCn_phi_pert, FCn_vrt, FCn_div,
				i_simulation_timestamp + dt
		);

		sweet::Data::Sphere2D::DataSpectral R0_phi_pert = io_phi_pert;
		sweet::Data::Sphere2D::DataSpectral R0_vrt = io_vrt;
		sweet::Data::Sphere2D::DataSpectral R0_div = io_div;

		sweet::Data::Sphere2D::DataSpectral &R1_phi_pert = FUn_phi_pert;
		sweet::Data::Sphere2D::DataSpectral &R1_vrt = FUn_vrt;
		sweet::Data::Sphere2D::DataSpectral &R1_div = FUn_div;

		sweet::Data::Sphere2D::DataSpectral R2_phi_pert = FAn_phi_pert + FBn_phi_pert;
		sweet::Data::Sphere2D::DataSpectral R2_vrt = FAn_vrt + FBn_vrt;
		sweet::Data::Sphere2D::DataSpectral R2_div = FAn_div + FBn_div;

		sweet::Data::Sphere2D::DataSpectral &R3_phi_pert = FCn_phi_pert;
		sweet::Data::Sphere2D::DataSpectral &R3_vrt = FCn_vrt;
		sweet::Data::Sphere2D::DataSpectral &R3_div = FCn_div;


		/*
		 * U_{n+1} =
		 * 		\psi_{0}(\Delta tL) R_{0}
		 * 			+ \Delta t
		 * 			(
		 * 				  \upsilon_{1}(\Delta tL) R_{1} +
		 * 				2*\upsilon_{2}(\Delta tL) R_{2} +
		 * 				  \upsilon_{3}(\Delta tL) R_{3}
		 * 			)
		 */
		timestepping_lg_exp_phi0.runTimestep(
				R0_phi_pert, R0_vrt, R0_div,
				dt,		i_simulation_timestamp
			);

		timestepping_lg_exp_ups1.runTimestep(
				R1_phi_pert, R1_vrt, R1_div,
				dt,		i_simulation_timestamp
			);

		timestepping_lg_exp_ups2.runTimestep(
				R2_phi_pert, R2_vrt, R2_div,
				dt,		i_simulation_timestamp
			);

		timestepping_lg_exp_ups3.runTimestep(
				R3_phi_pert, R3_vrt, R3_div,
				dt,		i_simulation_timestamp
			);

		io_phi_pert = R0_phi_pert + dt*(R1_phi_pert + 2.0*R2_phi_pert + R3_phi_pert);
		io_vrt = R0_vrt + dt*(R1_vrt + 2.0*R2_vrt + R3_vrt);
		io_div = R0_div + dt*(R1_div + 2.0*R2_div + R3_div);
		return;
	}

	SWEETErrorFatal("This time stepping order is not yet supported!");
}




PDESWESphere2DTS_l_exp_direct_special::~PDESWESphere2DTS_l_exp_direct_special()
{
}

