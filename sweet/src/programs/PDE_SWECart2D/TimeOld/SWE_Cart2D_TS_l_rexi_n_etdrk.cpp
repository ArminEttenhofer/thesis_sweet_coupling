/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_cart2d.cpp
 *					which was also written by Pedro Peixoto
 */

#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_rexi_n_etdrk.hpp>




/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Cart2D_TS_l_rexi_n_etdrk::euler_timestep_update_nonlinear(
		const sweet::Data::Cart2D::DataSpectral &i_h,	//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_u,	//!< prognostic variables
		const sweet::Data::Cart2D::DataSpectral &i_v,	//!< prognostic variables

		sweet::Data::Cart2D::DataSpectral &o_h_t,	//!< time updates
		sweet::Data::Cart2D::DataSpectral &o_u_t,	//!< time updates
		sweet::Data::Cart2D::DataSpectral &o_v_t,	//!< time updates

		double i_timestamp
)
{
	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */


	o_u_t = -i_u*ops->diff_c_x(i_u) - i_v*ops->diff_c_y(i_u);
	o_v_t = -i_u*ops->diff_c_x(i_v) - i_v*ops->diff_c_y(i_v);

	if (use_only_linear_divergence)
	{
		//only nonlinear advection left to solve
		o_h_t = - (i_u*ops->diff_c_x(i_h) + i_v*ops->diff_c_y(i_h));
	}
	else //full nonlinear equation on h
	{
		if (shackPDESWECart2D->use_nonlinear_only_visc != 0)
		{
			//solve nonlinear divergence
			o_h_t = - (i_h*ops->diff_c_x(i_u) + i_h*ops->diff_c_y(i_v));
			//filter
#if !SWEET_USE_CART2D_SPECTRAL_SPACE
			SWEETErrorFatal("Implicit diffusion only supported with spectral space activated");
#else
			o_h_t = ops->implicit_diffusion(o_h_t, shackTimestepControl->currentTimestepSize*shackPDESWECart2D->viscosity, shackPDESWECart2D->viscosity_order);
#endif
			//add nonlinear advection
			o_h_t = o_h_t - (i_u*ops->diff_c_x(i_h) + i_v*ops->diff_c_y(i_h));
		}
		else
		{
			o_h_t = -ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);
		}
	}

}



void SWE_Cart2D_TS_l_rexi_n_etdrk::runTimestep(
		sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETErrorFatal("SWE_Cart2D_TS_l_phi0_n_edt: Only constant time step size allowed");


	const sweet::Data::Cart2D::Config *cart2DDataConfig = io_h.cart2DDataConfig;

	if (timestepping_order == 1)
	{
		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}
		 * 			+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 */

		sweet::Data::Cart2D::DataSpectral phi0_Un_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi0_Un_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi0_Un_v(cart2DDataConfig);
		ts_phi0_rexi.runTimestep(
				io_h, io_u, io_v,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_dt,
				i_simulation_timestamp
			);

		sweet::Data::Cart2D::DataSpectral FUn_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral FUn_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral FUn_v(cart2DDataConfig);
		euler_timestep_update_nonlinear(
				io_h, io_u, io_v,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);

		sweet::Data::Cart2D::DataSpectral phi1_FUn_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi1_FUn_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi1_FUn_v(cart2DDataConfig);

		ts_phi1_rexi.runTimestep(
				FUn_h, FUn_u, FUn_v,
				phi1_FUn_h, phi1_FUn_u, phi1_FUn_v,
				i_dt,
				i_simulation_timestamp
			);

		io_h = phi0_Un_h + i_dt*phi1_FUn_h;
		io_u = phi0_Un_u + i_dt*phi1_FUn_u;
		io_v = phi0_Un_v + i_dt*phi1_FUn_v;
	}
	else if (timestepping_order == 2)
	{

		/*
		 * A_{n}=\psi_{0}(\Delta tL)U_{n}+\Delta t\psi_{1}(\Delta tL)F(U_{n})
		 */

		sweet::Data::Cart2D::DataSpectral phi0_Un_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi0_Un_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi0_Un_v(cart2DDataConfig);

		ts_phi0_rexi.runTimestep(
				io_h, io_u, io_v,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_dt,
				i_simulation_timestamp
			);

		sweet::Data::Cart2D::DataSpectral FUn_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral FUn_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral FUn_v(cart2DDataConfig);

		euler_timestep_update_nonlinear(
				io_h, io_u, io_v,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);

		sweet::Data::Cart2D::DataSpectral phi1_FUn_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi1_FUn_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi1_FUn_v(cart2DDataConfig);

		ts_phi1_rexi.runTimestep(
				FUn_h, FUn_u, FUn_v,
				phi1_FUn_h, phi1_FUn_u, phi1_FUn_v,
				i_dt,
				i_simulation_timestamp
			);

		sweet::Data::Cart2D::DataSpectral A_h = phi0_Un_h + i_dt*phi1_FUn_h;
		sweet::Data::Cart2D::DataSpectral A_u = phi0_Un_u + i_dt*phi1_FUn_u;
		sweet::Data::Cart2D::DataSpectral A_v = phi0_Un_v + i_dt*phi1_FUn_v;

		/*
		 * U_{n+1} = A_{n}+ \Delta t \psi_{2}(\Delta tL)
		 * 				\left(F(A_{n},t_{n}+\Delta t)-F(U_{n})\right)
		 */

		sweet::Data::Cart2D::DataSpectral FAn_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral FAn_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral FAn_v(cart2DDataConfig);

		euler_timestep_update_nonlinear(
				A_h, A_u, A_v,
				FAn_h, FAn_u, FAn_v,
				i_simulation_timestamp
		);


		sweet::Data::Cart2D::DataSpectral phi2_X_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi2_X_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi2_X_v(cart2DDataConfig);

		ts_phi2_rexi.runTimestep(
				FAn_h - FUn_h,
				FAn_u - FUn_u,
				FAn_v - FUn_v,

				phi2_X_h,
				phi2_X_u,
				phi2_X_v,

				i_dt,
				i_simulation_timestamp
			);

		io_h = A_h + i_dt*phi2_X_h;
		io_u = A_u + i_dt*phi2_X_u;
		io_v = A_v + i_dt*phi2_X_v;
	}
	else if (timestepping_order == 4)
	{
		double dt = i_dt;
		double dt_half = dt*0.5;



		/*
		 * Precompute commonly used terms
		 */
		sweet::Data::Cart2D::DataSpectral phi0_Un_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi0_Un_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi0_Un_v(cart2DDataConfig);

		ts_phi0_rexi.runTimestep(
				io_h, io_u, io_v,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				dt_half,
				i_simulation_timestamp
			);

		sweet::Data::Cart2D::DataSpectral FUn_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral FUn_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral FUn_v(cart2DDataConfig);

		euler_timestep_update_nonlinear(
				io_h, io_u, io_v,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);



		/*
		 * Some commonly shared buffers
		 */

		sweet::Data::Cart2D::DataSpectral phi1_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi1_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi1_v(cart2DDataConfig);


		/*
		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + \Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
		 */
		ts_phi1_rexi.runTimestep(
				FUn_h, FUn_u, FUn_v,
				phi1_h, phi1_u, phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		sweet::Data::Cart2D::DataSpectral A_h = phi0_Un_h + dt_half*phi1_h;
		sweet::Data::Cart2D::DataSpectral A_u = phi0_Un_u + dt_half*phi1_u;
		sweet::Data::Cart2D::DataSpectral A_v = phi0_Un_v + dt_half*phi1_v;



		/*
		 * B_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(A_{n}, t_{n} + 0.5*\Delta t)
		 */

		sweet::Data::Cart2D::DataSpectral FAn_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral FAn_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral FAn_v(cart2DDataConfig);

		euler_timestep_update_nonlinear(
				A_h, A_u, A_v,
				FAn_h, FAn_u, FAn_v,
				i_simulation_timestamp + dt_half
		);

		ts_phi1_rexi.runTimestep(
				FAn_h, FAn_u, FAn_v,
				phi1_h, phi1_u, phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		sweet::Data::Cart2D::DataSpectral B_h = phi0_Un_h + dt_half*phi1_h;
		sweet::Data::Cart2D::DataSpectral B_u = phi0_Un_u + dt_half*phi1_u;
		sweet::Data::Cart2D::DataSpectral B_v = phi0_Un_v + dt_half*phi1_v;



		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 */

		sweet::Data::Cart2D::DataSpectral phi0_An_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi0_An_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral phi0_An_v(cart2DDataConfig);

		ts_phi0_rexi.runTimestep(
				A_h, A_u, A_v,
				phi0_An_h, phi0_An_u, phi0_An_v,
				dt_half,
				i_simulation_timestamp
			);


		sweet::Data::Cart2D::DataSpectral FBn_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral FBn_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral FBn_v(cart2DDataConfig);

		euler_timestep_update_nonlinear(
				B_h, B_u, B_v,
				FBn_h, FBn_u, FBn_v,
				i_simulation_timestamp + dt_half
		);

		ts_phi1_rexi.runTimestep(
				2.0*FBn_h - FUn_h,
				2.0*FBn_u - FUn_u,
				2.0*FBn_v - FUn_v,
				phi1_h,	phi1_u,	phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		sweet::Data::Cart2D::DataSpectral C_h = phi0_An_h + dt_half*phi1_h;
		sweet::Data::Cart2D::DataSpectral C_u = phi0_An_u + dt_half*phi1_u;
		sweet::Data::Cart2D::DataSpectral C_v = phi0_An_v + dt_half*phi1_v;



		/*
		 * R0 - R3
		 */
		sweet::Data::Cart2D::DataSpectral FCn_h(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral FCn_u(cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral FCn_v(cart2DDataConfig);

		euler_timestep_update_nonlinear(
				C_h, C_u, C_v,
				FCn_h, FCn_u, FCn_v,
				i_simulation_timestamp + dt
		);

		sweet::Data::Cart2D::DataSpectral R0_h = io_h;
		sweet::Data::Cart2D::DataSpectral R0_u = io_u;
		sweet::Data::Cart2D::DataSpectral R0_v = io_v;

		sweet::Data::Cart2D::DataSpectral &R1_h = FUn_h;
		sweet::Data::Cart2D::DataSpectral &R1_u = FUn_u;
		sweet::Data::Cart2D::DataSpectral &R1_v = FUn_v;

		sweet::Data::Cart2D::DataSpectral R2_h = FAn_h + FBn_h;
		sweet::Data::Cart2D::DataSpectral R2_u = FAn_u + FBn_u;
		sweet::Data::Cart2D::DataSpectral R2_v = FAn_v + FBn_v;

		sweet::Data::Cart2D::DataSpectral &R3_h = FCn_h;
		sweet::Data::Cart2D::DataSpectral &R3_u = FCn_u;
		sweet::Data::Cart2D::DataSpectral &R3_v = FCn_v;


		/*
		 * U_{n+1} =
		 * 		\psi_{0}(\Delta tL)R_{0}
		 * 			+ \Delta t
		 * 			(
		 * 				  \upsilon_{1}(\Delta tL) R_{1} +
		 * 				2*\upsilon_{2}(\Delta tL) R_{2} +
		 * 				  \upsilon_{3}(\Delta tL) R_{3}
		 * 			)
		 */
		ts_phi0_rexi.runTimestep(
				R0_h, R0_u, R0_v,
				dt,		i_simulation_timestamp
			);

		ts_ups1_rexi.runTimestep(
				R1_h, R1_u, R1_v,
				dt,		i_simulation_timestamp
			);

		ts_ups2_rexi.runTimestep(
				R2_h, R2_u, R2_v,
				dt,		i_simulation_timestamp
			);

		ts_ups3_rexi.runTimestep(
				R3_h, R3_u, R3_v,
				dt,		i_simulation_timestamp
			);

		io_h = R0_h + dt*(R1_h + 2.0*R2_h + R3_h);
		io_u = R0_u + dt*(R1_u + 2.0*R2_u + R3_u);
		io_v = R0_v + dt*(R1_v + 2.0*R2_v + R3_v);
	}
	else
	{
		SWEETErrorFatal("TODO: This order is not implemented, yet!");
	}
}


bool SWE_Cart2D_TS_l_rexi_n_etdrk::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	PDESWECart2DTS_BaseInterface::shackRegistration(io_shackDict);

	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

	ts_phi0_rexi.shackRegistration(io_shackDict);
	ts_phi1_rexi.shackRegistration(io_shackDict);
	ts_phi2_rexi.shackRegistration(io_shackDict);

	ts_ups0_rexi.shackRegistration(io_shackDict);
	ts_ups1_rexi.shackRegistration(io_shackDict);
	ts_ups2_rexi.shackRegistration(io_shackDict);
	ts_ups3_rexi.shackRegistration(io_shackDict);

	return true;
}


bool SWE_Cart2D_TS_l_rexi_n_etdrk::setup(
		sweet::Data::Cart2D::Operators *io_ops
)
{
	PDESWECart2DTS_BaseInterface::setup(io_ops);

	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	use_only_linear_divergence = shackPDESWECart2D->use_only_linear_divergence;

	if (timestepping_order == 1)
	{
		ts_phi0_rexi.setup(ops, "phi0");
		ts_phi1_rexi.setup(ops, "phi1");
	}
	else if (timestepping_order == 2)
	{
		ts_phi0_rexi.setup(ops, "phi0");
		ts_phi1_rexi.setup(ops, "phi1");
		ts_phi2_rexi.setup(ops, "phi2");
	}
	else if (timestepping_order == 4)
	{
		ts_phi0_rexi.setup(ops, "phi0");
		ts_phi1_rexi.setup(ops, "phi1");
		ts_phi2_rexi.setup(ops, "phi2");

		ts_ups0_rexi.setup(ops, "phi0");
		ts_ups1_rexi.setup(ops, "ups1");
		ts_ups2_rexi.setup(ops, "ups2");
		ts_ups3_rexi.setup(ops, "ups3");
	}

	return true;
}



SWE_Cart2D_TS_l_rexi_n_etdrk::~SWE_Cart2D_TS_l_rexi_n_etdrk()
{
}

