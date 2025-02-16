/*
 * PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv
 *
 * Created on: 24 Mar 2022
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 */

#include "PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv.hpp"


bool PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{

	if (	i_timestepping_method == "lg_exp_na_sl_lc_nr_etdrk_uv"	||
		i_timestepping_method == "lg_exp_na_sl_lc_etdrk_uv"	||
			false
	)
		return true;

	return false;
}


std::string PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv::getIDString()
{
	return "lg_exp_na_sl_lc_nr_etdrk_uv";
}


bool PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	ops = io_ops;

	timestepping_method = i_timestepping_method;

	if (shackPDESWESphere2D->sphere2d_use_fsphere2D)
		SWEETErrorFatal("TODO: Not yet supported");

	NLRemainderTreatment_enum nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR;

	if (shackPDESWETimeDisc->timestepping_method == "lg_exp_na_sl_lc_nr_etdrk_uv")
	{
		nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR;
	}
	else if (shackPDESWETimeDisc->timestepping_method == "lg_exp_na_sl_lc_etdrk_uv")
	{
		nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_IGNORE;
	}
	else
	{
		SWEETErrorFatal(std::string("Timestepping method '")+shackPDESWETimeDisc->timestepping_method+std::string("' not known"));
	}

	return setup(
			ops,
			shackExpIntegration,
			shackPDESWETimeDisc->timestepping_order,
			shackPDESWETimeDisc->timestepping_order2,
			shackTimestepControl->currentTimestepSize,

			nonlinear_remainder_treatment
		);
}



/*
 * Setup
 */
bool PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv::setup(
		const sweet::Data::Sphere2D::Operators *io_ops,
		sweet::ExpIntegration::Shack *i_shackExpIntegration,
		int i_timestepping_order,
		int i_timestepping_order2,
		double i_timestepSize,

		NLRemainderTreatment_enum i_nonlinear_remainder_treatment
)
{
	ops = io_ops;

	timestepping_order = i_timestepping_order;
	timestepping_order2 = i_timestepping_order2;

	nonlinear_remainder_treatment = i_nonlinear_remainder_treatment;

	ts_ln_erk_split_uv.setup_main(ops, i_timestepping_order, true, true, true, true, false);

	// Setup semi-lag
	semiLagrangian.setup(ops->sphere2DDataConfig, shackTimesteppingSemiLagrangianSphere2DData, timestepping_order);

	if (timestepping_order != timestepping_order2)
		SWEETErrorFatal("Mismatch of orders, should be equal");

	if (timestepping_order == 0 || timestepping_order == 1 || timestepping_order == -2)
	{
		ts_phi0_exp.setup_variant_50(ops, i_shackExpIntegration, "phi0", i_timestepSize, false, true, timestepping_order);	/* NO Coriolis */
		ts_psi1_exp.setup_variant_50(ops, i_shackExpIntegration, "psi1", i_timestepSize, false, true, timestepping_order);
	}
	else if (timestepping_order == 2)
	{
		ts_phi0_exp.setup_variant_50(ops, i_shackExpIntegration, "phi0", i_timestepSize, false, true, timestepping_order);	/* NO Coriolis */
		ts_psi1_exp.setup_variant_50(ops, i_shackExpIntegration, "psi1", i_timestepSize, false, true, timestepping_order);
		ts_psi2_exp.setup_variant_50(ops, i_shackExpIntegration, "psi2", i_timestepSize, false, true, timestepping_order);
	}
	else if  (timestepping_order == 4)
	{
		SWEETErrorFatal("4th order method not (yet) supported");
	}
	else
	{
		SWEETErrorFatal("TODO: This order is not implemented, yet!");
	}

	return true;
}

void PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv::printHelp()
{
	std::cout << "	Exponential SL ETD:" << std::endl;
	std::cout << "		+ lg_exp_na_sl_lc_nr_etdrk_uv" << std::endl;
}



void PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_U_phi,
		sweet::Data::Sphere2D::DataSpectral &io_U_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_U_div,

		double i_dt,
		double i_simulation_timestamp
)
{
	const sweet::Data::Sphere2D::Config *sphere2DDataConfig = io_U_phi.sphere2DDataConfig;

	// Keep io data unchanged
	sweet::Data::Sphere2D::DataSpectral U_phi = io_U_phi;
	sweet::Data::Sphere2D::DataSpectral U_vrt = io_U_vrt;
	sweet::Data::Sphere2D::DataSpectral U_div = io_U_div;

	if (i_simulation_timestamp == 0)
	{
#if !SWEET_PARAREAL && !SWEET_XBRAID
		/*
		 * First time step:
		 * Simply backup existing fields for multi-step parts of this algorithm.
		 */
		U_phi_prev = U_phi;
		U_vrt_prev = U_vrt;
		U_div_prev = U_div;
#endif
	}


	/*************************************************************************************************
	 * Step 1) Compute departure points
	 *************************************************************************************************/
	sweet::Data::Sphere2D::DataGrid U_u_lon_prev, U_v_lat_prev;
	ops->vrtdiv_2_uv(U_vrt_prev, U_div_prev, U_u_lon_prev, U_v_lat_prev);

	sweet::Data::Sphere2D::DataGrid U_u_lon, U_v_lat;
	ops->vrtdiv_2_uv(U_vrt, U_div, U_u_lon, U_v_lat);

	double dt_div_radius = i_dt / shackSphere2DDataOps->sphere2d_radius;

	// Calculate departure points
	sweet::Data::Vector::Vector<double> pos_lon_d, pos_lat_d;
	semiLagrangian.semi_lag_departure_points_settls_specialized(
			dt_div_radius*U_u_lon_prev, dt_div_radius*U_v_lat_prev,
			dt_div_radius*U_u_lon, dt_div_radius*U_v_lat,
			pos_lon_d, pos_lat_d		// OUTPUT
		);


	// N(Un) = Nonlinear term applied to Un
	sweet::Data::Sphere2D::DataSpectral FUn_phi(sphere2DDataConfig, 0);
	sweet::Data::Sphere2D::DataSpectral FUn_vrt(sphere2DDataConfig, 0);
	sweet::Data::Sphere2D::DataSpectral FUn_div(sphere2DDataConfig, 0);

	// psi1 applied to N(Un)
	sweet::Data::Sphere2D::DataSpectral psi1_FUn_phi(sphere2DDataConfig, 0);
	sweet::Data::Sphere2D::DataSpectral psi1_FUn_vrt(sphere2DDataConfig, 0);
	sweet::Data::Sphere2D::DataSpectral psi1_FUn_div(sphere2DDataConfig, 0);

	// Interpolated (.5 * psi1NUn + U) to departure points ( )_* //
	sweet::Data::Sphere2D::DataSpectral psi1_FUn_phi_D(sphere2DDataConfig, 0);
	sweet::Data::Sphere2D::DataSpectral psi1_FUn_vrt_D(sphere2DDataConfig, 0);
	sweet::Data::Sphere2D::DataSpectral psi1_FUn_div_D(sphere2DDataConfig, 0);

	// All schemes use the result provided by ETD1RK
	// timestepping_order == -2 --> SL-ETD2RK-bis (double application of SL-ETD1RK; 2nd order scheme)
	if (timestepping_order == 1 || timestepping_order == 2 || timestepping_order == -2)
	{

		////////////////////
		// Compute N(U_n) //
		////////////////////
		ts_ln_erk_split_uv.euler_timestep_update_lc(
				U_phi, U_vrt, U_div,
				FUn_phi, FUn_vrt, FUn_div,
				i_simulation_timestamp
		);

		if (nonlinear_remainder_treatment == NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					U_phi, U_vrt, U_div,
					FUn_phi, FUn_vrt, FUn_div,
					i_simulation_timestamp
			);
		}

		///////////////////////////
		// Apply psi_1 to N(U_n) //
		///////////////////////////
		ts_psi1_exp.runTimestep(
				FUn_phi, FUn_vrt, FUn_div,
				psi1_FUn_phi, psi1_FUn_vrt, psi1_FUn_div,
				i_dt,
				i_simulation_timestamp
			);


		////////////////////////////////////////////////////////
		// Add this to Un and interpolate to departure points //
		////////////////////////////////////////////////////////
		U_phi = U_phi + i_dt * psi1_FUn_phi;
		U_vrt = U_vrt + i_dt * psi1_FUn_vrt;
		U_div = U_div + i_dt * psi1_FUn_div;

		sweet::Data::Sphere2D::DataSpectral U_phi_D, U_vrt_D, U_div_D;
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				U_phi, U_vrt, U_div,
				pos_lon_d, pos_lat_d,
				U_phi_D, U_vrt_D, U_div_D
			);

		//////////////////////////////////////////
		// Apply phi_0 to interpolated solution //
		//////////////////////////////////////////
		sweet::Data::Sphere2D::DataSpectral phi0_Un_phi(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral phi0_Un_vrt(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral phi0_Un_div(sphere2DDataConfig, 0);

		ts_phi0_exp.runTimestep(
				U_phi_D, U_vrt_D, U_div_D,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_dt,
				i_simulation_timestamp
			);

		U_phi = phi0_Un_phi;
		U_vrt = phi0_Un_vrt;
		U_div = phi0_Un_div;

	}
	// SL-ETD2RK-bis (2nd order)
	if (timestepping_order == -2)
	{

		////////////////////////////////////////////////////////
		//Calculate order 1 SL-ETDRK (as above) : {U}_1^{n+1} //
		////////////////////////////////////////////////////////
		// Save SL-ETD1RK from above
		sweet::Data::Sphere2D::DataSpectral U_phi1(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral U_vrt1(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral U_div1(sphere2DDataConfig, 0);
		U_phi1 = U_phi;
		U_vrt1 = U_vrt;
		U_div1 = U_div;

		/////////////////////////////
		// Apply psi_1 to N(Un) //
		/////////////////////////////
		// already computed in the previous if block

		//////////////////////////////////////////////////////////////////////////
		// Add half of this to (original) U and interpolate to departure points //
		//////////////////////////////////////////////////////////////////////////
		U_phi = io_U_phi + .5 * i_dt * psi1_FUn_phi;
		U_vrt = io_U_vrt + .5 * i_dt * psi1_FUn_vrt;
		U_div = io_U_div + .5 * i_dt * psi1_FUn_div;

		/////////////////////////////////////////////////////////////
		//Interpolate (.5 * psi1NUn + U) to departure points ( )_* //
		/////////////////////////////////////////////////////////////
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				U_phi, U_vrt, U_div,
				pos_lon_d, pos_lat_d,
				psi1_FUn_phi_D, psi1_FUn_vrt_D, psi1_FUn_div_D
			);

		////////////////////////
		// Compute N(U_{n+1}) //
		////////////////////////
		sweet::Data::Sphere2D::DataSpectral FU1_phi(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral FU1_vrt(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral FU1_div(sphere2DDataConfig, 0);
		ts_ln_erk_split_uv.euler_timestep_update_lc(
				U_phi1, U_vrt1, U_div1,
				FU1_phi, FU1_vrt, FU1_div,
				i_simulation_timestamp
		);

		if (nonlinear_remainder_treatment == NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					U_phi1, U_vrt1, U_div1,
					FU1_phi, FU1_vrt, FU1_div,
					i_simulation_timestamp
			);
		}

		//////////////////////////////
		// Apply psi1 to N(U_{n+1}) //
		//////////////////////////////
		sweet::Data::Sphere2D::DataSpectral psi1_FU1_phi(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral psi1_FU1_vrt(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral psi1_FU1_div(sphere2DDataConfig, 0);
		ts_psi1_exp.runTimestep(
				FU1_phi, FU1_vrt, FU1_div,
				psi1_FU1_phi, psi1_FU1_vrt, psi1_FU1_div,
				i_dt,
				i_simulation_timestamp
			);

		/////////////////////////////////////////////////////
		// Add half of this to already interpolated values //
		/////////////////////////////////////////////////////
		U_phi = psi1_FUn_phi_D + .5 * i_dt * psi1_FU1_phi;
		U_vrt = psi1_FUn_vrt_D + .5 * i_dt * psi1_FU1_vrt;
		U_div = psi1_FUn_div_D + .5 * i_dt * psi1_FU1_div;

		////////////////////////////////////////////////////////////////////////
		// Calculate phi_0 of [interpolated (.5 * psi1NUN + U) + .5 * ps1NU1] //
		////////////////////////////////////////////////////////////////////////
		sweet::Data::Sphere2D::DataSpectral phi0_Un_phi(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral phi0_Un_vrt(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral phi0_Un_div(sphere2DDataConfig, 0);
		ts_phi0_exp.runTimestep(
				U_phi, U_vrt, U_div,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_dt,
				i_simulation_timestamp
		);

		U_phi = phi0_Un_phi;
		U_vrt = phi0_Un_vrt;
		U_div = phi0_Un_div;

	}
	else if (timestepping_order == 2)
	{

		////////////////////////////////////////////////////////
		//Calculate order 1 SL-ETDRK (as above) : {U}_1^{n+1} //
		////////////////////////////////////////////////////////
		// Save SL-ETD1RK from above
		sweet::Data::Sphere2D::DataSpectral U_phi1(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral U_vrt1(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral U_div1(sphere2DDataConfig, 0);
		U_phi1 = U_phi;
		U_vrt1 = U_vrt;
		U_div1 = U_div;

		////////////////////////
		// Compute N(U_{n+1}) //
		////////////////////////
		sweet::Data::Sphere2D::DataSpectral FU1_phi(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral FU1_vrt(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral FU1_div(sphere2DDataConfig, 0);
		ts_ln_erk_split_uv.euler_timestep_update_lc(
				U_phi1, U_vrt1, U_div1,
				FU1_phi, FU1_vrt, FU1_div,
				i_simulation_timestamp
		);

		if (nonlinear_remainder_treatment == NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR)
		{
			ts_ln_erk_split_uv.euler_timestep_update_nr(
					U_phi1, U_vrt1, U_div1,
					FU1_phi, FU1_vrt, FU1_div,
					i_simulation_timestamp
			);
		}

		//////////////////////////////
		// Apply psi2 to N(U_{n+1}) //
		//////////////////////////////
		sweet::Data::Sphere2D::DataSpectral psi2_FU1_phi(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral psi2_FU1_vrt(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral psi2_FU1_div(sphere2DDataConfig, 0);
		ts_psi2_exp.runTimestep(
				FU1_phi, FU1_vrt, FU1_div,
				psi2_FU1_phi, psi2_FU1_vrt, psi2_FU1_div,
				i_dt,
				i_simulation_timestamp
		);

		//////////////////////////////
		// Apply psi2 to N(U_{n}) //
		//////////////////////////////
		sweet::Data::Sphere2D::DataSpectral psi2_FUn_phi(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral psi2_FUn_vrt(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral psi2_FUn_div(sphere2DDataConfig, 0);
		ts_psi2_exp.runTimestep(
				FUn_phi, FUn_vrt, FUn_div,
				psi2_FUn_phi, psi2_FUn_vrt, psi2_FUn_div,
				i_dt,
				i_simulation_timestamp
		);

		////////////////////////////////////////////////////
		// Interpolate psi2NUn to departure points ( )_* //
		////////////////////////////////////////////////////
		sweet::Data::Sphere2D::DataSpectral psi2_FUn_phi_D(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral psi2_FUn_vrt_D(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral psi2_FUn_div_D(sphere2DDataConfig, 0);
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				psi2_FUn_phi, psi2_FUn_vrt, psi2_FUn_div,
				pos_lon_d, pos_lat_d,
				psi2_FUn_phi_D, psi2_FUn_vrt_D, psi2_FUn_div_D
			);

		////////////////////////////////////////
		// Apply phi0 to psi2NU1 - (psi2NUn)_* //
		////////////////////////////////////////
		sweet::Data::Sphere2D::DataSpectral phi0_dif2_phi(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral phi0_dif2_vrt(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral phi0_dif2_div(sphere2DDataConfig, 0);
		ts_phi0_exp.runTimestep(
				psi2_FU1_phi - psi2_FUn_phi_D, psi2_FU1_vrt - psi2_FUn_vrt_D, psi2_FU1_div - psi2_FUn_div_D,
				phi0_dif2_phi, phi0_dif2_vrt, phi0_dif2_div,
				i_dt,
				i_simulation_timestamp
		);

		U_phi = U_phi1 + i_dt*phi0_dif2_phi;
		U_vrt = U_vrt1 + i_dt*phi0_dif2_vrt;
		U_div = U_div1 + i_dt*phi0_dif2_div;

	}

	// Save current time step for next step
	U_phi_prev = io_U_phi;
	U_vrt_prev = io_U_vrt;
	U_div_prev = io_U_div;

	io_U_phi = U_phi;
	io_U_vrt = U_vrt;
	io_U_div = U_div;

}



PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv::PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv()
{
}



PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv::~PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv()
{
}

