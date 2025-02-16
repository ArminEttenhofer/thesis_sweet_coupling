/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etd_uv.hpp"


bool PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etd_uv::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;

	if (	i_timestepping_method == "lg_exp_na_sl_lc_nr_etd_uv"	||
			i_timestepping_method == "lg_exp_na_sl_lc_etd_uv"	||
			false
	)
		return true;

	return false;
}



bool PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etd_uv::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	if (shackPDESWESphere2D->sphere2d_use_fsphere2D)
		SWEETErrorFatal("TODO: Not yet supported");

	NLRemainderTreatment_enum _nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR;

	if (timestepping_method == "lg_exp_na_sl_lc_nr_etd_uv")
	{
		_nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR;
	}
	else if (timestepping_method == "lg_exp_na_sl_lc_etd_uv")
	{
		_nonlinear_remainder_treatment = NLRemainderTreatment_enum::NL_REMAINDER_IGNORE;
	}
	else
	{
		SWEETErrorFatal(std::string("Timestepping method '")+timestepping_method+std::string("' not known"));
	}

	return setup(
			io_ops,
			shackExpIntegration,
			timestepping_order,
			timestepping_order2,
			shackTimestepControl->currentTimestepSize,

			_nonlinear_remainder_treatment
		);
}


bool PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etd_uv::setup(
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


std::string PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etd_uv::getIDString()
{
	return "lg_exp_na_sl_lc_nr_etd_uv";
}

void PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etd_uv::printHelp()
{
	std::cout << "	Exponential SL ETD:" << std::endl;
	//std::cout << "		+ lg_exp_na_sl_etd_uv" << std::endl;
	//std::cout << "		+ lg_exp_na_sl_nr_etd_uv" << std::endl;
	//std::cout << "		+ lg_exp_na_sl_lc_etd_uv" << std::endl;
	std::cout << "		+ lg_exp_na_sl_lc_nr_etd_uv" << std::endl;
#if 0
	std::cout << "		+ l_exp_na_sl_etd_uv" << std::endl;
	std::cout << "		+ l_exp_na_sl_nr_etd_uv" << std::endl;
#endif
}



void PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etd_uv::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_U_phi,
		sweet::Data::Sphere2D::DataSpectral &io_U_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_U_div,

		double i_dt,
		double i_simulation_timestamp
)
{
	const sweet::Data::Sphere2D::Config *sphere2DDataConfig = io_U_phi.sphere2DDataConfig;


	if (timestepping_order != 2)
	{
		SWEETErrorFatal("Only 2nd order supported");
	}


	const sweet::Data::Sphere2D::DataSpectral &U_phi = io_U_phi;
	const sweet::Data::Sphere2D::DataSpectral &U_vrt = io_U_vrt;
	const sweet::Data::Sphere2D::DataSpectral &U_div = io_U_div;

	if (i_simulation_timestamp == 0)
	{
#if !SWEET_PARAREAL
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


	/*************************************************************************************************
	 * Step 2) Compute U_D at departure points
	 *************************************************************************************************/

	sweet::Data::Sphere2D::DataSpectral U_phi_D, U_vrt_D, U_div_D;
	semiLagrangian.apply_sl_timeintegration_uv(
			ops,
			U_phi, U_vrt, U_div,
			pos_lon_d, pos_lat_d,
			U_phi_D, U_vrt_D, U_div_D
		);



	/*
	 * u^{n+1} = \varphi_{0} (\Delta t L )u^{n}_D + \Delta t (K + K_D)
	 *
	 * K = \left[
	 * 		\varphi_{1}(\Delta t L)F(u^{n}) +
	 * 		\varphi_{2}(\Delta t L) (F(u^{n})-F(u^{n-1}))
	 * 	\right]
	 */


	/*
	 * phi0
	 */
	sweet::Data::Sphere2D::DataSpectral phi0_phi(sphere2DDataConfig);
	sweet::Data::Sphere2D::DataSpectral phi0_vrt(sphere2DDataConfig);
	sweet::Data::Sphere2D::DataSpectral phi0_div(sphere2DDataConfig);

	ts_phi0_exp.runTimestep(
			//U_phi, U_vrt, U_div,
			U_phi_D, U_vrt_D, U_div_D,
			phi0_phi, phi0_vrt, phi0_div,
			i_dt,
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

	if (nonlinear_remainder_treatment == NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR)
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
			i_dt,
			i_simulation_timestamp
		);


	/*
	 * phi2
	 */
	sweet::Data::Sphere2D::DataSpectral F_phi_prev(sphere2DDataConfig, 0);
	sweet::Data::Sphere2D::DataSpectral F_vrt_prev(sphere2DDataConfig, 0);
	sweet::Data::Sphere2D::DataSpectral F_div_prev(sphere2DDataConfig, 0);

	ts_ln_erk_split_uv.euler_timestep_update_lc(
			U_phi_prev, U_vrt_prev, U_div_prev,
			F_phi_prev, F_vrt_prev, F_div_prev,
			i_simulation_timestamp
	);



	if (nonlinear_remainder_treatment == NLRemainderTreatment_enum::NL_REMAINDER_NONLINEAR)
	{
		ts_ln_erk_split_uv.euler_timestep_update_nr(
				U_phi_prev, U_vrt_prev, U_div_prev,
				F_phi_prev, F_vrt_prev, F_div_prev,
				i_simulation_timestamp
		);
	}

	sweet::Data::Sphere2D::DataSpectral phi2_phi(sphere2DDataConfig);
	sweet::Data::Sphere2D::DataSpectral phi2_vrt(sphere2DDataConfig);
	sweet::Data::Sphere2D::DataSpectral phi2_div(sphere2DDataConfig);

	ts_phi2_exp.runTimestep(
			F_phi - F_phi_prev,
			F_vrt - F_vrt_prev,
			F_div - F_div_prev,

			phi2_phi,
			phi2_vrt,
			phi2_div,

			i_dt,
			i_simulation_timestamp
		);

	sweet::Data::Sphere2D::DataSpectral K_phi = phi1_phi + phi2_phi;
	sweet::Data::Sphere2D::DataSpectral K_vrt = phi1_vrt + phi2_vrt;
	sweet::Data::Sphere2D::DataSpectral K_div = phi1_div + phi2_div;


	/*
	 * Compute K at departure points
	 */

	sweet::Data::Sphere2D::DataSpectral K_phi_D, K_vrt_D, K_div_D;
	semiLagrangian.apply_sl_timeintegration_uv(
			ops,
			K_phi, K_vrt, K_div,
			pos_lon_d, pos_lat_d,
			K_phi_D, K_vrt_D, K_div_D
		);

	io_U_phi = phi0_phi + (i_dt*0.5)*(K_phi + K_phi_D);
	io_U_vrt = phi0_vrt + (i_dt*0.5)*(K_vrt + K_vrt_D);
	io_U_div = phi0_div + (i_dt*0.5)*(K_div + K_div_D);


	{
		U_phi_prev = U_phi;
		U_vrt_prev = U_vrt;
		U_div_prev = U_div;
	}
}




PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etd_uv::PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etd_uv()
{
}


PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etd_uv::~PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etd_uv()
{
}

