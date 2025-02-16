/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2019-10-24: Partly based on cart2d version
 */

#include "PDESWESphere2DTS_ln_settls_uv.hpp"




bool PDESWESphere2DTS_ln_settls_uv::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	PDESWESphere2DTS_ln_settls_uv::LinearCoriolisTreatment_enum _linear_coriolis_treatment = PDESWESphere2DTS_ln_settls_uv::CORIOLIS_IGNORE;
	PDESWESphere2DTS_ln_settls_uv::NLRemainderTreatment_enum _nonlinear_divergence_treatment = PDESWESphere2DTS_ln_settls_uv::NL_REMAINDER_IGNORE;
	bool _original_linear_operator_sl_treatment = true;

	if (timestepping_method == "ln_settls")
	{
		_linear_coriolis_treatment = PDESWESphere2DTS_ln_settls_uv::CORIOLIS_LINEAR;
		_nonlinear_divergence_treatment = PDESWESphere2DTS_ln_settls_uv::NL_REMAINDER_NONLINEAR;
		_original_linear_operator_sl_treatment = true;
	}
	else
	{
		// Search for Coriolis
		if (timestepping_method.find("l_irk") != std::string::npos || timestepping_method.find("l_exp") != std::string::npos)
			_linear_coriolis_treatment = PDESWESphere2DTS_ln_settls_uv::CORIOLIS_LINEAR;
		else if (timestepping_method.find("lc_na_sl") != std::string::npos)
			_linear_coriolis_treatment = PDESWESphere2DTS_ln_settls_uv::CORIOLIS_SEMILAGRANGIAN;
		else if (timestepping_method.find("lc_") != std::string::npos)
			_linear_coriolis_treatment = PDESWESphere2DTS_ln_settls_uv::CORIOLIS_NONLINEAR;

		// Search for Nonlinear divergence
		if (timestepping_method.find("_nr_") != std::string::npos)
			_nonlinear_divergence_treatment = PDESWESphere2DTS_ln_settls_uv::NL_REMAINDER_NONLINEAR;

		if (timestepping_method.find("_ver1") != std::string::npos)
			_original_linear_operator_sl_treatment = false;
		else if (timestepping_method.find("_ver0") != std::string::npos)
			_original_linear_operator_sl_treatment = true;
		else
			_original_linear_operator_sl_treatment = true;

#if 1
		string_id_storage = "";

		if (_linear_coriolis_treatment == PDESWESphere2DTS_ln_settls_uv::CORIOLIS_LINEAR)
			string_id_storage += "l";
		else
			string_id_storage += "lg";

		string_id_storage += "_irk";

		if (_linear_coriolis_treatment == PDESWESphere2DTS_ln_settls_uv::CORIOLIS_SEMILAGRANGIAN)
			string_id_storage += "_lc";

		string_id_storage += "_na";

		string_id_storage += "_sl";

		if (_linear_coriolis_treatment == PDESWESphere2DTS_ln_settls_uv::CORIOLIS_NONLINEAR)
			string_id_storage += "_lc";

		if (_nonlinear_divergence_treatment == PDESWESphere2DTS_ln_settls_uv::NL_REMAINDER_NONLINEAR)
			string_id_storage += "_nr";

		string_id_storage += "_settls";

		if (!_original_linear_operator_sl_treatment)
			string_id_storage += "_ver1";

		std::string string_id_storage_ = string_id_storage + "_uv";

		if (timestepping_method == string_id_storage_)
		{
			string_id_storage = string_id_storage_;
		}
		else
		{
			if (!_original_linear_operator_sl_treatment)
			{
				// there must be a ver1 which is likely missing
				std::cerr << "Detected time stepping method: "+string_id_storage_ << std::endl;
				std::cerr << "Provided time stepping method: "+timestepping_method << std::endl;
				SWEETErrorFatal("Autodetection of parts of time stepping methods failed!");
			}

			std::string string_id_storage2 = string_id_storage+"_ver0"+"_uv";
			if (timestepping_method == string_id_storage2)
			{
				string_id_storage = string_id_storage2;
			}
			else
			{
				std::cerr << "Detected time stepping method: "+string_id_storage_ << std::endl;
				std::cerr << "Provided time stepping method: "+timestepping_method << std::endl;
				std::cerr << "Detected alternative time stepping method: "+string_id_storage2 << std::endl;
				SWEETErrorFatal("Autodetection of parts of time stepping methods failed!");
			}
		}
	}
#endif

	return setup_main(
			io_ops,
			shackPDESWETimeDisc->timestepping_order,
			_linear_coriolis_treatment,				// Coriolis treatment
			_nonlinear_divergence_treatment,			// Nonlinear divergence treatment
			_original_linear_operator_sl_treatment	// original SL linear operator treatment
		);
}




bool PDESWESphere2DTS_ln_settls_uv::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_timestepping_order,
		LinearCoriolisTreatment_enum i_coriolis_treatment,
		NLRemainderTreatment_enum i_nonlinear_divergence_treatment,
		bool i_original_linear_operator_sl_treatment
)
{
	ops = io_ops;

	coriolis_treatment = i_coriolis_treatment;
	nonlinear_remainder_treatment = i_nonlinear_divergence_treatment;
	timestepping_order = i_timestepping_order;
	original_linear_operator_sl_treatment = i_original_linear_operator_sl_treatment;


	if (timestepping_order != 1 && timestepping_order != 2)
		SWEETErrorFatal("Invalid time stepping order, must be 1 or 2");


	// Setup sampler for future interpolations
	//sphere2DSampler.setup(ops->sphere2DDataConfig);

	// Setup semi-Lagrangian
	semiLagrangian.setup(ops->sphere2DDataConfig, shackTimesteppingSemiLagrangianSphere2DData, timestepping_order);

	swe_sphere2d_ts_ln_erk_split_uv.setup_main(ops, 1, true, true, true, true, false);


	if (coriolis_treatment == CORIOLIS_SEMILAGRANGIAN)
	{
		setupFG();
		coriolis_arrival_spectral = fg;
	}

	if (coriolis_treatment == CORIOLIS_LINEAR)
	{
		if (timestepping_order == 1)
		{
			swe_sphere2d_ts_l_irk.setup(
					ops,
					1,
					shackTimestepControl->currentTimestepSize
				);
		}
		else
		{
			// initialize with 1st order and half time step size
			swe_sphere2d_ts_l_irk.setup(
					ops,
					1,
					0.5 * shackTimestepControl->currentTimestepSize
				);
		}
	}
	else
	{
		if (timestepping_order == 1)
		{
			swe_sphere2d_ts_lg_irk.setup(
					ops,
					1,
					shackTimestepControl->currentTimestepSize
				);
		}
		else
		{
			// initialize with 1st order and half time step size
			swe_sphere2d_ts_lg_irk.setup(
					ops,
					1,
					0.5 * shackTimestepControl->currentTimestepSize
				);
		}
	}

	return true;
}



bool PDESWESphere2DTS_ln_settls_uv::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	/*
	 * Should contain _exp and _settls
	 */
	timestepping_method = i_timestepping_method;
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	return (
		!(i_timestepping_method.find("_exp") != std::string::npos)		&&
		(i_timestepping_method.find("_settls") != std::string::npos)	&&
		(i_timestepping_method.find("_uv") != std::string::npos)		&&
		!(i_timestepping_method.find("_only") != std::string::npos)		&&
		true
	);
}


std::string PDESWESphere2DTS_ln_settls_uv::getIDString()
{
	return string_id_storage;
}


void PDESWESphere2DTS_ln_settls_uv::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 2)
	{
		run_timestep_2nd_order(io_phi_pert, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
		return;
	}
}



void PDESWESphere2DTS_ln_settls_uv::run_timestep_2nd_order(
	sweet::Data::Sphere2D::DataSpectral &io_U_phi_pert,
	sweet::Data::Sphere2D::DataSpectral &io_U_vrt,
	sweet::Data::Sphere2D::DataSpectral &io_U_div,

	double i_dt,		
	double i_simulation_timestamp
)
{
	const sweet::Data::Sphere2D::DataSpectral &U_phi = io_U_phi_pert;
	const sweet::Data::Sphere2D::DataSpectral &U_vrt = io_U_vrt;
	const sweet::Data::Sphere2D::DataSpectral &U_div = io_U_div;

	if (i_simulation_timestamp == 0)
	{
		/*
		 * First time step:
		 * Simply backup existing fields for multi-step parts of this algorithm.
		 */
		U_phi_prev = U_phi;
		U_vrt_prev = U_vrt;
		U_div_prev = U_div;
	}


	/*
	 * Step 1) SL
	 * Compute Lagrangian trajectories based on SETTLS.
	 * This depends on V(t-\Delta t) and V(t).
	 *
	 * See Hortal's paper for equation.
	 */
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


	/*
	 * Step 2) Midpoint rule
	 * Put everything together with midpoint rule and solve resulting Helmholtz problem
	 */

	/*
	 * Step 2a) Compute RHS
	 * R = X_D + 1/2 dt L_D + dt N*
	 */

	/*
	 * Compute X_D
	 */
	sweet::Data::Sphere2D::DataSpectral U_phi_D, U_vrt_D, U_div_D;
	semiLagrangian.apply_sl_timeintegration_uv(
			ops,
			U_phi, U_vrt, U_div,
			pos_lon_d, pos_lat_d,
			U_phi_D, U_vrt_D, U_div_D
		);


	sweet::Data::Sphere2D::DataSpectral coriolis_departure_spectral;
	if (coriolis_treatment == CORIOLIS_SEMILAGRANGIAN)
	{
		/*
		 * Compute Coriolis effect at the departure points
		 */

		SWEET_ASSERT(!shackPDESWESphere2D->sphere2d_use_fsphere2D);

		sweet::Data::Sphere2D::DataGrid coriolis_depature_lat = sweet::Data::Vector::Convert::Vector_2_Sphere2DDataGrid::convert(pos_lon_d, ops->sphere2DDataConfig);

		sweet::Data::Sphere2D::DataGrid coriolis_departure_physical(ops->sphere2DDataConfig);
		coriolis_departure_physical.grid_update_lambda_array_idx(
			[&](int i_index, double &o_data)
			{
				o_data = 2.0*shackPDESWESphere2D->sphere2d_rotating_coriolis_omega*std::sin(coriolis_depature_lat.grid_space_data[i_index]);
			}
		);

		coriolis_departure_spectral = coriolis_departure_physical;
	}

	/*
	 * Compute L_D
	 */
	sweet::Data::Sphere2D::DataSpectral L_U_phi_D, L_U_vrt_D, L_U_div_D;
	if (original_linear_operator_sl_treatment)
	{
		/*
		 * Method 1) First evaluate L, then sample result at departure point
		 */
		const sweet::Data::Sphere2D::Config *sphere2DDataConfig = io_U_phi_pert.sphere2DDataConfig;
		sweet::Data::Sphere2D::DataSpectral L_U_phi(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral L_U_vrt(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral L_U_div(sphere2DDataConfig, 0);


		if (coriolis_treatment == CORIOLIS_SEMILAGRANGIAN)
		{
			/*
			 * L_g(U): Linear gravity modes
			 */
			swe_sphere2d_ts_ln_erk_split_uv.euler_timestep_update_lg(
					U_phi, U_vrt - coriolis_arrival_spectral, U_div,		// SL treatment of Coriolis!!!
					L_U_phi, L_U_vrt, L_U_div,
					i_simulation_timestamp
				);
		}
		else if (coriolis_treatment == CORIOLIS_LINEAR)
		{
		
			/*
			 * L_g(U): Linear gravity modes
			 */
			swe_sphere2d_ts_ln_erk_split_uv.euler_timestep_update_lg(
					U_phi, U_vrt, U_div,
					L_U_phi, L_U_vrt, L_U_div,
					i_simulation_timestamp
				);

			/*
			 * L_c(U): Linear Coriolis effect
			 */
			swe_sphere2d_ts_ln_erk_split_uv.euler_timestep_update_lc(
					U_phi, U_vrt, U_div,
					L_U_phi, L_U_vrt, L_U_div,
					i_simulation_timestamp
				);
		}
		else
		{
			/*
			 * L_g(U): Linear gravity modes
			 */
			swe_sphere2d_ts_ln_erk_split_uv.euler_timestep_update_lg(
					U_phi, U_vrt, U_div,
					L_U_phi, L_U_vrt, L_U_div,
					i_simulation_timestamp
				);
		}

		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				L_U_phi, L_U_vrt, L_U_div,
				pos_lon_d, pos_lat_d,
				L_U_phi_D, L_U_vrt_D, L_U_div_D
			);

		if (coriolis_treatment == CORIOLIS_SEMILAGRANGIAN)
		{
			L_U_vrt_D += coriolis_departure_spectral;							// SL treatment of Coriolis!!!
		}
	}
	else
	{
		/*
		 * Method 2) First get variables on departure points, then evaluate L
		 */

		sweet::Data::Sphere2D::DataSpectral U_phi_D, U_vrt_D, U_div_D;


		if (coriolis_treatment == CORIOLIS_SEMILAGRANGIAN)
		{
			semiLagrangian.apply_sl_timeintegration_uv(
					ops,
					U_phi, U_vrt - coriolis_arrival_spectral, U_div,	// SL treatment of Coriolis!!!
					pos_lon_d, pos_lat_d,
					U_phi_D, U_vrt_D, U_div_D
				);

			U_vrt_D += coriolis_departure_spectral;						// SL treatment of Coriolis!!!
		}
		else
		{
			semiLagrangian.apply_sl_timeintegration_uv(
					ops,
					U_phi, U_vrt, U_div,
					pos_lon_d, pos_lat_d,
					U_phi_D, U_vrt_D, U_div_D
				);
		}

		const sweet::Data::Sphere2D::Config *sphere2DDataConfig = io_U_phi_pert.sphere2DDataConfig;
		L_U_phi_D.setup(sphere2DDataConfig, 0);
		L_U_vrt_D.setup(sphere2DDataConfig, 0);
		L_U_div_D.setup(sphere2DDataConfig, 0);

		/*
		 * L_g(U): Linear gravity modes
		 */
		swe_sphere2d_ts_ln_erk_split_uv.euler_timestep_update_lg(
				U_phi_D, U_vrt_D, U_div_D,
				L_U_phi_D, L_U_vrt_D, L_U_div_D,
				i_simulation_timestamp
			);

		if (coriolis_treatment == CORIOLIS_LINEAR)
		{
			/*
			 * L_c(U): Linear Coriolis effect
			 */
			swe_sphere2d_ts_ln_erk_split_uv.euler_timestep_update_lc(
					U_phi_D, U_vrt_D, U_div_D,
					L_U_phi_D, L_U_vrt_D, L_U_div_D,
					i_simulation_timestamp
				);
		}
	}

	/*
	 * Compute R = X_D + 1/2 dt L_D
	 */
	sweet::Data::Sphere2D::DataSpectral R_phi = U_phi_D + (0.5 * i_dt) * L_U_phi_D;
	sweet::Data::Sphere2D::DataSpectral R_vrt = U_vrt_D + (0.5 * i_dt) * L_U_vrt_D;
	sweet::Data::Sphere2D::DataSpectral R_div = U_div_D + (0.5 * i_dt) * L_U_div_D;

	/*
	 * Nonlinear remainder term starts here
	 */
	if (nonlinear_remainder_treatment == NL_REMAINDER_NONLINEAR || coriolis_treatment == CORIOLIS_NONLINEAR)
	{
		/*
		 * N*(t+0.5dt) = 1/2 ([ 2*N(t) - N(t-dt) ]_D + N(t))
		 *
		 * R += dt*N*(t+0.5dt)
		 */

		/*
		 * Compute
		 * [ 2*N(t) - N(t-dt) ]_D
		 */

		/*
		 * N(t-dt)
		 */
		const sweet::Data::Sphere2D::Config *sphere2DDataConfig = io_U_phi_pert.sphere2DDataConfig;
		sweet::Data::Sphere2D::DataSpectral N_U_phi_prev_nr(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral N_U_vrt_prev_nr(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral N_U_div_prev_nr(sphere2DDataConfig, 0);

		if (nonlinear_remainder_treatment == NL_REMAINDER_NONLINEAR)
		{
			swe_sphere2d_ts_ln_erk_split_uv.euler_timestep_update_nr(
					U_phi_prev, U_vrt_prev, U_div_prev,
					N_U_phi_prev_nr, N_U_vrt_prev_nr, N_U_div_prev_nr,
					i_simulation_timestamp-i_dt
			);
		}

		if (coriolis_treatment == CORIOLIS_NONLINEAR)
		{
			swe_sphere2d_ts_ln_erk_split_uv.euler_timestep_update_lc(
					U_phi_prev, U_vrt_prev, U_div_prev,
					N_U_phi_prev_nr, N_U_vrt_prev_nr, N_U_div_prev_nr,
					i_simulation_timestamp-i_dt
				);
		}

		/*
		 * N(t)
		 */
		sweet::Data::Sphere2D::DataSpectral N_U_phi_nr(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral N_U_vrt_nr(sphere2DDataConfig, 0);
		sweet::Data::Sphere2D::DataSpectral N_U_div_nr(sphere2DDataConfig, 0);

		if (nonlinear_remainder_treatment == NL_REMAINDER_NONLINEAR)
		{
			swe_sphere2d_ts_ln_erk_split_uv.euler_timestep_update_nr(
					U_phi, U_vrt, U_div,
					N_U_phi_nr, N_U_vrt_nr, N_U_div_nr,
					i_simulation_timestamp
			);
		}

		if (coriolis_treatment == CORIOLIS_NONLINEAR)
		{
			swe_sphere2d_ts_ln_erk_split_uv.euler_timestep_update_lc(
					U_phi, U_vrt, U_div,
					N_U_phi_nr, N_U_vrt_nr, N_U_div_nr,
					i_simulation_timestamp-i_dt
				);
		}

		/*
		 * N(t+dt)_D = [ 2*N(t) - N(t-dt) ]_D
		 */
		sweet::Data::Sphere2D::DataSpectral N_U_phi_next_D, N_U_vrt_next_D, N_U_div_next_D;
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				2.0 * N_U_phi_nr - N_U_phi_prev_nr,
				2.0 * N_U_vrt_nr - N_U_vrt_prev_nr,
				2.0 * N_U_div_nr - N_U_div_prev_nr,

				pos_lon_d, pos_lat_d,
				N_U_phi_next_D, N_U_vrt_next_D, N_U_div_next_D
			);

		/*
		 * Compute N*(t+0.5dt) = 1/2 ([ 2*N(t) - N(t-dt) ]_D + N(t))
		 * and add to R terms
		 */
		R_phi += (i_dt * 0.5) * (N_U_phi_next_D + N_U_phi_nr);
		R_vrt += (i_dt * 0.5) * (N_U_vrt_next_D + N_U_vrt_nr);
		R_div += (i_dt * 0.5) * (N_U_div_next_D + N_U_div_nr);
	}

	/*
	 * Step 2b) Solve Helmholtz problem
	 * X - 1/2 dt LX = R
	 */

	if (coriolis_treatment == CORIOLIS_LINEAR)
	{
		swe_sphere2d_ts_l_irk.runTimestep(
				R_phi, R_vrt, R_div,
				0.5 * i_dt,
				i_simulation_timestamp
			);
	}
	else
	{
		swe_sphere2d_ts_lg_irk.runTimestep(
				R_phi, R_vrt, R_div,
				0.5 * i_dt,
				i_simulation_timestamp
			);
	}

	/*
	 * Backup previous variables for multi-step SL method
	 */
	U_phi_prev = U_phi;
	U_vrt_prev = U_vrt;
	U_div_prev = U_div;

	/*
	 * Return new fields stored in R_*
	 */
	io_U_phi_pert = R_phi;
	io_U_vrt = R_vrt;
	io_U_div = R_div;
}





PDESWESphere2DTS_ln_settls_uv::PDESWESphere2DTS_ln_settls_uv()
{
}

PDESWESphere2DTS_ln_settls_uv::~PDESWESphere2DTS_ln_settls_uv()
{
}

