/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere2DTS_l_irk_na_sl_nr_settls_vd_only.hpp"



bool PDESWESphere2DTS_l_irk_na_sl_nr_settls_vd_only::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	if (i_timestepping_method == "l_irk_na_sl_nr_settls_vd_only")
		return true;

	return false;
}

std::string PDESWESphere2DTS_l_irk_na_sl_nr_settls_vd_only::getIDString()
{
	return "l_irk_na_sl_nr_settls_vd_only";
}


bool PDESWESphere2DTS_l_irk_na_sl_nr_settls_vd_only::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup_main(
		io_ops,
		shackPDESWETimeDisc->timestepping_order
	);
}

bool PDESWESphere2DTS_l_irk_na_sl_nr_settls_vd_only::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_timestepping_order
)
{
	ops = io_ops;

	timestepping_order = i_timestepping_order;

	if (timestepping_order != 1 && timestepping_order != 2)
		SWEETErrorFatal("Invalid time stepping order, must be 1 or 2");

	// Setup semi-lag
	semiLagrangian.setup(ops->sphere2DDataConfig, shackTimesteppingSemiLagrangianSphere2DData, timestepping_order);

	// Initialize with 1st order
	swe_sphere2d_ts_ln_erk_split_vd__l_erk_1st_order.setup_main(
			ops,
			1, true, true, false, false, false
		);

	// Initialize with 1st order and half time step size
	swe_sphere2d_ts_l_irk.setup_main(
			ops,
			1,
			0.5 * shackTimestepControl->currentTimestepSize,
			shackPDESWETimeDisc->timestepping_crank_nicolson_filter,
			false
		);

	return true;
}




void PDESWESphere2DTS_l_irk_na_sl_nr_settls_vd_only::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		SWEETErrorFatal("TODO run_timestep_1st_order_pert");
	}
	else if (timestepping_order == 2)
	{
		run_timestep_2nd_order(io_phi_pert, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
	}
	else
	{
		SWEETErrorFatal("Only orders 1/2 supported (ERRID 098nd89eje)");
	}
}




void PDESWESphere2DTS_l_irk_na_sl_nr_settls_vd_only::run_timestep_2nd_order(
	sweet::Data::Sphere2D::DataSpectral &io_U_phi,
	sweet::Data::Sphere2D::DataSpectral &io_U_vrt,
	sweet::Data::Sphere2D::DataSpectral &io_U_div,

	double i_dt,		
	double i_simulation_timestamp
)
{
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
	semiLagrangian.apply_sl_timeintegration_vd(
			ops,
			U_phi, U_vrt, U_div,
			pos_lon_d, pos_lat_d,
			U_phi_D, U_vrt_D, U_div_D
		);

	/*
	 * Compute L_D
	 */
	const sweet::Data::Sphere2D::Config *sphere2DDataConfig = io_U_phi.sphere2DDataConfig;

	/*
	 * Method 1) First evaluate L, then sample result at departure point
	 */
	sweet::Data::Sphere2D::DataSpectral L_U_phi(sphere2DDataConfig, 0), L_U_vrt(sphere2DDataConfig, 0), L_U_div(sphere2DDataConfig, 0);

	/*
	 * L_g(U): Linear gravity modes
	 */
	swe_sphere2d_ts_ln_erk_split_vd__l_erk_1st_order.euler_timestep_update_lg(
			U_phi, U_vrt, U_div,
			L_U_phi, L_U_vrt, L_U_div,
			i_simulation_timestamp
		);

	/*
	 * L_c(U): Linear Coriolis effect
	 */
	swe_sphere2d_ts_ln_erk_split_vd__l_erk_1st_order.euler_timestep_update_lc(
			U_phi, U_vrt, U_div,
			L_U_phi, L_U_vrt, L_U_div,
			i_simulation_timestamp
		);


	sweet::Data::Sphere2D::DataSpectral L_U_phi_D, L_U_vrt_D, L_U_div_D;
	semiLagrangian.apply_sl_timeintegration_vd(
			ops,
			L_U_phi, L_U_vrt, L_U_div,
			pos_lon_d, pos_lat_d,
			L_U_phi_D, L_U_vrt_D, L_U_div_D
		);

	/*
	 * Compute R = X_D + 1/2 dt L_D
	 */
	sweet::Data::Sphere2D::DataSpectral R_phi = U_phi_D + (0.5 * i_dt) * L_U_phi_D;
	sweet::Data::Sphere2D::DataSpectral R_vrt = U_vrt_D + (0.5 * i_dt) * L_U_vrt_D;
	sweet::Data::Sphere2D::DataSpectral R_div = U_div_D + (0.5 * i_dt) * L_U_div_D;

	/*
	 * Nonlinear remainder term starts here
	 */
	if (1)
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
		sweet::Data::Sphere2D::DataSpectral N_U_phi_prev_nr(sphere2DDataConfig, 0), N_U_vrt_prev_nr(sphere2DDataConfig, 0), N_U_div_prev_nr(sphere2DDataConfig, 0);

		swe_sphere2d_ts_ln_erk_split_vd__l_erk_1st_order.euler_timestep_update_nr(
				U_phi_prev, U_vrt_prev, U_div_prev,
				N_U_phi_prev_nr, N_U_vrt_prev_nr, N_U_div_prev_nr,
				i_simulation_timestamp-i_dt
		);

		/*
		 * N(t)
		 */
		sweet::Data::Sphere2D::DataSpectral N_U_phi(sphere2DDataConfig, 0), N_U_vrt(sphere2DDataConfig, 0), N_U_div(sphere2DDataConfig, 0);

		swe_sphere2d_ts_ln_erk_split_vd__l_erk_1st_order.euler_timestep_update_nr(
				U_phi, U_vrt, U_div,
				N_U_phi, N_U_vrt, N_U_div,
				i_simulation_timestamp
		);

		/*
		 * N(t+dt)_D = [ 2*N(t) - N(t-dt) ]_D
		 */
		sweet::Data::Sphere2D::DataSpectral N_U_phi_next_D, N_U_vrt_next_D, N_U_div_next_D;
		semiLagrangian.apply_sl_timeintegration_vd(
				ops,
				2.0 * N_U_phi - N_U_phi_prev_nr,
				2.0 * N_U_vrt - N_U_vrt_prev_nr,
				2.0 * N_U_div - N_U_div_prev_nr,

				pos_lon_d, pos_lat_d,
				N_U_phi_next_D, N_U_vrt_next_D, N_U_div_next_D
			);

		/*
		 * Compute midpoint
		 * N(t+dt/2)_m = dt * 1/2 ([ 2*N(t) - N(t-dt) ]_D + N(t)_A)
		 */
		R_phi += (i_dt * 0.5) * (N_U_phi_next_D + N_U_phi);
		R_vrt += (i_dt * 0.5) * (N_U_vrt_next_D + N_U_vrt);
		R_div += (i_dt * 0.5) * (N_U_div_next_D + N_U_div);
	}

	/*
	 * Step 2b) Solve Helmholtz problem
	 * X - 1/2 dt LX = R
	 */
	swe_sphere2d_ts_l_irk.runTimestep(
			R_phi, R_vrt, R_div,
			0.5 * i_dt,
			i_simulation_timestamp
		);

	/*
	 * Backup previous variables for multi-step SL method
	 */
	U_phi_prev = U_phi;
	U_vrt_prev = U_vrt;
	U_div_prev = U_div;

	/*
	 * Return new fields stored in R_*
	 */
	io_U_phi = R_phi;
	io_U_vrt = R_vrt;
	io_U_div = R_div;
}





PDESWESphere2DTS_l_irk_na_sl_nr_settls_vd_only::PDESWESphere2DTS_l_irk_na_sl_nr_settls_vd_only()
{
}



PDESWESphere2DTS_l_irk_na_sl_nr_settls_vd_only::~PDESWESphere2DTS_l_irk_na_sl_nr_settls_vd_only()
{
}

