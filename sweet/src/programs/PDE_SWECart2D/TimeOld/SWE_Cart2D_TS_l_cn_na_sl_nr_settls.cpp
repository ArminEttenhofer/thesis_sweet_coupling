/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_cart2d.cpp
 *					which was also written by Pedro Peixoto
 */

#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_cn_na_sl_nr_settls.hpp>

/**
 * Solve  SWE with Crank-Nicolson implicit time stepping
 *  (spectral formulation for Helmholtz eq) with semi-Lagrangian
 *   SL-SI-SP
 *
 * U_t = L U(0)
 *
 * Fully implicit version:
 *
 * (U(tau) - U(0)) / tau = 0.5*(L U(tau) + L U(0))
 *
 * <=> U(tau) - U(0) =  tau * 0.5*(L U(tau) + L U(0))
 *
 * <=> U(tau) - 0.5* L tau U(tau) = U(0) + tau * 0.5*L U(0)
 *
 * <=> (1 - 0.5 L tau) U(tau) = (1 + tau*0.5*L) U(0)
 *
 * <=> (2/tau - L) U(tau) = (2/tau + L) U(0)
 *
 * <=> U(tau) = (2/tau - L)^{-1} (2/tau + L) U(0)
 *
 * Semi-implicit has Coriolis term as totally explicit
 *
 * Semi-Lagrangian:
 *   U(tau) is on arrival points
 *   U(0) is on departure points
 *
 * Nonlinear term is added following Hortal (2002)
 * http://onlinelibrary.wiley.com/doi/10.1002/qj.200212858314/pdf
 *
 */
void SWE_Cart2D_TS_l_cn_na_sl_nr_settls::runTimestep(
		sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETErrorFatal("SWE_Cart2D_TS_l_cn_na_sl_nd_settls: Only constant time step size allowed (Please set --dt)");

	if (i_simulation_timestamp == 0)
	{
#if (!SWEET_PARAREAL) && (!SWEET_XBRAID)
		/*
		 * First time step
		 */
		h_prev = io_h;
		u_prev = io_u;
		v_prev = io_v;
#endif
	}

	// Out vars
	sweet::Data::Cart2D::DataSpectral h(io_h.cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral u(io_h.cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral v(io_h.cart2DDataConfig);

	// Departure points and arrival points
	sweet::Data::Vector::Vector<double> posx_d = posx_a;
	sweet::Data::Vector::Vector<double> posy_d = posy_a;

	// Parameters
	double h_bar = shackPDESWECart2D->h0;
	double g = shackPDESWECart2D->gravitation;
	double f0 = shackPDESWECart2D->cart2d_rotating_f0;
	double dt = i_dt;
	double alpha = 2.0/dt;
	double kappa = alpha*alpha;
	kappa += f0*f0;

	sweet::Data::Cart2D::Staggering staggering;
	SWEET_ASSERT(staggering.staggering_type == 'a');

	// Calculate departure points
	semiLagrangian.semi_lag_departure_points_settls(
			u_prev.toGrid(),	v_prev.toGrid(),
			io_u.toGrid(),		io_v.toGrid(),
			posx_a,	posy_a,
			dt,
			posx_d,	posy_d,
			shackCart2DDataOps->cart2d_domain_size,
			&staggering,
			shackPDESWETimeDisc->timestepping_order,

			shackSemiLagrangian->semi_lagrangian_max_iterations,
			shackSemiLagrangian->semi_lagrangian_convergence_threshold

	);


	// Calculate Divergence and vorticity spectrally
	sweet::Data::Cart2D::DataSpectral div = ops->diff_c_x(io_u) + ops->diff_c_y(io_v);

	// This could be pre-stored
	sweet::Data::Cart2D::DataSpectral div_prev = ops->diff_c_x(u_prev) + ops->diff_c_y(v_prev);

	/**
	 * Calculate the RHS
	 * 
	 * The terms related to alpha include the current solution.
	 * 
	 * In this implementation, the original formulation is rearranged to
	 * 
	 *    U + 1/2 dt L U + dt N(U)
	 * 
	 *    = 1/2 * dt (2.0/dt U + L U + 2.0 * N(U))
	 */
	sweet::Data::Cart2D::DataSpectral rhs_u = alpha * io_u + f0 * io_v    - g * ops->diff_c_x(io_h);
	sweet::Data::Cart2D::DataSpectral rhs_v =  - f0 * io_u + alpha * io_v - g * ops->diff_c_y(io_h);
	sweet::Data::Cart2D::DataSpectral rhs_h = alpha * io_h - h_bar * div;

	// All the RHS are to be evaluated at the departure points
	sweet::Data::Cart2D::DataGrid rhs_u_phys = rhs_u.toGrid();
	sweet::Data::Cart2D::DataGrid rhs_v_phys = rhs_v.toGrid();
	sweet::Data::Cart2D::DataGrid rhs_h_phys = rhs_h.toGrid();
	rhs_u = sampler2D.bicubic_scalar(rhs_u_phys, posx_d, posy_d, -0.5, -0.5);
	rhs_v = sampler2D.bicubic_scalar(rhs_v_phys, posx_d, posy_d, -0.5, -0.5);
	rhs_h = sampler2D.bicubic_scalar(rhs_h_phys, posx_d, posy_d, -0.5, -0.5);

	// Calculate nonlinear term at half timestep and add to RHS of h eq.

	if (!use_only_linear_divergence) //full nonlinear case
	{
		// Extrapolation
		sweet::Data::Cart2D::DataSpectral hdiv = 2.0 * io_h * div - h_prev * div_prev;
		sweet::Data::Cart2D::DataSpectral nonlin(io_h.cart2DDataConfig);
		if (shackPDESWECart2D->use_nonlinear_only_visc != 0)
		{
#if !SWEET_USE_CART2D_SPECTRAL_SPACE
			SWEETErrorFatal("Implicit diffusion only supported with spectral space activated");
#else
			// Add diffusion (stabilisation)
			hdiv = ops->implicit_diffusion(hdiv, shackTimestepControl->currentTimestepSize*shackPDESWECart2D->viscosity, shackPDESWECart2D->viscosity_order);
#endif
		}
		// Average
		sweet::Data::Cart2D::DataGrid hdiv_phys = hdiv.toGrid();
		nonlin = 0.5*(io_h*div) + 0.5*sampler2D.bicubic_scalar(hdiv_phys, posx_d, posy_d, -0.5, -0.5);

		// Add to RHS h (TODO (2020-03-16): No clue why there's a -2.0)
		rhs_h = rhs_h - 2.0*nonlin;

	}


	// Build Helmholtz eq.
	sweet::Data::Cart2D::DataSpectral rhs_div = ops->diff_c_x(rhs_u)+ops->diff_c_y(rhs_v);
	sweet::Data::Cart2D::DataSpectral rhs_vort = ops->diff_c_x(rhs_v)-ops->diff_c_y(rhs_u);
	sweet::Data::Cart2D::DataSpectral rhs     = kappa* rhs_h / alpha - h_bar * rhs_div - f0 * h_bar * rhs_vort / alpha;

	// Helmholtz solver
	helmholtz_spectral_solver(kappa, g*h_bar, rhs, h);

	// Update u and v
	u = (1/kappa)*
			( alpha *rhs_u + f0 * rhs_v
					- g * alpha * ops->diff_c_x(h)
					- g * f0 * ops->diff_c_y(h))
					;

	v = (1/kappa)*
			( alpha *rhs_v - f0 * rhs_u
					+ g * f0 * ops->diff_c_x(h)
					- g * alpha * ops->diff_c_y(h))
					;

	// Set time (n) as time (n-1)
	h_prev = io_h;
	u_prev = io_u;
	v_prev = io_v;

	// output data
	io_h = h;
	io_u = u;
	io_v = v;
}



bool SWE_Cart2D_TS_l_cn_na_sl_nr_settls::setup(
		sweet::Data::Cart2D::Operators *io_ops
)
{
	PDESWECart2DTS_BaseInterface::setup(io_ops);

	use_only_linear_divergence = shackPDESWECart2D->use_only_linear_divergence;

	h_prev.setup(io_ops->cart2DDataConfig);
	u_prev.setup(io_ops->cart2DDataConfig);
	v_prev.setup(io_ops->cart2DDataConfig);

	posx_a.setup(io_ops->cart2DDataConfig->grid_number_elements);
	posy_a.setup(io_ops->cart2DDataConfig->grid_number_elements);

	posx_d.setup(io_ops->cart2DDataConfig->grid_number_elements);
	posy_d.setup(io_ops->cart2DDataConfig->grid_number_elements);


	if (shackCart2DDataOps->space_grid_use_c_staggering)
		SWEETErrorFatal("SWE_Cart2D_TS_l_cn_na_sl_nd_settls: Staggering not supported for l_cn_na_sl_nd_settls");

	// Setup sampler for future interpolations
	sampler2D.setup(shackCart2DDataOps->cart2d_domain_size, ops->cart2DDataConfig);

	// Setup semi-lag
	semiLagrangian.setup(shackCart2DDataOps->cart2d_domain_size, ops->cart2DDataConfig);


	sweet::Data::Cart2D::DataGrid tmp_x(ops->cart2DDataConfig);
	tmp_x.grid_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
		io_data = ((double)i)*shackCart2DDataOps->cart2d_domain_size[0]/(double)shackCart2DDataOps->space_res_physical[0];
			},
			false
	);
	sweet::Data::Cart2D::DataGrid tmp_y(ops->cart2DDataConfig);
	tmp_y.grid_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
		io_data = ((double)j)*shackCart2DDataOps->cart2d_domain_size[1]/(double)shackCart2DDataOps->space_res_physical[1];
			},
			false
	);

	//Initialize arrival points with h position
	sweet::Data::Vector::Vector<double> pos_x = sweet::Data::Cart2D::Convert::DataGrid_2_Vector::convert(tmp_x);
	sweet::Data::Vector::Vector<double> pos_y = sweet::Data::Cart2D::Convert::DataGrid_2_Vector::convert(tmp_y);

	double cell_size_x = shackCart2DDataOps->cart2d_domain_size[0]/(double)shackCart2DDataOps->space_res_physical[0];
	double cell_size_y = shackCart2DDataOps->cart2d_domain_size[1]/(double)shackCart2DDataOps->space_res_physical[1];

	// Initialize arrival points with h position
	posx_a = pos_x+0.5*cell_size_x;
	posy_a = pos_y+0.5*cell_size_y;

	return true;
}



SWE_Cart2D_TS_l_cn_na_sl_nr_settls::~SWE_Cart2D_TS_l_cn_na_sl_nr_settls()
{
}

