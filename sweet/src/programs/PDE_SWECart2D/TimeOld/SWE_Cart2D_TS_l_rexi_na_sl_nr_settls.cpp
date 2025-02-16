/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_cart2d.cpp
 *					which was also written by Pedro Peixoto
 */


#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_rexi_na_sl_nr_settls.hpp>


void SWE_Cart2D_TS_l_rexi_na_sl_nr_settls::runTimestep(
		sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables - perturbation of height!
		sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables - zonal velocity
		sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables - meridional velocity

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETErrorFatal("SWE_Cart2D_TS_l_rexi_na_sl_nd_settls: Only constant time step size allowed (Please set --dt)");

	const sweet::Data::Cart2D::Config *cart2DDataConfig = io_h.cart2DDataConfig;

	//Out vars
	sweet::Data::Cart2D::DataSpectral h(io_h.cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral u(io_h.cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral v(io_h.cart2DDataConfig);

	//Temporary
	sweet::Data::Cart2D::DataSpectral N_h(io_h.cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral N_u(io_h.cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral N_v(io_h.cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral hdiv(io_h.cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral N_h_prev(io_h.cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral N_h_ext(io_h.cart2DDataConfig);

	// Departure points
	sweet::Data::Vector::Vector<double> posx_d(io_h.cart2DDataConfig->grid_number_elements);
	sweet::Data::Vector::Vector<double> posy_d(io_h.cart2DDataConfig->grid_number_elements);

	// Parameters
	double dt = i_dt;

	sweet::Data::Cart2D::Staggering staggering;
	SWEET_ASSERT(staggering.staggering_type == 'a');

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

	//Preserve io unmodified
	u = io_u;
	v = io_v;
	h = io_h;

	// Calculate departure points
	//Calculate departure points - always force to be second order accurate!
	semiLagrangian.semi_lag_departure_points_settls(
			u_prev.toGrid(),	v_prev.toGrid(),
			u.toGrid(),		v.toGrid(),
			posx_a,		posy_a,
			dt,
			posx_d,	posy_d,			// output
			shackCart2DDataOps->cart2d_domain_size,
			&staggering,
			2, //shackDict.disc.timestepping_order

			shackSemiLagrangian->semi_lagrangian_max_iterations,
			shackSemiLagrangian->semi_lagrangian_convergence_threshold
	);

	N_u.spectral_setZero();
	N_v.spectral_setZero();
	N_h.spectral_setZero();
	N_h_prev.spectral_setZero();
	N_h_ext.spectral_setZero();
	hdiv.spectral_setZero();

	//Original more stable scheme

	//Calculate nonlinear terms - not done in case of only linear divergence (linear div is already in linear part)
	if (!use_only_linear_divergence) // Full nonlinear case
	{
		// Calculate nonlinear term for the previous time step
		N_h = -h_prev * (ops->diff_c_x(u_prev) + ops->diff_c_y(v_prev));

		if (shackPDESWECart2D->use_nonlinear_only_visc != 0)
		{
#if !SWEET_USE_CART2D_SPECTRAL_SPACE
			SWEETErrorFatal("Implicit diffusion only supported with spectral space activated");
#else
			N_h = ops->implicit_diffusion(N_h, shackTimestepControl->currentTimestepSize*shackPDESWECart2D->viscosity, shackPDESWECart2D->viscosity_order);
#endif
		}

		//Calculate exp(Ldt)N(n-1), relative to previous timestep
		ts_l_rexi.runTimestep(N_h, N_u, N_v, i_dt, i_simulation_timestamp);

		//Use N_h to store now the nonlinearity of the current time (prev will not be required anymore)
		//Update the nonlinear terms with the constants relative to dt
		// N=dtN^n-0.5dt exp(dtL)N^n-1 from paper
		// N=-h*div is calculate in cartesian space (pseudo-spectrally)
		hdiv =  - h * (ops->diff_c_x(io_u) + ops->diff_c_y(io_v));
		if (shackPDESWECart2D->use_nonlinear_only_visc != 0)
		{
#if !SWEET_USE_CART2D_SPECTRAL_SPACE
			SWEETErrorFatal("Implicit diffusion only supported with spectral space activated");
#else
			hdiv = ops->implicit_diffusion(hdiv, shackTimestepControl->currentTimestepSize*shackPDESWECart2D->viscosity, shackPDESWECart2D->viscosity_order);
#endif
		}

		N_u = -0.5 * dt * N_u; // N^n of u term is zero
		N_v = -0.5 * dt * N_v; // N^n of v term is zero
		N_h = dt * hdiv - 0.5 * dt * N_h ; //N^n of h has the nonlin term

		//Build variables to be interpolated to dep. points
		// This is the W^n term in documentation
		u = u + N_u;
		v = v + N_v;
		h = h + N_h;
	}

	// Interpolate W to departure points
	h = sampler2D.bicubic_scalar(h, posx_d, posy_d, -0.5, -0.5);
	u = sampler2D.bicubic_scalar(u, posx_d, posy_d, -0.5, -0.5);
	v = sampler2D.bicubic_scalar(v, posx_d, posy_d, -0.5, -0.5);


	// Add nonlinearity in h
	if (!use_only_linear_divergence) // Full nonlinear case
	{
		h = h + 0.5 * dt * hdiv;
	}

	/*
	 * Calculate the exp(Ldt) of the resulting u,v,h
	 */
	//Calculate phi_0 of interpolated U
	sweet::Data::Cart2D::DataSpectral phi0_Un_h(cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral phi0_Un_u(cart2DDataConfig);
	sweet::Data::Cart2D::DataSpectral phi0_Un_v(cart2DDataConfig);
	//ts_l_rexi.runTimestep(h, u, v, i_dt, i_simulation_timestamp);
	ts_l_rexi.runTimestep(
			h, u, v,
			phi0_Un_h, phi0_Un_u, phi0_Un_v,
			i_dt,
			i_simulation_timestamp
	);

	h = phi0_Un_h;
	u = phi0_Un_u;
	v = phi0_Un_v;

	// Set time (n) as time (n-1)
	h_prev = io_h;
	u_prev = io_u;
	v_prev = io_v;

	// output data
	io_h = h;
	io_u = u;
	io_v = v;
}



/*
 * Setup
 */
bool SWE_Cart2D_TS_l_rexi_na_sl_nr_settls::setup(
		sweet::Data::Cart2D::Operators *io_ops
)
{
	PDESWECart2DTS_BaseInterface::setup(io_ops);

	h_prev.setup(ops->cart2DDataConfig);
	u_prev.setup(ops->cart2DDataConfig);
	v_prev.setup(ops->cart2DDataConfig);

	posx_a.setup(ops->cart2DDataConfig->grid_number_elements);
	posy_a.setup(ops->cart2DDataConfig->grid_number_elements);

	posx_d.setup(ops->cart2DDataConfig->grid_number_elements);
	posy_d.setup(ops->cart2DDataConfig->grid_number_elements);


	use_only_linear_divergence = shackPDESWECart2D->use_only_linear_divergence;

	ts_l_rexi.setup(io_ops, "phi0");

	if (shackCart2DDataOps->space_grid_use_c_staggering)
		SWEETErrorFatal("SWE_Cart2D_TS_l_rexi_na_sl_nd_settls: Staggering not supported for l_rexi_na_sl_nd_settls");

	//with_linear_div_only = i_use_linear_div;

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

	// Initialize arrival points with h position
	sweet::Data::Vector::Vector<double> pos_x = sweet::Data::Cart2D::Convert::DataGrid_2_Vector::convert(tmp_x);
	sweet::Data::Vector::Vector<double> pos_y = sweet::Data::Cart2D::Convert::DataGrid_2_Vector::convert(tmp_y);


	double cell_size_x = shackCart2DDataOps->cart2d_domain_size[0]/(double)shackCart2DDataOps->space_res_physical[0];
	double cell_size_y = shackCart2DDataOps->cart2d_domain_size[1]/(double)shackCart2DDataOps->space_res_physical[1];

	// Initialize arrival points with h position
	posx_a = pos_x+0.5*cell_size_x;
	posy_a = pos_y+0.5*cell_size_y;

	return true;
}



bool SWE_Cart2D_TS_l_rexi_na_sl_nr_settls::shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	PDESWECart2DTS_BaseInterface::shackRegistration(io_shackDict);

	ts_l_rexi.shackRegistration(io_shackDict);
	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ts_l_rexi);
	return true;
}
