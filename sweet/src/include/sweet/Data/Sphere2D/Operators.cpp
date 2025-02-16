#include <sweet/Data/Sphere2D/Operators.hpp>


namespace sweet {
namespace Data {
namespace Sphere2D {


Operators::Operators(
	const Sphere2D::Config *i_sphere2DDataConfig,
	const Sphere2D::Shack *i_shackSphere2DDataOps
)
{
	setup(i_sphere2DDataConfig, i_shackSphere2DDataOps);
}



Operators::Operators()	:
	sphere2DDataConfig(nullptr),
	r(-1),
	r2(-1),
	ir(-1),
	ir2(-1)
{
}


Sphere2D::DataGrid Operators::getFG_fSphere2D(
		double i_fsphere2d_f0
)	const
{
	SWEET_ASSERT(sphere2DDataConfig != nullptr);

	Sphere2D::DataGrid fg;
	fg.setup(sphere2DDataConfig);

	fg.grid_update_lambda_gaussian_grid(
		[&](double lon, double mu, double &o_data)
		{
			o_data = i_fsphere2d_f0;
		}
	);

	return fg;
}

void Operators::getFG_fSphere2D(
		double i_fsphere2d_f0,
		Sphere2D::DataGrid o_fg
)	const
{
	SWEET_ASSERT(sphere2DDataConfig != nullptr);

	o_fg.grid_update_lambda_gaussian_grid(
		[&](double lon, double mu, double &o_data)
		{
			o_data = i_fsphere2d_f0;
		}
	);
}


Sphere2D::DataGrid Operators::getFG_rotatingSphere2D(
		double i_sphere2d_rotating_coriolis_omega
)	const
{
	SWEET_ASSERT(sphere2DDataConfig != nullptr);

	Sphere2D::DataGrid fg;
	fg.setup(sphere2DDataConfig);

	fg.grid_update_lambda_gaussian_grid(
		[&](double lon, double mu, double &o_data)
		{
			o_data = mu*2.0*i_sphere2d_rotating_coriolis_omega;
		}
	);
	return fg;
}


void Operators::getFG_rotatingSphere2D(
		double i_sphere2d_rotating_coriolis_omega,
		Sphere2D::DataGrid &o_fg
)	const
{
	SWEET_ASSERT(sphere2DDataConfig != nullptr);

	o_fg.grid_update_lambda_gaussian_grid(
		[&](double lon, double mu, double &o_data)
		{
			o_data = mu*2.0*i_sphere2d_rotating_coriolis_omega;
		}
	);
}



void Operators::setup(
	const Sphere2D::Config *i_sphere2DDataConfig,
	const Sphere2D::Shack *i_shackSphere2DDataOps
)
{
	sphere2DDataConfig = i_sphere2DDataConfig;

	r = i_shackSphere2DDataOps->sphere2d_radius;
	r2 = r*r;

	ir = 1.0/r;
	ir2 = ir*ir;

#if 1
	std::vector<double> mx;
	mx.resize(2*sphere2DDataConfig->shtns->nlm);

	st_dt_matrix(sphere2DDataConfig->shtns, mx.data());

	for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		int idx = sphere2DDataConfig->getArrayIndexByModes(m, m);

		for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			double a = (-n+1.0)*R(n-1,m);
			double b = (n+2.0)*S(n+1,m);

			if (n+1 > sphere2DDataConfig->spectral_modes_n_max)
				b = 0;

			double errora = std::abs(a+mx[idx*2+0]);
			double errorb = std::abs(b+mx[idx*2+1]);

			if (errora > 1e-12 || errorb > 1e-12)
			{
				std::cout << idx << ": n=" << n << ", m=" << m << " | " << errora << "\t" << errorb << std::endl;
				SWEETErrorFatal("SAFETY CHECK NOT SUCCESSFUL");
			}

			idx++;
		}
	}
#endif

}


void Operators::clear()
{
	if (sphere2DDataConfig == nullptr)
		return;

	sphere2DDataConfig = nullptr;
}


void Operators::grad_2_vec(
		const Sphere2D::DataSpectral &i_phi,
		Sphere2D::DataGrid &o_u,
		Sphere2D::DataGrid &o_v,
		double i_radius

)	const
{
	double ir = 1.0/i_radius;

	Sphere2D::DataSpectral psi(sphere2DDataConfig);
	psi.spectral_setZero();

	o_u.setup_if_required(i_phi.sphere2DDataConfig);
	o_v.setup_if_required(i_phi.sphere2DDataConfig);
	SHsphtor_to_spat(
					sphere2DDataConfig->shtns,
					psi.spectral_space_data,
					i_phi.spectral_space_data,
					o_u.grid_space_data,
					o_v.grid_space_data
			);

	o_u *= ir;
	o_v *= ir;
}



void Operators::vrtdiv_2_uv(
		const Sphere2D::DataSpectral &i_vrt,
		const Sphere2D::DataSpectral &i_div,
		Sphere2D::DataGrid &o_u,
		Sphere2D::DataGrid &o_v

)	const
{
	/* Calculate stream function and velocity potential (multiplied by earth radius) */
	Sphere2D::DataSpectral psi = inv_laplace(i_vrt)*ir;
	Sphere2D::DataSpectral chi = inv_laplace(i_div)*ir;

	#if SWEET_DEBUG && 0
		#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
			if (omp_in_parallel())
				SWEETErrorFatal("IN PARALLEL REGION!!!");
		#endif
	#endif


	o_u.setup_if_required(i_vrt.sphere2DDataConfig);
	o_v.setup_if_required(i_vrt.sphere2DDataConfig);

	SHsphtor_to_spat(
			sphere2DDataConfig->shtns,
			psi.spectral_space_data,
			chi.spectral_space_data,
			o_u.grid_space_data,
			o_v.grid_space_data
	);
}


void Operators::scalar_spectral_2_physical(
		const Sphere2D::DataSpectral &i_spectral,
		Sphere2D::DataGrid &o_physical

)	const
{
	o_physical.setup_if_required(i_spectral.sphere2DDataConfig);

	SH_to_spat(
			sphere2DDataConfig->shtns,
			i_spectral.spectral_space_data,
			o_physical.grid_space_data
	);
}


Sphere2D::DataGrid Operators::scalar_spectral_2_physical(
		const Sphere2D::DataSpectral &i_spectral

)	const
{
	Sphere2D::DataGrid retdata(i_spectral.sphere2DDataConfig);

	SH_to_spat(
			sphere2DDataConfig->shtns,
			i_spectral.spectral_space_data,
			retdata.grid_space_data
	);

	return retdata;
}


void Operators::scalar_grid_2_spectral(
		const Sphere2D::DataGrid &i_physical,
		Sphere2D::DataSpectral &o_spectral
)	const
{
	o_spectral.setup_if_required(i_physical.sphere2DDataConfig);

	spat_to_SH(
			sphere2DDataConfig->shtns,
			i_physical.grid_space_data,
			o_spectral.spectral_space_data
	);
}


DataSpectral Operators::scalar_grid_2_spectral(
		const Sphere2D::DataGrid &i_physical
)	const
{
	Sphere2D::DataSpectral retdata(i_physical.sphere2DDataConfig);

	spat_to_SH(
			sphere2DDataConfig->shtns,
			i_physical.grid_space_data,
			retdata.spectral_space_data
	);

	return retdata;
}


DataSpectral Operators::uv_2_vort(
		const Sphere2D::DataGrid &i_u,
		const Sphere2D::DataGrid &i_v

)	const
{
	Sphere2D::DataSpectral tmp(sphere2DDataConfig);
	Sphere2D::DataSpectral vort(sphere2DDataConfig);

	#if SWEET_DEBUG && 0
		#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
			if (omp_in_parallel())
				SWEETErrorFatal("IN PARALLEL REGION!!!");
		#endif
	#endif

	spat_to_SHsphtor(
			sphere2DDataConfig->shtns,
			i_u.grid_space_data,
			i_v.grid_space_data,
			vort.spectral_space_data,
			tmp.spectral_space_data
	);

	return laplace(vort)*r;
}


DataSpectral Operators::V_dot_grad_scalar(
		const Sphere2D::DataGrid &i_u_phys,		//!< u velocity
		const Sphere2D::DataGrid &i_v_phys,		//!< v velocity
		const Sphere2D::DataGrid &i_div_phys,		//!< divergence in physical space to avoid transformation
		const Sphere2D::DataGrid &i_scalar_phys	//!< scalar
)	const
{
	return uv_2_div(
			i_u_phys*i_scalar_phys,
			i_v_phys*i_scalar_phys
		)
			- i_div_phys*i_scalar_phys;
}



DataSpectral Operators::uv_2_div(
		const Sphere2D::DataGrid &i_u,
		const Sphere2D::DataGrid &i_v
)	const
{
	Sphere2D::DataSpectral tmp(sphere2DDataConfig);
	Sphere2D::DataSpectral div(sphere2DDataConfig);

	#if SWEET_DEBUG && 0
		#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
			if (omp_in_parallel())
				SWEETErrorFatal("IN PARALLEL REGION!!!");
		#endif
	#endif

	spat_to_SHsphtor(
			sphere2DDataConfig->shtns,
			i_u.grid_space_data,
			i_v.grid_space_data,
			tmp.spectral_space_data,
			div.spectral_space_data
	);

	return laplace(div)*r;
}



void Operators::uv_2_vrtdiv(
		const Sphere2D::DataGrid &i_u,
		const Sphere2D::DataGrid &i_v,
		Sphere2D::DataSpectral &o_vrt,
		Sphere2D::DataSpectral &o_div

)	const
{
	o_vrt.setup_if_required(i_u.sphere2DDataConfig);
	o_div.setup_if_required(i_u.sphere2DDataConfig);

	#if SWEET_DEBUG && 0
		#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
			if (omp_in_parallel())
				SWEETErrorFatal("IN PARALLEL REGION!!!");
		#endif
	#endif

	spat_to_SHsphtor(
			sphere2DDataConfig->shtns,
			i_u.grid_space_data,
			i_v.grid_space_data,
			o_vrt.spectral_space_data,
			o_div.spectral_space_data
	);

	o_vrt = laplace(o_vrt)*r;
	o_div = laplace(o_div)*r;
}



DataSpectral Operators::spectral_one_minus_sinphi_squared_diff_lat_mu(
		const Sphere2D::DataSpectral &i_sph_data
)	const
{
	return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
}



DataSpectral Operators::spectral_cosphi2_diff_lat_mu(
		const Sphere2D::DataSpectral &i_sph_data
)	const
{
	return spectral_one_minus_mu_squared_diff_lat_mu(i_sph_data);
}


DataSpectral Operators::spectral_one_minus_mu_squared_diff_lat_mu(
		const Sphere2D::DataSpectral &i_sph_data
)	const
{
	const Sphere2D::Config *sphere2DDataConfig = i_sph_data.sphere2DDataConfig;

	Sphere2D::DataSpectral out_sph_data(sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= i_sph_data.sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		int idx = i_sph_data.sphere2DDataConfig->getArrayIndexByModes(m, m);

		for (int n = m; n <= i_sph_data.sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			out_sph_data.spectral_space_data[idx] =
					((double)(-n+1.0)*R(n-1,m))*i_sph_data.spectral_get_DEPRECATED(n-1, m) +
					((double)(n+2.0)*S(n+1,m))*i_sph_data.spectral_get_DEPRECATED(n+1, m);

			idx++;
		}
	}

	return out_sph_data;
}


DataSpectral Operators::mu(
		const Sphere2D::DataSpectral &i_sphere2d_data
)	const
{
	Sphere2D::DataSpectral out_sph_data(i_sphere2d_data.sphere2DDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		int idx = i_sphere2d_data.sphere2DDataConfig->getArrayIndexByModes(m, m);

		for (int n = m; n <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			out_sph_data.spectral_space_data[idx] = 0;

			if (n-1 >= m)
				out_sph_data.spectral_space_data[idx] +=
						R_real(n-1,m)*i_sphere2d_data.spectral_get_(n-1, m);

			if (n+1 <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max)
				out_sph_data.spectral_space_data[idx] +=
						S_real(n+1,m)*i_sphere2d_data.spectral_get_(n+1, m);

			idx++;
		}
	}

	return out_sph_data;
}


DataSpectral Operators::implicit_helmholtz(
		const Sphere2D::DataSpectral &i_sphere2d_data,
		const double &i_b,
		double i_sphere2d_radius
)	const
{
	Sphere2D::DataSpectral out(i_sphere2d_data);

	const double b = i_b/(i_sphere2d_radius*i_sphere2d_radius);

	out.spectral_update_lambda(
		[&](
			int n, int m,
			std::complex<double> &io_data
		)
		{
			/*
			 * Laplace operator in SPH space is given by
			 * 	-n*(n+1.0))/r^2
			 */
			io_data /= (1 + (-b*(double)n*((double)n+1.0)));
		}
	);

	return out;
}



DataSpectral Operators::implicit_J(
		const Sphere2D::DataSpectral &i_sphere2d_data,
		double i_dt_two_omega
)	const
{
	Sphere2D::DataSpectral out_sph_data(i_sphere2d_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		int idx = i_sphere2d_data.sphere2DDataConfig->getArrayIndexByModes(m, m);

		for (int n = m; n <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			out_sph_data[idx] = i_sphere2d_data[idx] * implicit_J_scalar(n, m, i_dt_two_omega);
			idx++;
		}
	}

	return out_sph_data;
}


DataSpectral Operators::implicit_Jinv(
		const Sphere2D::DataSpectral &i_sphere2d_data,
		double i_dt_two_omega
)	const
{
	Sphere2D::DataSpectral out_sph_data(i_sphere2d_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		int idx = i_sphere2d_data.sphere2DDataConfig->getArrayIndexByModes(m, m);
		for (int n = m; n <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			out_sph_data[idx] =	i_sphere2d_data.spectral_space_data[idx] / implicit_J_scalar(n, m, i_dt_two_omega);
			idx++;
		}
	}

	return out_sph_data;
}


DataSpectral Operators::implicit_F(
		const Sphere2D::DataSpectral &i_sphere2d_data,
		double i_dt_two_omega
)	const
{
	Sphere2D::DataSpectral out_sph_data(i_sphere2d_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		int idx = i_sphere2d_data.sphere2DDataConfig->getArrayIndexByModes(m, m);

		for (int n = m; n <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			out_sph_data.spectral_space_data[idx] = 0;

			// out of boundary check for P(n-1, m)
			if (n-1 >= m)
			{
				out_sph_data.spectral_space_data[idx] +=
						i_dt_two_omega
						* implicit_f_minus(n, m)
						* i_sphere2d_data.spectral_get_(n-1, m);
			}

			// out of boundary check for P(n+1, m)
			if (n+1 <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max)
			{
				out_sph_data.spectral_space_data[idx] +=
						i_dt_two_omega
						* implicit_f_plus(n, m)
						* i_sphere2d_data.spectral_get_(n+1, m);
			}

			idx++;
		}
	}

	return out_sph_data;
}


DataSpectral Operators::implicit_FJinv(
		const Sphere2D::DataSpectral &i_sphere2d_data,
		double i_dt_two_omega
)	const
{
	Sphere2D::DataSpectral out_sph_data(i_sphere2d_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		int idx = i_sphere2d_data.sphere2DDataConfig->getArrayIndexByModes(m, m);

		for (int n = m; n <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			out_sph_data.spectral_space_data[idx] = 0;

			// Out of boundary check for P(n-1, m)
			if (n-1 >= m)
			{
				out_sph_data.spectral_space_data[idx] +=
						i_dt_two_omega
						* implicit_f_minus(n, m)
						/ implicit_J_scalar(n-1, m, i_dt_two_omega)
						* i_sphere2d_data.spectral_get_(n-1, m);
			}

			// Out of boundary check for P(n+1, m)
			if (n+1 <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max)
			{
				out_sph_data.spectral_space_data[idx] +=
						i_dt_two_omega
						* implicit_f_plus(n, m)
						/ implicit_J_scalar(n+1, m, i_dt_two_omega)
						* i_sphere2d_data.spectral_get_(n+1, m);
			}

			idx++;
		}
	}

	return out_sph_data;
}


DataSpectral Operators::implicit_L(
		const Sphere2D::DataSpectral &i_sphere2d_data,
		double i_dt
)	const
{
	Sphere2D::DataSpectral out_sph_data(i_sphere2d_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		int idx = i_sphere2d_data.sphere2DDataConfig->getArrayIndexByModes(m, m);

		for (int n = m; n <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			/*
			 * Note, the Laplace operator in SPH space is given by
			 * 	-(double)n*((double)n+1.0))/r^2
			 */
			out_sph_data[idx] = i_dt*ir2*(n*(n+1)) * i_sphere2d_data[idx];
			idx++;
		}
	}

	return out_sph_data;
}


DataSpectral Operators::implicit_Linv(
		const Sphere2D::DataSpectral &i_sphere2d_data,
		double i_dt
)	const
{
	Sphere2D::DataSpectral out_sph_data(i_sphere2d_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		int idx = i_sphere2d_data.sphere2DDataConfig->getArrayIndexByModes(m, m);

		for (int n = m; n <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			if (n == 0)
				out_sph_data[idx] = 0.0;
			else
				out_sph_data[idx] = 1.0/(i_dt*ir2*(n*(n+1))) * i_sphere2d_data[idx];

			idx++;
		}
	}

	return out_sph_data;
}


DataSpectral Operators::mu2(
		const Sphere2D::DataSpectral &i_sph_data
)	const
{
	const Sphere2D::Config *sphere2DDataConfig = i_sph_data.sphere2DDataConfig;

	Sphere2D::DataSpectral out_sph_data = DataSpectral(sphere2DDataConfig);


	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= i_sph_data.sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		int idx = i_sph_data.sphere2DDataConfig->getArrayIndexByModes(m, m);

		for (int n = m; n <= i_sph_data.sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			out_sph_data.spectral_space_data[idx] =
					+A(n-2,m)*i_sph_data.spectral_get_DEPRECATED(n-2, m)
					+B(n+0,m)*i_sph_data.spectral_get_DEPRECATED(n+0, m)
					+C(n+2,m)*i_sph_data.spectral_get_DEPRECATED(n+2, m)
					;
			idx++;
		}
	}

	return out_sph_data;
}


DataSpectral Operators::laplace(
		const Sphere2D::DataSpectral &i_sph_data
)	const
{
	Sphere2D::DataSpectral out_sph_data(i_sph_data);

	out_sph_data.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &o_data)
			{
				o_data *= -(double)n*((double)n+1.0)*(ir*ir);
			}
		);

	return out_sph_data;
}


DataSpectral Operators::root_laplace(
		const Sphere2D::DataSpectral &i_sph_data
)	const
{
	Sphere2D::DataSpectral out_sph_data(i_sph_data);

	out_sph_data.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &o_data)
			{
				o_data *= std::sqrt((double)n*((double)n+1.0))*(ir);
			}
		);

	return out_sph_data;
}


DataSpectral Operators::inv_laplace(
		const Sphere2D::DataSpectral &i_sph_data
)	const
{
	Sphere2D::DataSpectral out(i_sph_data);

	out.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &o_data)
			{
				if (n != 0)
					o_data /= -(double)n*((double)n+1.0)*ir*ir;
				else
					o_data = 0;
			}
		);

	return out;
}


DataSpectral Operators::inv_root_laplace(
		const Sphere2D::DataSpectral &i_sph_data
)	const
{
	Sphere2D::DataSpectral out(i_sph_data);

	out.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &o_data)
			{
				if (n != 0)
					o_data /= std::sqrt((double)n*((double)n+1.0))*ir;
				else
					o_data = 0;
			}
		);

	return out;
}


DataSpectral Operators::implicit_diffusion(
		const Sphere2D::DataSpectral &i_data,
		double i_coef,
		/*int i_order*/
		double i_r
)
{
	Sphere2D::DataSpectral out = i_data;

	const double scalar = i_coef;
	const double r      = i_r;

	out  = out.spectral_solve_helmholtz(1.0,  -scalar, r);

	return out;
}


DataSpectral Operators::implicit_hyperdiffusion(
		const Sphere2D::DataSpectral &i_data,
		double i_coef,
		int i_order,
		double i_r
)
{
	Sphere2D::DataSpectral out = i_data;

	const double r      = i_r;

	std::array<double, 4> visc_factors;
	if (i_order == 2)
		visc_factors = {-i_coef, 0, 0, 0};
	else if (i_order == 4)
		visc_factors = {0, -i_coef, 0, 0};
	else if (i_order == 6)
		visc_factors = {0, 0, -i_coef, 0};
	else if (i_order == 8)
		visc_factors = {0, 0, 0, -i_coef};
	else
		SWEETErrorFatal("This viscosity order is not supported: " + std::to_string(i_order));

	out  = out.spectral_solve_helmholtz_higher_order(1.0, visc_factors, r);

	return out;
}

}}}
