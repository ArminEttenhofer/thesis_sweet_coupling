#include <sweet/Data/Sphere2DComplex/Operators.hpp>

namespace sweet {
namespace Data {
namespace Sphere2DComplex {


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
	ir(-1)
{

}



void Operators::setup(
	const Sphere2D::Config *i_sphere2DDataConfig,
	const Sphere2D::Shack *i_shackSphere2DDataOps
)
{
	sphere2DDataConfig = i_sphere2DDataConfig;

	r = i_shackSphere2DDataOps->sphere2d_radius;
	ir = 1.0/r;
}





DataSpectral Operators::implicit_helmholtz(
		const Sphere2DComplex::DataSpectral &i_sphere2d_data,
		const std::complex<double> &i_b,
		double i_radius
)	const
{
	Sphere2DComplex::DataSpectral out(i_sphere2d_data);

	const std::complex<double> b = i_b/(i_radius*i_radius);

	out.spectral_update_lambda(
		[&](
			int n, int m,
			std::complex<double> &io_data
		)
		{
			io_data /= (1.0 + (-b*(double)n*((double)n+1.0)));
		}
	);

	return out;
}


DataSpectral Operators::implicit_J(
		const Sphere2DComplex::DataSpectral &i_sphere2d_data,
		const std::complex<double>& i_dt_two_omega
)	const
{
	Sphere2DComplex::DataSpectral out_sph_data(i_sphere2d_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int n = 0; n <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max; n++)
	{
		int idx = i_sphere2d_data.sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);
		for (int m = -n; m <= n; m++)
		{
			out_sph_data[idx] = i_sphere2d_data[idx] * implicit_J_scalar(n, m, i_dt_two_omega);
			idx++;
		}
	}

	return out_sph_data;
}


DataSpectral Operators::implicit_Jinv(
		const Sphere2DComplex::DataSpectral &i_sphere2d_data,
		const std::complex<double>& i_dt_two_omega
)	const
{
	Sphere2DComplex::DataSpectral out_sph_data(i_sphere2d_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int n = 0; n <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max; n++)
	{
		int idx = i_sphere2d_data.sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);
		for (int m = -n; m <= n; m++)
		{
			out_sph_data[idx] =	i_sphere2d_data.spectral_space_data[idx] / implicit_J_scalar(n, m, i_dt_two_omega);
			idx++;
		}
	}

	return out_sph_data;
}


DataSpectral Operators::implicit_F(
		const Sphere2DComplex::DataSpectral &i_sphere2d_data,
		const std::complex<double>& i_dt_two_omega
)	const
{
	Sphere2DComplex::DataSpectral out_sph_data(i_sphere2d_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int n = 0; n <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max; n++)
	{
		int idx = i_sphere2d_data.sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);
		for (int m = -n; m <= n; m++)
		{
			out_sph_data.spectral_space_data[idx] = 0;

			// out of boundary check for P(n-1, m)
			if (n-1 >= std::abs(m))
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
		const Sphere2DComplex::DataSpectral &i_sphere2d_data,
		const std::complex<double>& i_dt_two_omega
)	const
{
	Sphere2DComplex::DataSpectral out_sph_data(i_sphere2d_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int n = 0; n <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max; n++)
	{
		int idx = i_sphere2d_data.sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);
		for (int m = -n; m <= n; m++)
		{
			out_sph_data.spectral_space_data[idx] = 0;

			// Out of boundary check for P(n-1, m)
			if (n-1 >= std::abs(m))
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
		const Sphere2DComplex::DataSpectral &i_sphere2d_data,
		const std::complex<double>& i_dt
)	const
{
	Sphere2DComplex::DataSpectral out_sph_data(i_sphere2d_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int n = 0; n <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max; n++)
	{
		int idx = i_sphere2d_data.sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);
		for (int m = -n; m <= n; m++)
		{
			out_sph_data[idx] = (i_dt*ir*ir*(double)(n*(n+1))) * i_sphere2d_data[idx];
			idx++;
		}
	}

	return out_sph_data;
}



DataSpectral Operators::implicit_Linv(
		const Sphere2DComplex::DataSpectral &i_sphere2d_data,
		const std::complex<double>& i_dt
)	const
{
	Sphere2DComplex::DataSpectral out_sph_data(i_sphere2d_data.sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int n = 0; n <= i_sphere2d_data.sphere2DDataConfig->spectral_modes_n_max; n++)
	{
		int idx = i_sphere2d_data.sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);
		for (int m = -n; m <= n; m++)
		{
			if (n == 0)
				out_sph_data[idx] = 0.0;
			else
				out_sph_data[idx] = 1.0/(i_dt*ir*ir*Treal(n*(n+1))) * i_sphere2d_data[idx];

			idx++;
		}
	}

	return out_sph_data;
}


DataSpectral Operators::diff_lon(
		const Sphere2DComplex::DataSpectral &i_sph_data
)	const
{
	Sphere2DComplex::DataSpectral out(i_sph_data.sphere2DDataConfig);

	// compute d/dlambda in spectral space
	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int n = 0; n <= i_sph_data.sphere2DDataConfig->spectral_modes_n_max; n++)
	{
		int idx = i_sph_data.sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);
		for (int m = -n; m <= n; m++)
		{
			out.spectral_space_data[idx] = i_sph_data.spectral_space_data[idx]*std::complex<double>(0, m);
			idx++;
		}
	}

	return out;
}


DataSpectral Operators::mu(
		const Sphere2DComplex::DataSpectral &i_sph_data
)
{
	const Sphere2D::Config *sphere2DDataConfig = i_sph_data.sphere2DDataConfig;

	Sphere2DComplex::DataSpectral out = DataSpectral(sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int n = 0; n <= i_sph_data.sphere2DDataConfig->spectral_modes_n_max; n++)
	{
		int idx = i_sph_data.sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);
		for (int m = -n; m <= n; m++)
		{
			out.spectral_space_data[idx] = 0;

			if (n-1 >= std::abs(m))
				out.spectral_space_data[idx] += R(n-1,m)*i_sph_data.spectral_get_(n-1, m);

			if (n+1 <= i_sph_data.sphere2DDataConfig->spectral_modes_n_max)
				out.spectral_space_data[idx] += S(n+1,m)*i_sph_data.spectral_get_(n+1, m);

			idx++;
		}
	}

	return out;
}


DataSpectral Operators::mu2(
		const Sphere2DComplex::DataSpectral &i_sph_data
)	const
{
	const Sphere2D::Config *sphere2DDataConfig = i_sph_data.sphere2DDataConfig;

	Sphere2DComplex::DataSpectral out = DataSpectral(sphere2DDataConfig);

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int n = 0; n <= i_sph_data.sphere2DDataConfig->spectral_modes_n_max; n++)
	{
		int idx = i_sph_data.sphere2DDataConfig->getArrayIndexByModes_Complex(n, -n);
		for (int m = -n; m <= n; m++)
		{
			out.spectral_space_data[idx] = 0;

			if (n-2 >= std::abs(m))
				out.spectral_space_data[idx] += A(n-2,m)*i_sph_data.spectral_get_(n-2, m);

			out.spectral_space_data[idx] += B(n+0,m)*i_sph_data.spectral_get_(n+0, m);

			if (n+2 <= i_sph_data.sphere2DDataConfig->spectral_modes_n_max)
				out.spectral_space_data[idx] += C(n+2,m)*i_sph_data.spectral_get_(n+2, m);

			idx++;
		}
	}

	return out;
}


DataSpectral Operators::laplace(
		const Sphere2DComplex::DataSpectral &i_sph_data
)	const
{
	Sphere2DComplex::DataSpectral out(i_sph_data);

	out.spectral_update_lambda(
			[&](int n, int m, std::complex<double> &o_data)
			{
				o_data *= -(double)n*((double)n+1.0)*ir*ir;
			}
		);

	return out;
}


void Operators::vrtdiv_2_uv(
		const Sphere2DComplex::DataSpectral &i_vrt,
		const Sphere2DComplex::DataSpectral &i_div,
		Sphere2DComplex::DataGrid &o_u,
		Sphere2DComplex::DataGrid &o_v

)	const
{
	Sphere2DComplex::DataSpectral psi = inv_laplace(i_vrt)*ir;
	Sphere2DComplex::DataSpectral chi = inv_laplace(i_div)*ir;

	#if SWEET_DEBUG
		#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
			if (omp_in_parallel())
				SWEETErrorFatal("IN PARALLEL REGION!!!");
		#endif
	#endif


	o_u.setup_if_required(i_vrt.sphere2DDataConfig);
	o_v.setup_if_required(i_vrt.sphere2DDataConfig);

	SHsphtor_to_spat_cplx(
			sphere2DDataConfig->shtns,
			psi.spectral_space_data,
			chi.spectral_space_data,
			o_u.grid_space_data,
			o_v.grid_space_data
	);
}


void Operators::uv_2_vrtdiv(
		const Sphere2DComplex::DataGrid &i_u,
		const Sphere2DComplex::DataGrid &i_v,
		Sphere2DComplex::DataSpectral &o_vrt,
		Sphere2DComplex::DataSpectral &o_div

)	const
{
	o_vrt.setup_if_required(i_u.sphere2DDataConfig);
	o_div.setup_if_required(i_u.sphere2DDataConfig);

	#if SWEET_DEBUG
		#if SWEET_THREADING_SPACE || SWEET_THREADING_TIME_REXI
			if (omp_in_parallel())
				SWEETErrorFatal("IN PARALLEL REGION!!!");
		#endif
	#endif

	spat_cplx_to_SHsphtor(
			sphere2DDataConfig->shtns,
			i_u.grid_space_data,
			i_v.grid_space_data,
			o_vrt.spectral_space_data,
			o_div.spectral_space_data
	);

	o_vrt = laplace(o_vrt)*r;
	o_div = laplace(o_div)*r;
}


DataSpectral Operators::inv_laplace(
		const Sphere2DComplex::DataSpectral &i_sph_data,
		double i_radius
)	const
{
	if (i_radius == -1)
		i_radius = r;

	double ir = 1.0/i_radius;

	Sphere2DComplex::DataSpectral out(i_sph_data);

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

}}}

