/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_BENCHMARKS_PDEADVECTIONSPHERE2DBENCHMARK_ADVECTION_VECTOR_UV_GAUSS_BUMPS_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_BENCHMARKS_PDEADVECTIONSPHERE2DBENCHMARK_ADVECTION_VECTOR_UV_GAUSS_BUMPS_HPP

#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmark_williamson_1_advection_gauss_bump.hpp>
#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmarks_BaseInterface.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <ostream>
#include <sweet/LibMath/VectorMath.hpp>
#include <sweet/Shacks/Dictionary.hpp>



class PDEAdvectionSphere2DBenchmark_advection_vector_uv_gauss_bumps	:
	public PDEAdvectionSphere2DBenchmarks_BaseInterface
{
	PDEAdvectionSphere2DBenchmark_williamson_1_advection_gauss_bump benchmark;

public:
	PDEAdvectionSphere2DBenchmark_advection_vector_uv_gauss_bumps()
	{
	}



	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		return (
				i_benchmark_name == "advection_vector_uv_gauss_bumps"	||
				i_benchmark_name == "advection_vector_uv_gaussian_bumps"	||
				false
			);
	}


	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		PDEAdvectionSphere2DBenchmarks_BaseInterface::shackRegistration(io_shackDict);
		benchmark.shackRegistration(io_shackDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(benchmark);
		return true;
	}

	void setup_1_shackData()
	{
		benchmark.setup_1_shackData();
	}

	void setup_2_withOps(
			sweet::Data::Sphere2D::Operators *io_ops
	)
	{
		ops = io_ops;

		benchmark.setup_2_withOps(ops);
	}



	std::string printHelp()
	{
		std::ostringstream stream;

		stream << " * Advection test case with 2d vector in lat-lon space, setup with gaussian bumps:" << std::endl;
		stream << "    + 'advection_vector_uv_gauss_bumps'" << std::endl;

		return stream.str();
	}



	/*
	 * Return number of prognostic fields to be used
	 */
	int getNumPrognosticFields()
	{
		return 2;
	}



	void getInitialState(
		std::vector<sweet::Data::Sphere2D::DataSpectral> &o_prog_vec,
		sweet::Data::Sphere2D::DataGrid &o_u,
		sweet::Data::Sphere2D::DataGrid &o_v
	)
	{
		sweet::Data::Sphere2D::DataSpectral tmp(ops->sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral vrt(ops->sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral div(ops->sphere2DDataConfig);
		benchmark.getInitialState_Spectral(tmp, vrt, div);

		// Convert vrt/div to velocity field
		ops->vrtdiv_2_uv(vrt, div, o_u, o_v);

		/*
		 * IMPORTANT INFORMATION:
		 * The prognostic fields here are the vrt/div of the velocities!
		 * So we first convert the bumps in u/v field to vrt/div fields.
		 * Be careful about pole problems! The bump needs to be numerically 0 at the poles!!!
		 */


		double lambda_c = 3.0*M_PI/2.0;
		double theta_c = 0.0;
		double i_exp_fac = 20.0;

		sweet::Data::Sphere2D::DataGrid u_phys(ops->sphere2DDataConfig);
		u_phys.grid_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
		{
				double d = std::acos(
						std::sin(theta_c)*std::sin(i_theta) +
						std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
				);

				io_data = std::exp(-d*d*i_exp_fac);
			}
		);

		sweet::Data::Sphere2D::DataGrid v_phys(ops->sphere2DDataConfig);
		v_phys.grid_update_lambda(
			[&](double i_lambda, double i_theta, double &io_data)
		{
				double d = std::acos(
						std::sin(theta_c)*std::sin(i_theta) +
						std::cos(theta_c)*std::cos(i_theta)*std::cos(i_lambda-lambda_c)
				);

				io_data = std::exp(-d*d*i_exp_fac);

			}
		);

		ops->uv_2_vrtdiv(
				u_phys, v_phys,
				o_prog_vec[0], o_prog_vec[1]
			);
	}
};

#endif
