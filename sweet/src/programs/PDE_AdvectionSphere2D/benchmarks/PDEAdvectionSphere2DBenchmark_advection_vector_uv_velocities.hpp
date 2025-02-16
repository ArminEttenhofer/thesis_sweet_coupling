/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_BENCHMARKS_PDEADVECTIONSPHERE2DBENCHMARK_ADVECTION_VECTOR_UV_VELOCITIES_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_BENCHMARKS_PDEADVECTIONSPHERE2DBENCHMARK_ADVECTION_VECTOR_UV_VELOCITIES_HPP

#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmark_williamson_1_advection_gauss_bump.hpp>
#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmarks_BaseInterface.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <ostream>
#include <sweet/LibMath/VectorMath.hpp>
#include <sweet/Shacks/Dictionary.hpp>



class PDEAdvectionSphere2DBenchmark_advection_vector_uv_velocities	:
	public PDEAdvectionSphere2DBenchmarks_BaseInterface
{
	PDEAdvectionSphere2DBenchmark_williamson_1_advection_gauss_bump benchmark;

public:
	PDEAdvectionSphere2DBenchmark_advection_vector_uv_velocities()
	{
	}

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		return (
				i_benchmark_name == "advection_vector_uv_velocities"	||
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

		stream << " * Advection test case with 2d vector in lat-lon space:" << std::endl;
		stream << "    + 'advection_vector_uv_velocities'" << std::endl;

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
		std::vector<sweet::Data::Sphere2D::DataSpectral> &o_prognostic_fields,
		sweet::Data::Sphere2D::DataGrid &o_u,
		sweet::Data::Sphere2D::DataGrid &o_v
	)
	{
		SWEET_ASSERT_MSG(o_prognostic_fields.size() == 2, "Only a vectorial field (3 elements) supported for this benchmark!");

		const sweet::Data::Sphere2D::Config *sphere2DDataConfig = o_prognostic_fields[0].sphere2DDataConfig;

		sweet::Data::Sphere2D::DataSpectral tmp(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral vrt(sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral div(sphere2DDataConfig);
		benchmark.getInitialState_Spectral(tmp, vrt, div);

		/*
		 * Setup velocity field
		 */
		ops->vrtdiv_2_uv(vrt, div, o_u, o_v);

		/*
		 * IMPORTANT INFORMATION:
		 * Here, we like to use the velocities as prognostic fields.
		 *
		 * The prognostic fields here are the vrt/div of the velocities!!!
		 */
		o_prognostic_fields[0] = vrt;
		o_prognostic_fields[1] = div;
	}
};

#endif
