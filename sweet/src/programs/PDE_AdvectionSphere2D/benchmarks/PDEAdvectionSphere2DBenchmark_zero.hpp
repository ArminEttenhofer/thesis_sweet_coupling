/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_BENCHMARKS_PDEADVECTIONSPHERE2DBENCHMARK_ZERO_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_BENCHMARKS_PDEADVECTIONSPHERE2DBENCHMARK_ZERO_HPP


#include <programs/PDE_AdvectionSphere2D/benchmarks/PDEAdvectionSphere2DBenchmarks_BaseInterface.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>



class PDEAdvectionSphere2DBenchmark_zero	:
	public PDEAdvectionSphere2DBenchmarks_BaseInterface
{
public:
	PDEAdvectionSphere2DBenchmark_zero()
	{
	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
			i_benchmark_name == "zero"
		;
	}

	void setup_1_shackData()
	{
	}

	void setup_2_withOps(
			sweet::Data::Sphere2D::Operators *io_ops
	)
	{
		ops = io_ops;
	}


	std::string printHelp()
	{
		std::ostringstream stream;
		stream << " * ZERO FIELD TEST CASES:" << std::endl;
		stream << "    - 'zero':	Initialize every field with 0" << std::endl;
		return stream.str();
	}


	void getInitialState(
		std::vector<sweet::Data::Sphere2D::DataSpectral> &o_prognostic_fields,
		sweet::Data::Sphere2D::DataGrid &o_u,
		sweet::Data::Sphere2D::DataGrid &o_v
	)
	{
		SWEET_ASSERT_MSG(o_prognostic_fields.size() == 1, "Only scalar field supported for this benchmark!");

		getInitialState(o_prognostic_fields[0], o_u, o_v);
	}


	void getInitialState(
		sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
		sweet::Data::Sphere2D::DataGrid &o_u,
		sweet::Data::Sphere2D::DataGrid &o_v
	)
	{
		o_phi_pert.spectral_setZero();
		o_u.grid_setZero();
		o_v.grid_setZero();
	}
};

#endif
