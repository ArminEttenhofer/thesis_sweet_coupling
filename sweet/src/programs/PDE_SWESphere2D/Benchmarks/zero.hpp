/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_ZERO_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_ZERO_HPP


#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "BaseInterface.hpp"


namespace PDE_SWESphere2D {
namespace Benchmarks {

class zero	: public BaseInterface
{
public:
	zero()
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


	void clear()
	{
	}


	std::string getHelp()
	{
		std::ostringstream stream;
		stream << "  'zero':	Initialize every field with 0" << std::endl;
		return stream.str();
	}

	void getInitialState(
		sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &o_vrt,
		sweet::Data::Sphere2D::DataSpectral &o_div
	)
	{
		o_phi_pert.spectral_setZero();
		o_vrt.spectral_setZero();
		o_div.spectral_setZero();
	}



	void setup_topography() override
	{
		// initialize the topography
		shackPDESWEBenchmarks->topography.setup(ops->sphere2DDataConfig);
	}

};

}}

#endif
