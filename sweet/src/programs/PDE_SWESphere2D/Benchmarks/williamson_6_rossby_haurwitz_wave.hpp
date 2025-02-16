/*
 * Author: Francois Hamon, Martin Schreiber <SchreiberX@Gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_WILLIAMSON_6_ROSSBY_HAURWITZ_WAVE_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_WILLIAMSON_6_ROSSBY_HAURWITZ_WAVE_HPP


#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "HelperGeostropicBalance.hpp"
#include "BaseInterface.hpp"

namespace PDE_SWESphere2D {
namespace Benchmarks {


class williamson_6_rossby_haurwitz_wave	:
		public BaseInterface
{
	HelperGeostropicBalance helperGeostropicBalance;

public:
	williamson_6_rossby_haurwitz_wave()
	{
	}

	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		BaseInterface::shackRegistration(io_shackDict);
		helperGeostropicBalance.shackRegistration(io_shackDict);
		return true;
	}


	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
				benchmark_name == "williamson6"	||
				benchmark_name == "rossby_haurwitz_wave" ||
				false
		;
	}


	void setup_1_shackData()
	{
		if (shackParallelization->isMPIRoot)
		{
			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;
		}

		//! Setup Williamson's parameters
		shackPDESWESphere2D->sphere2d_rotating_coriolis_omega = 7.292e-5;
		shackPDESWESphere2D->gravitation = 9.80616;
		shackSphere2DDataOps->sphere2d_radius = 6.37122e6;
		shackPDESWESphere2D->h0 = 8000;
	}

	void setup_2_withOps(
			sweet::Data::Sphere2D::Operators *io_ops
	)
	{
		ops = io_ops;

		helperGeostropicBalance.setup(ops);

	}


	void clear()
	{
	}



	std::string getHelp()
	{
		std::ostringstream stream;
		stream << "  WILLIAMSON #6:" << std::endl;
		stream << "     'williamson6'/'rossby_haurwitz_wave': Rossby Haurwitz wave" << std::endl;
		return stream.str();
	}



	void getInitialState(
		sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &o_vrt,
		sweet::Data::Sphere2D::DataSpectral &o_div
	)
	{
		double gh0 = shackPDESWESphere2D->gravitation*shackPDESWESphere2D->h0;

		const double omega = 7.484e-6;
		const double K = omega;
		const int R	= 4;
		const double a = shackSphere2DDataOps->sphere2d_radius;

		/*
		 * Setup U=...
		 */
		sweet::Data::Sphere2D::DataGrid ug(o_phi_pert.sphere2DDataConfig);
		ug.grid_update_lambda(
						[&](double lon, double phi, double &o_data)
							{
								o_data = a * omega * cos(phi)
										+ a * K * pow(cos(phi), R-1) * (R * sin(phi)*sin(phi) - cos(phi)*cos(phi)) * cos(R*lon);
							}
						);


		/*
		 * Setup V=...
		 */
		sweet::Data::Sphere2D::DataGrid vg(o_phi_pert.sphere2DDataConfig);
		vg.grid_update_lambda(
						  [&](double lon, double phi, double &o_data)
								{
									o_data = - a * K * R * pow(cos(phi), R-1) * sin(phi) * sin(R*lon);
								}
						  );

		ops->uv_2_vrtdiv(ug, vg, o_vrt, o_div);

		sweet::Data::Sphere2D::DataGrid hg(o_phi_pert.sphere2DDataConfig);

		helperGeostropicBalance.computeGeostrophicBalance_nonlinear(
				o_vrt,
				o_div,
				o_phi_pert
		);

		o_phi_pert -= gh0;
	}



	void setup_topography() override
	{
		// initialize the topography
		shackPDESWEBenchmarks->topography.setup(ops->sphere2DDataConfig);
	}

};

}}

#endif
