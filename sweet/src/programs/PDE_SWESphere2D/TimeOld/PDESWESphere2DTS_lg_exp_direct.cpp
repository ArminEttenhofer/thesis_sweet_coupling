/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Sphere2DComplex_DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/Convert/DataSpectral_2_Sphere2D_DataSpectral.hpp>
#include <sweet/ExpIntegration/REXI/REXI.hpp>
#include "PDESWESphere2DTS_lg_exp_direct.hpp"

#include <iostream>
#include <utility>
#include <sweet/Tools/StopwatchBox.hpp>



bool PDESWESphere2DTS_lg_exp_direct::implementsTimesteppingMethod(
		const std::string &i_timestepping_method
)
{
	timestepping_method = i_timestepping_method;

	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	if (i_timestepping_method == "lg_exp_direct")
		return true;

	return false;
}


bool PDESWESphere2DTS_lg_exp_direct::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup_main(ops, "phi0");
}


bool PDESWESphere2DTS_lg_exp_direct::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		const std::string &i_function_name
)
{
	ops = io_ops;
	sphere2DDataConfig = ops->sphere2DDataConfig;

	function_name = i_function_name;
	expFunctions.setup(i_function_name);

	if (function_name == "phi2")
	{
		SWEETErrorFatal("This phi function has a buggy implementation. Checkout new TimeTree to use it.");
	}

	return true;
}


std::string PDESWESphere2DTS_lg_exp_direct::getIDString()
{
	return "l_exp_direct";
}


void PDESWESphere2DTS_lg_exp_direct::runTimestep(
	const sweet::Data::Sphere2D::DataSpectral &i_prog_phi0,
	const sweet::Data::Sphere2D::DataSpectral &i_prog_vrt0,
	const sweet::Data::Sphere2D::DataSpectral &i_prog_div0,

	sweet::Data::Sphere2D::DataSpectral &o_prog_phi0,
	sweet::Data::Sphere2D::DataSpectral &o_prog_vrt0,
	sweet::Data::Sphere2D::DataSpectral &o_prog_div0,

	double i_fixed_dt,
	double i_simulation_timestamp
)
{
	o_prog_phi0 = i_prog_phi0;
	o_prog_vrt0 = i_prog_vrt0;
	o_prog_div0 = i_prog_div0;

	runTimestep(o_prog_phi0, o_prog_vrt0, o_prog_div0, i_fixed_dt, i_simulation_timestamp);
}



PDESWESphere2DTS_lg_exp_direct::PDESWESphere2DTS_lg_exp_direct()
{
}


PDESWESphere2DTS_lg_exp_direct::~PDESWESphere2DTS_lg_exp_direct()
{
}


void PDESWESphere2DTS_lg_exp_direct::runTimestep(
	sweet::Data::Sphere2D::DataSpectral &io_prog_phi,
	sweet::Data::Sphere2D::DataSpectral &io_prog_vrt,
	sweet::Data::Sphere2D::DataSpectral &io_prog_div,

	double i_dt,
	double i_simulation_timestamp
)
{
	SWEET_ASSERT(shackSphere2DDataOps != nullptr);

	#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().rexi.start();
		sweet::Tools::StopwatchBox::getInstance().rexi_timestepping.start();
	#endif


	/*
	 * Using exponential integrators, we must compute an
	 */

	double ir = 1.0/shackSphere2DDataOps->sphere2d_radius;
	// avg. geopotential

	double G = -shackPDESWESphere2D->h0*shackPDESWESphere2D->gravitation;

	/*
	 * See doc/rexi/rexi_for_swe_on_nonrotating_sphere2D.pdf
	 */

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
	{
		std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
		for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
		{
			double D = (double)n*((double)n+1.0)*ir*ir;

			if (D == 0)
			{
				idx++;
				continue;
			}

			std::complex<double> &phi0 = io_prog_phi.spectral_space_data[idx];
			std::complex<double> &div0 = io_prog_div.spectral_space_data[idx];

			// result will be imaginary only!
			std::complex<double> sqrt_DG = std::sqrt(std::complex<double>(D*G));

			// Multiply with Q^{-1}
			std::complex<double> l0 = -sqrt_DG/(2*G) * phi0 + 0.5*div0;
			std::complex<double> l1 = +sqrt_DG/(2*G) * phi0 + 0.5*div0;

			std::complex<double> tmp;
			expFunctions.eval(i_dt*(-sqrt_DG), tmp);
			l0 = tmp*l0;
			expFunctions.eval(i_dt*sqrt_DG, tmp);
			l1 = tmp*l1;

			phi0 = -G/sqrt_DG * l0 + G/sqrt_DG* l1;
			div0 = l0 + l1;

			idx++;
		}
	}

	#if SWEET_BENCHMARK_TIMINGS
		sweet::Tools::StopwatchBox::getInstance().rexi_timestepping.stop();
		sweet::Tools::StopwatchBox::getInstance().rexi.stop();
	#endif
}

