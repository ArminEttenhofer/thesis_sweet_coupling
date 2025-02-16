/*
 * Sphere2DDiagnostics.hpp
 *
 *  Created on: 25 Feb 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_SPHERE2D_HELPERS_HELPERS_INTEGRAL_HPP
#define INCLUDE_SWEET_DATA_SPHERE2D_HELPERS_HELPERS_INTEGRAL_HPP

#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>


namespace sweet {
namespace Data {
namespace Sphere2D {


class Helpers_Integral
{
	const Sphere2D::Config *sphere2DDataConfig;

	/*
	 * Gaussian quadrature weights
	 */
	std::vector<double> gauss_weights;


public:
	Helpers_Integral()	:
		sphere2DDataConfig(nullptr)
	{

	}

	void setup(
			const Sphere2D::Config *i_sphere2DDataConfig,
			int i_verbose = 1
	)
	{
		sphere2DDataConfig = i_sphere2DDataConfig;

		gauss_weights.resize(sphere2DDataConfig->grid_num_lat);

		Sphere2D::DataSpectral modeSelector(sphere2DDataConfig);
		modeSelector.spectral_setZero();

		if (i_verbose > 0)
			std::cout << "Setting up Sphere2DDiagnostics" << std::endl;

		int n = shtns_gauss_wts(sphere2DDataConfig->shtns, gauss_weights.data());

		if (n*2 != sphere2DDataConfig->grid_num_lat)
		{
			std::cerr << "Returned " << n << " number of Gaussian quadrature point (halved)" << std::endl;
			SWEETErrorFatal("Wrong number of Gaussian quadrature points given!");
		}

		for (int i = 0; i < sphere2DDataConfig->grid_num_lat/2; i++)
			gauss_weights[sphere2DDataConfig->grid_num_lat-i-1] = gauss_weights[i];

		Sphere2D::DataSpectral modeIntegralValues(sphere2DDataConfig);

#if 0
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			if (i_verbose > 0)
				std::cout << "Mode m = " << m << " / " << sphere2DDataConfig->spectral_modes_m_max << std::endl;

			std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				double integral = 0;

#if !SWEET_DEBUG
				if (m != 0)
#endif
				{
					modeSelector.spectral_setZero();

					// activate mode
					modeSelector.spectral_space_data[idx] = 1.0;

					// conversion to physical space stores Gaussian quadrature weights
					modeSelector.request_data_physical();


					integral = modeSelector.grid_reduce_sum_quad();

					if (m != 0 && integral > 10e-12)
					{
						std::cout << n << ", " << m << ": Integral value expected to be close to zero, but is " << integral << std::endl;
						SWEETErrorFatal("Integral value not close to zero for m != 0");
					}
				}

				modeIntegralValues.spectral_space_data[idx] = integral;
				idx++;
			}
		}
#endif

		/*
		 * Test Gaussian quadrature
		 *
		 * Accurate for order (2n-1)
		 *
		 * test with y(x) = x^(2n-1)
		 */
		{
			int n = sphere2DDataConfig->grid_num_lat;
			double sum = 0;
			for (int jlat = 0; jlat < sphere2DDataConfig->grid_num_lat; jlat++)
			{
				double x = sphere2DDataConfig->lat_gaussian[jlat];
				double value = std::pow(x, 2.0*n-1.0) + std::pow(x, 2.0*n-2.0);

				sum += value*gauss_weights[jlat];
			}

			double accurate = 0;
			accurate += 1.0/(2.0*n)*(std::pow(1.0, 2.0*n)-std::pow(-1.0, 2.0*n));
			accurate += 1.0/(2.0*n-1.0)*(std::pow(1.0, 2.0*n-1.0)-std::pow(-1.0, 2.0*n-1.0));

			if (std::abs(accurate-sum) > 1e-10)
				SWEETErrorFatal("Error in quadrature");
		}
	}



#if 0

public:
	double compute_sphere2d_integral(
			const Sphere2D::DataGrid &i_data
	)	const
	{
		double sum = 0;

#if SPHERE2D_DATA_GRID_LAYOUT	== SPHERE2D_DATA_LAT_CONTIGUOUS
#error "TODO"
#else

		for (int jlat = 0; jlat < sphere2DDataConfig->grid_num_lat; jlat++)
		{
			for (int ilon = 0; ilon < sphere2DDataConfig->grid_num_lon; ilon++)
			{
				double value = i_data.grid_space_data[jlat*sphere2DDataConfig->grid_num_lon + ilon];

				value *= gauss_weights[jlat]*sphere2DDataConfig->lat_cogaussian[jlat];

				sum += value;
			}
		}
#endif

		sum /= (double)sphere2DDataConfig->grid_num_lon;

		sum *= (4.0/M_PI);

		return sum;
	}

#endif


public:
	/*
	 * The integral similar to the zylinder is used because of the lat-related scaling factor.
	 */
	double compute_zylinder_integral(
			const Sphere2D::DataSpectral &i_data
	)	const
	{
		Sphere2D::DataGrid data = i_data.toGrid();

		double sum = 0;

#if SPHERE2D_DATA_GRID_LAYOUT	== SPHERE2D_DATA_LAT_CONTIGUOUS
#error "TODO"
#else
		for (int jlat = 0; jlat < sphere2DDataConfig->grid_num_lat; jlat++)
		{
			for (int ilon = 0; ilon < sphere2DDataConfig->grid_num_lon; ilon++)
			{
				double value = data.grid_space_data[jlat*sphere2DDataConfig->grid_num_lon + ilon];

				sum += value*gauss_weights[jlat];
			}
		}
#endif
		sum /= (double)sphere2DDataConfig->grid_num_lon;

		sum *= 2.0*M_PI;

		return sum;
	}

};

}}}

#endif
