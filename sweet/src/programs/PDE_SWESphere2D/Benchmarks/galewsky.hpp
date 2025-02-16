/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_GALEWSKY_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_GALEWSKY_HPP

#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/LibMath/GaussQuadrature.hpp>
#include <sweet/Shacks/Dictionary.hpp>

#include "HelperGeostropicBalance.hpp"
#include "BaseInterface.hpp"

namespace PDE_SWESphere2D {
namespace Benchmarks {

class galewsky	:
		public BaseInterface
{
	HelperGeostropicBalance helperGeostropicBalance;

public:
	galewsky()
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
			i_benchmark_name == "galewsky" ||			//!< Standard Galewsky benchmark
			i_benchmark_name == "galewsky_linearbalance" ||	//!< Standard Galewsky benchmark with linear balanced initial conditions
			i_benchmark_name == "galewsky_nobump" ||	//!< Galewsky benchmark without bumps
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

		// Setup Galewski parameters
		shackPDESWESphere2D->sphere2d_rotating_coriolis_omega = 7.292e-5;
		shackPDESWESphere2D->gravitation = 9.80616;
		shackSphere2DDataOps->sphere2d_radius = 6.37122e6;

		/*
		 * This is just the h0 parameter for the PDE solver,
		 * not the h0 parameter in the Galwesky paper which
		 * is required to set the average height to 10000.
		 */
		shackPDESWESphere2D->h0 = 10000;
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
		helperGeostropicBalance.clear();
	}

	std::string getHelp()
	{
		std::ostringstream stream;
		stream << "  'galewsky': Galwesky benchmark" << std::endl;
		stream << "  'galewsky_nobump': Galwesky benchmark without any bump" << std::endl;
		return stream.str();
	}

	void getInitialState(
		sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &o_vrt,
		sweet::Data::Sphere2D::DataSpectral &o_div
	)
	{
		const sweet::Data::Sphere2D::Config *sphere2DDataConfig = o_phi_pert.sphere2DDataConfig;

		// Search for substrings
		bool benchmark_nobump = benchmark_name.find("nobump") != std::string::npos;
		bool benchmark_linearbalance = benchmark_name.find("linearbalance") != std::string::npos;

		bool use_analytical_geostrophic_setup;
		if (shackPDESWEBenchmarks->benchmark_galewsky_geostrophic_setup == "analytical")
		{
			use_analytical_geostrophic_setup = true;
		}
		else if (shackPDESWEBenchmarks->benchmark_galewsky_geostrophic_setup == "numerical")
		{
			use_analytical_geostrophic_setup = false;
		}
		else
		{
			SWEETErrorFatal("Invalid geostropic setup choosen");
			use_analytical_geostrophic_setup = false;	// Make compiler happy
		}



		/*
		 * Parameters from Galewsky paper setup
		 */
		double a = shackSphere2DDataOps->sphere2d_radius;
		double omega = shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;
		double umax = 80.;
		double phi0 = M_PI/7.;
		double phi1 = 0.5*M_PI - phi0;
		double phi2 = 0.25*M_PI;		//! latitude placement of gaussian bump
		double en = std::exp(-4.0/std::pow((phi1-phi0), 2.0));
		double alpha = 1./3.;
		double beta = 1./15.;
		double hamp = 120.;



		if (shackPDESWEBenchmarks->benchmark_galewsky_umax >= 0)
			umax = shackPDESWEBenchmarks->benchmark_galewsky_umax;

		if (shackPDESWEBenchmarks->benchmark_galewsky_hamp >= 0)
			hamp = shackPDESWEBenchmarks->benchmark_galewsky_hamp;

		if (shackPDESWEBenchmarks->benchmark_galewsky_phi2 >= 0)
			phi2 = shackPDESWEBenchmarks->benchmark_galewsky_phi2;

		/*
		 * Setup V=0
		 */
		sweet::Data::Sphere2D::DataGrid vg(o_phi_pert.sphere2DDataConfig);
		vg.grid_setZero();

		auto lambda_u = [&](double phi) -> double
		{
			if (phi >= phi1-1e-5 || phi <= phi0+1e-5)
				return 0.0;
			else
				return umax/en*std::exp(1.0/((phi-phi0)*(phi-phi1)));
		};

		auto lambda_f = [&](double phi) -> double
		{
			return a*lambda_u(phi)*(2.0*omega*std::sin(phi)+(std::tan(phi)/a)*lambda_u(phi));
		};

		/*
		 * Setup U=...
		 * initial velocity along longitude
		 */
		sweet::Data::Sphere2D::DataGrid ug(o_phi_pert.sphere2DDataConfig);
		ug.grid_update_lambda(
			[&](double lon, double phi, double &o_data)
			{
				o_data = lambda_u(phi);
			}
		);

		ops->uv_2_vrtdiv(ug, vg, o_vrt, o_div);

		if (use_analytical_geostrophic_setup)
		{
			if (shackParallelization->isMPIRoot)
				std::cout << "[MULE] use_analytical_geostrophic_setup: 1" << std::endl;

			if (!benchmark_linearbalance)
			{
				// use nonlinear balanced initial conditions
				helperGeostropicBalance.computeGeostrophicBalance_nonlinear(
						o_vrt,
						o_div,
						o_phi_pert
				);
			}
			else
			{
				// use linearly balanced initial conditions
				helperGeostropicBalance.computeGeostrophicBalance_linear(
						o_vrt,
						o_div,
						o_phi_pert
				);
			}

#if 0
			double h0_ = 10e3;
			o_phi_pert = shackPDESWESphere2D->gravitation * h0_ + o_phi_pert;
			o_phi_pert -= shackPDESWESphere2D->gravitation*shackPDESWESphere2D->h0;
#endif
		}
		else
		{
			if (benchmark_linearbalance)
			{
				SWEETErrorFatal("Not supported");
			}

			if (shackParallelization->isMPIRoot)
				std::cout << "[MULE] use_analytical_geostrophic_setup: 0" << std::endl;

			/*
			 * Initialization of SWE height
			 *
			 * Metric correction terms based on John Thuburn's code
			 */
#if 1
			const unsigned short nlat = sphere2DDataConfig->grid_num_lat;
			std::vector<double> hg_cached;
			hg_cached.resize(nlat);

			double h_metric_area = 0;
			//double hg_sum = 0;
			double int_start, int_end, int_delta;

			int j = sphere2DDataConfig->grid_num_lat-1;


			// start/end of first integration interval
			{
				SWEET_ASSERT(sphere2DDataConfig->lat[j] < 0);

				// start at the south pole
				int_start = -M_PI*0.5;

				// first latitude gaussian point
				int_end = sphere2DDataConfig->lat[j];

				// 1d area of integration
				int_delta = int_end - int_start;

				SWEET_ASSERT(int_delta > 0);
				SWEET_ASSERT(int_delta < 1);

				double hg = sweet::LibMath::GaussQuadrature::integrate5_intervals<double>(int_start, int_end, lambda_f, 20);
				//hg = (int_end+int_start)*0.5;
				hg_cached[j] = hg;

				/*
				 * cos scaling is required for 2D sphere2D coverage at this latitude
				 *
				 * metric term which computes the area coverage of each point
				 */
				// use integrated average as below instead of the following formulation
				// double mterm = cos((int_start+int_end)*0.5);
				//double mterm = (std::sin(int_end)-std::sin(int_start))*2.0*M_PI;
				double mterm = std::cos(sphere2DDataConfig->lat[j])*2.0*M_PI;
				SWEET_ASSERT(mterm > 0);

				h_metric_area += mterm;

				int_start = int_end;
			}
			j--;

			for (; j >= 0; j--)
			{
				double int_end = sphere2DDataConfig->lat[j];
				int_delta = int_end - int_start;
				SWEET_ASSERT(int_delta > 0);

				double hg = hg_cached[j+1] + sweet::LibMath::GaussQuadrature::integrate5_intervals<double>(int_start, int_end, lambda_f, 20);

				//hg = (int_end+int_start)*0.5;
				hg_cached[j] = hg;

				// metric term which computes the area coverage of each point
				//double mterm = (std::sin(int_end)-std::sin(int_start))*2.0*M_PI;
				double mterm = std::cos(sphere2DDataConfig->lat[j])*2.0*M_PI;

				//hg_sum += hg*mterm;
				h_metric_area += mterm;

				// continue at the end of the last integration interval
				int_start = int_end;
			}

			// last integration interval
			{
				SWEET_ASSERT(int_start > 0);
				int_end = M_PI*0.5;

				int_delta = int_end - int_start;
				SWEET_ASSERT(int_delta > 0);

				// metric term which computes the area coverage of each point
				//double mterm = (std::sin(int_end)-std::sin(int_start))*2.0*M_PI;
				double mterm = std::cos(sphere2DDataConfig->lat[0])*2.0*M_PI;

				//double hg = hg_cached[0] + sweet::LibMath::GaussQuadrature::integrate5_intervals<double>(int_start, int_end, lambda_f, 20);
				//hg = (int_end+int_start)*0.5;
				//hg_sum += hg*mterm;
				h_metric_area += mterm;
			}

			SWEET_ASSERT(h_metric_area > 0);

#else

			std::vector<double> hg_cached;
			hg_cached.resize(sphere2DDataConfig->grid_num_lat);

			double int_start = -M_PI*0.5;
			for (int j = sphere2DDataConfig->grid_num_lat-1; j >= 0; j--)
			{
				double int_end = sphere2DDataConfig->lat[j];
				double quad = sweet::LibMath::GaussQuadrature::integrate5_intervals<double>(int_start, int_end, lambda_f, 5);

				if (j == sphere2DDataConfig->grid_num_lat-1)
					hg_cached[j] = quad;
				else
					hg_cached[j] = hg_cached[j+1] + quad;

				int_start = int_end;


				std::cout << sphere2DDataConfig->lat[j] << ": " << hg_cached[j] << std::endl;
			}

#endif

			// update data
			sweet::Data::Sphere2D::DataGrid phig(sphere2DDataConfig);
			phig.grid_update_lambda_array(
				[&](int i, int j, double &o_data)
				{
					o_data = hg_cached[j];
				}
			);

			o_phi_pert.loadSphere2DDataGrid(phig);

			o_phi_pert = -o_phi_pert;
		}

		/*
		 * Now change global mean layer depth to 10km
		 *
		 * From Galewsky et al. paper:
		 * "and the constant h0 is chosen so that the global mean layer depth is equal to 10 km"
		 */

		// The average is given by the very first mode
		// magic constant which accomplishes that
		//double avg_h = 10158.186170454619;
		//double avg_h = 10000;
		//o_phi_pert += (avg_h - shackPDESWESphere2D->h0)*shackPDESWESphere2D->gravitation;

		// Using h0=10000 we get h0 in average as well if we set the first mode to 0
		o_phi_pert.spectral_set(0, 0, 0);

		sweet::Data::Sphere2D::DataGrid hbump(o_phi_pert.sphere2DDataConfig);
		if (!benchmark_nobump)
		{
			hbump.grid_update_lambda(
				[&](double lon, double phi, double &o_data)
				{
					o_data = hamp*std::cos(phi)*std::exp(-std::pow((lon-M_PI)/alpha, 2.0))*std::exp(-std::pow((phi2-phi)/beta, 2.0));
				}
			);
			o_phi_pert += hbump*shackPDESWESphere2D->gravitation;
		}


		if (shackParallelization->isMPIRoot)
		{
			std::cout << "phi min: " << o_phi_pert.toGrid().grid_reduce_min() << std::endl;
			std::cout << "phi max: " << o_phi_pert.toGrid().grid_reduce_max() << std::endl;
		}

	}

	void setup_topography() override
	{
		// initialize the topography
		shackPDESWEBenchmarks->topography.setup(ops->sphere2DDataConfig);
	}

};

}}

#endif
