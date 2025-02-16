/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWECART2D_PDESWECART2D_BENCHMARKSCOMBINED_HPP
#define PROGRAMS_PDE_SWECART2D_PDESWECART2D_BENCHMARKSCOMBINED_HPP

#include <programs/PDE_SWECart2D/Benchmarks/Shack_Polvani.hpp>
#include <programs/PDE_SWECart2D/Benchmarks/Shack.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2D/Operators.hpp>
#include <iostream>
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>


#if SWEET_USE_CART2D_SPECTRAL_SPACE
	#include <programs/PDE_SWECart2D/Benchmarks/polvani.hpp>
	#include <programs/PDE_SWECart2D/Benchmarks/merge_vortex.hpp>
	#include <programs/PDE_SWECart2D/Benchmarks/normal_modes.hpp>
#endif

#include <programs/PDE_SWECart2D/Benchmarks/unstable_jet.hpp>
#include <programs/PDE_SWECart2D/Benchmarks/unstable_jet_fast.hpp>
#include <programs/PDE_SWECart2D/Benchmarks/unstable_jet_adv.hpp>
#include <programs/PDE_SWECart2D/Benchmarks/gaussian_bump.hpp>
#include <programs/PDE_SWECart2D/Benchmarks/column.hpp>
#include <programs/PDE_SWECart2D/Benchmarks/impulse.hpp>
#include <programs/PDE_SWECart2D/Benchmarks/empty.hpp>


#include <programs/PDE_SWECart2D/Shack.hpp>


#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/TimeTree/Shack.hpp>

namespace PDE_SWECart2D {
namespace Benchmarks {

class BenchmarksCombined
{
public:
	sweet::Error::Base error;

	// cart2d or sphere2D data config
	const void* ext_forces_data_config;

	sweet::Shacks::Dictionary *shackDict;
	sweet::Data::Cart2D::Shack *shackCart2DDataOps;
	sweet::TimeTree::Shack *shackTimestepControl;
	PDE_SWECart2D::Shack *shackPDESWECart2D;
	PDE_SWECart2D::Benchmarks::Shack *shackPDESWECart2DBenchmarks;
	PDE_SWECart2D::Benchmarks::Shack_Polvani *shackPDESWECart2DBench_polvaniBenchmark;

	BenchmarksCombined()	:
		shackDict(nullptr),
		shackCart2DDataOps(nullptr),
		shackTimestepControl(nullptr),
		shackPDESWECart2D(nullptr),
		shackPDESWECart2DBenchmarks(nullptr)
	{
	}


	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackCart2DDataOps = shackDict->getAutoRegistration<sweet::Data::Cart2D::Shack>();
		shackTimestepControl = shackDict->getAutoRegistration<sweet::TimeTree::Shack>();
		shackPDESWECart2D = shackDict->getAutoRegistration<PDE_SWECart2D::Shack>();
		shackPDESWECart2DBenchmarks = shackDict->getAutoRegistration<PDE_SWECart2D::Benchmarks::Shack>();
		shackPDESWECart2DBench_polvaniBenchmark = shackDict->getAutoRegistration<PDE_SWECart2D::Benchmarks::Shack_Polvani>();

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}


	bool shackRegistration(
			sweet::Shacks::Dictionary &io_shackDict
	)
	{
		return shackRegistration(&io_shackDict);
	}


	bool clear()
	{
		// TODO
		return true;
	}


public:
	bool setupInitialConditions(
			sweet::Data::Cart2D::DataSpectral &o_h_pert,
			sweet::Data::Cart2D::DataSpectral &o_u,
			sweet::Data::Cart2D::DataSpectral &o_v,
			sweet::Data::Cart2D::Operators *io_ops,				//!< Make this IO, since changes in the simulation parameters might require to also update the operators
			sweet::Data::Cart2D::Config *io_cart2DDataConfig
	)
	{
		SWEET_ASSERT(io_ops != nullptr);
		SWEET_ASSERT(io_cart2DDataConfig != nullptr);

		auto callback_gaussian_bump =
				[&](
						double i_center_x, double i_center_y,
						double i_x, double i_y,
						double i_exp_fac
				)
				{
					double sx = shackCart2DDataOps->cart2d_domain_size[0];
					double sy = shackCart2DDataOps->cart2d_domain_size[1];

					// Gaussian
					double dx = i_x-i_center_x*sx;
					double dy = i_y-i_center_y*sy;

					if (dx > 0.5*shackCart2DDataOps->cart2d_domain_size[0])
						dx -= shackCart2DDataOps->cart2d_domain_size[0];
					else if (dx < -0.5*shackCart2DDataOps->cart2d_domain_size[0])
						dx += shackCart2DDataOps->cart2d_domain_size[0];

					if (dy > 0.5*shackCart2DDataOps->cart2d_domain_size[1])
						dy -= shackCart2DDataOps->cart2d_domain_size[1];
					else if (dy < -0.5*shackCart2DDataOps->cart2d_domain_size[1])
						dy += shackCart2DDataOps->cart2d_domain_size[1];

					dx /= sx*shackPDESWECart2DBenchmarks->object_scale;
					dy /= sy*shackPDESWECart2DBenchmarks->object_scale;

					return std::exp(-i_exp_fac*(dx*dx + dy*dy));
				};

		if (shackPDESWECart2DBenchmarks->benchmark_name == "")
			return error.set("SWECart2DBenchmarksCombined: Benchmark name not given, use --benchmark-name=[name]");

#if SWEET_USE_CART2D_SPECTRAL_SPACE
		if (shackPDESWECart2DBenchmarks->benchmark_name == "polvani")
		{
			PDE_SWECart2D::Benchmarks::polvani swe_polvani;
			swe_polvani.shackRegistration(shackDict);
			swe_polvani.setup(io_ops, io_cart2DDataConfig);
			swe_polvani.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
				);

			return true;
		}
		else if (shackPDESWECart2DBenchmarks->benchmark_name == "mergevortex")
		{
			PDE_SWECart2D::Benchmarks::merge_vortex swe_mergevortex;
			swe_mergevortex.shackRegistration(shackDict);
			swe_mergevortex.setup(io_ops, io_cart2DDataConfig);
			swe_mergevortex.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
				);

			return true;
		}
		else if (shackPDESWECart2DBenchmarks->benchmark_name == "unstablejet")
		{
			shackCart2DDataOps->cart2d_domain_size[0] = 10e6;
			shackCart2DDataOps->cart2d_domain_size[1] = 10e6;
			shackPDESWECart2D->h0 = 1000;

//			double r = 6.37122e6;
//			shackCart2DDataOps->cart2d_domain_size[0] = 2.0*M_PI*r;
//			shackCart2DDataOps->cart2d_domain_size[1] = 2.0*M_PI*r;
//            shackPDESWECart2D->h0 = 10000;
            shackPDESWECart2D->gravitation = 9.80616;
            shackPDESWECart2D->cart2d_rotating_f0 = 0.00014584;

			std::cout << "WARNING: OVERWRITING SIMULATION PARAMETERS FOR THIS BENCHMARK!" << std::endl;
			io_ops->clear();
			io_ops->setup(io_cart2DDataConfig, shackCart2DDataOps);


			PDE_SWECart2D::Benchmarks::unstable_jet swe_unstablejet;
			swe_unstablejet.shackRegistration(shackDict);
			swe_unstablejet.setup(io_ops, io_cart2DDataConfig);
			swe_unstablejet.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (shackPDESWECart2DBenchmarks->benchmark_name == "unstablejet_nobump")
		{
            shackCart2DDataOps->cart2d_domain_size[0] = 1e7;
            shackCart2DDataOps->cart2d_domain_size[1] = 1e7;
            shackPDESWECart2D->h0 = 1000;

//			double r = 6.37122e6;
//			shackCart2DDataOps->cart2d_domain_size[0] = 2.0*M_PI*r;
//			shackCart2DDataOps->cart2d_domain_size[1] = 2.0*M_PI*r;
//            shackPDESWECart2D->h0 = 10000;
            shackPDESWECart2D->gravitation = 9.80616;
            shackPDESWECart2D->cart2d_rotating_f0 = 0.00014584;


            std::cout << "WARNING: OVERWRITING SIMULATION PARAMETERS FOR THIS BENCHMARK!" << std::endl;
            io_ops->clear();
            io_ops->setup(io_cart2DDataConfig, shackCart2DDataOps);


			PDE_SWECart2D::Benchmarks::unstable_jet swe_unstablejet(false);
			swe_unstablejet.shackRegistration(shackDict);
			swe_unstablejet.setup(io_ops, io_cart2DDataConfig);
			swe_unstablejet.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
        else if (shackPDESWECart2DBenchmarks->benchmark_name == "column")
        {
            shackCart2DDataOps->cart2d_domain_size[0] = 1e5;
            shackCart2DDataOps->cart2d_domain_size[1] = 1e5;
            shackPDESWECart2D->h0 = 1000;

            shackPDESWECart2D->gravitation = 9.80616;
            shackPDESWECart2D->cart2d_rotating_f0 = 0;

            std::cout << "WARNING: OVERWRITING SIMULATION PARAMETERS FOR THIS BENCHMARK!" << std::endl;
            io_ops->clear();
            io_ops->setup(io_cart2DDataConfig, shackCart2DDataOps);


            PDE_SWECart2D::Benchmarks::column swe_column{};
            swe_column.shackRegistration(shackDict);
            swe_column.setup(io_ops, io_cart2DDataConfig);
            swe_column.setupBenchmark(
                    o_h_pert,
                    o_u,
                    o_v
            );

            return true;
        }
        else if (shackPDESWECart2DBenchmarks->benchmark_name == "impulse")
        {
            shackCart2DDataOps->cart2d_domain_size[0] = 1e7;
            shackCart2DDataOps->cart2d_domain_size[1] = 1e7;
            shackPDESWECart2D->h0 = 1000;

            shackPDESWECart2D->gravitation = 9.80616;
            shackPDESWECart2D->cart2d_rotating_f0 = 0;

            std::cout << "WARNING: OVERWRITING SIMULATION PARAMETERS FOR THIS BENCHMARK!" << std::endl;
            io_ops->clear();
            io_ops->setup(io_cart2DDataConfig, shackCart2DDataOps);


            PDE_SWECart2D::Benchmarks::impulse swe_impulse{};
            swe_impulse.shackRegistration(shackDict);
            swe_impulse.setup(io_ops, io_cart2DDataConfig);
            swe_impulse.setupBenchmark(
                    o_h_pert,
                    o_u,
                    o_v
            );

            return true;
        }
        else if (shackPDESWECart2DBenchmarks->benchmark_name == "empty")
        {
            shackCart2DDataOps->cart2d_domain_size[0] = 1e5;
            shackCart2DDataOps->cart2d_domain_size[1] = 1e5;
            shackPDESWECart2D->h0 = 5000;

            shackPDESWECart2D->gravitation = 9.80616;
            shackPDESWECart2D->cart2d_rotating_f0 = 0;

            std::cout << "WARNING: OVERWRITING SIMULATION PARAMETERS FOR THIS BENCHMARK!" << std::endl;
            io_ops->clear();
            io_ops->setup(io_cart2DDataConfig, shackCart2DDataOps);


            PDE_SWECart2D::Benchmarks::empty swe_empty{};
            swe_empty.shackRegistration(shackDict);
            swe_empty.setup(io_ops, io_cart2DDataConfig);
            swe_empty.setupBenchmark(
                    o_h_pert,
                    o_u,
                    o_v
            );

            return true;
        }
		else if (shackPDESWECart2DBenchmarks->benchmark_name == "unstablejetfast")
		{
			PDE_SWECart2D::Benchmarks::unstable_jet_fast swe_unstablejetfast;
			swe_unstablejetfast.shackRegistration(shackDict);
			swe_unstablejetfast.setup(io_ops, io_cart2DDataConfig);
			swe_unstablejetfast.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (shackPDESWECart2DBenchmarks->benchmark_name == "unstablejetadv")
		{
			PDE_SWECart2D::Benchmarks::unstable_jet_adv swe_unstablejetadv;
			swe_unstablejetadv.shackRegistration(shackDict);
			swe_unstablejetadv.setup(io_ops, io_cart2DDataConfig);
			swe_unstablejetadv.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (shackPDESWECart2DBenchmarks->benchmark_name == "normalmodes")
		{
			PDE_SWECart2D::Benchmarks::normal_modes swe_normalmodes;
			swe_normalmodes.shackRegistration(shackDict);
			swe_normalmodes.setup(io_ops, io_cart2DDataConfig);
			swe_normalmodes.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
#endif
		else if (shackPDESWECart2DBenchmarks->benchmark_name == "gaussian_bump" || shackPDESWECart2DBenchmarks->benchmark_name == "gaussian_bump_phi_pint")
		{
			PDE_SWECart2D::Benchmarks::gaussian_bump swe_gaussian_bump;
			swe_gaussian_bump.shackRegistration(shackDict);
			swe_gaussian_bump.setup(io_ops, io_cart2DDataConfig);
			swe_gaussian_bump.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (shackPDESWECart2DBenchmarks->benchmark_name == "gaussian_bump_advection")
		{

			auto callback_external_forces_advection_field =
					[](
							int i_field_id,
							double i_simulation_timestamp,
							sweet::Data::Cart2D::DataSpectral* io_data,			//! cart2ddata or sphere2Ddata
							PDE_SWECart2D::Benchmarks::Shack* i_shackBenchmark		//! user data (pointer to this class)
			)
			{
				sweet::Data::Cart2D::DataGrid cart2d_data_phys(io_data->cart2DDataConfig);

				if (i_field_id >= 1 && i_field_id <= 2)
				{
					double u = i_shackBenchmark->advection_velocity[0];
					double v = i_shackBenchmark->advection_velocity[1];

					double r;
					if (i_shackBenchmark->advection_velocity[2] == 0)
						r = 0;
					else
						r = i_simulation_timestamp/i_shackBenchmark->advection_velocity[2]*2.0*M_PI;

					if (i_field_id == 1)
					{
						// u-velocity
						//*o_cart2d_data = std::cos(r)*u - std::sin(r)*v;
						cart2d_data_phys = u*(1.0+std::sin(r));
					}
					else if (i_field_id == 2)
					{
						// v-velocity
						//*o_cart2d_data = std::sin(r)*u + std::cos(r)*v;
						cart2d_data_phys = v*(1.0+std::cos(r));
					}

					io_data->loadCart2DDataGrid(cart2d_data_phys);

					return;
				}

				SWEETErrorFatal("Non-existing external field requested!");
				return;
			};

			if (shackPDESWECart2DBenchmarks->advection_velocity[2] != 0)
			{
				// backup data config
				ext_forces_data_config = o_h_pert.cart2DDataConfig;

				// set callback
				shackPDESWECart2DBenchmarks->getExternalForcesCallback = callback_external_forces_advection_field;

				// set user data to this class
				shackPDESWECart2DBenchmarks->getExternalForcesUserData = shackPDESWECart2DBenchmarks;

				// setup velocities with initial time stamp
				callback_external_forces_advection_field(1, shackTimestepControl->currentSimulationTime, &o_u, shackPDESWECart2DBenchmarks->getExternalForcesUserData);
				callback_external_forces_advection_field(2, shackTimestepControl->currentSimulationTime, &o_v, shackPDESWECart2DBenchmarks->getExternalForcesUserData);
			}
			else
			{
				sweet::Data::Cart2D::DataGrid u_phys(o_u.cart2DDataConfig);
				sweet::Data::Cart2D::DataGrid v_phys(o_u.cart2DDataConfig);

				u_phys = shackPDESWECart2DBenchmarks->advection_velocity[0];
				v_phys = shackPDESWECart2DBenchmarks->advection_velocity[1];

				o_u.loadCart2DDataGrid(u_phys);
				o_v.loadCart2DDataGrid(v_phys);
			}

			double center_x = 0.5;
			double center_y = 0.5;
			double exp_fac = 50.0;

			sweet::Data::Cart2D::DataGrid h_pert_phys(o_h_pert.cart2DDataConfig);
			h_pert_phys.grid_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackCart2DDataOps->cart2d_domain_size[0]/(double)shackCart2DDataOps->space_res_physical[0]);
					double y = (double)j*(shackCart2DDataOps->cart2d_domain_size[1]/(double)shackCart2DDataOps->space_res_physical[1]);

					io_data = callback_gaussian_bump(center_x, center_y, x, y, exp_fac);
				}
			);
			o_h_pert.loadCart2DDataGrid(h_pert_phys);

			return true;
		}
		else if (
				shackPDESWECart2DBenchmarks->benchmark_name == "benchmark_id_0" ||
				shackPDESWECart2DBenchmarks->benchmark_name == "cylinder"
		)
		{
			double sx = shackCart2DDataOps->cart2d_domain_size[0];
			double sy = shackCart2DDataOps->cart2d_domain_size[1];


			sweet::Data::Cart2D::DataGrid h_pert_phys(o_h_pert.cart2DDataConfig);
			h_pert_phys.grid_setZero();
			h_pert_phys.grid_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackCart2DDataOps->cart2d_domain_size[0]/(double)shackCart2DDataOps->space_res_physical[0]);
					double y = (double)j*(shackCart2DDataOps->cart2d_domain_size[1]/(double)shackCart2DDataOps->space_res_physical[1]);

					// radial dam break
					double dx = x-shackPDESWECart2DBenchmarks->object_coord_x*sx;
					double dy = y-shackPDESWECart2DBenchmarks->object_coord_y*sy;

					double radius = shackPDESWECart2DBenchmarks->object_scale*sqrt(sx*sx+sy*sy);
					if (dx*dx+dy*dy < radius*radius)
						io_data = 1.0;
					else
						io_data = 0.0;
				}
			);

			o_h_pert.loadCart2DDataGrid(h_pert_phys);
			o_u.spectral_setZero();
			o_v.spectral_setZero();

			return true;
		}
		else if (
				shackPDESWECart2DBenchmarks->benchmark_name == "benchmark_id_1" ||
				shackPDESWECart2DBenchmarks->benchmark_name == "radial_gaussian_bump"
		)
		{
			double sx = shackCart2DDataOps->cart2d_domain_size[0];
			double sy = shackCart2DDataOps->cart2d_domain_size[1];


			sweet::Data::Cart2D::DataGrid h_pert_phys(o_h_pert.cart2DDataConfig);
			h_pert_phys.grid_setZero();
			h_pert_phys.grid_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackCart2DDataOps->cart2d_domain_size[0]/(double)shackCart2DDataOps->space_res_physical[0]);
					double y = (double)j*(shackCart2DDataOps->cart2d_domain_size[1]/(double)shackCart2DDataOps->space_res_physical[1]);

					// radial dam break
					double dx = x-shackPDESWECart2DBenchmarks->object_coord_x*sx;
					double dy = y-shackPDESWECart2DBenchmarks->object_coord_y*sy;

					double radius = shackPDESWECart2DBenchmarks->object_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
					dx /= radius;
					dy /= radius;

					io_data = std::exp(-50.0*(dx*dx + dy*dy));
				}
			);

			o_h_pert.loadCart2DDataGrid(h_pert_phys);
			o_u.spectral_setZero();
			o_v.spectral_setZero();

			return true;
		}
		else if (
				shackPDESWECart2DBenchmarks->benchmark_name == "benchmark_id_2" ||
				shackPDESWECart2DBenchmarks->benchmark_name == "steady_state_meridional_flow"
		)
		{
			double f = shackPDESWECart2D->cart2d_rotating_f0;
			double sx = shackCart2DDataOps->cart2d_domain_size[0];
			//double sy = shackSim->domain_size[1];

			if (shackPDESWECart2D->cart2d_rotating_f0 == 0)
				SWEETErrorFatal("Coriolis = 0!");

			sweet::Data::Cart2D::DataGrid h_pert_phys(o_h_pert.cart2DDataConfig);
			h_pert_phys.grid_setZero();
			h_pert_phys.grid_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i/(double)shackCart2DDataOps->space_res_physical[0];
					//double y = (double)j*(shackSim->domain_size[1]/(double)disc->res_physical[1]);

					io_data = std::sin(2.0*M_PI*x);
				}
			);

			o_h_pert.loadCart2DDataGrid(h_pert_phys);

			o_u.spectral_setZero();

			sweet::Data::Cart2D::DataGrid v_phys(o_v.cart2DDataConfig);
			v_phys.grid_setZero();
			v_phys.grid_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i/(double)shackCart2DDataOps->space_res_physical[0];
					//double y = (double)j*(shackSim->domain_size[1]/(double)disc->res_physical[1]);

					io_data = shackPDESWECart2D->gravitation/f*2.0*M_PI*std::cos(2.0*M_PI*x)/sx;
				}
			);

			o_v.loadCart2DDataGrid(v_phys);

			return true;
		}
		else if (
				shackPDESWECart2DBenchmarks->benchmark_name == "benchmark_id_3" ||
				shackPDESWECart2DBenchmarks->benchmark_name == "steady_state_zonal_flow"
		)
		{
			double f = shackPDESWECart2D->cart2d_rotating_f0;
			//double sx = shackSim->domain_size[0];
			double sy = shackCart2DDataOps->cart2d_domain_size[1];

			if (shackPDESWECart2D->cart2d_rotating_f0 == 0)
				SWEETErrorFatal("Coriolis = 0!");

			sweet::Data::Cart2D::DataGrid h_pert_phys(o_h_pert.cart2DDataConfig);
			h_pert_phys.grid_setZero();
			h_pert_phys.grid_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double y = (double)j*(shackCart2DDataOps->cart2d_domain_size[1]/(double)shackCart2DDataOps->space_res_physical[1]);

					io_data = std::sin(2.0*M_PI*y/sy);
				}
			);

			o_h_pert.loadCart2DDataGrid(h_pert_phys);

			sweet::Data::Cart2D::DataGrid u_phys(o_u.cart2DDataConfig);
			u_phys.grid_setZero();
			u_phys.grid_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double y = (double)j*(shackCart2DDataOps->cart2d_domain_size[1]/(double)shackCart2DDataOps->space_res_physical[1]);

					io_data = -shackPDESWECart2D->gravitation*2.0*M_PI*std::cos(2.0*M_PI*y/sy)/(f*sy);
				}
			);

			o_u.loadCart2DDataGrid(u_phys);

			o_v.spectral_setZero();

			return true;
		}
		else if (
				shackPDESWECart2DBenchmarks->benchmark_name == "benchmark_id_4" ||
				shackPDESWECart2DBenchmarks->benchmark_name == "yadda_yadda_whatever_this_is"
		)
		{
			double sx = shackCart2DDataOps->cart2d_domain_size[0];
			double sy = shackCart2DDataOps->cart2d_domain_size[1];

			if (shackPDESWECart2D->cart2d_rotating_f0 == 0)
				SWEETErrorFatal("Coriolis = 0!");

			sweet::Data::Cart2D::DataGrid h_pert_phys(o_h_pert.cart2DDataConfig);
			h_pert_phys.grid_setZero();
			h_pert_phys.grid_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackCart2DDataOps->cart2d_domain_size[0]/(double)shackCart2DDataOps->space_res_physical[0]);
					double y = (double)j*(shackCart2DDataOps->cart2d_domain_size[1]/(double)shackCart2DDataOps->space_res_physical[1]);

					// radial dam break
					double dx = x-shackPDESWECart2DBenchmarks->object_coord_x*sx;
					double dy = y-shackPDESWECart2DBenchmarks->object_coord_y*sy;

					io_data = (std::abs(dx-0.5) < 0.3)*(std::abs(dy-0.5) < 0.1);
				}
			);

			o_h_pert.loadCart2DDataGrid(h_pert_phys);
			o_u.spectral_setZero();
			o_v.spectral_setZero();

			return true;
		}


		else if (
				shackPDESWECart2DBenchmarks->benchmark_name == "benchmark_id_14" ||
				shackPDESWECart2DBenchmarks->benchmark_name == "rotated_steady_state"
		)
		{
			double freq = 10.0;

			double sx = shackCart2DDataOps->cart2d_domain_size[0];
			double sy = shackCart2DDataOps->cart2d_domain_size[1];

			sweet::Data::Cart2D::DataGrid h_pert_phys(o_h_pert.cart2DDataConfig);
			h_pert_phys.grid_setZero();
			h_pert_phys.grid_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackCart2DDataOps->cart2d_domain_size[0]/(double)shackCart2DDataOps->space_res_physical[0]);
					double y = (double)j*(shackCart2DDataOps->cart2d_domain_size[1]/(double)shackCart2DDataOps->space_res_physical[1]);

					io_data = std::cos(2.0*M_PI*freq*(x/sx+y/sy));
				}
			);

			sweet::Data::Cart2D::DataGrid u_phys(o_u.cart2DDataConfig);
			u_phys.grid_setZero();
			u_phys.grid_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackCart2DDataOps->cart2d_domain_size[0]/(double)shackCart2DDataOps->space_res_physical[0]);
					double y = (double)j*(shackCart2DDataOps->cart2d_domain_size[1]/(double)shackCart2DDataOps->space_res_physical[1]);

					double factor = shackPDESWECart2D->gravitation*2.0*M_PI*freq/(shackPDESWECart2D->cart2d_rotating_f0*sy);
					io_data = factor*std::sin(2.0*M_PI*freq*(x/sx+y/sy));
				}
			);

			sweet::Data::Cart2D::DataGrid v_phys(o_v.cart2DDataConfig);
			v_phys.grid_setZero();
			v_phys.grid_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackCart2DDataOps->cart2d_domain_size[0]/(double)shackCart2DDataOps->space_res_physical[0]);
					double y = (double)j*(shackCart2DDataOps->cart2d_domain_size[1]/(double)shackCart2DDataOps->space_res_physical[1]);

					double factor = -shackPDESWECart2D->gravitation*2.0*M_PI*freq/(shackPDESWECart2D->cart2d_rotating_f0*sx);
					io_data = factor*std::sin(2.0*M_PI*freq*(x/sx+y/sy));
				}
			);

			o_h_pert.loadCart2DDataGrid(h_pert_phys);
			o_u.loadCart2DDataGrid(u_phys);
			o_v.loadCart2DDataGrid(v_phys);
			return true;
		}

		printBenchmarkInformation();
		SWEETErrorFatal(std::string("Benchmark ")+shackPDESWECart2DBenchmarks->benchmark_name+ " not found (or not available)");


		return false;
	}

	void printBenchmarkInformation()
	{
		std::cout << "Some available benchmark scenarios (--benchmark-name):" << std::endl;
		std::cout << "		polvani : Polvani et al (1994) initial condition" << std::endl;
		std::cout << "		mergevortex : Vortex merging initial conditions from McRae QJ 2014 paper" << std::endl;
		std::cout << "		unstablejet : Unstable Jet test case" << std::endl;
		std::cout << "		gaussian_bump : Gaussian bump" << std::endl;
		std::cout << "		normalmodes : Normal mode initialization" << std::endl;
		std::cout << "Check out more options in src/include/benchmarks_swe_cart2d/SWECart2DBenchmarksCombined.hpp" << std::endl;
	}
};



}}

#endif
