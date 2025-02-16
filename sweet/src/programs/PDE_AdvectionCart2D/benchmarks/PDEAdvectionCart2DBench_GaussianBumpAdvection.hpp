/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef PROGRAMS_PDE_ADVECTIONCART2D_BENCHMARKS_PDEADVECTIONCART2DBENCH_GAUSSIANBUMPADVECTION_HPP
#define PROGRAMS_PDE_ADVECTIONCART2D_BENCHMARKS_PDEADVECTIONCART2DBENCH_GAUSSIANBUMPADVECTION_HPP


#include <programs/PDE_AdvectionCart2D/benchmarks/PDEAdvectionCart2DBench_BaseInterface.hpp>
#include <programs/PDE_AdvectionCart2D/benchmarks/ShackPDEAdvectionCart2DBenchmarks.hpp>
#include <stdlib.h>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <cmath>

/**
 * Setup Gaussian Bump with advection
 */
class PDEAdvectionCart2DBenchGaussianBumpAdvection	:
		public PDEAdvectionCart2DBench_BaseInterface
{


public:
	bool setupBenchmark(
			sweet::Data::Cart2D::DataSpectral &o_h_pert,
			sweet::Data::Cart2D::DataSpectral &o_u,
			sweet::Data::Cart2D::DataSpectral &o_v
	)
	{

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

					dx /= sx*shackBenchmarks->object_scale;
					dy /= sy*shackBenchmarks->object_scale;

					return std::exp(-i_exp_fac*(dx*dx + dy*dy));
				};


		auto callback_external_forces_advection_field =
				[](
						int i_field_id,
						double i_simulation_timestamp,
						sweet::Data::Cart2D::DataSpectral* o_cart2d_data,			//! cart2ddata or sphere2Ddata
						ShackPDEAdvectionCart2DBenchmarks* o_data_user_void		//! user data (pointer to this class)
		)
		{
			sweet::Data::Cart2D::DataGrid cart2d_data_phys(o_cart2d_data->cart2DDataConfig);
			ShackPDEAdvectionCart2DBenchmarks* shackBenchmarks = (ShackPDEAdvectionCart2DBenchmarks*)o_data_user_void;

			if (i_field_id >= 1 && i_field_id <= 2)
			{
				double u = shackBenchmarks->advection_velocity[0];
				double v = shackBenchmarks->advection_velocity[1];

				double r;
				if (shackBenchmarks->advection_velocity[2] == 0)
					r = 0;
				else
					r = i_simulation_timestamp/shackBenchmarks->advection_velocity[2]*2.0*M_PI;

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

				o_cart2d_data->loadCart2DDataGrid(cart2d_data_phys);

				return;
			}

			SWEETErrorFatal("Non-existing external field requested!");
			return;
		};

		if (shackBenchmarks->advection_velocity[2] != 0)
		{
			// set callback
			shackBenchmarks->getExternalForcesCallback = callback_external_forces_advection_field;

			// set user data to this class
			shackBenchmarks->getExternalForcesUserData = this;

			// setup velocities with initial time stamp
			callback_external_forces_advection_field(
					1,
					shackTimestepControl->currentSimulationTime,
					&o_u,
					shackBenchmarks
				);

			callback_external_forces_advection_field(
					2,
					shackTimestepControl->currentSimulationTime,
					&o_v,
					shackBenchmarks
				);
		}
		else
		{
			sweet::Data::Cart2D::DataGrid u_phys(o_u.cart2DDataConfig);
			sweet::Data::Cart2D::DataGrid v_phys(o_u.cart2DDataConfig);

			u_phys = shackBenchmarks->advection_velocity[0];
			v_phys = shackBenchmarks->advection_velocity[1];

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
};


#endif
