/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef PROGRAMS_PDE_ADVECTIONCART2D_BENCHMARKS_PDEADVECTIONCART2DBENCH_GAUSSIANBUMP_HPP
#define PROGRAMS_PDE_ADVECTIONCART2D_BENCHMARKS_PDEADVECTIONCART2DBENCH_GAUSSIANBUMP_HPP


#include <programs/PDE_AdvectionCart2D/benchmarks/PDEAdvectionCart2DBench_BaseInterface.hpp>
#include <programs/PDE_AdvectionCart2D/benchmarks/ShackPDEAdvectionCart2DBenchmarks.hpp>
#include <stdlib.h>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <cmath>

/**
 * Setup Gaussian Bump
 */
class PDEAdvectionCart2DBenchGaussianBump	:
		public PDEAdvectionCart2DBench_BaseInterface
{

public:
	bool setupBenchmark(
			sweet::Data::Cart2D::DataSpectral &o_h_pert,
			sweet::Data::Cart2D::DataSpectral &o_u,
			sweet::Data::Cart2D::DataSpectral &o_v
	)
	{
		sweet::Data::Cart2D::DataGrid h_pert_phys(o_h_pert.cart2DDataConfig);
		sweet::Data::Cart2D::DataGrid u_phys(o_u.cart2DDataConfig);
		sweet::Data::Cart2D::DataGrid v_phys(o_v.cart2DDataConfig);

		double sx = shackCart2DDataOps->cart2d_domain_size[0];
		double sy = shackCart2DDataOps->cart2d_domain_size[1];

		h_pert_phys.grid_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
			{
				double x = (double)i*(shackCart2DDataOps->cart2d_domain_size[0]/(double)shackCart2DDataOps->space_res_physical[0]);
				double y = (double)j*(shackCart2DDataOps->cart2d_domain_size[1]/(double)shackCart2DDataOps->space_res_physical[1]);

				// Gaussian
				double dx = x-shackBenchmarks->object_coord_x*sx;
				double dy = y-shackBenchmarks->object_coord_y*sy;


				double radius = shackBenchmarks->object_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
				dx /= radius;
				dy /= radius;

				io_data = std::exp(-50.0*(dx*dx + dy*dy));
			}
		);

		u_phys.grid_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				io_data = shackBenchmarks->advection_velocity[0];
			}
		);

		v_phys.grid_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
			{
				io_data = shackBenchmarks->advection_velocity[1];
			}
		);

		o_h_pert.loadCart2DDataGrid(h_pert_phys);
		o_u.loadCart2DDataGrid(u_phys);
		o_v.loadCart2DDataGrid(v_phys);

		return true;
	}
};


#endif
