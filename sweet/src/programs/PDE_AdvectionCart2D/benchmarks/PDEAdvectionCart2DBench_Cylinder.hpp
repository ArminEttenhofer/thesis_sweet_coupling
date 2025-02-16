/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef PROGRAMS_PDE_ADVECTIONCART2D_BENCHMARKS_PDEADVECTIONCART2DBENCH_CYLINDER_HPP
#define PROGRAMS_PDE_ADVECTIONCART2D_BENCHMARKS_PDEADVECTIONCART2DBENCH_CYLINDER_HPP


#include <programs/PDE_AdvectionCart2D/benchmarks/PDEAdvectionCart2DBench_BaseInterface.hpp>
#include <programs/PDE_AdvectionCart2D/benchmarks/ShackPDEAdvectionCart2DBenchmarks.hpp>
#include <stdlib.h>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <cmath>


/**
 * Setup Cylinder
 */
class PDEAdvectionCart2DBenchCylinder	:
		public PDEAdvectionCart2DBench_BaseInterface
{

public:
	bool setupBenchmark(
			sweet::Data::Cart2D::DataSpectral &o_h_pert,
			sweet::Data::Cart2D::DataSpectral &o_u,
			sweet::Data::Cart2D::DataSpectral &o_v
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
				double dx = x-shackBenchmarks->object_coord_x*sx;
				double dy = y-shackBenchmarks->object_coord_y*sy;

				double radius = shackBenchmarks->object_scale*sqrt(sx*sx+sy*sy);
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
};


#endif
