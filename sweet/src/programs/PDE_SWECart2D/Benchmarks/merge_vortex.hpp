/*
 * SWEMergeVortex.hpp
 *
 *  Created on: 01 Nov 2017
 * Author: Pedro Peixoto <pedrosp@ime.usp.br>
 */
#ifndef PROGRAMS_PDE_SWECART2D_BENCHMARKS_PDESWECART2DBENCH_MERGEVORTEX_HPP
#define PROGRAMS_PDE_SWECART2D_BENCHMARKS_PDESWECART2DBENCH_MERGEVORTEX_HPP


#include <programs/PDE_SWECart2D/Benchmarks/BaseInterface.hpp>
#include <stdlib.h>
#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <cmath>


/**
 * Implement merging vortex initial conditions
 * See Energy- and enstrophy-conserving schemes for the shallow-water
 * equations, based on mimetic finite elements
 * Andrew T. T. McRae and Colin J. Cotter
 *
 * IMPORTANT: TO BE USED IN [0,1]x[0,1] domain
 *
 *  To match paper use:
 * f = 1
 * g = 1
 * h0 = 1
 * [0,1]x[0,1
 **/

namespace PDE_SWECart2D {
namespace Benchmarks {

class merge_vortex	:
		public PDE_SWECart2D::Benchmarks::BaseInterface
{
	double f;
	double g;
	double sx;
	double sy;

	double stream(
			double x,
			double y
			)
	{

		//double radius = shackDict.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
		double radius = 1; //shackDict.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
		double factor = 500.0;

		//Andrew's parameters
		double k = 0.03;
		double a = 10.0;
		//my parameters
		//a=20.0;
		factor=(a*a);

		// Gaussian Vortice 1
		//double dx = x-0.45*sx;
		//double dy = y-0.5*sy;
		double dx = x-0.425*sx;
		double dy = y-0.5*sy;

		dx /= radius;
		dy /= radius;

		double exp1 = std::exp(-factor*(dx*dx + dy*dy));

		// Gaussian Vortice 2
		//dx = x-0.55*sx;
		//dy = y-0.5*sy;
		dx = x-0.575*sx;
		dy = y-0.5*sy;

		dx /= radius;
		dy /= radius;

		double exp2 = std::exp(-factor*(dx*dx + dy*dy));

		return (k/a)*(exp1+exp2);

	}

	void setup_stream(
			sweet::Data::Cart2D::DataSpectral &o_psi
	)
	{

		sweet::Data::Cart2D::DataGrid psi_phys(o_psi.cart2DDataConfig);

		for (int j = 0; j < shackCart2DDataOps->space_res_physical[1]; j++)
		{
			for (int i = 0; i < shackCart2DDataOps->space_res_physical[0]; i++)
			{
				// h - lives in the center of the cell
				double x = (((double)i+0.5)/(double)shackCart2DDataOps->space_res_physical[0])*shackCart2DDataOps->cart2d_domain_size[0];
				double y = (((double)j+0.5)/(double)shackCart2DDataOps->space_res_physical[1])*shackCart2DDataOps->cart2d_domain_size[1];

				psi_phys.grid_setValue(j, i, stream(x, y));
			}
		}

		o_psi.loadCart2DDataGrid(psi_phys);

	}


public:
	bool setupBenchmark(
			sweet::Data::Cart2D::DataSpectral &o_h,
			sweet::Data::Cart2D::DataSpectral &o_u,
			sweet::Data::Cart2D::DataSpectral &o_v
	)
	{
		f = shackPDESWECart2D->cart2d_rotating_f0;
		g = shackPDESWECart2D->gravitation;
		sx = shackCart2DDataOps->cart2d_domain_size[0];
		sy = shackCart2DDataOps->cart2d_domain_size[1];

		sweet::Data::Cart2D::DataSpectral psi(o_h.cart2DDataConfig);

		/*
		 * Prepare laplace operator
		 */
		sweet::Data::Cart2D::DataSpectral laplace = ops->diff2_c_x + ops->diff2_c_y;


		/*
		 * Setup stream function
		 */
		setup_stream(psi);
		//psi.file_grid_saveData_ascii("ouput_stream");

		//Calculate Velocities
		o_u = ops->diff_c_y(psi);
		o_v = -ops->diff_c_x(psi);

		//Calculate vorticity
		sweet::Data::Cart2D::DataSpectral vort = ops->vort(o_u, o_v);

		//Solve Poisson equation for height to get balance initial condition
		sweet::Data::Cart2D::DataSpectral lap_h = (f/g)*vort;
		o_h = lap_h.spectral_div_element_wise(laplace);

		return true;
	}


};

}}

#endif
