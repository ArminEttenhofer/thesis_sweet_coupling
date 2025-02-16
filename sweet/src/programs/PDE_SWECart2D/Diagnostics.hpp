/*
 * ShackDiagnostics.hpp
 *
 *  Created on: Feb 21, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWECART2D_DIAGNOSTICS_HPP
#define PROGRAMS_PDE_SWECART2D_DIAGNOSTICS_HPP

#include <programs/PDE_SWECart2D/Shack.hpp>
#include <sweet/Shacks/Base.hpp>
#include <string>
#include <iostream>
#include <sweet/Tools/ProgramArguments.hpp>



/**
 * Diagnostic variables
 */
namespace PDE_SWECart2D {

class Diagnostics
///	public sweet::Shacks::Base
{
	sweet::Data::Cart2D::Operators *cart2DOperators;
	Shack *shackPDESWECart2D;

	//sweet::Data::Cart2D::Helpers_Integral sphere2DHelpers_integral;
	//sweet::Data::Cart2D::DataGrid fg;

public:
	//! total mass
	double total_mass = 0;

	//! kinetic energy
	double kinetic_energy = 0;

	//! potential energy
	double potential_energy = 0;

	//! total energy
	double total_energy = 0;

	//! total potential enstropy
	double total_potential_enstrophy = 0;

	double ref_total_mass = -1;
	double ref_kinetic_energy = -1;
	double ref_potential_energy = -1;
	double ref_total_energy = -1;
	double ref_total_potential_enstrophy = -1;

	int last_update_timestep_nr;

public:
	Diagnostics()	:
		cart2DOperators(nullptr),
		shackPDESWECart2D(nullptr),
		total_mass(0),
		potential_energy(0),
		kinetic_energy(0),
		total_potential_enstrophy(0),
		total_energy(0),

		ref_total_mass(0),
		ref_potential_energy(0),
		ref_kinetic_energy(0)
	{
	}

	void setup(
			sweet::Data::Cart2D::Operators *i_cart2DOperators,
			Shack *i_shackPDESWECart2D,
			int i_verbose = 1
	)
	{
		cart2DOperators = i_cart2DOperators;
		shackPDESWECart2D = i_shackPDESWECart2D;

		//sphere2DHelpers_integral.setup(sphere2DOperators->sphere2DDataConfig, i_verbose);

		//setupFG();

		last_update_timestep_nr = -1;
	}


	void backup_reference()
	{
		ref_total_mass = total_mass;
		ref_kinetic_energy = kinetic_energy;
		ref_potential_energy = potential_energy;
		ref_total_energy = total_energy;
		ref_total_potential_enstrophy = total_potential_enstrophy;
	}

	void print(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "DIAGNOSTICS:" << std::endl;
		std::cout << " + total_mass: " << total_mass << std::endl;
		std::cout << " + total_energy: " << total_energy << std::endl;
		std::cout << " + kinetic_energy: " << kinetic_energy << std::endl;
		std::cout << " + potential_energy: " << potential_energy << std::endl;
		std::cout << " + total_potential_enstrophy: " << total_potential_enstrophy << std::endl;
		std::cout << std::endl;
	}


	void update_nonstaggered_huv_2_mass_energy_enstrophy(
			sweet::Data::Cart2D::Operators &ops,
			sweet::Data::Cart2D::Shack *shackCart2DDataOps,
			PDE_SWECart2D::Shack *shackPDESWECart2D,
			sweet::Data::Cart2D::DataSpectral &i_prog_h, //h perturbation
			sweet::Data::Cart2D::DataSpectral &i_prog_u,
			sweet::Data::Cart2D::DataSpectral &i_prog_v
	)
	{
		double normalization = (shackCart2DDataOps->cart2d_domain_size[0]*shackCart2DDataOps->cart2d_domain_size[1]) /
								((double)shackCart2DDataOps->space_res_physical[0]*(double)shackCart2DDataOps->space_res_physical[1]);

		//std::cout << "Size x, sixe y" << (shackCart2DDataOps->domain_size[0]) << (shackCart2DDataOps->domain_size[1]) << std::endl;
		//std::cout << "resphysx, resphysy" << (double)shackCart2DDataOps->res_physical[0] << (double)shackCart2DDataOps->res_physical[1] << std::endl;
		//std::cout << "normal" << normalization << std::endl;

		sweet::Data::Cart2D::DataGrid h_phys = i_prog_h.toGrid();
		sweet::Data::Cart2D::DataGrid u_phys = i_prog_u.toGrid();
		sweet::Data::Cart2D::DataGrid v_phys = i_prog_v.toGrid();

		// mass (mean depth needs to be added)
		total_mass = (h_phys+ shackPDESWECart2D->h0).grid_reduce_sum_quad() * normalization;

		// energy
		sweet::Data::Cart2D::DataGrid pot_energy = (h_phys + shackPDESWECart2D->h0)*(shackPDESWECart2D->gravitation*normalization);
		sweet::Data::Cart2D::DataGrid kin_energy = (h_phys + shackPDESWECart2D->h0)*(u_phys*u_phys + v_phys*v_phys)*(0.5*normalization);

		potential_energy = pot_energy.grid_reduce_sum_quad();
		kinetic_energy = kin_energy.grid_reduce_sum_quad();

		total_energy = kinetic_energy + potential_energy;

		// absolute vorticity
		sweet::Data::Cart2D::DataSpectral eta = (ops.diff_c_x(i_prog_v) - ops.diff_c_y(i_prog_u) + shackPDESWECart2D->cart2d_rotating_f0);

		// enstrophy
		total_potential_enstrophy = 0.5*(eta*eta).toGrid().grid_reduce_sum_quad() * normalization;
	}



public:
	void update_staggered_huv_2_mass_energy_enstrophy(
			sweet::Data::Cart2D::Operators &op,
			sweet::Data::Cart2D::Shack *shackCart2DDataOps,
			PDE_SWECart2D::Shack *shackPDESWECart2D,
			sweet::Data::Cart2D::DataSpectral &i_prog_h,
			sweet::Data::Cart2D::DataSpectral &i_prog_u,
			sweet::Data::Cart2D::DataSpectral &i_prog_v
	)
	{
		double normalization = (shackCart2DDataOps->cart2d_domain_size[0]*shackCart2DDataOps->cart2d_domain_size[1]) /
								((double)shackCart2DDataOps->space_res_physical[0]*(double)shackCart2DDataOps->space_res_physical[1]);

		sweet::Data::Cart2D::DataGrid h_phys = i_prog_h.toGrid();

		// mass
		total_mass = (h_phys + shackPDESWECart2D->h0).grid_reduce_sum_quad() * normalization;

		sweet::Data::Cart2D::DataGrid u_phys = op.avg_b_x(i_prog_u.toGrid());
		sweet::Data::Cart2D::DataGrid v_phys = op.avg_b_y(i_prog_v.toGrid());

		sweet::Data::Cart2D::DataSpectral u(u_phys.cart2DDataConfig);
		sweet::Data::Cart2D::DataSpectral v(v_phys.cart2DDataConfig);
		u.loadCart2DDataGrid(u_phys);
		v.loadCart2DDataGrid(v_phys);

		// energy
		sweet::Data::Cart2D::DataGrid pot_energy = (h_phys + shackPDESWECart2D->h0)*(shackPDESWECart2D->gravitation*normalization);
		sweet::Data::Cart2D::DataGrid kin_energy = (h_phys + shackPDESWECart2D->h0)*(u_phys*u_phys+v_phys*v_phys)*(0.5*normalization);

		total_energy = (pot_energy + kin_energy).grid_reduce_sum_quad();

		// total vorticity
		sweet::Data::Cart2D::DataSpectral eta = (op.diff_c_x(v) - op.diff_c_y(u) + shackPDESWECart2D->cart2d_rotating_f0);

		// enstrophy
		total_potential_enstrophy = 0.5*(eta*eta).toGrid().grid_reduce_sum_quad() * normalization;
	}
};

}

#endif
