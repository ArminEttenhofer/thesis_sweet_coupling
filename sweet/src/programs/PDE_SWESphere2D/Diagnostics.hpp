/*
 * PDESWESphere2DDiagnostics.hpp
 *
 *  Created on: Mar 10, 2023
 * Author: martin
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_DIAGNOSTICS_HPP
#define PROGRAMS_PDE_SWESPHERE2D_DIAGNOSTICS_HPP


#include <sweet/Data/Sphere2D/Helpers/Helpers_Integral.hpp>
#include "Shack.hpp"


namespace PDE_SWESphere2D {

class Diagnostics
{
	sweet::Data::Sphere2D::Operators *sphere2DOperators;
	Shack *shackPDESWESphere2D;

	sweet::Data::Sphere2D::Helpers_Integral sphere2DHelpers_integral;
	sweet::Data::Sphere2D::DataGrid fg;

public:
	double total_mass;
	double potential_energy;
	double kinetic_energy;
	double total_potential_enstrophy;
	double total_energy;

	double ref_total_mass;
	double ref_potential_energy;
	double ref_kinetic_energy;

	int last_update_timestep_nr;

public:
	Diagnostics()	:
		sphere2DOperators(nullptr),
		shackPDESWESphere2D(nullptr),
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
			sweet::Data::Sphere2D::Operators *i_sphere2DOperators,
			Shack *i_shackPDESWESphere2D,
			int i_verbose = 1
	)
	{
		sphere2DOperators = i_sphere2DOperators;
		shackPDESWESphere2D = i_shackPDESWESphere2D;

		sphere2DHelpers_integral.setup(sphere2DOperators->sphere2DDataConfig, i_verbose);

		setupFG();

		last_update_timestep_nr = -1;
	}


	bool setupFG()
	{
		if (shackPDESWESphere2D->sphere2d_use_fsphere2D)
			fg = sphere2DOperators->getFG_fSphere2D(shackPDESWESphere2D->sphere2d_fsphere2d_f0);
		else
			fg = sphere2DOperators->getFG_rotatingSphere2D(shackPDESWESphere2D->sphere2d_rotating_coriolis_omega);

		return true;
	}


public:
	void update_phi_vrt_div_2_mass_energy_enstrophy(
			const sweet::Data::Sphere2D::Operators *i_ops,
			const sweet::Data::Sphere2D::DataSpectral &i_prog_phi,
			const sweet::Data::Sphere2D::DataSpectral &i_prog_vort,
			const sweet::Data::Sphere2D::DataSpectral &i_prog_div,

			double i_sphere2d_radius,
			double i_gravitation
	)
	{
		SWEET_ASSERT(sphere2DOperators != nullptr);

		sweet::Data::Sphere2D::DataGrid h(sphere2DOperators->sphere2DDataConfig);
		sweet::Data::Sphere2D::DataGrid u(sphere2DOperators->sphere2DDataConfig);
		sweet::Data::Sphere2D::DataGrid v(sphere2DOperators->sphere2DDataConfig);

		h = i_prog_phi.toGrid()*(1.0/i_gravitation);
		i_ops->vrtdiv_2_uv(i_prog_vort, i_prog_div, u, v);

		double normalization = i_sphere2d_radius*i_sphere2d_radius;

		// mass
		total_mass = sphere2DHelpers_integral.compute_zylinder_integral(h) * normalization;

		// energy
		//Sphere2DDataGrid pot_energy = h*(io_simVars.sim.gravitation*normalization);
		sweet::Data::Sphere2D::DataGrid pot_energy = h*h*0.5*normalization;
		sweet::Data::Sphere2D::DataGrid kin_energy = h*(u*u+v*v)*(0.5*normalization);

		potential_energy = sphere2DHelpers_integral.compute_zylinder_integral(pot_energy);
		kinetic_energy = sphere2DHelpers_integral.compute_zylinder_integral(kin_energy);

		/*
		 * We follow the Williamson et al. equation (137) here
		 */
//		double dummy_energy = compute_zylinder_integral(h*h*(0.5*normalization));
//		io_simVars.diag.total_energy = io_simVars.diag.kinetic_energy + dummy_energy;//io_simVars.diag.potential_energy;

		/*
		 * We follow pot/kin energy here
		 */
		total_energy = kinetic_energy + potential_energy;

		// total vorticity
		// TODO: maybe replace this with the i_vort parameter
		sweet::Data::Sphere2D::DataGrid eta(h.sphere2DDataConfig);
		eta = i_ops->uv_2_vort(u, v).toGrid();

		eta += fg;

		// enstrophy (Williamson paper, equation 138)
		total_potential_enstrophy = 0.5*sphere2DHelpers_integral.compute_zylinder_integral(eta*eta/h) * normalization;
	}


	void print(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << i_prefix << "DIAGNOSTICS:" << std::endl;
		std::cout << i_prefix << " + total_mass: " << total_mass << std::endl;
		std::cout << i_prefix << " + total_energy: " << total_energy << std::endl;
		std::cout << i_prefix << " + kinetic_energy: " << kinetic_energy << std::endl;
		std::cout << i_prefix << " + potential_energy: " << potential_energy << std::endl;
		std::cout << i_prefix << " + total_potential_enstrophy: " << total_potential_enstrophy << std::endl;
		std::cout << std::endl;
	}

	void printTabularHeader(
		const std::string& i_prefix = ""
	)
	{
		std::cout << i_prefix << "T\tTOTAL_MASS\tPOT_ENERGY\tKIN_ENERGY\tTOT_ENERGY\tPOT_ENSTROPHY\tREL_TOTAL_MASS\tREL_POT_ENERGY\tREL_KIN_ENERGY\tREL_TOT_ENERGY\tREL_POT_ENSTROPHY";
	}

	void printTabularRow(
		double i_current_simulation_time,
		const std::string& i_prefix = ""
	)
	{
		std::cout << i_prefix;

		// Print simulation time, energy and pot enstrophy
		std::cout << i_current_simulation_time << "\t";
		std::cout << total_mass << "\t";
		std::cout << potential_energy << "\t";
		std::cout << kinetic_energy << "\t";
		std::cout << total_energy << "\t";
		std::cout << total_potential_enstrophy << "\t";

		std::cout << (total_mass-ref_total_mass)/total_mass << "\t";
		std::cout << (potential_energy-ref_potential_energy)/potential_energy << "\t";
		std::cout << (kinetic_energy-ref_kinetic_energy)/kinetic_energy << "\t";
		std::cout << (total_energy-total_energy)/total_energy << "\t";
		std::cout << (total_potential_enstrophy-total_potential_enstrophy)/total_potential_enstrophy;
		std::cout << std::endl;
	}
};

}

#endif
