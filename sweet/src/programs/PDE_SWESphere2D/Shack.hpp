/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_SHACK_HPP
#define PROGRAMS_PDE_SWESPHERE2D_SHACK_HPP

#include <sweet/Tools/ProgramArguments.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/Base.hpp>

namespace PDE_SWESphere2D {

class Shack	:
		public sweet::Shacks::Base
{
public:
	/**
	 * Average height for perturbed formulation
	 *
	 * We use a default value similar to the Galewski benchmark
	 */
	double h0 = 10000.0;


	/*
	 * For more information on viscosity,
	 *
	 * see 13.3.1 "Generic Form of the Explicit Diffusion Mechanism"
	 * in "Numerical Techniques for Global Atmospheric Models"
	 */

	/*
	 * viscosity-term on velocities with 2nd order diff operator
	 */
	double viscosity = 1e-6;

	/**
	 * Order of viscosity
	 */
	int viscosity_order = 2;


	/**
	 * Gravitational constant
	 */
	double gravitation = 9.80616;


	/**
	 * Simulation on f-sphere2D? (constant f0 term over entire sphere2D)
	 */
	bool sphere2d_use_fsphere2D = false;


	/**
	 * Coriolis effect
	 * 7.2921 x 10^{-5}
	 */
	double sphere2d_rotating_coriolis_omega = 0.000072921;

	double sphere2d_fsphere2d_f0 = 0.00007292*2; //Sphere2D


	/*
	 * Do a normal mode analysis, see
	 * Hillary Weller, John Thuburn, Collin J. Cotter,
	 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
	 */
	int normal_mode_analysis_generation = 0;

	/*
	 * Compute errors compared to analytical solution
	 */
	bool compute_errors = false;


	/*
	 * Compute diagnostics
	 */
	bool compute_diagnostics = false;


	/*
	 * Check for instabilities and stop
	 */
	bool instability_checks = false;


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << "PDESWESphere2D parameters:" << std::endl;
		std::cout << i_prefix << "	--pde-h0 [float]	average (initial) height of water" << std::endl;
		std::cout << i_prefix << "	--pde-viscosity [visc]	Viscosity, default=0" << std::endl;
		std::cout << i_prefix << "	--pde-viscosity-order [visc]	Viscosity order, default=2" << std::endl;
		std::cout << i_prefix << "	--pde-gravitation [float]	Gravitation" << std::endl;
		std::cout << i_prefix << "	--compute-errors [bool]	Compute errors to analytical solution (if available)" << std::endl;
		std::cout << i_prefix << "	--compute-diagnostics [bool]	Compute diagnostics" << std::endl;
		std::cout << i_prefix << "	--normal-mode-analysis-generation=[int]			Control generation of normal mode analysis (default: 0)" << std::endl;
		std::cout << i_prefix << "	--instability-checks=[bool]			Check for instabilities (default: 0)" << std::endl;
	}


	void printShack(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << "PDESWESphere2D parameters:" << std::endl;
		std::cout << i_prefix << "[MULE] shack.PDE_SWESphere2D.h0: " << h0 << std::endl;
		std::cout << i_prefix << "[MULE] shack.PDE_SWESphere2D.viscosity: " << viscosity << std::endl;
		std::cout << i_prefix << "[MULE] shack.PDE_SWESphere2D.viscosity_order: " << viscosity_order << std::endl;
		std::cout << i_prefix << "[MULE] shack.PDE_SWESphere2D.gravitation: " << gravitation << std::endl;
		std::cout << i_prefix << "[MULE] shack.PDE_SWESphere2D.sphere2d_rotating_coriolis_omega: " << sphere2d_rotating_coriolis_omega << std::endl;
		std::cout << i_prefix << "[MULE] shack.PDE_SWESphere2D.sphere2d_use_fsphere2D: " << sphere2d_use_fsphere2D << std::endl;
		std::cout << i_prefix << "[MULE] shack.PDE_SWESphere2D.sphere2d_fsphere2d_f0: " << sphere2d_fsphere2d_f0 << std::endl;
		std::cout << i_prefix << "[MULE] shack.PDE_SWESphere2D.compute_errors: " << compute_errors << std::endl;
		std::cout << i_prefix << "[MULE] shack.PDE_SWESphere2D.compute_diagnostics: " << compute_diagnostics << std::endl;
		std::cout << i_prefix << "[MULE] shack.PDE_SWESphere2D.instability_checks: " << instability_checks << std::endl;
	}



	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueBy3Keys("--pde-h0", "-H", "--h0", h0);
		i_pa.getArgumentValueBy2Keys("--pde-viscosity", "-u", viscosity);
		i_pa.getArgumentValueBy2Keys("--pde-viscosity-order", "-U", viscosity_order);
		i_pa.getArgumentValueBy3Keys("--pde-g", "-g", "--pde-gravitation", gravitation);
		i_pa.getArgumentValueByKey("-F", sphere2d_use_fsphere2D);

		if (i_pa.getArgumentValueByKey("-f", sphere2d_rotating_coriolis_omega))
			sphere2d_fsphere2d_f0 = sphere2d_rotating_coriolis_omega;

		i_pa.getArgumentValueByKey("--normal-mode-analysis-generation", normal_mode_analysis_generation);
		i_pa.getArgumentValueByKey("--compute-errors", compute_errors);
		i_pa.getArgumentValueByKey("--compute-diagnostics", compute_diagnostics);
		i_pa.getArgumentValueByKey("--instability-checks", instability_checks);
		
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_pa);
		return true;
	}
};

}

#endif
