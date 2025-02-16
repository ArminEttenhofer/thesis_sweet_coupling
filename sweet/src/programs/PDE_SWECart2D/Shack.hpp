/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWECART2D_SHACKPDESWECART2D_HPP
#define PROGRAMS_PDE_SWECART2D_SHACKPDESWECART2D_HPP

#include <sweet/Shacks/Base.hpp>
#include <string>
#include <iostream>
#include <sweet/Tools/ProgramArguments.hpp>


namespace PDE_SWECart2D {

/**
 * simulation coefficients
 */
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
	 * Cart2D with f-Coriolis rotation
	 */
	double cart2d_rotating_f0 = 1.0; //Cart2D

	/**
	 * Avoid nonlinear divergence and only solve linear one
	 */
	bool use_only_linear_divergence = false;


	/**
	 * Diffusion applied only on nonlinear divergence
	 */
	int use_nonlinear_only_visc = 0;

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
	 * Check for instabilities and stop
	 */
	bool instability_checks = false;

	/*
	 * Compute diagnostics
	 */
	bool compute_diagnostics = false;

	void printProgramArguments(const std::string& i_prefix = "") override
	{
		std::cout << i_prefix << "PDESWECart2D parameters:" << std::endl;
		std::cout << i_prefix << "	--pde-h0 [float]	average (initial) height of water" << std::endl;
		std::cout << i_prefix << "	--pde-viscosity [visc]	viscosity, , default=0" << std::endl;
		std::cout << i_prefix << "	--pde-viscosity-order [visc]	viscosity order, default=2" << std::endl;
		std::cout << i_prefix << "	--pde-gravitation [float]	gravity" << std::endl;
		std::cout << i_prefix << "	-f [float]	f-parameter for f-cart2d, default=0" << std::endl;
		std::cout << i_prefix << "	--use-only-linear-divergence [bool]	Use only linear divergence" << std::endl;
		std::cout << i_prefix << "	--use-nonlinear-only-visc [bool]	Use only viscosity on nonlinear part" << std::endl;
		std::cout << i_prefix << "	--compute-errors [bool]	Compute errors to analytical solution (if available)" << std::endl;
		std::cout << i_prefix << "	--normal-mode-analysis-generation=[int]			Control generation of normal mode analysis (default: 0)" << std::endl;
		std::cout << i_prefix << "	--instability-checks=[bool]			Check for instabilities (default: 0)" << std::endl;
		std::cout << i_prefix << "	--compute-diagnostics [bool]	Compute diagnostics" << std::endl;
		std::cout << i_prefix << "" << std::endl;
	}


	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa) override
	{
		i_pa.getArgumentValueBy3Keys("--pde-h0", "-H", "--h0", h0);
		i_pa.getArgumentValueBy2Keys("--pde-viscosity", "-u", viscosity);
		i_pa.getArgumentValueBy2Keys("--pde-viscosity-order", "-U", viscosity_order);
		i_pa.getArgumentValueBy3Keys("--pde-g", "-g", "--pde-gravitation", gravitation);

		i_pa.getArgumentValueByKey("-f", cart2d_rotating_f0);

		i_pa.getArgumentValueByKey("--use-only-linear-divergence", use_only_linear_divergence);
		i_pa.getArgumentValueByKey("--use-nonlinear-only-visc", use_nonlinear_only_visc);

		i_pa.getArgumentValueByKey("--normal-mode-analysis-generation", normal_mode_analysis_generation);
		i_pa.getArgumentValueByKey("--compute-errors", compute_errors);

		i_pa.getArgumentValueByKey("--instability-checks", instability_checks);

		i_pa.getArgumentValueByKey("--compute-diagnostics", compute_diagnostics);

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_pa);
		return true;
	}


	void printShack(
		const std::string& i_prefix = ""
	) override
	{
		std::cout << std::endl;
		std::cout << i_prefix << "SIMULATION SWE PLANE:" << std::endl;
		std::cout << i_prefix << " + h0: " << h0 << std::endl;
		std::cout << i_prefix << " + viscosity: " << viscosity << std::endl;
		std::cout << i_prefix << " + viscosity_order: " << viscosity_order << std::endl;
		std::cout << i_prefix << " + gravitation: " << gravitation << std::endl;
		std::cout << i_prefix << " + cart2d_rotating_f0: " << cart2d_rotating_f0 << std::endl;
		std::cout << i_prefix << " + use_only_linear_divergence: " << use_only_linear_divergence << std::endl;
		std::cout << i_prefix << " + use_nonlinear_only_visc: " << use_nonlinear_only_visc << std::endl;
		std::cout << i_prefix << " + compute_errors: " << compute_errors << std::endl;
		std::cout << i_prefix << " + normal_mode_analysis_generation: " << normal_mode_analysis_generation << std::endl;
		std::cout << i_prefix << " + instability_checks: " << instability_checks << std::endl;
		std::cout << i_prefix << " + compute_diagnostics: " << compute_diagnostics << std::endl;
		std::cout << std::endl;
	}
};

}

#endif
