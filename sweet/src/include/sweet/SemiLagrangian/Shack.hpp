/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_SEMILAGRANGIAN_SHACK_HPP
#define INCLUDE_SWEET_SEMILAGRANGIAN_SHACK_HPP


#include <sweet/Shacks/Base.hpp>
#include <iostream>


namespace sweet {
namespace SemiLagrangian {

/*!
 * \brief Shack for semi-Lagrangian methods
 */
class Shack	:
		public sweet::Shacks::Base
{
public:
	//! Number of iterations for semi-Lagrangian methods
	int semi_lagrangian_max_iterations = 2;

	//! Convergence threshold for semi-Lagrangian methods (set to -1 to ignore error)
	double semi_lagrangian_convergence_threshold = -1;

	//! Method to use for computing depature points
	std::string semi_lagrangian_departure_point_method = "settls";

	//! Use limiter for higher-order SL interpolation (cubic)
	bool semi_lagrangian_interpolation_limiter = false;

	//! Create pseudo points at poles for interpolation
	bool semi_lagrangian_sampler_use_pole_pseudo_points = false;

	//! Use accurate spherical geometry (???) or approximation (Ritchie 1995)
	bool semi_lagrangian_approximate_sphere2d_geometry = true;



	void printProgramArguments(const std::string& i_prefix = "") override
	{
		std::cout << "Semi-Lagrangian time Discretization:" << std::endl;
		std::cout << "	--semi-lagrangian-max-iterations [int]		Number of max. iterations during semi-Lagrangian time integration" << std::endl;
		std::cout << "	--semi-lagrangian-convergence-threshold [float]	Threshold to stop iterating, Use -1 to disable" << std::endl;

		std::cout << "	--semi-lagrangian-departure-point-method [str]		'settls' (default), 'midpoint', 'std'" << std::endl;
		std::cout << "	--semi-lagrangian-sampler-use-pole-pseudo-points [bool]" << std::endl;
		std::cout << "								Use pseudo poles for sampling" << std::endl;
		std::cout << "									false (default)" << std::endl;
		std::cout << "	--semi-lagrangian-interpolation-limiter [bool]	Use limiter for cubic interpolation" << std::endl;
		std::cout << "	--semi-lagrangian-convergence-threshold [float]	Threshold to stop iterating, Use -1 to disable" << std::endl;
		std::cout << "	--semi-lagrangian-approximate-sphere2d-geometry [int]	0: no approximation, 1: Richies approximation, default: 0" << std::endl;
	}


	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa) override
	{
		i_pa.getArgumentValueByKey("--semi-lagrangian-max-iterations", semi_lagrangian_max_iterations);
		i_pa.getArgumentValueByKey("--semi-lagrangian-convergence-threshold", semi_lagrangian_convergence_threshold);

		i_pa.getArgumentValueByKey("--semi-lagrangian-departure-point-method", semi_lagrangian_departure_point_method);
		i_pa.getArgumentValueByKey("--semi-lagrangian-sampler-use-pole-pseudo-points", semi_lagrangian_sampler_use_pole_pseudo_points);
		i_pa.getArgumentValueByKey("--semi-lagrangian-interpolation-limiter", semi_lagrangian_interpolation_limiter);
		i_pa.getArgumentValueByKey("--semi-lagrangian-approximate-sphere2d-geometry", semi_lagrangian_approximate_sphere2d_geometry);

		if (i_pa.error.exists())
			return error.forwardWithPositiveReturn(i_pa.error);

		return true;
	}


	void printShack(
		const std::string& i_prefix = ""
	) override
	{
		std::cout << std::endl;
		std::cout << "SEMI-LAGRANGIAN TIME DISCRETIZATION:" << std::endl;
		std::cout << " + semi_lagrangian_max_iterations: " << semi_lagrangian_max_iterations << std::endl;
		std::cout << " + semi_lagrangian_convergence_threshold: " << semi_lagrangian_convergence_threshold << std::endl;
		std::cout << " + semi_lagrangian_departure_point_method: " << semi_lagrangian_departure_point_method << std::endl;
		std::cout << " + semi_lagrangian_interpolation_limiter: " << semi_lagrangian_interpolation_limiter << std::endl;
		std::cout << " + semi_lagrangian_sampler_use_pole_pseudo_points: " << semi_lagrangian_sampler_use_pole_pseudo_points << std::endl;
		std::cout << " + semi_lagrangian_approximate_sphere2d_geometry: " << semi_lagrangian_approximate_sphere2d_geometry << std::endl;
		std::cout << std::endl;
	}
};

}}

#endif
