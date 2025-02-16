/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_CART2D_SHACK_HPP
#define INCLUDE_SWEET_DATA_CART2D_SHACK_HPP


#include <sweet/Shacks/Base.hpp>
#include <string>
#include <iostream>
#include <sweet/Tools/ProgramArguments.hpp>
#include <sweet/Tools/TransformationPlans.hpp>
#include <sweet/Tools/StringSplit.hpp>


namespace sweet {
namespace Data {
namespace Cart2D {

/*!
 * \brief Shack for Cart2D
 */
class Shack	:
		public sweet::Shacks::Base
{
public:
	/*!
	 * resolution in physical space (grid cells)
	 */
	int space_res_physical[2] = {0, 0};


	/*!
	 * resolution in spectral space (number of modes)
	 */
	int space_res_spectral[2] = {0, 0};


	/*!
	 * domain size if running simulation on the cart2d
	 */
	double cart2d_domain_size[2] = {1.0, 1.0};


	/*!
	 * use spectral differential operators
	 */
	bool space_use_spectral_basis_diffs =
#if SWEET_USE_CART2D_SPECTRAL_SPACE || SWEET_USE_SPHERE2D_SPECTRAL_SPACE
			true;
#else
			false;
#endif

	/*!
	 * Use C-grid staggering
	 */
	bool space_grid_use_c_staggering = false;


	/*!
	 * Load / Save plans for SHTNS (useful for reproducibility)
	 *
	 * Can also exist in the Sphere2DDataOps
	 */
	sweet::Tools::TransformationPlans::TRANSFORMATION_PLAN_CACHE reuse_spectral_transformation_plans = sweet::Tools::TransformationPlans::QUICK;

	void printProgramArguments(const std::string& i_prefix = "") override
	{
		std::cout << i_prefix << "Discretization:" << std::endl;
		std::cout << i_prefix << "  >Space:" << std::endl;
		std::cout << i_prefix << "	--space-grid-use-c-staggering [0/1]	Use staggering" << std::endl;
		std::cout << i_prefix << "	-N [res]		resolution in x and y direction, default=0" << std::endl;
		std::cout << i_prefix << "	-n [resx]		resolution in x direction, default=0" << std::endl;
		std::cout << i_prefix << "	-m [resy]		resolution in y direction, default=0" << std::endl;
		std::cout << i_prefix << "	-M [modes]		modes in x/y, default=0" << std::endl;
		std::cout << i_prefix << "	-S [0/1]		Control Operator discretization for Cart2DData" << std::endl;
		std::cout << i_prefix << "					0: FD, 1: spectral derivatives, default: ";
		std::cout << i_prefix << "	-X [length]	length of simulation domain in x direction, default=1" << std::endl;
		std::cout << i_prefix << "	-Y [width]	width of simulation domain in y direction, default=1" << std::endl;
		std::cout << i_prefix << "  --reuse-plans [str]" << std::endl;
		std::cout << i_prefix << "                  'quick': Quick setup" << std::endl;
		std::cout << i_prefix << "                  'save': Save plans" << std::endl;
		std::cout << i_prefix << "                  'load': Load plans if available" << std::endl;
		std::cout << i_prefix << "                  'require_load': Require to load plan and terminate program if plan doesn't exist" << std::endl;
		std::cout << i_prefix << "                  'load_save': Load plan and save it back (E.g., to accumulate plans)" << std::endl;
	}

	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa) override
	{
		i_pa.getArgumentValueByKey("--space-grid-use-c-staggering", space_grid_use_c_staggering);
		i_pa.getArgumentValueByKey("-S", space_use_spectral_basis_diffs);

		std::string tmp_N;
		if (i_pa.getArgumentValueByKey("-N", tmp_N))
		{

			int c = sweet::Tools::StringSplit::split2int(tmp_N, &space_res_physical[0], &space_res_physical[1]);
			if (c == 1)
				space_res_physical[1] = space_res_physical[0];
		}

		std::string tmp_M;
		if (i_pa.getArgumentValueByKey("-M", tmp_M))
		{

			int c = sweet::Tools::StringSplit::split2int(tmp_M, &space_res_spectral[0], &space_res_spectral[1]);
			if (c == 1)
				space_res_spectral[1] = space_res_spectral[0];
		}

		i_pa.getArgumentValueByKey("-X", cart2d_domain_size[0]);
		i_pa.getArgumentValueByKey("-Y", cart2d_domain_size[1]);

		/*
		 * We allow multiple --reuse-plans here since this can be also processed by sphere2D
		 */
		std::string tmp;
		if (i_pa.getArgumentValueByKey("--reuse-plans", tmp, false, false))
		{
			sweet::Tools::TransformationPlans tp;
			tp.getEnumFromString(tmp, reuse_spectral_transformation_plans);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(tp);
		}

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(i_pa);
	}

	bool validateResolution()
	{
		if (
				(space_res_physical[0] == 0 || space_res_physical[1] == 0)	&&
				(space_res_spectral[0] == 0 || space_res_spectral[1] == 0)
		)
			return error.set("Select physical resolution or spectral modes (use -N (or -n, -m) for physical and -M for spectral");

		return true;
	}

	void printShack(
		const std::string& i_prefix = ""
	) override
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "DISCRETIZATION:" << std::endl;
		std::cout << i_prefix << " + space_res_physical: " << space_res_physical[0] << " x " << space_res_physical[1] << std::endl;
		std::cout << i_prefix << " + space_res_spectral: " << space_res_spectral[0] << " x " << space_res_spectral[1] << std::endl;
		std::cout << i_prefix << " + space_use_spectral_basis_diffs: " << space_use_spectral_basis_diffs << std::endl;
		std::cout << i_prefix << " + space_grid_use_c_staggering: " << space_grid_use_c_staggering << std::endl;
		std::cout << i_prefix << " + cart2d_dealiasing (compile time): " <<
#if SWEET_USE_CART2D_SPECTRAL_DEALIASING
		1
#else
		0
#endif
		<< std::endl;
		std::cout << i_prefix << " + domain_size: " << cart2d_domain_size[0] << " x " << cart2d_domain_size[1] << std::endl;
		std::string tmp;
		sweet::Tools::TransformationPlans tp;
		tp.getStringFromEnum(reuse_spectral_transformation_plans, tmp);
		std::cout << i_prefix << " + reuse_spectral_transformation_plans: " << tmp << std::endl;
		std::cout << i_prefix << std::endl;
	}
};

}}}

#endif
