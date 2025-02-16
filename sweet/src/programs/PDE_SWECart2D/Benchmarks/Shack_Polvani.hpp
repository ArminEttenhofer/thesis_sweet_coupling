/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef PROGRAMS_PDE_SWECART2D_BENCHMARKS_SHACKPDESWECART2DBENCH_POLVANIBENCH_HPP
#define PROGRAMS_PDE_SWECART2D_BENCHMARKS_SHACKPDESWECART2DBENCH_POLVANIBENCH_HPP


#include <sweet/Error/Fatal.hpp>
#include <sweet/Shacks/Base.hpp>


namespace PDE_SWECart2D {
namespace Benchmarks {

struct Shack_Polvani	:
	public sweet::Shacks::Base
{
	/**
	 * Rossby number
	 */
	double r = 0.01;

	/**
	 * Froude number
	 */
	double f = 0.04;


	void printProgramArguments(const std::string& i_prefix = "") override
	{
		std::cout << "" << std::endl;
		std::cout << "Polvani benchmark settings (on the cart2d):" << std::endl;
		std::cout << "	--polvani-rossby [float]	Choose Rossby number, default:0" << std::endl;
		std::cout << "	--polvani-froude [float]	Choose Froude number, default:0" << std::endl;
		std::cout << "" << std::endl;
	}

	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa) override
	{
		i_pa.getArgumentValueByKey("--polvani-rossby", r);
		i_pa.getArgumentValueByKey("--polvani-froude", f);

		return error.forwardWithPositiveReturn(i_pa.error);
	}

	void printShack(
		const std::string& i_prefix = ""
	) override
	{
		std::cout << std::endl;
		std::cout << "SWEPolvani:" << std::endl;
		std::cout << " + Rossby number: " << r << std::endl;
		std::cout << " + Froude number: " << f << std::endl;
		std::cout << std::endl;
	}
};


}}


#endif
