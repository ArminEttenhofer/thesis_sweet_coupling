/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef INCLUDE_SWEET_DATA_SCALAR_SHACK_HPP
#define INCLUDE_SWEET_DATA_SCALAR_SHACK_HPP


#include <sweet/Shacks/Base.hpp>
#include <string>
#include <iostream>
#include <sweet/Tools/ProgramArguments.hpp>
#include <sweet/Tools/TransformationPlans.hpp>
#include <sweet/Tools/StringSplit.hpp>


namespace sweet {
namespace Data {
namespace Scalar {

/*!
 * \brief Dummy Shack for Scalar data
 */
class Shack	:
		public sweet::Shacks::Base
{
public:
	/**
	 * resolution in physical space (grid cells)
	 */
	int space_res_physical[1] = {1};

	/**
	 * resolution in spectral space (number of modes)
	 */
	int space_res_spectral[1] = {1};

	void printProgramArguments(const std::string& i_prefix = "") override
	{
	}

	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa) override
	{
		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(i_pa);
	}

	bool validateResolution()
	{
		return true;
	}

	void printShack(
		const std::string& i_prefix = ""
	) override
	{
	}
};

}}}

#endif
