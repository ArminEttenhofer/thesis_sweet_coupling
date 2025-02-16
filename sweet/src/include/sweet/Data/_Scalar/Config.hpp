/*
 * SPHSetup.hpp
 *
 *  Created on: 11 Jul 2023
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef INCLUDE_SWEET_DATA_SCALAR_CONFIG_HPP
#define INCLUDE_SWEET_DATA_SCALAR_CONFIG_HPP

#include <iostream>
#include <cmath>
#include <sweet/IO/FileOperations.hpp>
#include <sweet/Tools/TransformationPlans.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Error/Fatal.hpp>
#include <sweet/LibMath/shtns_inc.hpp>
#include <stdexcept>

#include "Shack.hpp"

#if SWEET_MPI
#	include <mpi.h>
#endif

namespace sweet {
namespace Data {
namespace Scalar {

#if SWEET_XBRAID

/*!
 * \brief Dummy Config data class for Generic_ODE
 */
class Config
{
	friend class Operators;
	///friend class Sphere2DData;

public:
	Error::Base error;

	int ndim = 1;
	std::size_t grid_res[1] = {1};
	std::size_t spectral_data_size[1] = {1};
	int spectral_array_data_number_of_elements = 1;

public:
	Config()
	{
	}

public:
	bool setupAuto(
		sweet::Data::Scalar::Shack* i_shackDataOps
	)
	{
		return true;
	}


	void clear()
	{
	}


	~Config()
	{
		clear();
	}
};

#endif

}}}

#endif
