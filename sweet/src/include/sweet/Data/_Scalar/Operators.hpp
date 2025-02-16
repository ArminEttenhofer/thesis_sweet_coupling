/*
 * Operators.hpp
 *
 *  Created on: 11 Jul 2023
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef INCLUDE_SWEET_DATA_SCALAR_OPERATORS_HPP
#define INCLUDE_SWEET_DATA_SCALAR_OPERATORS_HPP

#include "Config.hpp"
#include <sweet/Memory/MemBlockAlloc.hpp>
#include <sweet/Error/Base.hpp>


namespace sweet {
namespace Data {
namespace Scalar {


#if SWEET_XBRAID

class Operators
{
	friend Config;

public:
	Error::Base error;

	const Scalar::Config* scalarDataConfig;
	const Scalar::Shack* shackScalarDataOps;

public:
	Operators(
		const Scalar::Config* i_scalarDataConfig,
		const Scalar::Shack* i_shackScalarDataOps
	)
	{
		setup(i_scalarDataConfig, i_shackScalarDataOps);
	}

public:
	Operators()	:
		scalarDataConfig(nullptr)
	{
	}


public:
	void setup(
		const Scalar::Config* i_scalarDataConfig,
		const Scalar::Shack* i_shackScalarDataOps
	)
	{
		scalarDataConfig = i_scalarDataConfig;
		shackScalarDataOps = i_shackScalarDataOps;

	}


	void clear()
	{
		if (scalarDataConfig == nullptr)
			return;

		scalarDataConfig = nullptr;
	}

};

#endif

}}}

#endif
