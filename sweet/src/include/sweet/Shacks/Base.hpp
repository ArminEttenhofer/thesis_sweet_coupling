/*
 * ClassDictionaryInterface.hpp
 *
 *  Created on: Feb 19, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_SHACKS_BASE_HPP
#define INCLUDE_SWEET_SHACKS_BASE_HPP


#include <typeinfo>
#include <sweet/Tools/ProgramArguments.hpp>
#include <sweet/Error/Base.hpp>
#include <string>


namespace sweet {
namespace Shacks {

/*!
 * \brief Base class for Shacks
 *
 * Use this class if you want to implement new program parameters, etc.
 */
class Base
{
public:
	Error::Base error;

	/*!
	 * True if arguments processed of this Shack.
	 *
	 * This is helpful to skip processing if this should be done partially.
	 */
	bool argumentsProcessed;

	Base()	:
		argumentsProcessed(false)
	{
	}

	/*!
	 * Print out all program arguments supported by the shack
	 */
	virtual void printProgramArguments(
			const std::string &i_prefix = ""	//!< Prefix to be used before each line.
	)
	= 0;

	/*!
	 * Process program arguments
	 */
	virtual bool processProgramArguments(
			sweet::Tools::ProgramArguments &i_pa	//!< Program arguments to be used for parsing arguments of this Shack
	)
	= 0;

	/*!
	 * Print information about this shack
	 */
	virtual void printShack(
			const std::string &i_prefix = ""	//!< Prefix to be used before each line.
	)
	= 0;


	virtual ~Base()
	{}
};

}}

#endif
