/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef PROGRAMS_ODE_GENERIC_FILEOUTPUT_BASE_HPP
#define PROGRAMS_ODE_GENERIC_FILEOUTPUT_BASE_HPP


// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>

///#include "Shack.hpp"

namespace ODE_Generic {
namespace FileOutput {

class Base
{
public:
	typedef std::complex<double> T;

	sweet::Error::Base error;

	sweet::IO::Shack *shackIOData;
	sweet::TimeTree::Shack *shackTimestepControl;

	virtual
	~Base(){};

	void setup(
			sweet::IO::Shack *i_shackIOData,
			sweet::TimeTree::Shack *i_shackTimestepControl
	)
	{
		shackIOData = i_shackIOData;
		shackTimestepControl = i_shackTimestepControl;
	}

	virtual
	bool fileSave(
			const T &i_U,
			///const std::string &i_filename
			const std::string &i_name
	)	const = 0;

};

}}

#endif
