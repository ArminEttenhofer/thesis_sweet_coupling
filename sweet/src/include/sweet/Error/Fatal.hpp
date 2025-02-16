#ifndef INCLUDE_SWEET_ERROR_FATAL_HPP
#define INCLUDE_SWEET_ERROR_FATAL_HPP


#include <cassert>
#include <string>
#include <iostream>
#include <signal.h>
#include <stdlib.h>

#include <sweet/Error/StackBacktrace.hpp>


namespace sweet{
namespace Error{

/************************************************************
 * SWEETError in release and debug mode
 ************************************************************/
class _Fatal
{
public:
//	[[ noreturn ]]
	_Fatal(
			const std::string &i_error_type,
			const std::string &i_error_message,
			const char* i_filename,
			int i_line_no,
			const char* i_func,
			bool stop_after_error = true
		)
	{
		std::cerr << std::flush << std::endl;
		std::cerr << "********************************************" << std::endl;

		if (i_error_message != "")
			std::cerr << " " << i_error_type << ": " << i_error_message << std::endl;
		else
			std::cerr << " " << i_error_type << std::endl;

		std::cerr << "********************************************" << std::endl;
		std::cerr << " +     File: " << i_filename << std::endl;
		std::cerr << " + Line Nr.: " << i_line_no << std::endl;
		std::cerr << " + Function: " << i_func << std::endl;
		std::cerr << "********************************************" << std::endl;
		std::cerr << std::endl;
		if (stop_after_error)
		{
			std::string gdb = StackBacktrace::getGDBBacktrace();

			std::cout << gdb << std::endl;
			assert(false);	// Trigger every possible error we can trigger to suport debuggers
			raise(SIGABRT);
			exit(-1);
		}
	}
};

}}


// Errors which should never happen
//#define SWEETErrorInternal(msg)	sweet::core::_Fatal("Internal SWEET ERROR", msg, __FILE__, __LINE__, __func__)

// Errors which should never happen
//#define SWEETErrorTODO(msg)		sweet::core::_Fatal("TODO", msg, __FILE__, __LINE__, __func__)

//! Regular errors such as wrong time integration method, negative resolution, etc.
#define SWEETErrorFatal(msg)			sweet::Error::_Fatal("ERROR", msg, __FILE__, __LINE__, __func__)

// Regular errors such as wrong time integration method, negative resolution, etc.
//#define _SWEETErrornostop(msg)			sweet::core::_Fatal("ERROR", msg, __FILE__, __LINE__, __func__, false)


#endif
