/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_ERROR_BASE_HPP
#define INCLUDE_SWEET_ERROR_BASE_HPP

#include <sweet/Error/Fatal.hpp>
#include <sweet/Error/StackBacktrace.hpp>
#include <string>
#include <iostream>
#include <ostream>

/*!
 * Do an error check.
 * If there's an error: forward error and return false
 */
#define ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(classWithError)	\
	{ if ((classWithError).error.exists()) return error.forwardWithPositiveReturn((classWithError).error); }

/*!
 * Do an error check.
 * If there's an error: print error and return with EXIT_FAILURE
 */
#define ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(classWithError)	\
	{ if ((classWithError).error.exists()) { (classWithError).error.print(); return EXIT_FAILURE; } }

/*!
 * Do an error check.
 * If there's an error: forward error and return with EXIT_FAILURE
 */
#define ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_BOOLEAN(classWithError)	\
	{ if ((classWithError).error.exists()) { (classWithError).error.print(); return false; } }

/*!
 * Do an error check and return in case of an error
 */
#define ERROR_CHECK_COND_RETURN_BOOLEAN(classWithError)	\
	{ if ((classWithError).error.exists()) { return false; } }


/*!
 * Do an error check.
 * If there's an error: forward error and return
 * If there's no error: continue
 */
#define ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN(classWithError) \
	{ if ((classWithError).error.exists()) { error.forward((classWithError).error); return; } }


/*!
 * If there's an error: forward error and return false
 * If there's no error: return true
 */
#define ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(classWithError)	\
	{ return error.forwardWithPositiveReturn((classWithError).error); }

/*!
 * This is useful for contructors
 * If there's an error: simply forward
 */
#define ERROR_FORWARD(classWithError)	\
	{ error.forward((classWithError).error); }


namespace sweet {
namespace Error {

class Base
{
	bool _hasError;
	std::string _errorMessage;


private:
	/*
	 * Check for environment variable SWEET_ERROR_WITH_STACKTRACE and if it exists,
	 * add stack trace to error information
	 */
	bool _withStacktrace()
	{
		static bool envLoaded = false;
		static bool envExists = false;

		if (!envLoaded)
		{
			char *val = std::getenv("SWEET_ERROR_WITH_STACKTRACE");
			envExists = (val != nullptr);
			envLoaded = true;
		}
		return envExists;
	}


public:
	Base()	:
		_hasError(false)
	{
	}

	/*
	 * Set an error
	 *
	 * \return always *false* to be able to use a one-liner error message
	 */
	bool set(const std::string &i_errorMessage)
	{
#if SWEET_THREADING_TIME_REXI
#pragma omp critical
#endif
		{
			if (_hasError)
			{
				_errorMessage += " | " + i_errorMessage;
			}
			else
			{
				_hasError = true;
				_errorMessage = i_errorMessage;
			}


			if (_withStacktrace())
			{
				std::string gdbBacktrace = StackBacktrace::getGDBBacktrace();

				if (gdbBacktrace != "")
				{
					_errorMessage += "\n";
					_errorMessage += "Stacktrace (from GDB):\n";
					_errorMessage += gdbBacktrace;
				}
			}
		}
		return false;
	}

	/*!
	 * Forward errors
	 *
	 * \return **false** if there's no error
	 */
	bool forward(Base &i_error)
	{
		bool retval;

#if SWEET_THREADING
#pragma omp critical
#endif
		{
			if (!i_error._hasError)
			{
				retval = false;
			}
			else
			{
				if (this == &i_error)
				{
					retval = true;
				}
				else
				{
					_hasError = i_error._hasError;
					_errorMessage = i_error._errorMessage;
					i_error.clear();
					retval = true;
				}
			}
		}

		return retval;
	}

	/*!
	 * Forward errors with positive return
	 *
	 * \return **true** if there's no error
	 */
	bool forwardWithPositiveReturn(Base &i_error)
	{
		bool retval;
#if SWEET_THREADING_TIME_REXI
#pragma omp critical
#endif
		{
			if (!i_error._hasError)
			{
				retval = true;
			}
			else
			{

				if (this == &i_error)
				{
					retval = false;
				}
				else
				{
					_hasError = i_error._hasError;
					_errorMessage = i_error._errorMessage;

					i_error.clear();
					retval = false;
				}
			}
		}

		return retval;
	}

	bool exists() const
	{
		return _hasError;
	}

	void assertNoError() const
	{
		if (_hasError)
		{
			std::cerr << "Internal problem detected" << std::endl;
			SWEETErrorFatal("Use backtrace to debug this problem!");
			exit(1);
		}
	}

	void clear()
	{
		_hasError = false;
		_errorMessage = "";
	}

	std::string get()
	{
		if (!_hasError)
		{
			std::cerr << "Error message requested, but has no error!" << std::endl;
			SWEETErrorFatal("Use backtrace to debug this problem!");
			exit(1);
		}

		std::string tmp = _errorMessage;
		clear();
		return tmp;
	}

public:
	void print(std::ostream &io_os = std::cerr)
	{
		io_os << "ERROR: " << get() << std::endl;
	}

	~Base()
	{
		if (_hasError)
		{
			std::cerr << "ERROR was not processed!" << std::endl;
			std::cerr << "************************************************************" << std::endl;
			std::cerr << _errorMessage << std::endl;
			std::cerr << "************************************************************" << std::endl;
			exit(1);
		}
	}
};

}}

#endif
