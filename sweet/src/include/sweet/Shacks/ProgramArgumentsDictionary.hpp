/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_SHACKS_PROGRAMARGUMENTSDICTIONARY_HPP
#define INCLUDE_SWEET_SHACKS_PROGRAMARGUMENTSDICTIONARY_HPP


#include <list>
#include <memory>
#include <typeinfo>
#include <sweet/Tools/ProgramArguments.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Error/Assert.hpp>
#include <sweet/Shacks/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>


namespace sweet {
namespace Shacks {

/*!
 * \brief Extends Shack 'Dictionary' to process program arguments
 */
class ProgramArgumentsDictionary	:
		public Dictionary
{
private:
	sweet::Tools::ProgramArguments _programArguments;

	int _argc;
	char *const *_argv;

public:
	ProgramArgumentsDictionary(
			int i_argc,
			char *const *i_argv
	)	:
		Dictionary(),
		_argc(i_argc),
		_argv(i_argv)
	{
	}

public:
	ProgramArgumentsDictionary()	:
		Dictionary(),
		_argc(-1),
		_argv(nullptr)
	{
	}

	/*!
	 * Setup with arguments given in constructor
	 */
public:
	bool setup()
	{
		_programArguments.setup(_argc, _argv);

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(_programArguments);
		return true;
	}


	/*!
	 * Setup program arguments
	 */
public:
	void setupArgs(
		int i_argc,	//!< main's argc argument
		char *const *i_argv	//!< main's argv argument
	)
	{
		_argc = i_argc;
		_argv = i_argv;
	}

	/*!
	 * Clear everything
	 */
public:
	void clear(bool i_withArgReset = false)
	{
		_programArguments.clear();
		Dictionary::clear();

		if (i_withArgReset)
		{
			_argc = -1;
			_argv = nullptr;
		}
	}

	/*!
	 * Process all program arguments using all registered shacks
	 */
public:
	bool processProgramArguments(
			bool i_skipProcessedShacks = true 	//!< Skip already processed shacks
	)
	{
		return Dictionary::processProgramArguments(_programArguments, i_skipProcessedShacks);
	}

	/*!
	 * Process program arguments only for a particular Shack.
	 *
	 * This is helpful if we need some early evaluation of program arguments.
	 */
	bool processProgramArgumentsForShack(
			Base *io_shackInstance	//!< Shack to process arguments for
	)
	{
		bool retval = io_shackInstance->processProgramArguments(_programArguments);
		io_shackInstance->argumentsProcessed = true;
		return retval;
	}


	/*!
	 * Check for help arguments -h or --help.
	 */
	bool processHelpArguments(
			bool i_withError = true	//<! Trigger an error if help was triggered
	)
	{
		/*
		 * First, check for --help or -h
		 */
		if (_programArguments.argumentWithKeyExists("-h") || _programArguments.argumentWithKeyExists("--help"))
		{
			std::cout << "Printing help:" << std::endl;
			printProgramArguments();
			if (i_withError)
				error.set("Help requested");

			return false;
		}

		return true;
	}

	/*!
	 * Check if all program arguments have been processed.
	 *
	 * This ensures that no program argument is left out, e.g., because of a typo
	 */
	bool checkAllArgumentsProcessed(
			bool i_createError = true		//!< Trigger an error if program argument wasn't processed
	)
	{
		_programArguments.checkAllArgumentsProcessed(i_createError);
		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(_programArguments);
	}

public:
	~ProgramArgumentsDictionary()
	{
	}
};

}}

#endif
