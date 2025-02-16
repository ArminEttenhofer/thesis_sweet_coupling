/*
 * ProgramArguments.hpp
 *
 *  Created on: Feb 18, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_TOOLS_PROGRAMARGUMENTS_HPP
#define INCLUDE_SWEET_TOOLS_PROGRAMARGUMENTS_HPP

#include <vector>
#include <list>
#include <string>
#include <complex>
#include <ostream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <sweet/Error/Base.hpp>
#include <sweet/Error/Assert.hpp>
#include "StringSplit.hpp"


namespace sweet {
namespace Tools {

/*!
 * \brief Program argument parser
 *
 * Get arguments of the form
 *
 * -T 123
 * --foo-bar=123
 *
 * and split them into a list of "key => value"
 *
 * Supported, but not recommended:
 * --foo-foo message
 */
class ProgramArguments
{
public:
	Error::Base error;

private:
	/*!
	 * Single program argument
	 */
	class _ProgramArgument
	{
	public:
		std::string key;
		std::string value;
		bool argumentParsedAndAccessed;

		_ProgramArgument(
				const std::string &i_key,
				const std::string &i_value
		)	:
			key(i_key),
			value(i_value),
			argumentParsedAndAccessed(false)
		{
		}


		_ProgramArgument()	:
			argumentParsedAndAccessed(false)
		{
		}

		void reset()
		{
			key = "";
			value = "";
			argumentParsedAndAccessed = false;
		}
	};

	/*!
	 * Vector of all program arguments
	 */
	std::vector<_ProgramArgument> _arguments;

	/*!
	 * Provided argument string to be parsed
	 */
	std::string _arg0;

	/*!
	 * True if setup is finished.
	 */
	bool _setupFinished;

	/*!
	 * Trigger an error if key is parsed twice
	 */
	bool _errorForDoubleParsing;

	/*!
	 * List to keep track of keys already parsed.
	 */
	std::list<std::string> _keysInParsing;

	/*!
	 * Error if duplicated keys were provided by user.
	 * This typically means that one argument is parsed at two different locations.
	 */
	bool _errorForDuplicateKeysInParsing;

	/*!
	 * Strip the dashed of a key? E.g. "--help" becomes "help"
	 */
	bool _stripKeyDashes;


public:
	ProgramArguments(
			bool i_errorForDoubleProgramArgumentParsing = true,	//!< Trigger an error if an argument is parsed twice
			bool i_errorForDuplicateKeysInParsing = true,		//!< Trigger an error if there are duplicate keys
			bool i_stripKeyDashes = false			//!< Strip dashes at the keys (e.g., will convert "--help" to "help" so that only "help" needs to be provided as a key	)
	)	:
		_setupFinished(false),
		_errorForDoubleParsing(i_errorForDoubleProgramArgumentParsing),
		_errorForDuplicateKeysInParsing(i_errorForDuplicateKeysInParsing),
		_stripKeyDashes(i_stripKeyDashes)
	{
	}

	/*!
	 * Clear everything to reuse class from scratch
	 */
	void clear()
	{
		_arg0 = "";
		_arguments.clear();

		_keysInParsing.clear();

		_setupFinished = false;
	}


	/*!
	 * Add a new program argument
	 */
private:
	void _addArguments(
			const std::string i_key,	//!< Key, e.g., from key=value
			const std::string i_value	//!< Value, e.g., from key=value
	)
	{
		if (argumentWithKeyExists(i_key))
		{
			error.set("Argument with key '"+i_key+"' already exists - did you specify this twice?");
			return;
		}

		_arguments.push_back(
			_ProgramArgument(
					i_key,
					i_value
				)
			);
	}

	/*!
	 * Initialize this class given the program's arguments
	 */
public:
	bool setup(
			int i_argc,				//!< argc from main()
			char *const *i_argv		//!< argv from main()
	)
	{
		if (i_argc == 0)
		{
			error.set("No argument at all available!");
			return false;
		}

		SWEET_ASSERT(i_argv[0] != nullptr);

		_arg0 = i_argv[0];

		/*
		 * 0: no previous processed data
		 * 1: key processed
		 */
		int state = 0;

		std::string tmp_key;

		for (int i = 1; i < i_argc; i++)
		{
			std::string arg = i_argv[i];

			if (state == 0)
			{
				/*
				 * Basic sanity check
				 */

				if (arg.size() <= 1)
				{
					std::stringstream ss;
					ss << "Error parsing argument " << i << ": '" << arg << "' (too short)" << std::endl;
					return error.set(ss.str());
				}

				/*
				 * Check for single or double dash
				 */
				int num_dashes = 0;

				if (arg[0] == '-')
				{
					if (arg[1] == '-')
					{
						num_dashes = 2;
					}
					else
					{
						num_dashes = 1;
					}
				}

				if (num_dashes == 0)
				{
					std::stringstream ss;
					ss << "Error parsing argument " << i << ": '" << arg << "' (missing dashes)" << std::endl;
					return error.set(ss.str());
				}

				if (!_stripKeyDashes)
					num_dashes = 0;

				/*
				 * Search for "=" separator
				 */
				std::size_t pos = arg.find('=');

				if (pos == std::string::npos)
				{
					/*
					 * Special handling for
					 * "--help", "-h"
					 * which will be simply stored with an empty key
					 */
					if (arg == "--help" || arg == "-h")
					{
						_addArguments(
								arg,
								""
							);
						continue;
					}
					/*
					 * Not found => continue
					 */
					state = 1;
					tmp_key = arg.substr(num_dashes, std::string::npos);
					continue;
				}

				_addArguments(
						arg.substr(num_dashes, pos-num_dashes),
						arg.substr(pos+1, std::string::npos)
					);
				continue;
			}

			if (state == 1)
			{
				_addArguments(tmp_key, arg);

				state = 0;
				continue;
			}
		}

		if (state == 1)
			return error.set("Invalid format of program arguments (last one could not be parsed)");

		_setupFinished = true;

		return true;
	}

	/*!
	 * Return whether an argument with the given key exists
	 *
	 * \return True if key exists
	 */
public:
	bool argumentWithKeyExists(
			const std::string& i_key	//!< Key to search for
	)
	{
		for (std::size_t i = 0; i < _arguments.size(); i++)
		{
			_ProgramArgument &a = _arguments[i];
			if (a.key == i_key)
			{
				return true;
			}
		}

		return false;
	}



	/*!
	 * Check if all program arguments have been processed.
	 *
	 * This ensures that no program argument is left out, e.g., because of a typo
	 */
public:
	bool checkAllArgumentsProcessed(
			bool i_createError = true	//!< True to create an error message in case that not all program arguments were processed
	)
	{
		if (!_setupFinished)
			error.set("You need to call setup() before searching for arguments!");

		for (std::size_t i = 0; i < _arguments.size(); i++)
		{
			_ProgramArgument &a = _arguments[i];
			if (!a.argumentParsedAndAccessed)
			{
				if (i_createError)
					error.set("Argument '"+a.key+"' not processed!");

				return false;
			}
		}

		return true;
	}

	/*!
	 * Check whether the key already exists and trigger an error if it exists.
	 *
	 * Otherwise, insert into the list to check for this.
	 */
	bool checkAndAddDuplicateKeys(const std::string& i_key)
	{
		for (std::list<std::string>::iterator i = _keysInParsing.begin(); i != _keysInParsing.end(); i++)
		{
			if (i_key == *i)
			{
				error.set("Key '"+i_key+"' parsed twice");
				return false;
			}
		}

		_keysInParsing.push_back(i_key);

		return true;
	}


	/*!
	 * Get the full argument by providing the key.
	 *
	 * This is a private function with specializations
	 * to particular output types provided by other functions.
	 */
private:
	bool _getFullProgramArgumentByKey(
			const std::string& i_key,	//!< Key to search for
			_ProgramArgument** o_pa,	//!< Program argument
			bool i_errorIfKeyNotFound,	//!< Trigger an error if key is not found
			bool i_argErrorForDuplicatedKeysInParsing	//!< Trigger an error if this key was already searched for.
	)
	{
		if (!_setupFinished)
			error.set("You need to call setup() before searching for arguments!");

		/*
		 * Check whether key has been already processed
		 */
		if (_errorForDuplicateKeysInParsing && i_argErrorForDuplicatedKeysInParsing)
			if (!checkAndAddDuplicateKeys(i_key))
				return false;

		for (std::size_t i = 0; i < _arguments.size(); i++)
		{
			_ProgramArgument &a = _arguments[i];
			if (a.key == i_key)
			{
				*o_pa = &a;

				/*
				 * Check whether argument was already parsed
				 */
				if (_errorForDoubleParsing && a.argumentParsedAndAccessed)
				{
					error.set("Argument with key '"+i_key+"' already parsed");
					return false;
				}

				return true;
			}
		}

		if (i_errorIfKeyNotFound)
			error.set(std::string("")+"Key '"+i_key+"' not found");

		return false;
	}

	/*!
	 * Get the value of type string related to a given key.
	 */
public:
	bool getArgumentValueByKey(
			const std::string& i_key,
			std::string &o_value,
			bool i_errorIfKeyNotFound = false,
			bool i_argErrorForDuplicatedKeysInParsing = true
	)
	{
		_ProgramArgument *pa;
		if (!_getFullProgramArgumentByKey(i_key, &pa, i_errorIfKeyNotFound, i_argErrorForDuplicatedKeysInParsing))
			return false;

		o_value = pa->value;
		pa->argumentParsedAndAccessed = true;
		return true;
	}


	/*!
	 * Get the value of type double precision related to a given key.
	 */
public:
	bool getArgumentValueByKey(
			const std::string& i_key,
			std::complex<double> &o_value,
			bool i_error_if_key_is_missing = false,
			bool i_argErrorForDuplicatedKeysInParsing = true
	)
	{
		_ProgramArgument *pa;
		if (!_getFullProgramArgumentByKey(i_key, &pa, i_error_if_key_is_missing, i_argErrorForDuplicatedKeysInParsing))
			return false;

		try
		{
			std::vector<std::string> c = StringSplit::split(pa->value, ",");

			o_value = 0;
			if (c.size() == 1)
			{
				o_value.real(std::stod(c[0]));
			}
			else if (c.size() == 2)
			{
				o_value.real(std::stod(c[0]));
				o_value.imag(std::stod(c[1]));
			}
			else
				throw std::invalid_argument("Neither 1, nor 2 numbers detected for complex number");
		}
		catch (const std::exception &e)
		{
			error.set("Exception caught during conversion of value '"+pa->value+"' to double: "+e.what());
			return false;
		}

		pa->argumentParsedAndAccessed = true;

		return true;
	}

	/*!
	 * Get the value of type double precision related to a given key.
	 */
public:
	bool getArgumentValueByKey(
			const std::string& i_key,
			double &o_value,
			bool i_error_if_key_is_missing = false,
			bool i_argErrorForDuplicatedKeysInParsing = true
	)
	{
		_ProgramArgument *pa;
		if (!_getFullProgramArgumentByKey(i_key, &pa, i_error_if_key_is_missing, i_argErrorForDuplicatedKeysInParsing))
			return false;

		try
		{
			o_value = std::stod(pa->value);
		}
		catch (const std::exception &e)
		{
			error.set("Exception caught during conversion of value '"+pa->value+"' to double: "+e.what());
			return false;
		}

		pa->argumentParsedAndAccessed = true;

		return true;
	}


	/*!
	 * Get the value of type integer related to a given key.
	 */
	bool getArgumentValueByKey(
			const std::string& i_key,
			int &o_value,
			bool i_errorIfKeyNotFound = false,
			bool i_argErrorForDuplicatedKeysInParsing = true
	)
	{
		_ProgramArgument *pa;
		if (!_getFullProgramArgumentByKey(i_key, &pa, i_errorIfKeyNotFound, i_argErrorForDuplicatedKeysInParsing))
			return false;


		try
		{
			o_value = std::stoi(pa->value);
		}
		catch (const std::exception &e)
		{
			error.set("Exception caught during conversion of value '"+pa->value+"' to integer: "+e.what());
			return false;
		}

		pa->argumentParsedAndAccessed = true;

		return true;
	}

	/*!
	 * Get the value of type boolean related to a given key.
	 */
public:
	bool getArgumentValueByKey(
			const std::string& i_key,
			bool &o_value,
			bool i_errorIfKeyNotFound = false,
			bool i_argErrorForDuplicatedKeysInParsing = true
	)
	{
		_ProgramArgument* pa;
		if (!_getFullProgramArgumentByKey(i_key, &pa, i_errorIfKeyNotFound, i_argErrorForDuplicatedKeysInParsing))
			return false;

		std::string val = pa->value;

		for (std::size_t i = 0; i < val.length(); i++)
		{
			if (val[i] >= 'A' && val[i] <= 'Z')
				val[i] -= 'A'-'a';
		}

		if (val == "1" || val == "true")
		{
			o_value = true;
			pa->argumentParsedAndAccessed = true;
			return true;
		}

		if (val == "0" || val == "false")
		{
			o_value = false;
			pa->argumentParsedAndAccessed = true;
			return true;
		}

		error.set("Cannot parse value '" + val +"' as boolean type");

		return false;
	}

	/*!
	 * Get the value of given type related to a given key.
	 *
	 * Here, 2 keys are handed over and the first match is used.
	 *
	 * This is helpful in case that there are multiple keys (e.g. alternatives) available.
	 */
	template <typename T>
	bool getArgumentValueBy2Keys(
			const std::string& i_key1,
			const std::string& i_key2,
			T &o_value,
			bool i_errorIfKeyNotFound = false,
			bool i_argErrorForDuplicatedKeysInParsing = true
	)
	{
		_ProgramArgument pa;

		error.assertNoError();

		if (getArgumentValueByKey(i_key1, o_value, i_errorIfKeyNotFound, i_argErrorForDuplicatedKeysInParsing))
			return true;

		if (error.exists())
			return false;

		if (getArgumentValueByKey(i_key2, o_value, i_errorIfKeyNotFound, i_argErrorForDuplicatedKeysInParsing))
			return true;

		return false;
	}

	/*!
	 * Get the value of given type related to a given key.
	 *
	 * Here, 3 keys are handed over and the first match is used.
	 *
	 * This is helpful in case that there are multiple keys (e.g. alternatives) available.
	 */
	template <typename T>
	bool getArgumentValueBy3Keys(
			const std::string& i_key1,
			const std::string& i_key2,
			const std::string& i_key3,
			T &o_value,
			bool i_errorIfKeyNotFound = false,
			bool i_argErrorForDuplicatedKeysInParsing = true
	)
	{
		_ProgramArgument pa;

		error.assertNoError();
		if (getArgumentValueByKey(i_key1, o_value, i_errorIfKeyNotFound, i_argErrorForDuplicatedKeysInParsing))
			return true;

		if (error.exists())
			return false;

		if (getArgumentValueByKey(i_key2, o_value, i_errorIfKeyNotFound, i_argErrorForDuplicatedKeysInParsing))
			return true;

		if (error.exists())
			return false;

		if (getArgumentValueByKey(i_key3, o_value, i_errorIfKeyNotFound, i_argErrorForDuplicatedKeysInParsing))
			return true;

		return false;
	}


	/*!
	 * output stream helper
	 */
	friend
	std::ostream&
	operator<<(std::ostream &io_os, const ProgramArguments &i_pa)
	{
		for (std::size_t i = 0; i < i_pa._arguments.size(); i++)
		{
			const _ProgramArgument &a = i_pa._arguments[i];
			std::cout << "'" << a.key << "' => '" << a.value << "' (processed=" << a.argumentParsedAndAccessed << ")" << std::endl;
		}
		return io_os;
	}
};

}}

#endif
