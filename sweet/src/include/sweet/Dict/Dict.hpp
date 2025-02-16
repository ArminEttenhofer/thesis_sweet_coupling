/*
 * Dict.hpp
 *
 *  Created on: Feb 13, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DICT_DICT_HPP
#define INCLUDE_SWEET_DICT_DICT_HPP

#include <string>
#include <vector>
#include <complex>
#include <fstream>
#include <sweet/Dict/_DictElement.hpp>
#include <sweet/Dict/_DictFileWrite.hpp>
#include <sweet/Dict/_TypesBasic.hpp>
#include <sweet/Dict/TypesArrayND.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Error/Fatal.hpp>

namespace sweet {
namespace Dict {

/*!
 * \brief Simple dictionary allowing to store scalar and array data
 *
 * It also supports reading/writing to/from files in C++ and Python
 */
class Dict:
		public _TypesBasic
{
public:
	sweet::Error::Base error;

private:
	std::vector<_DictElement> _dict;

	/*!
	 * Return the number of elements in the dictionary
	 */
public:
	int getNumElements()	const
	{
		return _dict.size();
	}

	/*!
	 * Return element at given index
	 */
public:
	const _DictElement& operator[](int i_index)	const
	{
		SWEET_ASSERT(i_index < getNumElements());

		return _dict[i_index];
	}
	_DictElement& operator[](int i_index)
	{
		SWEET_ASSERT(i_index < getNumElements());

		return _dict[i_index];
	}

private:
	bool _debug;

	const char* SWEET_MAGIC_FILE_CODE = "SWEETFileDict";


private:
	void _basicTest()
	{
		if (sizeof(int64) != 8)
			SWEETErrorFatal("Something is weird with int64!");

		if (sizeof(float64) != 8)
			SWEETErrorFatal("Something is weird with float64!");
	}


	/*!
	 * Setup empty dictionary
	 */
public:
	Dict(bool i_debug = false)	:
		_debug(i_debug)
	{
		_basicTest();
	}

	/*!
	 * Load dictionary from given file
	 */
public:
	Dict(const std::string &i_filename, bool i_debug = false)
	{
		_debug = i_debug;

		_basicTest();

		fileLoad(i_filename);
	}

	/*!
	 * Get element's value given a key
	 */
public:
	template <typename T>
	void get(
			const std::string &i_key,
			T &o_value
	)	const
	{
		for (std::size_t i = 0; i < _dict.size(); i++)
		{
			if (_dict[i].getKey() == i_key)
			{
				_dict[i].get(o_value);
				return;
			}
		}

		std::ostringstream ss;
		ss << "Key '" << i_key << "' not found!";
		SWEETErrorFatal(ss.str());
	}



	/*!
	 * Return whether a key exists
	 */
	bool keyExists(
			const std::string &i_key
	)	const
	{
		for (std::size_t i = 0; i < _dict.size(); i++)
		{
			if (_dict[i].getKey() == i_key)
				return true;
		}

		return false;
	}

	/*!
	 * Return the index of a key
	 */
	int keyIndex(const std::string &i_key)	const
	{
		for (std::size_t i = 0; i < _dict.size(); i++)
		{
			if (_dict[i].getKey() == i_key)
				return i;
		}

		return -1;
	}


	/*!
	 * Set (or update) a value given a key
	 */
	template <typename T>
	bool set(
			const std::string &i_key,		//!< Key to set or update
			const T &i_value,				//!< Value to set
			bool i_overwriteExistingValue = true	//<! Overwrite value if key already exists, otherwise trigger an error.
	)
	{
		int idx = keyIndex(i_key);

		if (idx < 0)
		{
			_dict.push_back(_DictElement(i_key, i_value));
			return true;
		}

		if (!i_overwriteExistingValue)
			return error.set("Key already exists");

		_dict[idx].set(i_key, i_value);
		return true;
	}


public:
	friend
	std::ostream& operator<<(
			std::ostream& io_os,
			const Dict &i_dict
	)
	{
		for (std::size_t i = 0; i < i_dict._dict.size(); i++)
		{
			const _DictElement &e = i_dict._dict[i];

			std::cout << " + " << e << std::endl;
		}

		return io_os;
	}


	/*!
	 * Load dictionary from a file
	 *
	 * See corresponding Python file for description of file format
	 */
public:
	bool fileLoad(const std::string &i_filename)
	{
		_DictFileRead f(i_filename);

		std::string magic_start = f.loadStr0();

		if (magic_start != SWEET_MAGIC_FILE_CODE)
		{
			std::ostringstream ss;
			ss << "Invalid Magic code '" << magic_start << "' at beginning!";
			SWEETErrorFatal(ss.str());
		}

		std::size_t num_entries = f.loadData<int64>();


		if (_debug)
			std::cout << "Dict: Found " << num_entries << " dictionary entries" << std::endl;

		_dict.resize(num_entries);
		for (std::size_t i = 0; i < num_entries; i++)
		{
			_DictElement &e = _dict[i];

			e.fileLoadKeyTypeValue(f);

			if (_debug)
				std::cout << " + Loaded '" << e.getKey() << "' => '" << e.getValueAsString() << "'" << std::endl;
		}

		std::string magic_end = f.loadStr0();

		if (magic_start != SWEET_MAGIC_FILE_CODE)
		{
			std::ostringstream ss;
			ss << "Invalid Magic code '" << magic_end << "' at end!";
			SWEETErrorFatal(ss.str());
		}

		return false;
	}


	/*!
	 * Save dictionary data to a file
	 */
public:
	bool fileSave(
			const std::string &i_filename,
			int i_verbosity = 0
	)
	{
		_DictFileWrite f(i_filename);

		f.writeStr0(SWEET_MAGIC_FILE_CODE);

		f.writeData<int64>(_dict.size());

		for (std::size_t i = 0; i < _dict.size(); i++)
		{
			_DictElement &e = _dict[i];

			e.fileWriteKeyTypeValue(f);

			if (_debug)
				std::cout << " + Loaded '" << e.getKey() << "' => '" << e.getValueAsString() << "'" << std::endl;
		}

		f.writeStr0(SWEET_MAGIC_FILE_CODE);

		return false;
	}

};

}}

#endif
