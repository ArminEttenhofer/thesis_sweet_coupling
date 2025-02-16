/*
 * Dict.hpp
 *
 *  Created on: Feb 18, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DICT__DICTFILEREAD_HPP
#define INCLUDE_SWEET_DICT__DICTFILEREAD_HPP

#include <sweet/Dict/_TypesBasic.hpp>
#include <sweet/Dict/TypesArrayND.hpp>
#include <sweet/Error/Fatal.hpp>
#include <string>
#include <vector>
#include <complex>
#include <fstream>

namespace sweet {
namespace Dict {

/*!
 * A helper class for Dict to read from files
 */
class _DictFileRead
{
	std::ifstream is;

public:
	_DictFileRead(const std::string &i_filename)
	{
		is = std::ifstream(i_filename, std::ios::in | std::ios::binary);

		if (!is.is_open())
			SWEETErrorFatal(std::string("Unable to open file ")+i_filename);
	}


public:
	/*!
	 * Read 0 terminated string from file
	 */
	std::string loadStr0()
	{
		std::vector<char> buffer;
		buffer.resize(1024);

		bool found = false;
		for (std::size_t i = 0; i < buffer.size()-1; i++)
		{
			is.read((char*)&(buffer[i]), sizeof(char));
			if (buffer[i] == 0)
			{
				found = true;
				break;
			}
		}

		if (!found)
			SWEETErrorFatal("0 termination character not found, stopping here");

		// convert to std::string
		return (std::string)(char*)buffer.data();
	}

	template <typename T>
	T loadData()
	{
		T retval;
		is.read((char*)&retval, sizeof(retval));
		return retval;
	}

	template <int D, typename T>
	void loadArray(
			TypesArrayND<D,T> &array
	)
	{
		for (std::size_t i = 0; i < array.size(); i++)
			array.data()[i] = loadData<T>();
	}
};


}}

#endif
