/*
 * Dict.hpp
 *
 *  Created on: Feb 18, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DICT__DICTFILEWRITE_HPP
#define INCLUDE_SWEET_DICT__DICTFILEWRITE_HPP

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
 * Helper class for Dict to write to file
 */
class _DictFileWrite
{
	std::ofstream os;

public:
	_DictFileWrite(const std::string &i_filename)
	{
		os = std::ofstream(i_filename, std::ios::out | std::ios::binary);

		if (!os.is_open())
			SWEETErrorFatal(std::string("Unable to open file ")+i_filename);
	}


public:
	/*!
	 * Write 0 terminated string from file
	 */
	void writeStr0(const std::string& i_value)
	{
		os.write(i_value.c_str(), i_value.length()+1);
	}

	template <typename T>
	void writeData(const T& i_value)
	{
		os.write((const char*)&i_value, sizeof(i_value));
	}

	template <int D, typename T>
	void writeDictArrayRawData(
			const TypesArrayND<D,T> &i_array
	)
	{
		for (std::size_t i = 0; i < i_array.size(); i++)
			writeData<T>(i_array.data()[i]);
	}

	template <int D, typename T>
	void writeDictArray(
			const TypesArrayND<D,T> &i_array
	)
	{
		// Write out shape
		for (int d = 0; d < D; d++)
			writeData<_TypesBasic::int64>(i_array.shape()[d]);

		writeDictArrayRawData(i_array);
	}
};

}}

#endif
