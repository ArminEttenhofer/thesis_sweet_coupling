/*
 * Dict.hpp
 *
 *  Created on: Feb 18, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DICT__TYPESBASIC_HPP
#define INCLUDE_SWEET_DICT__TYPESBASIC_HPP

#include <sweet/Dict/TypesArrayND.hpp>
#include <sweet/Error/Fatal.hpp>
#include <string>
#include <vector>
#include <complex>
#include <fstream>

namespace sweet {
namespace Dict {

/*!
 * \brief A class with type information only available to Dict implementations
 */
class _TypesBasic
{
public:
	typedef long long int64;
	typedef double float64;
	typedef std::complex<double> complex128;
};

}}

#endif
