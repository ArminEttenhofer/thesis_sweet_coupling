/*
 * parmemcpy.hpp
 *
 *  Created on: 11 Oct 2018
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_MEMORY_PARMEMCPY_HPP
#define INCLUDE_SWEET_MEMORY_PARMEMCPY_HPP



#include <complex>
#include <sweet/Parallelization/openmp_helper.hpp>

namespace sweet {
namespace Memory {

/*
 * Just put this function here right now and find a better place later on
 */
inline void parmemcpy(
		void *i_dst,
		const void *i_src,
		std::size_t i_size
)
{
	if ((i_size & 0xf) == 0)
	{
		/*
		 * 8 byte blocks
		 */
		std::complex<double> *dst = (std::complex<double>*)i_dst;
		std::complex<double> *src = (std::complex<double>*)i_src;
		i_size >>= 4;

		SWEET_OMP_PARALLEL_FOR
		for (std::size_t i = 0; i < i_size; i++)
		{
			dst[i] = src[i];
		}
		return;
	}


	if ((i_size & 0x7) == 0)
	{
		/*
		 * 8 byte blocks
		 */
		double *dst = (double*)i_dst;
		double *src = (double*)i_src;
		i_size >>= 3;

		SWEET_OMP_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < i_size; i++)
		{
			dst[i] = src[i];
		}
		return;
	}


	if ((i_size & 0x3) == 0)
	{
		/*
		 * 4 byte blocks
		 */
		float *dst = (float*)i_dst;
		float *src = (float*)i_src;
		i_size >>= 2;

		SWEET_OMP_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < i_size; i++)
		{
			dst[i] = src[i];
		}
		return;
	}

	/*
	 * Fallback
	 */
	char *dst = (char*)i_dst;
	char *src = (char*)i_src;

	SWEET_OMP_PARALLEL_FOR
	for (std::size_t i = 0; i < i_size; i++)
	{
		dst[i] = src[i];
	}
}

}}

#endif
