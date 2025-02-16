/*
 * DiagBandedMatrix.hpp
 *
 *  Created on: 24 Aug 2016
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_LIBMATH_LAPACKBANDEDMATRIXSOLVERCOMPLEX_HPP
#define INCLUDE_SWEET_LIBMATH_LAPACKBANDEDMATRIXSOLVERCOMPLEX_HPP

#include <complex>
#include <string.h>
#include <stdlib.h>
#include <limits>
#include <sweet/Memory/parmemcpy.hpp>


namespace sweet {
namespace LibMath {


class BandedMatrixSolverComplex
{
	typedef std::complex<double> T;

	int max_N;
	int num_diagonals;
	int num_halo_size_diagonals;

	std::complex<double>* AB;


public:
	BandedMatrixSolverComplex()	:
		AB(nullptr)
	{

	}

	~BandedMatrixSolverComplex()
	{
		clear();
	}


	void setup(
			int i_max_N,				//!< size of the matrix
			int i_num_off_diagonals		//!< number of block diagonals
	)
	{
		clear();

		max_N = i_max_N;
		num_diagonals = 2*i_num_off_diagonals+1;
		num_halo_size_diagonals = i_num_off_diagonals;

		SWEET_ASSERT(2*num_halo_size_diagonals+1 == num_diagonals);

		AB = (std::complex<double>*)malloc(sizeof(std::complex<double>)*num_diagonals*i_max_N);
	}


	void clear()
	{
		if (AB != nullptr)
		{
			free(AB);
			AB = nullptr;
		}
	}

	/*!
	 * Output array
	 */
public:
	void print_array_c(
			const std::complex<double> *i_data,
			int i_cols,
			int i_rows
	)	const
	{
		// rows
		for (int j = 0; j < i_rows; j++)
		{
			std::cout << j << ": ";
			// cols
			for (int i = 0; i < i_cols; i++)
			{
				std::cout << i_data[j*i_cols+i] << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}


	/*!
	 * Solve for input matrix
	 *
	 * i_A: cols: num_diagonals
	 *      rows: i_num_rows
	 *
	 * The Fortran array will be transposed and with a size of (rows: LDAB, cols: i_num_rows)
	 *
	 * i_b: RHS of equation
	 *
	 * o_x: Solution
	 */
public:
	void solve_diagBandedInverse_Carray(
		const std::complex<double>* i_A,
		const std::complex<double>* i_b,
		std::complex<double>* o_x,
		int i_num_rows,
		int i_debug_block
	)	const
	{
		SWEET_ASSERT(max_N >= i_num_rows);

#if 1
		/*
		 * self-written solver
		 */
		int halo_size = num_halo_size_diagonals;
		int mat_size = i_num_rows;

		sweet::Memory::parmemcpy((void*)AB, (void*)i_A, sizeof(std::complex<double>)*num_diagonals*mat_size);
		sweet::Memory::parmemcpy((void*)o_x, (void*)i_b, sizeof(std::complex<double>)*mat_size);

		T* b = o_x;

#define aArray(i, j)	(AB[((i)*num_diagonals)+(j)])

#if 0
		std::cout << "mat_size: " << mat_size << std::endl;
		for (int j = 0; j < mat_size; j++)
		{
			for (int i = 0; i < num_diagonals; i++)
			{
				std::cout << j << ", " << i << " => " << aArray(j,i) << std::endl;
			}
			std::cout << std::endl;
		}
#endif

		/*
		 * Eliminate lower diagonal matrix
		 * #for ipiv in range(mat_size-1)
		 */
		for (int ipiv = 0; ipiv < mat_size-1; ipiv++)
		{
			// Pivot element given by A[ipiv,ipiv]

			/*
			 * Process rows below below A[i,i] to zero
		     * #for jrow in range(ipiv+1,min(ipiv+halo_size+1, mat_size)):
			 */
			for (int jrow = ipiv+1; jrow < std::min(ipiv+halo_size+1, mat_size); jrow++)
			{
	            SWEET_ASSERT(aArray(ipiv, halo_size).real() != 0 || aArray(ipiv, halo_size).imag() != 0);

	            T piva = aArray(jrow,ipiv-jrow+halo_size)/aArray(ipiv,halo_size);

	            /*
	             * Iterate over columns
		         * #for jcol in range(ipiv, min(ipiv+halo_size+1, mat_size)):
	             */
	            for (int jcol = ipiv; jcol < std::min(ipiv+halo_size+1, mat_size); jcol++)
	            {
	            	// #a[jrow,jcol-jrow+halo_size] -= a[ipiv,jcol-ipiv+halo_size]*piva
	            	aArray(jrow,jcol-jrow+halo_size) -= aArray(ipiv,jcol-ipiv+halo_size)*piva;
	            }

				b[jrow] -= b[ipiv]*piva;
			}
		}

		/*
		 * Eliminate upper diagonal matrix
		 * #for ipiv in range(mat_size-1, -1, -1):
		 */
		for (int ipiv = mat_size-1; ipiv >= 0; ipiv--)
		{
			/*
			 * Pivot element given by a[ipiv,ipiv]
			 */

			/*
			 * Set col above A[i,i] to zero
			 * => Iterate over rows below
			 *
			 * #for jrow in range(max(ipiv-halo_size, 0), ipiv):
			 */
			for (int jrow = std::max(ipiv-halo_size, 0); jrow < ipiv; jrow++)
			{
	            SWEET_ASSERT(aArray(ipiv, halo_size).real() != 0 || aArray(ipiv, halo_size).imag() != 0);

				T piva = aArray(jrow,ipiv-jrow+halo_size)/aArray(ipiv,halo_size);

				/*
				 * Iterate over columns
				 * # for jcol in range(max(ipiv-halo_size+1, 0), min(ipiv+1, mat_size)):
				 */
				for (int jcol = std::max(ipiv-halo_size+1, 0); jcol < std::min(ipiv+1, mat_size); jcol++)
					aArray(jrow,jcol-jrow+halo_size) -= aArray(ipiv,jcol-ipiv+halo_size)*piva;

				b[jrow] -= b[ipiv]*piva;
			}

			b[ipiv] /= aArray(ipiv,halo_size);
		}

#undef aArray
		return;
#endif
	}


	/*!
	 * Solve for banded matrix "A"
	 *
	 * i_A: cols: num_diagonals
	 *      rows: i_num_rows
	 *
	 * io_b: RHS of equation and result for output
	 */
public:
	static
	void solve(
		std::complex<double>* i_A,	///<! Storage for banded matrix A
		int i_rows,							///<! Rows of matrix A
		int i_halo_size,					///<! Number of off diagonal
		std::complex<double>* io_b	///<! RHS for input and solution for output
	)
	{
		int halo_size = i_halo_size;
		int num_diagonals = halo_size*2+1;
		int rows = i_rows;

		std::complex<double>* b = io_b;

#define aArray(i, j)	(i_A[((i)*num_diagonals)+(j)])

#if 0
		std::cout << "mat_size: " << rows << std::endl;
		for (int j = 0; j < rows; j++)
		{
			for (int i = 0; i < num_diagonals; i++)
			{
				std::cout << j << ", " << i << " => " << aArray(j,i) << std::endl;
			}
			std::cout << std::endl;
		}
#endif

		/*
		 * Eliminate lower diagonal matrix
		 * #for ipiv in range(mat_size-1)
		 */
		for (int ipiv = 0; ipiv < rows-1; ipiv++)
		{
			// Pivot element given by A[ipiv,ipiv]

			/*
			 * Process rows below below A[i,i] to zero
		     * #for jrow in range(ipiv+1,min(ipiv+halo_size+1, mat_size)):
			 */
			for (int jrow = ipiv+1; jrow < std::min(ipiv+halo_size+1, rows); jrow++)
			{
	            SWEET_ASSERT(aArray(ipiv, halo_size).real() != 0 || aArray(ipiv, halo_size).imag() != 0);

	            T piva = aArray(jrow,ipiv-jrow+halo_size)/aArray(ipiv,halo_size);

	            /*
	             * Iterate over columns
		         * #for jcol in range(ipiv, min(ipiv+halo_size+1, mat_size)):
	             */
	            for (int jcol = ipiv; jcol < std::min(ipiv+halo_size+1, rows); jcol++)
	            {
	            	// #a[jrow,jcol-jrow+halo_size] -= a[ipiv,jcol-ipiv+halo_size]*piva
	            	aArray(jrow,jcol-jrow+halo_size) -= aArray(ipiv,jcol-ipiv+halo_size)*piva;
	            }

				b[jrow] -= b[ipiv]*piva;
			}
		}

		/*
		 * Eliminate upper diagonal matrix
		 * #for ipiv in range(mat_size-1, -1, -1):
		 */
		for (int ipiv = rows-1; ipiv >= 0; ipiv--)
		{
			/*
			 * Pivot element given by a[ipiv,ipiv]
			 */

			/*
			 * Set col above A[i,i] to zero
			 * => Iterate over rows below
			 *
			 * #for jrow in range(max(ipiv-halo_size, 0), ipiv):
			 */
			for (int jrow = std::max(ipiv-halo_size, 0); jrow < ipiv; jrow++)
			{
	            SWEET_ASSERT(aArray(ipiv, halo_size).real() != 0 || aArray(ipiv, halo_size).imag() != 0);

				T piva = aArray(jrow,ipiv-jrow+halo_size)/aArray(ipiv,halo_size);

				/*
				 * Iterate over columns
				 * # for jcol in range(max(ipiv-halo_size+1, 0), min(ipiv+1, mat_size)):
				 */
				for (int jcol = std::max(ipiv-halo_size+1, 0); jcol < std::min(ipiv+1, rows); jcol++)
					aArray(jrow,jcol-jrow+halo_size) -= aArray(ipiv,jcol-ipiv+halo_size)*piva;

				b[jrow] -= b[ipiv]*piva;
			}

			b[ipiv] /= aArray(ipiv,halo_size);
		}

#undef aArray
		return;
	}
};

}}

#endif
