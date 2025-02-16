/*
 * SPHMatrix.hpp
 *
 *  Created on: 24 Aug 2016
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_LIBMATH_BANDEDMATRIXGRIDCOMPLEX_HPP
#define INCLUDE_SWEET_LIBMATH_BANDEDMATRIXGRIDCOMPLEX_HPP

#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Memory/MemBlockAlloc.hpp>


namespace sweet {
namespace LibMath {


/*!
 * Matrix to store coefficients related to Spherical Harmonics
 *
 * This version is for complex physical valued transformations
 */
template <typename T>
class BandedMatrixGridComplex
{
	/*!
	 * Data stores the diagonal and off-diagonal components for a matrix which is to be inverted.
	 * This matrix is partitioned by independent chunks for each mode m.
	 *
	 * The off-diagonal components connect the l-modes only!
	 *
	 * EXAMPLE:
	 *
	 *   Multiplying the data matrix with a vector P would require P to be in the following format:
	 *   P is given by P_n^m
	 *
	 * The the Vector P is given by
	 *   P = (P_0^0, P_1^0, P_2^0, P_3^0, P_1^1, P_2^1, P_3^1, P_2^2, P_3^2, P_3^3)^T
	 */
public:
	T *data;

public:
	/*!
	 * Data storage format for Fortran
	 */
	T *fortran_data;

public:
	int halosize_off_diagonal;
	int num_diagonals;

	const sweet::Data::Sphere2D::Config *sphere2DDataConfig;


public:
	BandedMatrixGridComplex()	:
		data(nullptr),
		fortran_data(nullptr),
		halosize_off_diagonal(-1),
		num_diagonals(-1),
		sphere2DDataConfig(nullptr)
	{
	}


	/*!
	 * Zero all matrix coefficients
	 */
	void zeroAll()
	{
		for (std::size_t i = 0; i < sphere2DDataConfig->spectral_complex_array_data_number_of_elements*num_diagonals; i++)
			data[i] = T(0);
	}



	/*!
	 * Setup data storage
	 */
	void setup(
			const sweet::Data::Sphere2D::Config *i_sphere2DDataConfig,				//!< Handler to sphere2DDataConfig
			int i_halosize_off_diagonal = 0		//!< Size of the halo around. A value of 2 allocates data for 5 diagonals.
	)
	{
		shutdown();
		SWEET_ASSERT(data == nullptr);

		sphere2DDataConfig = i_sphere2DDataConfig;

		halosize_off_diagonal = i_halosize_off_diagonal;
		num_diagonals = 2*halosize_off_diagonal+1;

		data = sweet::Memory::MemBlockAlloc::alloc<T>( sizeof(T)*sphere2DDataConfig->spectral_complex_array_data_number_of_elements*num_diagonals );

		zeroAll();
	}



	/*!
	 * Return matrix row  which is related to the specified modes
	 */
	T *getMatrixRow(
			int n,		//!< row related to P Legendre mode n
			int m		//!< row related to P Fourier mode n
	)
	{
		std::size_t idx = sphere2DDataConfig->getArrayIndexByModes_Complex_NCompact(n, m);
		return data+idx*num_diagonals;
	}





	/*!
	 * Return reference to an element in the row to the specified value
	 */
	const T &rowElement_getRef(
			T *io_row,		//!< pointer to current row
			int i_row_n,	//!< row related to P Legendre mode n
			int i_row_m,	//!< row related to P Fourier mode n
			int rel_n		//!< Relative Legendre mode n (e.g. -1 or +2)
	)
	{
		static T dummy = 0;

		if (i_row_n < i_row_m)
			return dummy;

//		SWEET_ASSERT(i_row_n >= i_row_m);
		//SWEET_ASSERT(i_row_m >= 0);
		SWEET_ASSERT(i_row_m <= sphere2DDataConfig->spectral_modes_m_max);

		int n = i_row_n+rel_n;

		if (n < 0 || n < std::abs(i_row_m) || n > sphere2DDataConfig->spectral_modes_n_max)
			return dummy;

		int idx = rel_n + halosize_off_diagonal;

		SWEET_ASSERT(idx >= 0 && idx < num_diagonals);

		return io_row[idx];
	}



	/*!
	 * Return reference to an element in the row to the specified value
	 */
	void rowElement_set(
			T *io_row,		//!< pointer to current row
			int i_row_n,	//!< row related to P Legendre mode n
			int i_row_m,	//!< row related to P Fourier mode n
			int rel_n,		//!< Relative Legendre mode n (e.g. -1 or +2)
			T i_value
	)
	{
		if (i_row_n < i_row_m)
			return;

//		SWEET_ASSERT(i_row_n >= i_row_m);
		//SWEET_ASSERT(i_row_m >= 0);
		SWEET_ASSERT(i_row_m <= sphere2DDataConfig->spectral_modes_m_max);

		int n = i_row_n+rel_n;

		if (n < 0 || n < i_row_m || n > sphere2DDataConfig->spectral_modes_n_max)
			return;

		int idx = rel_n + halosize_off_diagonal;

		SWEET_ASSERT(idx >= 0 && idx < num_diagonals);

		io_row[idx] = i_value;
	}


	/*!
	 * Return reference to an element in the row to the specified value
	 */
	void rowElement_add(
			T *io_row,		//!< pointer to current row
			int i_row_n,	//!< row related to P Legendre mode n
			int i_row_m,	//!< row related to P Fourier mode n
			int rel_n,		//!< Relative Legendre mode n (e.g. -1 or +2)
			const T &i_value
	)
	{
		if (i_row_n < i_row_m)
			return;

//		SWEET_ASSERT(i_row_n >= i_row_m);
		//SWEET_ASSERT(i_row_m >= 0);
		SWEET_ASSERT(i_row_m <= sphere2DDataConfig->spectral_modes_m_max);

		int n = i_row_n+rel_n;

		if (n < 0 || n < i_row_m || n > sphere2DDataConfig->spectral_modes_n_max)
			return;

		int idx = rel_n + halosize_off_diagonal;

		SWEET_ASSERT(idx >= 0 && idx < num_diagonals);

		io_row[idx] += i_value;
	}


	/*!
	 * Return reference to an element in the row to the specified value
	 */
	void rowElement_add_NEW(
			T *io_row,		//!< pointer to current row
			int i_row_n,	//!< row related to P Legendre mode n
			int i_row_m,	//!< row related to P Fourier mode n
			int rel_n,		//!< Relative Legendre mode n (e.g. -1 or +2)
			const T &i_value
	)
	{
		if (i_row_n < i_row_m)
			return;

		int n = i_row_n+rel_n;

		int idx = rel_n + halosize_off_diagonal;

		SWEET_ASSERT(i_row_m <= sphere2DDataConfig->spectral_modes_m_max);
		SWEET_ASSERT(n >= 0);
		SWEET_ASSERT(n >= i_row_m);
		SWEET_ASSERT(n <= sphere2DDataConfig->spectral_modes_n_max);
		SWEET_ASSERT(idx >= 0 && idx < num_diagonals);

		io_row[idx] += i_value;
	}


	void shutdown()
	{
		if (data != nullptr)
		{
			sweet::Memory::MemBlockAlloc::free(data, sizeof(T)*sphere2DDataConfig->spectral_complex_array_data_number_of_elements*num_diagonals);
			data = nullptr;
		}
	}


	~BandedMatrixGridComplex()
	{
		shutdown();
	}


	void print()
	{
		std::size_t idx = 0;
		for (int m = -sphere2DDataConfig->spectral_modes_m_max; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			std::cout << "Meridional block M=" << m << " with N=[" << m << ", " << sphere2DDataConfig->spectral_modes_n_max << "]" << std::endl;

			for (int n = std::abs(m); n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				if (n == m)
				{
					for (int hn = n-halosize_off_diagonal; hn < n+halosize_off_diagonal+1; hn++)
					{
						std::cout << hn;
						if (hn != n+halosize_off_diagonal)
							std::cout << "\t";
					}
					std::cout << std::endl;
					for (int hn = n-halosize_off_diagonal; hn < n+halosize_off_diagonal+1; hn++)
					{
						std::cout << "*******";
						if (hn != n+halosize_off_diagonal)
							std::cout << "\t";
					}
					std::cout << std::endl;
				}

				for (int i = 0; i < num_diagonals; i++)
				{
					std::cout << data[idx*num_diagonals+i];
					if (i != num_diagonals-1)
						std::cout << "\t";
				}
				std::cout << std::endl;

				idx++;
			}
		}
	}

	void print_mblock(int m)
	{
//		std::size_t idx = 0;
//		for (int m = 0; m <= sphere2DDataConfig->spec_m_max; m++)
		{
			std::size_t idx = sphere2DDataConfig->getArrayIndexByModes_Complex_NCompact(std::abs(m), m);
			std::cout << "Meridional block M=" << m << " with N=[" << m << ", " << sphere2DDataConfig->spectral_modes_n_max << "]" << std::endl;

			for (int n = std::abs(m); n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				if (n == m)
				{
					for (int hn = n-halosize_off_diagonal; hn < n+halosize_off_diagonal+1; hn++)
					{
						std::cout << hn;
						if (hn != n+halosize_off_diagonal)
							std::cout << "\t";
					}
					std::cout << std::endl;
					for (int hn = n-halosize_off_diagonal; hn < n+halosize_off_diagonal+1; hn++)
					{
						std::cout << "*******";
						if (hn != n+halosize_off_diagonal)
							std::cout << "\t";
					}
					std::cout << std::endl;
				}

				for (int i = 0; i < num_diagonals; i++)
				{
					std::cout << data[idx*num_diagonals+i];
					if (i != num_diagonals-1)
						std::cout << "\t";
				}
				std::cout << std::endl;

				idx++;
			}
		}
	}
};

}}

#endif
