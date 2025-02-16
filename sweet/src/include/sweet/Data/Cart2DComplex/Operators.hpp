/*
 * Cart2DOperatorsComplex.hpp
 *
 *  Created on: 30 Jun 2015
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef INCLUDE_SWEET_DATA_CART2DCOMPLEX_OPERATORS_HPP
#define INCLUDE_SWEET_DATA_CART2DCOMPLEX_OPERATORS_HPP

#include <sweet/Error/Base.hpp>
#include <sweet/Error/Fatal.hpp>
#include <sweet/Data/Cart2DComplex/DataSpectral.hpp>
#include <sweet/Data/Cart2D/Config.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>

namespace sweet {
namespace Data {
namespace Cart2DComplex {


/*!
 * \brief Operators which can be applied to data in grid or spectral space
 */
class Operators
{
public:
	Error::Base error;
	
	///! Handler to configuration
	Cart2D::Config *cart2DDataConfig = nullptr;

	///! differential operators coefficients in spectral space for derivative along x and y direction
	DataSpectral diff_c_x, diff_c_y;

	///! 2nd differential operators coefficients in spectral space for derivative along x and y direction
	DataSpectral diff2_c_x, diff2_c_y;

	/*!
	 * D2, e.g. for viscosity
	 */
	inline DataSpectral diff2(
			const DataSpectral &i_dataArray
	)
	{
		return diff2_c_x(i_dataArray) + diff2_c_y(i_dataArray);
	}



	/*!
	 *        __2
	 * apply  \/  operator (aka Laplace)
	 */
	inline DataSpectral laplace(
			const DataSpectral &i_a
	)
	{
		return diff2_c_x(i_a)+diff2_c_y(i_a);
	}


	/*!
	 *        __
	 * apply  \/ .  operator
	 */
	inline DataSpectral diff_dot(
			const DataSpectral &i_a
	)
	{
		return diff_c_x(i_a)+diff_c_y(i_a);
	}



	/*!
	 * Diff N operator for hyperviscosity, see
	 *
	 * "Numerical Techniques for Global Atmospheric Models", page 500
	 */
	inline DataSpectral diffN_x(
			const DataSpectral &io_u,
			int i_order
	)
	{
		if (i_order == 0)
			return io_u;

		DataSpectral tu = io_u;

		for (int i = 0; i < i_order/2; i++)
			tu = diff2_c_x(tu);

		if (i_order & 1)
			tu = diff_c_x(tu);

		return tu;
	}


	/*!
	 * Diff N operator for hyperviscosity, see
	 * "Numerical Techniques for Global Atmospheric Models", page 500
	 */
	DataSpectral diffN_y(
			const DataSpectral &io_v,
			int i_order
	)
	{
		if (i_order == 0)
			return io_v;

		DataSpectral tv = io_v;

		for (int i = 0; i < i_order/2; i++)
			tv = diff2_c_y(tv);

		if (i_order & 1)
			tv = diff_c_y(tv);

		return tv;
	}


#if SWEET_USE_CART2D_SPECTRAL_SPACE
	/*!
	 * Diffusion or hyperviscosity coefficients
	 * Simply calculates the spectral coefficients
	 *  see "Numerical Techniques for Global Atmospheric Models", page 500
	 *
	 * i_order (q) needs to be even!!! (second or forth order usually)
	 *
	 * Returns operator D^q
	 *
	 */
	inline DataSpectral diffusion_coefficient(
			int i_order
	)
	{
		//Check if even
		SWEET_ASSERT( i_order % 2 == 0);
		SWEET_ASSERT( i_order > 0);
		DataSpectral out = diff2_c_x+diff2_c_y;

		for (int i = 1; i < i_order/2; i++)
			out = std::pow(-1, i)*(diff2_c_x(out)+diff2_c_y(out));

		return out;
	}

	/*!
	 * Calculates implicit diffusion (applies 1/(1-mu*dt*D^q) to spectrum)
	 *  see "Numerical Techniques for Global Atmospheric Models", page 500
	 *
	 * i_order (q) needs to be even!!! (second or forth order usually)
	 * i_coef is mu*dt
	 *
	 * Only works in spectral space
	 *
	 */
	inline DataSpectral implicit_diffusion(
			const DataSpectral &i_data,
			double i_coef,
			int i_order
	)
	{
		DataSpectral out=i_data;
	
		// Get diffusion coefficients (these are the -mu*dt*D^q, where q is the order
		DataSpectral diff = -i_coef*diffusion_coefficient(i_order);

		// Add 1 to get denominator
		diff = diff.spectral_addScalarAll(1.0);

		// Invert
		diff = diff.spectral_invert();
		// apply to data
		out = diff(out);
		return out;
	}
#endif


	bool setup(
			Cart2D::Config *i_cart2DDataConfig,		//!< data config setup for spectral transformations
			Cart2D::Shack *i_shackCart2DDataOps
	)
	{
		return setup(
				i_cart2DDataConfig,
				i_shackCart2DDataOps->cart2d_domain_size
		);
	}


	bool setup(
			Cart2D::Config &i_cart2DDataConfig,		//!< data config setup for spectral transformations
			Cart2D::Shack &i_shackCart2DDataOps
	)
	{
		return setup(
				&i_cart2DDataConfig,
				i_shackCart2DDataOps.cart2d_domain_size
		);
	}

	bool setup(
			Cart2D::Config &i_cart2DDataConfig,		//!< data config setup for spectral transformations
			Cart2D::Shack *i_shackCart2DDataOps
	)
	{
		return setup(
				&i_cart2DDataConfig,
				i_shackCart2DDataOps->cart2d_domain_size
		);
	}

	bool setup(
			Cart2D::Config *i_cart2DDataConfig,
			const double i_domain_size[2]
	)
	{
		SWEET_ASSERT(cart2DDataConfig == nullptr);
		cart2DDataConfig = i_cart2DDataConfig;


		/*
		 * Setup spectral differential operators
		 * 		diff(e(ix), x)
		 */
#if !SWEET_USE_CART2D_SPECTRAL_SPACE
		std::cerr << "Activate spectral space during compile time to use spectral diffs. Otherwise, the convolution would be freakingly expensive" << std::endl;
		SWEET_ASSERT(false);
#endif

		/*
		 * DIFF X
		 */
		{
			diff_c_x.setup(cart2DDataConfig);
			diff_c_x.spectral_setZero();
			double scale_x = 2.0*M_PI/i_domain_size[0];

			/*
			 * left bottom
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = cart2DDataConfig->spectral_complex_ranges[0][1][0]; j < cart2DDataConfig->spectral_complex_ranges[0][1][1]; j++)
				for (std::size_t i = cart2DDataConfig->spectral_complex_ranges[0][0][0]; i < cart2DDataConfig->spectral_complex_ranges[0][0][1]; i++)
					diff_c_x.spectral_set(j, i, 0, (double)i*scale_x);

			/*
			 * left top
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = cart2DDataConfig->spectral_complex_ranges[1][1][0]; j < cart2DDataConfig->spectral_complex_ranges[1][1][1]; j++)
				for (std::size_t i = cart2DDataConfig->spectral_complex_ranges[1][0][0]; i < cart2DDataConfig->spectral_complex_ranges[1][0][1]; i++)
					diff_c_x.spectral_set(j, i, 0, (double)i*scale_x);

			/*
			 * right bottom
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = cart2DDataConfig->spectral_complex_ranges[2][1][0]; j < cart2DDataConfig->spectral_complex_ranges[2][1][1]; j++)
				for (std::size_t i = cart2DDataConfig->spectral_complex_ranges[2][0][0]; i < cart2DDataConfig->spectral_complex_ranges[2][0][1]; i++)
					diff_c_x.spectral_set(j, i, 0, -(double)(cart2DDataConfig->spectral_complex_ranges[2][0][1]-i)*scale_x);

			/*
			 * right top
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = cart2DDataConfig->spectral_complex_ranges[3][1][0]; j < cart2DDataConfig->spectral_complex_ranges[3][1][1]; j++)
				for (std::size_t i = cart2DDataConfig->spectral_complex_ranges[3][0][0]; i < cart2DDataConfig->spectral_complex_ranges[3][0][1]; i++)
					diff_c_x.spectral_set(j, i, 0, -(double)(cart2DDataConfig->spectral_complex_ranges[3][0][1]-i)*scale_x);
		}

		/*
		 * DIFF Y
		 */
		{
			diff_c_y.setup(cart2DDataConfig);

			diff_c_y.spectral_setZero();
			double scale_y = 2.0*M_PI/i_domain_size[1];


			/*
			 * left bottom
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = cart2DDataConfig->spectral_complex_ranges[0][1][0]; j < cart2DDataConfig->spectral_complex_ranges[0][1][1]; j++)
				for (std::size_t i = cart2DDataConfig->spectral_complex_ranges[0][0][0]; i < cart2DDataConfig->spectral_complex_ranges[0][0][1]; i++)
					diff_c_y.spectral_set(j, i, 0, (double)j*scale_y);


			/*
			 * left top
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = cart2DDataConfig->spectral_complex_ranges[1][1][0]; j < cart2DDataConfig->spectral_complex_ranges[1][1][1]; j++)
				for (std::size_t i = cart2DDataConfig->spectral_complex_ranges[1][0][0]; i < cart2DDataConfig->spectral_complex_ranges[1][0][1]; i++)
					diff_c_y.spectral_set(j, i, 0, -(double)(cart2DDataConfig->spectral_complex_ranges[1][1][1]-j)*scale_y);

			/*
			 * right bottom
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = cart2DDataConfig->spectral_complex_ranges[2][1][0]; j < cart2DDataConfig->spectral_complex_ranges[2][1][1]; j++)
				for (std::size_t i = cart2DDataConfig->spectral_complex_ranges[2][0][0]; i < cart2DDataConfig->spectral_complex_ranges[2][0][1]; i++)
					diff_c_y.spectral_set(j, i, 0, (double)(j)*scale_y);

			/*
			 * right top
			 */
			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t j = cart2DDataConfig->spectral_complex_ranges[3][1][0]; j < cart2DDataConfig->spectral_complex_ranges[3][1][1]; j++)
				for (std::size_t i = cart2DDataConfig->spectral_complex_ranges[3][0][0]; i < cart2DDataConfig->spectral_complex_ranges[3][0][1]; i++)
					diff_c_y.spectral_set(j, i, 0, -(double)(cart2DDataConfig->spectral_complex_ranges[3][1][1]-j)*scale_y);
		}

		/*
		 * 2nd order differential operators
		 */
		/*
		 * DIFF2 X
		 */
		{
			diff2_c_x.setup(cart2DDataConfig);

			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t i = 0; i < cart2DDataConfig->spectral_complex_array_data_number_of_elements; i++)
				diff2_c_x.spectral_space_data[i] = diff_c_x.spectral_space_data[i]*diff_c_x.spectral_space_data[i];
		}


		/*
		 * DIFF2 X
		 */
		{
			diff2_c_y.setup(cart2DDataConfig);

			SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
			for (std::size_t i = 0; i < cart2DDataConfig->spectral_complex_array_data_number_of_elements; i++)
				diff2_c_y.spectral_space_data[i] = diff_c_y.spectral_space_data[i]*diff_c_y.spectral_space_data[i];
		}

		return true;
	}



public:
	Operators()	:
		cart2DDataConfig(nullptr)
	{
	}



public:
	Operators(
		Cart2D::Config *i_cart2DDataConfig,
		const double i_domain_size[2]	//!< domain size
	)	:
		cart2DDataConfig(i_cart2DDataConfig),

		diff_c_x(i_cart2DDataConfig),
		diff_c_y(i_cart2DDataConfig),
		diff2_c_x(i_cart2DDataConfig),
		diff2_c_y(i_cart2DDataConfig)
	{
		setup(i_cart2DDataConfig, i_domain_size);
	}

public:
	Operators(
		Cart2D::Config *i_cart2DDataConfig,		//!< data config setup for spectral transformations
		Cart2D::Shack *i_cart2DDataOps
	)	//:
		//cart2DDataConfig(i_cart2DDataConfig),

		//diff_c_x(i_cart2DDataConfig),
		//diff_c_y(i_cart2DDataConfig),

		//diff2_c_x(i_cart2DDataConfig),
		//diff2_c_y(i_cart2DDataConfig)
	{
		setup(i_cart2DDataConfig, i_cart2DDataOps->cart2d_domain_size);
	}


};

}}}

#endif
