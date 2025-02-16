/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_CART2DCOMPLEX_CONVERT_DATASPECTRAL_2_CART2D_DATASPECTRAL_HPP
#define INCLUDE_SWEET_DATA_CART2DCOMPLEX_CONVERT_DATASPECTRAL_2_CART2D_DATASPECTRAL_HPP

#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/DataGrid.hpp>
#include <sweet/Data/Cart2DComplex/DataSpectral.hpp>
#include <sweet/Data/Vector/Vector.hpp>

namespace sweet {
namespace Data {
namespace Cart2DComplex {
namespace Convert {


/*!
 * \brief Convert from Cart2DComplex::DataSpectral to Cart2D::DataSpectral
 */
class DataSpectral_2_Cart2D_DataSpectral
{
public:
	static
	Cart2D::DataSpectral convert_real(
			const Cart2DComplex::DataSpectral &i_cart2DData
	)
	{
		Cart2DComplex::DataGrid tmp_cplx = i_cart2DData.toGrid();
		Cart2D::DataGrid tmp(i_cart2DData.cart2DDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < tmp_cplx.cart2DDataConfig->grid_number_elements; i++)
			tmp.grid_space_data[i] = tmp_cplx.grid_space_data[i].real();

		return Cart2D::DataSpectral(tmp);
	}

public:
	static
	Cart2D::DataSpectral convert_imag(
			const Cart2DComplex::DataSpectral &i_cart2DData
	)
	{
		Cart2DComplex::DataGrid tmp_cplx = i_cart2DData.toGrid();
		Cart2D::DataGrid tmp(i_cart2DData.cart2DDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < tmp_cplx.cart2DDataConfig->grid_number_elements; i++)
			tmp.grid_space_data[i] = tmp_cplx.grid_space_data[i].imag();

		return Cart2D::DataSpectral(tmp);
	}

public:
	static
	Cart2D::DataSpectral convert_grid_real_only(
			const Cart2DComplex::DataSpectral &i_cart2DData
	)
	{
		Cart2D::DataSpectral out(i_cart2DData.cart2DDataConfig);
		out.spectral_setZero();

		for (int r = 0; r < 2; r++)
		{
			for (	std::size_t j = out.cart2DDataConfig->spectral_data_iteration_ranges[r][1][0];
					j < out.cart2DDataConfig->spectral_data_iteration_ranges[r][1][1];
					j++
			) {
				for (	std::size_t i = out.cart2DDataConfig->spectral_data_iteration_ranges[r][0][0];
						i < out.cart2DDataConfig->spectral_data_iteration_ranges[r][0][1];
						i++
				) {
					const std::complex<double> &data = i_cart2DData.spectral_get(j, i);
					out.spectral_set(j, i, data);
				}
			}
		}

		return out;
	}

};

}}}}

#endif
