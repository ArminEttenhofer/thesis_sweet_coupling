/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_CART2D_CONVERT_DATASPECTRAL_2_CART2DCOMPLEX_DATASPECTRAL_HPP
#define INCLUDE_SWEET_DATA_CART2D_CONVERT_DATASPECTRAL_2_CART2DCOMPLEX_DATASPECTRAL_HPP

#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/DataGrid.hpp>
#include <sweet/Data/Cart2DComplex/DataSpectral.hpp>
#include <sweet/Data/Vector/Vector.hpp>

namespace sweet {
namespace Data {
namespace Cart2D {
namespace Convert {


/*!
 * \brief Convert from Cart2D::DataSpectral to Cart2D::DataSpectral
 */
class DataSpectral_2_Cart2DComplex_DataSpectral
{
public:
	static
	Cart2DComplex::DataSpectral grid_convert(
			const DataSpectral &i_cart2DData
	)
	{
		DataGrid tmp = i_cart2DData.toGrid();
		Cart2DComplex::DataGrid tmpc(i_cart2DData.cart2DDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < tmp.cart2DDataConfig->grid_number_elements; i++)
			tmpc.grid_space_data[i] = tmp.grid_space_data[i];


		Cart2DComplex::DataSpectral ret(tmpc);
		return ret;
	}


public:
	static
	Cart2DComplex::DataSpectral spectral_convert(
			const DataSpectral &i_cart2DData
	)
	{
		Cart2DComplex::DataSpectral out(i_cart2DData.cart2DDataConfig);
		out.spectral_setZero();

		for (int r = 0; r < 2; r++)
		{
			for (	std::size_t j = out.cart2DDataConfig->spectral_data_iteration_ranges[r][1][0];
					j < out.cart2DDataConfig->spectral_data_iteration_ranges[r][1][1];
					j++
			) {
				SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
				for (	std::size_t i = out.cart2DDataConfig->spectral_data_iteration_ranges[r][0][0];
						i < out.cart2DDataConfig->spectral_data_iteration_ranges[r][0][1];
						i++
				) {
					const std::complex<double> &data = i_cart2DData.spectral_get(j, i);
					out.spectral_set(j, i, data);
				}

				SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
				for (	std::size_t i = out.cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]+1;
						i < out.cart2DDataConfig->spectral_data_iteration_ranges[r][0][1];
						i++
				) {
					const std::complex<double> &data = i_cart2DData.spectral_get(j, i);
					std::complex<double> data2 = data;

					data2.imag(-data2.imag());
					if (j == 0)
						out.spectral_set(j, out.cart2DDataConfig->spectral_complex_data_size[0]-i, data2);
					else
						out.spectral_set(out.cart2DDataConfig->spectral_complex_data_size[1]-j, out.cart2DDataConfig->spectral_complex_data_size[0]-i, data2);
				}
			}
		}

		return out;
	}
};

}}}}

#endif
