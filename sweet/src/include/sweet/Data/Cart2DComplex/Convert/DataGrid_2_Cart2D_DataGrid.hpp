/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_CART2DCOMPLEX_CONVERT_DATAGRID_2_CART2D_DATAGRID_HPP
#define INCLUDE_SWEET_DATA_CART2DCOMPLEX_CONVERT_DATAGRID_2_CART2D_DATAGRID_HPP

#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Cart2DComplex/DataGrid.hpp>
#include <sweet/Data/Vector/Vector.hpp>

namespace sweet {
namespace Data {
namespace Cart2DComplex {
namespace Convert {


/*!
 * \brief Convert from Cart2DComplex::DataGrid to Cart2D::DataGrid
 */
class DataGrid_2_Cart2D_DataGrid
{
public:
	static
	Cart2D::DataGrid convert_real(
			const Cart2DComplex::DataGrid &i_cart2DData
	)
	{
		Cart2D::DataGrid out(i_cart2DData.cart2DDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < out.cart2DDataConfig->grid_number_elements; i++)
			out.grid_space_data[i] = i_cart2DData.grid_space_data[i].real();

		return out;
	}



public:
	static
	Cart2D::DataGrid convert_imag(
			const Cart2DComplex::DataGrid &i_cart2DData
	)
	{
		Cart2D::DataGrid out(i_cart2DData.cart2DDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < out.cart2DDataConfig->grid_number_elements; i++)
			out.grid_space_data[i] = i_cart2DData.grid_space_data[i].imag();


		return out;
	}

};

}}}}

#endif
