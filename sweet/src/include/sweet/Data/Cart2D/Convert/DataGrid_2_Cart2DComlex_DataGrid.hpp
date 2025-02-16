/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_CART2D_CONVERT_DATAGRID_2_CART2DCOMLEX_DATAGRID_HPP
#define INCLUDE_SWEET_DATA_CART2D_CONVERT_DATAGRID_2_CART2DCOMLEX_DATAGRID_HPP

#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Cart2DComplex/DataGrid.hpp>
#include <sweet/Data/Vector/Vector.hpp>
#include <sweet/Parallelization/openmp_helper.hpp>

namespace sweet {
namespace Data {
namespace Cart2D {
namespace Convert {

/*!
 * \brief Convert from Cart2D::DataGrid to Cart2DComplex::DataGrid
 */
class DataGrid_2_Cart2DComplex_DataGrid
{
public:
	static
	Cart2DComplex::DataGrid convert(
			const DataGrid &i_cart2DData
	)
	{
		Cart2DComplex::DataGrid out(i_cart2DData.cart2DDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < i_cart2DData.cart2DDataConfig->grid_number_elements; i++)
			out.grid_space_data[i] = i_cart2DData.grid_space_data[i];

		return out;
	}
};

}}}}

#endif
