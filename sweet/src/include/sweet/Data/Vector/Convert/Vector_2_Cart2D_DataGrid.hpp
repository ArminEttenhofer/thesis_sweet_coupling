/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_VECTOR_CONVERT_VECTOR_2_CART2D_DATAGRID_HPP
#define INCLUDE_SWEET_DATA_VECTOR_CONVERT_VECTOR_2_CART2D_DATAGRID_HPP

#include <sweet/Data/Cart2D/Config.hpp>
#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Vector/Vector.hpp>


namespace sweet {
namespace Data {
namespace Vector {
namespace Convert {


class Vector_2_Cart2D_DataGrid
{
public:
	static
	Cart2D::DataGrid convert(
			const Vector<double> &i_scalarDataArray,
			const Cart2D::Config *i_cart2DDataConfig
	)
	{
		Cart2D::DataGrid out(i_cart2DDataConfig);

		for (std::size_t i = 0; i < out.cart2DDataConfig->grid_number_elements; i++)
			out.grid_space_data[i] = i_scalarDataArray.data[i];

		return out;
	}
};

}}}}

#endif
