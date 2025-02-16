/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_CART2D_CONVERT_DATAGRID_2_VECTOR_HPP
#define INCLUDE_SWEET_DATA_CART2D_CONVERT_DATAGRID_2_VECTOR_HPP

#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Vector/Vector.hpp>

namespace sweet {
namespace Data {
namespace Cart2D {
namespace Convert {


/*!
 * \brief Convert from Cart2D::DataGrid to Vector
 */
class DataGrid_2_Vector
{
public:
	static
	Vector::Vector<double> convert(
			const Cart2D::DataGrid &i_cart2DData,
			bool i_raise_error_if_spectral = true
	)
	{
		Vector::Vector<double> out(i_cart2DData.cart2DDataConfig->grid_number_elements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < out.numberOfElements; i++)
			out.data[i] = i_cart2DData.grid_space_data[i];

		return out;
	}
};

}}}}

#endif
