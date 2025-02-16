/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_VECTOR_CONVERT_VECTOR_2_SPHERE2D_DATAGRID_HPP
#define INCLUDE_SWEET_DATA_VECTOR_CONVERT_VECTOR_2_SPHERE2D_DATAGRID_HPP

#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/DataSpectral.hpp>
#include <sweet/Data/Vector/Vector.hpp>


namespace sweet {
namespace Data {
namespace Vector {
namespace Convert {

class Vector_2_Sphere2DDataGrid
{
public:
	static
	Sphere2D::DataGrid convert(
			const Vector<double> &i_scalarDataArray,
			const Sphere2D::Config *i_sphere2DDataConfig
	)
	{
		Sphere2D::DataGrid out(i_sphere2DDataConfig);

		SWEET_ASSERT(out.sphere2DDataConfig->grid_number_elements == i_scalarDataArray.numberOfElements);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < (std::size_t)out.sphere2DDataConfig->grid_number_elements; i++)
			out.grid_space_data[i] = i_scalarDataArray.data[i];

		return out;
	}
};

}}}}

#endif
