/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_SPHERE2D_CONVERT_DATAGRID_2_CART2D_DATAGRID_HPP
#define INCLUDE_SWEET_DATA_SPHERE2D_CONVERT_DATAGRID_2_CART2D_DATAGRID_HPP

#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Sphere2D/DataGrid.hpp>

namespace sweet {
namespace Data {
namespace Sphere2D {
namespace Convert {

/*!
 * \brief Convert from Sphere2D::DataGrid to Cart2D::DataGrid
 */
class DataGrid_2_Cart2D_DataGrid
{
public:
	static
	Cart2D::DataGrid convert(
			const Sphere2D::DataGrid &i_sphere2DData,
			Cart2D::Config &i_cart2DDataConfig
	)
	{
		return convert(i_sphere2DData, &i_cart2DDataConfig);
	}

public:
	static
	Cart2D::DataGrid convert(
			const Sphere2D::DataGrid &i_sphere2DData,
			Cart2D::Config *i_cart2DDataConfig
	)
	{
		SWEET_ASSERT(i_sphere2DData.sphere2DDataConfig->grid_num_lon == (int)i_cart2DDataConfig->grid_res[0]);
		SWEET_ASSERT(i_sphere2DData.sphere2DDataConfig->grid_num_lat == (int)i_cart2DDataConfig->grid_res[1]);
		SWEET_ASSERT(i_cart2DDataConfig->grid_number_elements == i_sphere2DData.sphere2DDataConfig->grid_number_elements);

		Cart2D::DataGrid out(i_cart2DDataConfig);


#if SPHERE2D_DATA_GRID_LAYOUT	== SPHERE2D_DATA_LAT_CONTIGUOUS

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int i = 0; i < i_sphere2DData.sphere2DDataConfig->grid_num_lon; i++)
			for (int j = 0; j < i_sphere2DData.sphere2DDataConfig->grid_num_lat; j++)
				out.grid_space_data[(i_sphere2DData.sphere2DDataConfig->grid_num_lat-1-j)*i_sphere2DData.sphere2DDataConfig->grid_num_lon + i] = i_sphere2DData.grid_space_data[i*i_sphere2DData.sphere2DDataConfig->grid_num_lat + j];
#else

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int j = 0; j < i_sphere2DData.sphere2DDataConfig->grid_num_lat; j++)
			for (int i = 0; i < i_sphere2DData.sphere2DDataConfig->grid_num_lon; i++)
				out.grid_space_data[(i_sphere2DData.sphere2DDataConfig->grid_num_lat-1-j)*i_sphere2DData.sphere2DDataConfig->grid_num_lon + i] = i_sphere2DData.grid_space_data[j*i_sphere2DData.sphere2DDataConfig->grid_num_lon + i];
#endif

		return out;
	}
};

}}}}

#endif
