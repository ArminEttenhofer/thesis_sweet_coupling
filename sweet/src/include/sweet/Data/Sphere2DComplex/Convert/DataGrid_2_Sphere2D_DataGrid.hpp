/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_SPHERE2DCOMPLEX_CONVERT_DATAGRID_2_SPHERE2D_DATAGRID_HPP
#define INCLUDE_SWEET_DATA_SPHERE2DCOMPLEX_CONVERT_DATAGRID_2_SPHERE2D_DATAGRID_HPP

#include <sweet/Data/Sphere2DComplex/DataGrid.hpp>
#include <sweet/Data/Sphere2D/DataGrid.hpp>
#include <sweet/Data/Sphere2DComplex/DataSpectral.hpp>
#include <sweet/Data/Vector/Vector.hpp>


namespace sweet {
namespace Data {
namespace Sphere2DComplex {
namespace Convert {

class DataGrid_2_Sphere2D_DataGrid
{
public:
	static
	Sphere2D::DataGrid convert_real(
			const Sphere2DComplex::DataGrid &i_sphere2DData
	)
	{
		Sphere2D::DataGrid out(i_sphere2DData.sphere2DDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < out.sphere2DDataConfig->grid_number_elements; i++)
			out.grid_space_data[i] = i_sphere2DData.grid_space_data[i].real();

		return out;
	}



public:
	static
	Sphere2D::DataGrid convert_imag(
			const Sphere2DComplex::DataGrid &i_sphere2DData
	)
	{
		Sphere2D::DataGrid out(i_sphere2DData.sphere2DDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < out.sphere2DDataConfig->grid_number_elements; i++)
			out.grid_space_data[i] = i_sphere2DData.grid_space_data[i].imag();


		return out;
	}
};

}}}}

#endif
