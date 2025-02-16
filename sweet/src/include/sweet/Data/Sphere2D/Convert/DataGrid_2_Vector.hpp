/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_SPHERE2D_CONVERT_DATAGRID_2_VECTOR_HPP
#define INCLUDE_SWEET_DATA_SPHERE2D_CONVERT_DATAGRID_2_VECTOR_HPP

#include <sweet/Parallelization/openmp_helper.hpp>
#include <sweet/Data/Sphere2D/DataGrid.hpp>
#include <sweet/Data/Vector/Vector.hpp>



namespace sweet {
namespace Data {
namespace Sphere2D {
namespace Convert {


/*!
 * \brief Convert from Sphere2D::DataGrid to Vector
 */
class DataGrid_2_Vector
{
public:
	static
	Vector::Vector<double> convert(
			const Sphere2D::DataGrid &i_sphere2DData
	)
	{
		Vector::Vector<double> out(i_sphere2DData.sphere2DDataConfig->grid_number_elements);

		for (std::size_t i = 0; i < out.numberOfElements; i++)
			out.data[i] = i_sphere2DData.grid_space_data[i];

		return out;
	}

public:
	static
	Vector::Vector<double> convert(
			const Sphere2D::DataGrid &i_sphere2DData,
			Vector::Vector<double> &o_out
	)
	{
		for (std::size_t i = 0; i < o_out.numberOfElements; i++)
			o_out.data[i] = i_sphere2DData.grid_space_data[i];

		return o_out;
	}
};

}}}}

#endif
