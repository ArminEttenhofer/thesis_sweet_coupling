/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_SPHERE2D_CONVERT_DATASPECTRAL_2_SPHERE2DCOMPLEX_DATASPECTRAL_HPP
#define INCLUDE_SWEET_DATA_SPHERE2D_CONVERT_DATASPECTRAL_2_SPHERE2DCOMPLEX_DATASPECTRAL_HPP

#include <sweet/Data/Sphere2D/DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/DataSpectral.hpp>
#include <sweet/Data/Vector/Vector.hpp>


namespace sweet {
namespace Data {
namespace Sphere2D {
namespace Convert {


/*!
 * \brief Convert from Sphere2D::DataSpectral to Sphere2DComplex::DataSpectral
 */
class DataSpectral_2_Sphere2DComplex_DataSpectral
{
public:
	static
	Sphere2DComplex::DataSpectral convert(
			const Sphere2D::DataSpectral &i_sphere2DData
	)
	{
		Sphere2D::DataGrid tmp = i_sphere2DData.toGrid();
		Sphere2DComplex::DataGrid tmpc(i_sphere2DData.sphere2DDataConfig);

#if SWEET_THREADING_SPACE
#pragma omp parallel for
#endif
		for (std::size_t i = 0; i < tmp.sphere2DDataConfig->grid_number_elements; i++)
			tmpc.grid_space_data[i] = tmp.grid_space_data[i];


		Sphere2DComplex::DataSpectral ret(tmpc);
		return ret;
	}
};

}}}}

#endif
