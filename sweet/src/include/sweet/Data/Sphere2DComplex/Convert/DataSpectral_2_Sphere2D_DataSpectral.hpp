/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_SPHERE2DCOMPLEX_CONVERT_DATASPECTRAL_2_DATASPECTRAL_HPP
#define INCLUDE_SWEET_DATA_SPHERE2DCOMPLEX_CONVERT_DATASPECTRAL_2_DATASPECTRAL_HPP

#include <sweet/Data/Sphere2D/DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/DataGrid.hpp>
#include <sweet/Data/Vector/Vector.hpp>


namespace sweet {
namespace Data {
namespace Sphere2DComplex {
namespace Convert {

class DataSpectral_2_Sphere2D_DataSpectral
{
public:
	static
	Sphere2D::DataSpectral convert_real(
			const Sphere2DComplex::DataSpectral &i_sphere2DData
	)
	{
		Sphere2DComplex::DataGrid tmp_cplx = i_sphere2DData.toGrid();
		Sphere2D::DataGrid tmp(i_sphere2DData.sphere2DDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < tmp_cplx.sphere2DDataConfig->grid_number_elements; i++)
			tmp.grid_space_data[i] = tmp_cplx.grid_space_data[i].real();

		return Sphere2D::DataSpectral(tmp);
	}



public:
	static
	Sphere2D::DataSpectral convert_imag(
			const Sphere2DComplex::DataSpectral &i_sphere2DData
	)
	{
		Sphere2DComplex::DataGrid tmp_cplx = i_sphere2DData.toGrid();
		Sphere2D::DataGrid tmp(i_sphere2DData.sphere2DDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD
		for (std::size_t i = 0; i < tmp_cplx.sphere2DDataConfig->grid_number_elements; i++)
			tmp.grid_space_data[i] = tmp_cplx.grid_space_data[i].imag();

		return Sphere2D::DataSpectral(tmp);
	}
};

}}}}

#endif
