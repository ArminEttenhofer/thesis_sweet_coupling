/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_CART2D_CONVERT_DATASPECTRAL_2_CART2D_DATAGRID_HPP
#define INCLUDE_SWEET_DATA_CART2D_CONVERT_DATASPECTRAL_2_CART2D_DATAGRID_HPP

#include <sweet/Data/Cart2D/Config.hpp>
#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>

namespace sweet {
namespace Data {
namespace Cart2D {
namespace Convert {


/*!
 * \brief Convert from Cart2D::DataSpectral to Cart2D::DataGrid
 */
class DataSpectral_2_Cart2D_DataGrid
{
public:
	static
	void convert(
			const DataSpectral& i_cart2DDataSpectral,
			DataGrid& o_cart2DDataGrid
	)
	{
		const Config* cart2DDataConfig = i_cart2DDataSpectral.cart2DDataConfig;

		/*
		 * Warning: The FFTW functions are in-situ operations.
		 * Therefore, the data in the source array will be destroyed.
		 * Hence, we create a copy
		 */
		DataSpectral tmp(i_cart2DDataSpectral);
		cart2DDataConfig->fft_spectral_2_grid_INPLACE(tmp.spectral_space_data, o_cart2DDataGrid.grid_space_data);
	}
};

}}}}

#endif
