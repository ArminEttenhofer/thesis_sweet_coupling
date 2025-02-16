/*
 * Cart2DDataGridMapping.hpp
 *
 *  Created on: 18 Jul 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_CART2D_GRIDMAPPING_HPP
#define INCLUDE_SWEET_DATA_CART2D_GRIDMAPPING_HPP

#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Cart2D/DataSampler.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2D/Staggering.hpp>
#include <sweet/Data/Vector/Vector.hpp>
#include <sweet/Shacks/Dictionary.hpp>

namespace sweet {
namespace Data {
namespace Cart2D {


/*!
 * \brief Map one Cart2D data field to another one by using some interpolation
 */
class GridMapping
{
public:
	//(x,y) grid points, refers to lower left corner of cells
	Vector::Vector<double> pos_ll_x, pos_ll_y;

	// Interpolation stuff
	DataSampler sampler2D;

	// Staggering
	Staggering staggering;

public:
	GridMapping()
	{
	}


	void setup(
			Shack *io_cart2DDataOps,
			Config *i_cart2DDataConfig
	)
	{
		// ll  refers to lower left corner of the cell.
		pos_ll_x.setup(i_cart2DDataConfig->grid_number_elements);
		pos_ll_y.setup(i_cart2DDataConfig->grid_number_elements);

		std::size_t idx = 0;
		for (std::size_t j = 0; j < i_cart2DDataConfig->grid_res[1]; j++)
		{
			for (std::size_t i = 0; i < i_cart2DDataConfig->grid_res[0]; i++)
			{
				pos_ll_x.data[idx] = ((double)i)*io_cart2DDataOps->cart2d_domain_size[0]/(double)io_cart2DDataOps->space_res_physical[0];
				pos_ll_y.data[idx] = ((double)j)*io_cart2DDataOps->cart2d_domain_size[1]/(double)io_cart2DDataOps->space_res_physical[1];
				idx++;
			}
		}

		// Setup sampler for future interpolations
		sampler2D.setup(io_cart2DDataOps->cart2d_domain_size, i_cart2DDataConfig);

		if (io_cart2DDataOps->space_grid_use_c_staggering)
			staggering.setup_c_staggering();
		else
			staggering.setup_a_staggering();
	}


	void mapCtoA_u(
			const DataGrid &i_src,
			DataGrid &o_dst
	)
	{
		// remap solution to A grid
		sampler2D.bicubic_scalar(
				i_src,
				pos_ll_x,
				pos_ll_y,
				o_dst,
				staggering.u[0],
				staggering.u[1]
			);
	}


	void mapCtoA_v(
			const DataGrid &i_src,
			DataGrid &o_dst
	)
	{
		// remap solution to A grid
		sampler2D.bicubic_scalar(
				i_src,
				pos_ll_x,
				pos_ll_y,
				o_dst,
				staggering.v[0],
				staggering.v[1]
			);
	}




	void mapAtoC_u(
			const DataGrid &i_src,
			DataGrid &o_dst
	)
	{
		// remap solution to C grid
		sampler2D.bicubic_scalar(i_src, pos_ll_x, pos_ll_y, o_dst, -staggering.u[0], -staggering.u[1]);
	}


	void mapAtoC_v(
			const DataGrid &i_src,
			DataGrid &o_dst
	)
	{
		// remap solution to C grid
		sampler2D.bicubic_scalar(i_src, pos_ll_x, pos_ll_y, o_dst, -staggering.v[0], -staggering.v[1]);
	}




	void mapCtoA_u(
			const DataSpectral &i_src,
			DataSpectral &o_dst
	)
	{
		DataGrid i_src_phys = i_src.toGrid();
		DataGrid o_dst_phys = o_dst.toGrid();

		// remap solution to A grid
		sampler2D.bicubic_scalar(i_src_phys, pos_ll_x, pos_ll_y, o_dst_phys, staggering.u[0], staggering.u[1]);

		o_dst.loadCart2DDataGrid(o_dst_phys);
	}


	void mapCtoA_v(
			const DataSpectral &i_src,
			DataSpectral &o_dst
	)
	{
		DataGrid i_src_phys = i_src.toGrid();
		DataGrid o_dst_phys = o_dst.toGrid();

		// remap solution to A grid
		sampler2D.bicubic_scalar(i_src_phys, pos_ll_x, pos_ll_y, o_dst_phys, staggering.v[0], staggering.v[1]);

		o_dst.loadCart2DDataGrid(o_dst_phys);
	}




	void mapAtoC_u(
			const DataSpectral &i_src,
			DataSpectral &o_dst
	)
	{
		DataGrid i_src_phys = i_src.toGrid();
		DataGrid o_dst_phys = o_dst.toGrid();

		// remap solution to C grid
		sampler2D.bicubic_scalar(i_src_phys, pos_ll_x, pos_ll_y, o_dst_phys, -staggering.u[0], -staggering.u[1]);

		o_dst.loadCart2DDataGrid(o_dst_phys);
	}


	void mapAtoC_v(
			const DataSpectral &i_src,
			DataSpectral &o_dst
	)
	{
		DataGrid i_src_phys = i_src.toGrid();
		DataGrid o_dst_phys = o_dst.toGrid();

		// remap solution to C grid
		sampler2D.bicubic_scalar(i_src_phys, pos_ll_x, pos_ll_y, o_dst_phys, -staggering.v[0], -staggering.v[1]);

		o_dst.loadCart2DDataGrid(o_dst_phys);
	}


};

}}}

#endif
