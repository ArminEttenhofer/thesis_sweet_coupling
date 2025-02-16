#ifndef PROGRAMS_LIBPFASST_PDE_SWESPHERE2D_MLSDC_CTRANSFER_HPP
#define PROGRAMS_LIBPFASST_PDE_SWESPHERE2D_MLSDC_CTRANSFER_HPP

#include "../interface/Sphere2DDataVars.hpp"
#include "Sphere2DDataCtx.hpp"

extern "C"
{
  void c_sweet_data_restrict(
				 Sphere2DDataVars *io_Y_coarse, 
				 Sphere2DDataVars *i_Y_fine, 
				 int i_level_coarse,
				 int i_level_fine, 
				 Sphere2DDataCtx *i_ctx,
				 double i_t
				 );

  void c_sweet_data_interpolate(
				Sphere2DDataVars *io_Y_fine, 
				Sphere2DDataVars *i_Y_coarse, 
				int i_level_fine,
				int i_level_coarse,
				Sphere2DDataCtx *i_ctx,
				double i_t
				);

}

#endif
