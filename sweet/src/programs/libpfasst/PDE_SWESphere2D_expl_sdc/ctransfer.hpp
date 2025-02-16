#ifndef PROGRAMS_LIBPFASST_PDE_SWESPHERE2D_EXPL_SDC_CTRANSFER_HPP
#define PROGRAMS_LIBPFASST_PDE_SWESPHERE2D_EXPL_SDC_CTRANSFER_HPP

#include "../interface/Sphere2DDataVars.hpp"
#include "Sphere2DDataCtxSDC.hpp"

extern "C"
{
  void c_sweet_data_restrict(
				 Sphere2DDataVars *io_Y_coarse, 
				 Sphere2DDataVars *i_Y_fine, 
				 int i_level_coarse,
				 int i_level_fine, 
				 Sphere2DDataCtxSDC *i_ctx,
				 double i_t
				 );

  void c_sweet_data_interpolate(
				Sphere2DDataVars *io_Y_fine, 
				Sphere2DDataVars *i_Y_coarse, 
				int i_level_fine,
				int i_level_coarse,
				Sphere2DDataCtxSDC *i_ctx,
				double i_t
				);

}

#endif
