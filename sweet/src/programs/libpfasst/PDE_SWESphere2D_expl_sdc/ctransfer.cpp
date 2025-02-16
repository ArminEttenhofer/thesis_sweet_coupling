#include "ctransfer.hpp"

#include "../../PDE_SWESphere2D/BenchmarksCombined.hpp"

#include "cencap.hpp"


extern "C"
{
  void c_sweet_data_restrict(
				 Sphere2DDataVars *io_Y_coarse, 
				 Sphere2DDataVars *i_Y_fine, 
				 int i_level_coarse,
				 int i_level_fine, 
				 Sphere2DDataCtxSDC *i_ctx,
				 double i_t) 
  {
	const sweet::Data::Sphere2D::DataSpectral& phi_pert_Y_fine  = i_Y_fine->get_phi_pert();
	const sweet::Data::Sphere2D::DataSpectral& vrt_Y_fine = i_Y_fine->get_vrt();
	const sweet::Data::Sphere2D::DataSpectral& div_Y_fine  = i_Y_fine->get_div();

	sweet::Data::Sphere2D::DataSpectral& phi_pert_Y_coarse  = io_Y_coarse->get_phi_pert();
	sweet::Data::Sphere2D::DataSpectral& vrt_Y_coarse = io_Y_coarse->get_vrt();
	sweet::Data::Sphere2D::DataSpectral& div_Y_coarse  = io_Y_coarse->get_div();

	// restrict the fine variables to the coarse grid and copy into the coarse variables
	phi_pert_Y_coarse  = phi_pert_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_sphere2d_data_config());
	vrt_Y_coarse = vrt_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_sphere2d_data_config());
	div_Y_coarse  = div_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_sphere2d_data_config());
 
  }  

  void c_sweet_data_interpolate(
				Sphere2DDataVars *io_Y_fine, 
				Sphere2DDataVars *i_Y_coarse, 
				int i_level_fine,
				int i_level_coarse,
				Sphere2DDataCtxSDC *i_ctx,
				double i_t) 
  {
	const sweet::Data::Sphere2D::DataSpectral& phi_pert_Y_coarse  = i_Y_coarse->get_phi_pert();
	const sweet::Data::Sphere2D::DataSpectral& vrt_Y_coarse = i_Y_coarse->get_vrt();
	const sweet::Data::Sphere2D::DataSpectral& div_Y_coarse  = i_Y_coarse->get_div();
  
	sweet::Data::Sphere2D::DataSpectral& phi_pert_Y_fine  = io_Y_fine->get_phi_pert();
	sweet::Data::Sphere2D::DataSpectral& vrt_Y_fine = io_Y_fine->get_vrt();
	sweet::Data::Sphere2D::DataSpectral& div_Y_fine  = io_Y_fine->get_div();

	// interpolate the coarse variables on the fine grid and copy into the coarse variables
	phi_pert_Y_fine  = phi_pert_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_sphere2d_data_config());
	vrt_Y_fine = vrt_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_sphere2d_data_config());
	div_Y_fine  = div_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_sphere2d_data_config());

  }
}

