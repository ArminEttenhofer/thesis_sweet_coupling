#ifndef PROGRAMS_LIBPFASST_PDE_SWESPHERE2D_EXPL_SDC_CEVAL_HPP
#define PROGRAMS_LIBPFASST_PDE_SWESPHERE2D_EXPL_SDC_CEVAL_HPP

#include "../interface/Sphere2DDataVars.hpp"
#include "Sphere2DDataCtxSDC.hpp"

/**
 * Determine if output should be written for current time step & iteration
 */
bool timestep_check_output(Sphere2DDataCtxSDC *i_ctx,
						   int i_current_iter,
						   int i_niters);

/**
 * Write file to data and return string of file name
 */

std::string write_file(
	Sphere2DDataCtxSDC *i_ctx,
	const sweet::Data::Sphere2D::DataSpectral &i_sphere2DData,
	const char *i_name //!< name of output variable
);

/*
  Right-hand-side functions called from Fortran
*/

extern "C"
{
	// initialization of the variables (initial condition)
	void cinitial(
		Sphere2DDataCtxSDC *i_ctx,
		double i_t,
		double i_dt,
		Sphere2DDataVars *o_Y);

	// finalizes the time step when libpfasst is done
	void cfinal(
		Sphere2DDataCtxSDC *i_ctx,
		Sphere2DDataVars *i_Y,
		int i_nnodes,
		int i_niters);

	// evaluates the explicit piece
	void ceval(
		Sphere2DDataVars *i_Y,
		double i_t,
		Sphere2DDataCtxSDC *i_ctx,
		Sphere2DDataVars *o_F1);
}

#endif
