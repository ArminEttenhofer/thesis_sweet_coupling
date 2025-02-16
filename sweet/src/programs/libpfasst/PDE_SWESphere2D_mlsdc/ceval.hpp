#ifndef PROGRAMS_LIBPFASST_PDE_SWESPHERE2D_MLSDC_CEVAL_HPP
#define PROGRAMS_LIBPFASST_PDE_SWESPHERE2D_MLSDC_CEVAL_HPP

#include "../interface/Sphere2DDataVars.hpp"
#include "Sphere2DDataCtx.hpp"

/**
 * Write file to data and return string of file name
 */
std::string write_file(
		Sphere2DDataCtx *i_ctx,
		const sweet::Data::Sphere2D::DataSpectral &i_sphere2DData,
		const char* i_name	//!< name of output variable
);

/*
 * Right-hand-side functions called from Fortran
 */

extern "C"
{
// initialization of the variables (initial condition)
void cinitial(
		Sphere2DDataCtx *i_ctx,
		double i_t,
		double i_dtq, 
		Sphere2DDataVars *o_Y
);

// finalizes the time step when libpfasst is done
void cfinal(
		Sphere2DDataCtx *i_ctx,
		Sphere2DDataVars *i_Y,
		int i_nnodes,
		int i_niters
);

// evaluates the explicit piece
void ceval_f1(
		Sphere2DDataVars *i_Y, 
		double i_t, 
		Sphere2DDataCtx *i_ctx,
		Sphere2DDataVars *o_F1
);

// evaluates the first implicit piece
void ceval_f2 (
		Sphere2DDataVars *i_Y,
		double i_t,
		Sphere2DDataCtx *i_ctx,
		Sphere2DDataVars *o_F2
);

// solves the first implicit system
void ccomp_f2 (
		Sphere2DDataVars *io_Y,
		double i_t,
		double i_dtq,
		Sphere2DDataVars *i_Rhs,
		Sphere2DDataCtx *i_ctx,
		Sphere2DDataVars *o_F2
);

// evaluates the second implicit piece
void ceval_f3 (Sphere2DDataVars *i_Y,
		double i_t,
		int i_level,
		Sphere2DDataCtx *i_ctx,
		Sphere2DDataVars *o_F3
);

// solves the second implicit system
void ccomp_f3 (Sphere2DDataVars *io_Y,
		double i_t,
		double i_dtq,
		int i_level,
		Sphere2DDataVars *i_Rhs,
		Sphere2DDataCtx *i_ctx,
		Sphere2DDataVars *o_F3
);

}

#endif
