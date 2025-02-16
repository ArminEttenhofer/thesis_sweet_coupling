#ifndef PROGRAMS_LIBPFASST_PDE_SWESPHERE2D_EXPL_SDC_CENCAP_HPP
#define PROGRAMS_LIBPFASST_PDE_SWESPHERE2D_EXPL_SDC_CENCAP_HPP

#include <iomanip>
#include "../interface/Sphere2DDataVars.hpp"
#include "Sphere2DDataCtxSDC.hpp"

/*
  "Encap" functions called from Fortran to manipulate Sphere2DData
*/

extern "C"
{
	// instantiates and returns the sweet data encapsulated object
	void c_sweet_data_create(
		Sphere2DDataCtxSDC *i_ctx,
		int i_level,
		Sphere2DDataVars **o_Y,
		int *o_size);

	// calls the destructor of the sweet data encapsulated object
	void c_sweet_data_destroy(
		Sphere2DDataVars *i_Y);

	// sets the value of the sweet data encapsulated object
	void c_sweet_data_setval(
		Sphere2DDataVars *io_Y,
		double i_val);

	// copies i_src into o_dst
	void c_sweet_data_copy(
		Sphere2DDataVars *i_src,
		Sphere2DDataVars *o_dst);

	// computes the norm of the sweet data encapsulated object
	void c_sweet_data_norm(
		Sphere2DDataVars *io_Y,
		double *o_val);

	// packs all the values contained in the sweet data object into a flat array
	void c_sweet_data_pack(
		Sphere2DDataVars *io_Y,
		double **o_flat_data_ptr);
	// unpacks the flat array into the sweet data object
	void c_sweet_data_unpack(
		double **i_flat_data_ptr,
		Sphere2DDataVars *io_Y);

	// computes io_Y = i_a * i_X + io_Y
	void c_sweet_data_saxpy(
		double i_a,
		Sphere2DDataVars *i_X,
		Sphere2DDataVars *io_Y);

	// prints the data to the terminal
	void c_sweet_data_eprint(
		Sphere2DDataVars *i_Y);
}

#endif
