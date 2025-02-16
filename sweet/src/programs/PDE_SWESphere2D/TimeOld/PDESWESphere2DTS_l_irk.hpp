/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_IRK_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_IRK_HPP


#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <complex>

#include "../TimeHelpers/SphBandedMatrix_GridReal.hpp"
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_l_erk.hpp"
#include "PDESWESphere2DTS_lg_erk.hpp"



/**
 * Implicit solver
 */
class PDESWESphere2DTS_l_irk	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup(
			const sweet::Data::Sphere2D::Operators *io_ops,
			int i_timestep_order,
			double i_timestepSize
	);

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			int i_timestep_order,
			double i_timestepSize,
			double i_crank_nicolson_damping_factor,
			bool i_no_coriolis
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;

	std::string timestepping_method;

private:
	PDESWESphere2DTS_lg_erk swe_sphere2d_ts_lg_erk;
	PDESWESphere2DTS_l_erk swe_sphere2d_ts_l_erk;

	SphBandedMatrix_GridReal sphSolverDiv;

	double crank_nicolson_damping_factor;

	//! timestep size
	double timestep_size;

	//! individual time step size
	double dt_explicit = -1;
	double dt_implicit = -1;

	//! earth radius
	double sphere2d_radius;

	bool use_f_sphere2D;

	bool no_coriolis;

	//! f0
	double f0;

	//! Coriolis effect
	double two_coriolis;

	sweet::Data::Sphere2D::DataGrid mug;

public:
	PDESWESphere2DTS_l_irk();

public:
	void update_coefficients(double i_timestepSize);

public:
	void clear();

public:
	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;

	void solveImplicit(
		sweet::Data::Sphere2D::DataSpectral &io_phi,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double dt
	);


	virtual ~PDESWESphere2DTS_l_irk();
};


#endif
