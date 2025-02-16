/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_IRK_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_IRK_HPP



#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <complex>

#include "../TimeHelpers/SphBandedMatrix_GridReal.hpp"
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_lg_erk.hpp"


class PDESWESphere2DTS_lg_irk	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

public:
	bool setup(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_timestep_order,
		double i_timestepSize
	);

public:
	bool setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_timestep_order,
		double i_timestepSize,
		double i_crank_nicolson_damping_factor
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;


private:
	PDESWESphere2DTS_lg_erk lg_erk;

	//! alpha/beta (time step related component for implicit solver)
	double alpha;
	double beta;

	//! Crank-Nicolson damping factor
	double crank_nicolson_damping_factor = 0.5;

	//! timestep size
	double timestep_size;

	//! earth radius
	double r;

	//! inverse of earth radius
	double inv_r;

	//! Average geopotential
	double gh;


public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override
	{
		PDESWESphere2DTS_BaseInterface::shackRegistration(io_shackDict);

		lg_erk.shackRegistration(io_shackDict);
		return true;
	}


public:
	PDESWESphere2DTS_lg_irk();


public:
	void update_coefficients();


public:
	void runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi,
		sweet::Data::Sphere2D::DataSpectral &io_vort,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt = 0,
		double i_simulation_timestamp = -1
	) override;

	virtual ~PDESWESphere2DTS_lg_irk();
};


#endif
