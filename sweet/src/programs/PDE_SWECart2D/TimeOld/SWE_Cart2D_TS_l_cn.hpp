/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWECART2D_TIME_SWE_CART2D_TS_L_CN_HPP
#define PROGRAMS_PDE_SWECART2D_TIME_SWE_CART2D_TS_L_CN_HPP

#include <programs/PDE_SWECart2D/TimeOld/PDESWECart2DTS_BaseInterface.hpp>
#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_erk.hpp>
#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_erk.hpp>
#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_irk.hpp>
#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKCart2DData.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include "PDESWECart2DTS_BaseInterface.hpp"

class SWE_Cart2D_TS_l_cn	: public PDESWECart2DTS_BaseInterface
{
	double crank_nicolson_damping_factor = 0.5;

	SWE_Cart2D_TS_l_erk ts_l_erk;
	SWE_Cart2D_TS_l_irk ts_l_irk;

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	);

private:
	void backward_euler_timestep_linear(
			sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables
			double i_dt
	);



public:

	bool setup(sweet::Data::Cart2D::Operators *io_ops);

	void runTimestep(
			sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
			sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

	virtual ~SWE_Cart2D_TS_l_cn() {}
};

#endif
