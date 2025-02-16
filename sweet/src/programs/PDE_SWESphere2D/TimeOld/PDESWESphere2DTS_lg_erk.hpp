/*
 * PDESWESphere2DTS_lg_erk.hpp
 *
 *  Created on: 30 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */


#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_ERK_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_ERK_HPP

#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKSphere2DData.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include "PDESWESphere2DTS_BaseInterface.hpp"



class PDESWESphere2DTS_lg_erk	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			int i_timestepping_order
		);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackPDESWETimeDisc->timestepping_order;
		timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
		return i_timestepping_method == "lg_erk";
	}

	std::string getIDString() override
	{
		return "lg_erk";
	}


private:
	// Sampler
	sweet::DEPRECATED_TimesteppingExplicitRKSphere2DData timestepping_rk;

public:
	void euler_timestep_update(
			const sweet::Data::Sphere2D::DataSpectral &i_phi,
			const sweet::Data::Sphere2D::DataSpectral &i_vort,
			const sweet::Data::Sphere2D::DataSpectral &i_div,

			sweet::Data::Sphere2D::DataSpectral &o_phi_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_vort_t,	//!< time updates
			sweet::Data::Sphere2D::DataSpectral &o_div_t,	//!< time updates

			double i_simulation_timestamp = -1
	);

public:
	PDESWESphere2DTS_lg_erk();

	bool setup(
			int i_order	//!< order of RK time stepping method
	);

	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vort,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;

	virtual ~PDESWESphere2DTS_lg_erk();
};

#endif
