/*
 * PDESWESphere2DTS_l_irk_n_erk.hpp
 *
 *  Created on: 30 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_IRK_N_ERK_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_IRK_N_ERK_HPP

#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKSphere2DData.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_l_erk_n_erk.hpp"
#include "PDESWESphere2DTS_l_irk.hpp"



class PDESWESphere2DTS_l_irk_n_erk	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			int i_order,	//!< order of RK time stepping method for linear parts
			int i_order2,	//!< order of RK time stepping method for non-linear parts
			int i_version_id
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;

	double timestep_size;

	/*
	 * Linear time steppers
	 */
	PDESWESphere2DTS_l_irk timestepping_l_irk;

	/*
	 * Non-linear time steppers
	 */
	PDESWESphere2DTS_l_erk_n_erk timestepping_l_erk_n_erk;

	sweet::DEPRECATED_TimesteppingExplicitRKSphere2DData timestepping_rk_nonlinear;

	int version_id;

	int timestepping_order;
	int timestepping_order2;

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override
	{
		PDESWESphere2DTS_BaseInterface::shackRegistration(io_shackDict);

		timestepping_l_irk.shackRegistration(io_shackDict);
		timestepping_l_erk_n_erk.shackRegistration(io_shackDict);
		return true;
	}


public:
	PDESWESphere2DTS_l_irk_n_erk();

	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
			sweet::Data::Sphere2D::DataSpectral &io_vort,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;


	virtual ~PDESWESphere2DTS_l_irk_n_erk();
};

#endif
