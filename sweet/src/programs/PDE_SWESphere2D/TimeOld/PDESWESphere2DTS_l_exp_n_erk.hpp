/*
 * PDESWESphere2DTS_l_rexi_n_erk.hpp
 *
 *  Created on: 30 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_EXP_N_ERK_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_EXP_N_ERK_HPP

#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKSphere2DData.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_l_erk_n_erk.hpp"
#include "PDESWESphere2DTS_l_exp.hpp"



class PDESWESphere2DTS_l_exp_n_erk	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			sweet::ExpIntegration::Shack *i_shackExpIntegration,
			const std::string &i_exp_method,
			int i_order,	//!< order of RK time stepping method
			int i_order2,	//!< order of RK time stepping method of non-linear parts
			double i_timestepSize,
			bool i_use_f_sphere2D,
			int i_version_id,
			bool i_use_rexi_sphere2d_solver_preallocation
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override
	{
		timestepping_method = i_timestepping_method;

		if (
			i_timestepping_method == "l_exp_n_erk" || i_timestepping_method == "l_exp_n_erk_ver0" ||
			i_timestepping_method == "l_exp_n_erk_ver1"
		)
			return true;

		return false;
	}

public:
	std::string getIDString() override
	{
		std::string s = "l_exp_n_erk_ver";

		if (version_id == 0)
			s += "0";
		else if (version_id == 1)
			s += "1";
		else
			SWEETErrorFatal("Version ID");

		return s;
	}


	double timestep_size;

	/*
	 * Linear time steppers
	 */
	PDESWESphere2DTS_l_exp timestepping_l_rexi;

	/*
	 * Non-linear time steppers
	 */
	PDESWESphere2DTS_l_erk_n_erk timestepping_l_erk_n_erk;

	sweet::DEPRECATED_TimesteppingExplicitRKSphere2DData timestepping_rk_nonlinear;

	int version_id;


public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override
	{
		PDESWESphere2DTS_BaseInterface::shackRegistration(io_shackDict);

		timestepping_l_rexi.shackRegistration(io_shackDict);
		timestepping_l_erk_n_erk.shackRegistration(io_shackDict);
		return true;
	}


public:
	PDESWESphere2DTS_l_exp_n_erk();

	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vort,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;



	virtual ~PDESWESphere2DTS_l_exp_n_erk();
};

#endif
