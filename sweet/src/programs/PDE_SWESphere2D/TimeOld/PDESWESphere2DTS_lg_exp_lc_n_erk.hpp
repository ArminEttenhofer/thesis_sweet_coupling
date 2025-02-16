/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_EXP_LC_N_ERK_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_EXP_LC_N_ERK_HPP

#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKSphere2DData.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_l_exp.hpp"
#include "PDESWESphere2DTS_lg_erk_lc_n_erk.hpp"



class PDESWESphere2DTS_lg_exp_lc_n_erk	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			sweet::ExpIntegration::Shack *i_shackExpIntegration,
			int i_timestepping_order,
			int i_timestepping_order2,
			double i_timestepSize,
			int i_version_id
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackPDESWETimeDisc->timestepping_order;
		timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
		if (
			i_timestepping_method == "lg_exp_lc_n_erk" || i_timestepping_method == "lg_exp_lc_n_erk_ver0" ||
			i_timestepping_method == "lg_exp_lc_n_erk_ver1"
		)
			return true;

		return false;
	}

	std::string getIDString() override
	{
		std::string s = "lg_exp_lc_n_erk_ver";

		if (version_id == 0)
			s += "0";
		else if (version_id == 1)
			s += "1";
		else
			SWEETErrorFatal("Version ID");

		return s;
	}

private:
	int version_id;

	double timestep_size;

	/*
	 * Linear time steppers
	 */
	PDESWESphere2DTS_l_exp timestepping_lg_rexi;

	/*
	 * Non-linear time steppers
	 */
	PDESWESphere2DTS_lg_erk_lc_n_erk timestepping_lg_erk_lc_n_erk;

	sweet::DEPRECATED_TimesteppingExplicitRKSphere2DData timestepping_rk_nonlinear;

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override
	{
		PDESWESphere2DTS_BaseInterface::shackRegistration(io_shackDict);

		timestepping_lg_rexi.shackRegistration(io_shackDict);
		timestepping_lg_erk_lc_n_erk.shackRegistration(io_shackDict);
		return true;
	}

public:
	PDESWESphere2DTS_lg_exp_lc_n_erk();

	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;


	virtual ~PDESWESphere2DTS_lg_exp_lc_n_erk();
};

#endif
