#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_EXP_DIRECT_SPECIAL_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_L_EXP_DIRECT_SPECIAL_HPP

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_lg_exp_direct.hpp"
#include "PDESWESphere2DTS_ln_erk_split_vd.hpp"



class PDESWESphere2DTS_l_exp_direct_special	: public PDESWESphere2DTS_BaseInterface
{
	/*
	 * This class acts as a wrapper around the lg_exp_direct
	 * method and an ETDnRK version (lg + lc) to automatically
	 * choose the right one.
	 *
	 * Note, that the lc ETD method is actually not a direct method.
	 */

public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			int i_order,	//!< order of RK time stepping method
			bool i_use_coriolis,		//!< Include Coriolis term
			const std::string &i_function_name	//!< phi/ups function
	);

public:
	bool implementsTimesteppingMethod(
		const std::string &i_timestepping_method
	) override;

	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override;

public:

	std::string getIDString() override
	{
		return "l_exp_special";
	}

	bool use_coriolis;

	PDESWESphere2DTS_lg_exp_direct timestepping_lg_exp_phi0;
	PDESWESphere2DTS_lg_exp_direct timestepping_lg_exp_phi1;
	PDESWESphere2DTS_lg_exp_direct timestepping_lg_exp_phi2;

	PDESWESphere2DTS_lg_exp_direct timestepping_lg_exp_ups1;
	PDESWESphere2DTS_lg_exp_direct timestepping_lg_exp_ups2;
	PDESWESphere2DTS_lg_exp_direct timestepping_lg_exp_ups3;

	PDESWESphere2DTS_ln_erk_split_vd timestepping_lc_erk;

	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;

	void euler_timestep_store_update_lc(
			const sweet::Data::Sphere2D::DataSpectral &i_phi_pert,
			const sweet::Data::Sphere2D::DataSpectral &i_vrt,
			const sweet::Data::Sphere2D::DataSpectral &i_div,
			sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
			sweet::Data::Sphere2D::DataSpectral &o_vrt,
			sweet::Data::Sphere2D::DataSpectral &o_div,
			double i_simulation_timestamp
	);

	PDESWESphere2DTS_l_exp_direct_special();

	virtual ~PDESWESphere2DTS_l_exp_direct_special();
};

#endif

