#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_EXP_LC_TAYLOR_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LG_EXP_LC_TAYLOR_HPP

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_lg_exp_direct.hpp"
#include "PDESWESphere2DTS_ln_erk_split_vd.hpp"



class PDESWESphere2DTS_lg_exp_lc_taylor	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			int i_order	//!< order of RK time stepping method
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackPDESWETimeDisc->timestepping_order;
		//timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
		return i_timestepping_method == "lg_exp_lc_exp";
	}

public:
	std::string getIDString() override
	{
		return "lg_exp_lc_exp";
	}

	PDESWESphere2DTS_lg_exp_direct timestepping_lg_exp;
	PDESWESphere2DTS_ln_erk_split_vd timestepping_lc_erk;


public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override
	{
		PDESWESphere2DTS_BaseInterface::shackRegistration(io_shackDict);

		timestepping_lg_exp.shackRegistration(io_shackDict);
		timestepping_lc_erk.shackRegistration(io_shackDict);
		return true;
	}


public:
	PDESWESphere2DTS_lg_exp_lc_taylor();

	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;

	void run_timestep_lc(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

	void run_timestep_lg(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

	virtual ~PDESWESphere2DTS_lg_exp_lc_taylor();
};

#endif
