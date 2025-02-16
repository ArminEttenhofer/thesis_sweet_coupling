/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LN_SL_EXP_SETTLS_UV_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LN_SL_EXP_SETTLS_UV_HPP

#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKSphere2DData.hpp>
#include <sweet/Data/Sphere2D/DataGrid.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Sphere2D/Operators_Sampler_Sphere2DDataGrid.hpp>
#include <sweet/SemiLagrangian/Shack.hpp>
#include <sweet/SemiLagrangian/Sphere2D.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_l_exp.hpp"
#include "PDESWESphere2DTS_ln_erk_split_uv.hpp"



class PDESWESphere2DTS_ln_sl_exp_settls_uv	: public PDESWESphere2DTS_BaseInterface
{

public:
	enum LinearCoriolisTreatment_enum {
		CORIOLIS_IGNORE,
		CORIOLIS_LINEAR,
		CORIOLIS_NONLINEAR,
		CORIOLIS_SEMILAGRANGIAN,
	};

	enum NLRemainderTreatment_enum{
		NL_REMAINDER_IGNORE,
		NL_REMAINDER_NONLINEAR,
	};

public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			int i_timestepping_order,
			LinearCoriolisTreatment_enum i_coriolis_treatment,
			NLRemainderTreatment_enum i_nonlinear_divergence_treatment,
			bool original_linear_operator_sl_treatment
	);


public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;

	std::string string_id_storage;

private:
	LinearCoriolisTreatment_enum coriolis_treatment;
	NLRemainderTreatment_enum nonlinear_remainder_treatment;

	bool original_linear_operator_sl_treatment;

	sweet::SemiLagrangian::Sphere2D semiLagrangian;

	sweet::Data::Sphere2D::DataSpectral U_phi_prev, U_vrt_prev, U_div_prev;

	PDESWESphere2DTS_ln_erk_split_uv swe_sphere2d_ts_ln_erk_split_uv;
	PDESWESphere2DTS_l_exp swe_sphere2d_ts_l_rexi;


public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override
	{
		PDESWESphere2DTS_BaseInterface::shackRegistration(io_shackDict);

		swe_sphere2d_ts_ln_erk_split_uv.shackRegistration(io_shackDict);
		swe_sphere2d_ts_l_rexi.shackRegistration(io_shackDict);
		return true;
	}



public:
	PDESWESphere2DTS_ln_sl_exp_settls_uv();

	void printHelp() override;

	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_dt = 0,
			double i_simulation_timestamp = -1
	) override;


	void run_timestep_2nd_order(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vrt,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

#if (SWEET_PARAREAL && SWEET_PARAREAL_SPHERE2D) || (SWEET_XBRAID && SWEET_XBRAID_SPHERE2D)
	void set_previous_solution(
				sweet::Data::Sphere2D::DataSpectral &i_phi_prev,
				sweet::Data::Sphere2D::DataSpectral &i_vrt_prev,
				sweet::Data::Sphere2D::DataSpectral &i_div_prev
	) override
	{
		U_phi_prev = i_phi_prev;
		U_vrt_prev = i_vrt_prev;
		U_div_prev = i_div_prev;
	}
#endif

	virtual ~PDESWESphere2DTS_ln_sl_exp_settls_uv();
};

#endif
