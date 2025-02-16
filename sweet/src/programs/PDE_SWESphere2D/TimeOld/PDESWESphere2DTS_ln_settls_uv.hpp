/*
 * PDESWESphere2DTS_ln_settls_uv.hpp
 *
 *  Created on: 24 Sep 2019
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Based on cart2d code
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LN_SETTLS_UV_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2DTS_LN_SETTLS_UV_HPP

#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKSphere2DData.hpp>
#include <sweet/Data/Sphere2D/DataGrid.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Sphere2D/Operators_Sampler_Sphere2DDataGrid.hpp>
#include <sweet/SemiLagrangian/Shack.hpp>
#include <sweet/SemiLagrangian/Sphere2D.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <limits>
#include <string>
#include "PDESWESphere2DTS_BaseInterface.hpp"
#include "PDESWESphere2DTS_l_irk.hpp"
#include "PDESWESphere2DTS_lg_irk.hpp"
#include "PDESWESphere2DTS_ln_erk_split_uv.hpp"



class PDESWESphere2DTS_ln_settls_uv	: public PDESWESphere2DTS_BaseInterface
{
public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;

	std::string string_id_storage;

	std::string getIDString() override;

private:
	sweet::SemiLagrangian::Sphere2D semiLagrangian;
	//sweet::Data::Sphere2D::Sphere2DOperators_Sampler_Sphere2DDataGrid sphere2DSampler;

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

private:

	LinearCoriolisTreatment_enum coriolis_treatment;
	NLRemainderTreatment_enum nonlinear_remainder_treatment;

	int timestepping_order;
	bool original_linear_operator_sl_treatment;

	sweet::Data::Sphere2D::DataSpectral coriolis_arrival_spectral;
	sweet::Data::Sphere2D::DataSpectral U_phi_prev, U_vrt_prev, U_div_prev;

	PDESWESphere2DTS_ln_erk_split_uv swe_sphere2d_ts_ln_erk_split_uv;
	PDESWESphere2DTS_l_irk swe_sphere2d_ts_l_irk;
	PDESWESphere2DTS_lg_irk swe_sphere2d_ts_lg_irk;


public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override
	{
		PDESWESphere2DTS_BaseInterface::shackRegistration(io_shackDict);

		swe_sphere2d_ts_ln_erk_split_uv.shackRegistration(io_shackDict);
		swe_sphere2d_ts_l_irk.shackRegistration(io_shackDict);
		swe_sphere2d_ts_lg_irk.shackRegistration(io_shackDict);
		return true;
	}

public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::Data::Sphere2D::Operators *io_ops
		) override;

	bool setup_main(
			const sweet::Data::Sphere2D::Operators *io_ops,
			int i_timestepping_order,
			LinearCoriolisTreatment_enum i_coriolis_treatment,// = PDESWESphere2DTS_ln_settls::CORIOLIS_LINEAR,		// "ignore", "linear", "nonlinear", "semi-lagrangian"
			NLRemainderTreatment_enum i_nonlinear_divergence_treatment,// = PDESWESphere2DTS_ln_settls::NL_DIV_NONLINEAR,	// "ignore", "nonlinear"
			bool original_linear_operator_sl_treatment	// = true
	);



public:
	PDESWESphere2DTS_ln_settls_uv();

	void runTimestep(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vort,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_dt = 0,
			double i_simulation_timestamp = -1
	) override;

	void run_timestep_1st_order(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vort,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);


	void run_timestep_2nd_order(
			sweet::Data::Sphere2D::DataSpectral &io_phi,
			sweet::Data::Sphere2D::DataSpectral &io_vort,
			sweet::Data::Sphere2D::DataSpectral &io_div,

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

	virtual ~PDESWESphere2DTS_ln_settls_uv();
};

#endif
