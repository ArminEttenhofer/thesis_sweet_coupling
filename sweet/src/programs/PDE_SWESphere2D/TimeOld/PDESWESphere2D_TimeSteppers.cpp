/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere2D_TimeSteppers.hpp"

#include "PDESWESphere2DTS_l_erk.hpp"
#include "PDESWESphere2DTS_l_erk_n_erk.hpp"
#include "PDESWESphere2DTS_l_erk_na_erk_uv.hpp"
#include "PDESWESphere2DTS_l_erk_na_erk_vd.hpp"
#include "PDESWESphere2DTS_l_exp.hpp"
#include "PDESWESphere2DTS_l_exp_direct_special.hpp"
#include "PDESWESphere2DTS_l_exp_n_erk.hpp"
#include "PDESWESphere2DTS_l_exp_n_etdrk.hpp"
#include "PDESWESphere2DTS_l_irk.hpp"
#include "PDESWESphere2DTS_l_irk_n_erk.hpp"
#include "PDESWESphere2DTS_l_irk_na_erk_uv.hpp"
#include "PDESWESphere2DTS_l_irk_na_erk_vd.hpp"
#include "PDESWESphere2DTS_l_irk_na_sl_nr_settls_uv_only.hpp"
#include "PDESWESphere2DTS_l_irk_na_sl_nr_settls_vd_only.hpp"
#include "PDESWESphere2DTS_l_irk_na_sl_settls_uv_only.hpp"
#include "PDESWESphere2DTS_l_irk_na_sl_settls_vd_only.hpp"
#include "PDESWESphere2DTS_lg_erk.hpp"
#include "PDESWESphere2DTS_lg_erk_lc_erk.hpp"
#include "PDESWESphere2DTS_lg_erk_lc_n_erk.hpp"
#include "PDESWESphere2DTS_lg_exp_lc_erk.hpp"
#include "PDESWESphere2DTS_lg_exp_lc_n_erk.hpp"
#include "PDESWESphere2DTS_lg_exp_lc_n_etd_uv.hpp"
#include "PDESWESphere2DTS_lg_exp_lc_n_etd_vd.hpp"
#include "PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etd_uv.hpp"
#include "PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv.hpp"
#include "PDESWESphere2DTS_lg_exp_lc_n_etdrk.hpp"
#include "PDESWESphere2DTS_lg_irk.hpp"
#include "PDESWESphere2DTS_lg_irk_lc_erk.hpp"
#include "PDESWESphere2DTS_lg_irk_lc_n_erk_ver01.hpp"
#include "PDESWESphere2DTS_lg_irk_lc_na_erk_vd.hpp"
#include "PDESWESphere2DTS_ln_erk.hpp"
#include "PDESWESphere2DTS_ln_erk_split_uv.hpp"
#include "PDESWESphere2DTS_ln_erk_split_vd.hpp"
#include "PDESWESphere2DTS_ln_settls_uv.hpp"
#include "PDESWESphere2DTS_ln_settls_vd.hpp"
#include "PDESWESphere2DTS_ln_sl_exp_settls_uv.hpp"
#include "PDESWESphere2DTS_ln_sl_exp_settls_vd.hpp"
#include "PDESWESphere2DTS_lg_0_lc_n_erk_bv.hpp"
#include "PDESWESphere2DTS_lg_exp_direct.hpp"
#include "PDESWESphere2DTS_lg_exp_lc_taylor.hpp"





PDESWESphere2D_TimeSteppers::PDESWESphere2D_TimeSteppers()
{
}




void PDESWESphere2D_TimeSteppers::setup_1_registerAllTimesteppers()
{
	/*
	 * Register time integrators
	 */
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_erk));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_erk_lc_erk));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_exp_lc_taylor));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_irk_lc_erk));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_irk_n_erk));

	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_erk_n_erk));

	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_erk_na_erk_vd));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_erk_na_erk_uv));

	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_irk_na_erk_vd));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_irk_na_erk_uv));

	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_exp_n_erk));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_exp_direct));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_irk_lc_na_erk_vd));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_irk_lc_n_erk));

	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_exp_lc_erk));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_exp_lc_n_erk));

	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_erk_lc_n_erk));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_erk));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_ln_erk));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_ln_erk_split_uv));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_ln_erk_split_vd));

	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_exp_n_etdrk));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_exp_lc_n_etdrk));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_exp_lc_n_etd_uv));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_exp_lc_n_etd_vd));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etd_uv));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_exp_na_sl_lc_nr_etdrk_uv));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_irk));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_irk));

	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_exp));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_exp_direct_special));

	/*
	 * EXP SETTLS VERSION
	 */
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_ln_sl_exp_settls_vd));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_ln_sl_exp_settls_uv));

	/*
	 * ONLY SETTLS VERSION without any special variants
	 */
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_irk_na_sl_nr_settls_vd_only));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_irk_na_sl_nr_settls_uv_only));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_irk_na_sl_settls_vd_only));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_l_irk_na_sl_settls_uv_only));

	/*
	 * IRK SETTLS VERSION
	 */
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_ln_settls_vd));
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_ln_settls_uv));

	/*
	 * BAROTROPIC VORTICITY EQ
	 */
	_registered_integrators.push_back(static_cast<PDESWESphere2DTS_BaseInterface*>(new PDESWESphere2DTS_lg_0_lc_n_erk_bv));
}



bool PDESWESphere2D_TimeSteppers::setup_2_shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
)
{
	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		_registered_integrators[i]->shackRegistration(io_shackDict);
	}
	return true;
}



void PDESWESphere2D_TimeSteppers::printImplementedTimesteppingMethods(
	std::ostream &o_ostream,
	const std::string &i_prefix
)
{
	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Timestepping methods (START)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;

	std::string prefix = i_prefix+"  ";
	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		_registered_integrators[i]->printImplementedTimesteppingMethods(o_ostream, prefix);
	}

	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Timestepping methods (END)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << std::endl;
}


bool PDESWESphere2D_TimeSteppers::setup_3_timestepper(
		const std::string &i_timestepping_method,
		sweet::Shacks::Dictionary *io_shackDict,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	if (i_timestepping_method == "")
	{
		printImplementedTimesteppingMethods();
		return error.set("Please set time stepping method using --timestepping-method=...");
	}

	/*
	 * Find right one
	 */
	timestepper = nullptr;

	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		PDESWESphere2DTS_BaseInterface *ts = _registered_integrators[i];

		if (ts->implementsTimesteppingMethod(i_timestepping_method))
		{
			if (timestepper != nullptr)
			{
				return error.set("Duplicate implementation for method "+i_timestepping_method);
			}

			ts->setup_auto(i_timestepping_method, io_ops);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*ts);
			timestepper = ts;
		}
	}

	if (timestepper == nullptr)
		return error.set("No valid --timestepping-method '"+i_timestepping_method+"' provided");

	// Found integrator, freeing others
	_timesteppersFreeAll(timestepper);

	return true;
}


void PDESWESphere2D_TimeSteppers::_timesteppersFreeAll(
		PDESWESphere2DTS_BaseInterface *i_skip_this_timestepper
)
{

	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		PDESWESphere2DTS_BaseInterface *ts = _registered_integrators[i];

		if (ts == i_skip_this_timestepper)
			continue;

		delete ts;
	}

	_registered_integrators.clear();
}


void PDESWESphere2D_TimeSteppers::clear()
{
	delete timestepper;
	timestepper = nullptr;

	_timesteppersFreeAll();
}


PDESWESphere2D_TimeSteppers::~PDESWESphere2D_TimeSteppers()
{
	clear();
}
