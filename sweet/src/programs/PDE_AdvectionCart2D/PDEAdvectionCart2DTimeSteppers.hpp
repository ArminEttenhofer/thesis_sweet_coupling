/*
 * Adv_Cart2D_TimeSteppers.hpp
 *
 *  Created on: 4th April 2018
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONCART2D_PDEADVECTIONCART2DTIMESTEPPERS_HPP
#define PROGRAMS_PDE_ADVECTIONCART2D_PDEADVECTIONCART2DTIMESTEPPERS_HPP

#include <programs/PDE_AdvectionCart2D/time/PDEAdvectionCart2DTS_BaseInterface.hpp>
#include <programs/PDE_AdvectionCart2D/time/PDEAdvectionCart2DTS_na_erk.hpp>
#include <programs/PDE_AdvectionCart2D/time/PDEAdvectionCart2DTS_na_sl.hpp>
#include <programs/PDE_AdvectionCart2D/time/ShackPDEAdvectionCart2DTimeDiscretization.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>


/**
 * Advection for the cart2d time integrators
 */
class PDEAdvectionCart2DTimeSteppers
{
public:
	sweet::Error::Base error;

	PDEAdvectionCart2DTS_na_erk *na_erk;
	PDEAdvectionCart2DTS_na_sl *na_sl;
	PDEAdvectionCart2DTS_BaseInterface *master;

	ShackPDEAdvectionCart2DTimeDiscretization *shackTimeDisc;

	PDEAdvectionCart2DTimeSteppers()	:
		na_erk(nullptr),
		na_sl(nullptr),
		master(nullptr),
		shackTimeDisc(nullptr)
	{
	}

	void clear()
	{
		if (na_erk != nullptr)
		{
			delete na_erk;
			na_erk = nullptr;
		}

		if (na_sl != nullptr)
		{
			delete na_sl;
			na_sl = nullptr;
		}

		master = nullptr;

		shackTimeDisc = nullptr;
	}


	bool shackRegistration(
			sweet::Shacks::Dictionary &io_shackDict
	)
	{
		shackTimeDisc = io_shackDict.getAutoRegistration<ShackPDEAdvectionCart2DTimeDiscretization>();

		return error.forwardWithPositiveReturn(io_shackDict.error);
	}


	bool setup(
			sweet::Shacks::Dictionary &io_shackDict,
			sweet::Data::Cart2D::Operators &i_op
	)
	{
		SWEET_ASSERT(shackTimeDisc != nullptr);

		std::cout << "Setting up time stepping method '" << shackTimeDisc->timestepping_method << "'" << std::endl;

		if (shackTimeDisc->timestepping_method == "na_erk")
		{
			na_erk = new PDEAdvectionCart2DTS_na_erk;
			na_erk->shackRegistration(&io_shackDict);
			na_erk->setup(&i_op);

			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*na_erk);

			master = static_cast<PDEAdvectionCart2DTS_BaseInterface*>(na_erk);
			return true;
		}
		else if (shackTimeDisc->timestepping_method == "na_sl")
		{
			na_sl = new PDEAdvectionCart2DTS_na_sl;
			na_sl->shackRegistration(&io_shackDict);
			na_sl->setup(&i_op);

			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*na_sl);

			master = static_cast<PDEAdvectionCart2DTS_BaseInterface*>(na_sl);
			return true;
		}

		return error.set("No valid --timestepping-method=... provided ('"+shackTimeDisc->timestepping_method+"')");
	}


	~PDEAdvectionCart2DTimeSteppers()
	{
		clear();
	}
};




#endif
