/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2D_TIMESTEPPERS_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEOLD_PDESWESPHERE2D_TIMESTEPPERS_HPP


#include <sweet/Error/Base.hpp>

#include "PDESWESphere2DTS_BaseInterface.hpp"

/**
 * SWE Cart2D time steppers
 */
class PDESWESphere2D_TimeSteppers
{
public:
	sweet::Error::Base error;
	PDESWESphere2DTS_BaseInterface *timestepper = nullptr;

private:
	std::vector<PDESWESphere2DTS_BaseInterface*> _registered_integrators;

public:
	void setup_1_registerAllTimesteppers();

private:
	void _timesteppersFreeAll(
			PDESWESphere2DTS_BaseInterface *skip_this = nullptr
		);

public:
	PDESWESphere2D_TimeSteppers();

	void printImplementedTimesteppingMethods(
		std::ostream &o_ostream = std::cout,
		const std::string &i_prefix = ""
	);

	bool setup_2_shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	);

	bool setup_3_timestepper(
			const std::string &i_timestepping_method,
			sweet::Shacks::Dictionary *i_shackDict,
			sweet::Data::Sphere2D::Operators *io_ops
	);

	void clear();


	~PDESWESphere2D_TimeSteppers();
};




#endif
