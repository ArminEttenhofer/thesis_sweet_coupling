/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_PDEADVECTIONSPHERE2DTIMESTEPPERS_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_PDEADVECTIONSPHERE2DTIMESTEPPERS_HPP


#include <programs/PDE_AdvectionSphere2D/time/PDEAdvectionSphere2DTS_BaseInterface.hpp>
#include <sweet/Error/Base.hpp>



/**
 * SWE Cart2D time steppers
 */
class PDEAdvectionSphere2DTimeSteppers
{
public:
	sweet::Error::Base error;
	PDEAdvectionSphere2DTS_BaseInterface *timestepper = nullptr;

private:
	std::vector<PDEAdvectionSphere2DTS_BaseInterface*> _registered_integrators;


public:
	void setup_1_registerAllTimesteppers();

private:
	void _timesteppersFreeAll(
			PDEAdvectionSphere2DTS_BaseInterface *skip_this = nullptr
		);

public:
	PDEAdvectionSphere2DTimeSteppers();

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


	~PDEAdvectionSphere2DTimeSteppers();
};


#endif
