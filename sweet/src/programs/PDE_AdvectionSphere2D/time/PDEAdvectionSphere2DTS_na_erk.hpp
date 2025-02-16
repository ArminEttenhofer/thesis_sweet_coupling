/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_TIME_PDEADVECTIONSPHERE2DTS_NA_ERK_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_TIME_PDEADVECTIONSPHERE2DTS_NA_ERK_HPP

#include <programs/PDE_AdvectionSphere2D/PDEAdvectionSphere2DBenchmarksCombined.hpp>
#include <programs/PDE_AdvectionSphere2D/time/PDEAdvectionSphere2DTS_BaseInterface.hpp>
#include <sweet/_DEPRECATED/DEPRECATED_TimesteppingExplicitRKSphere2DData.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <limits>



class PDEAdvectionSphere2DTS_na_erk	:
		public PDEAdvectionSphere2DTS_BaseInterface
{
	int timestepping_order;

	sweet::DEPRECATED_TimesteppingExplicitRKSphere2DData timestepping_rk;

public:
	bool testImplementsTimesteppingMethod(
			const std::string &i_timestepping_method
	);

	std::string getStringId();

	void printImplementedTimesteppingMethods(
			std::ostream &o_ostream,
			const std::string &i_prefix
	);


private:
	void euler_timestep_update(
			const sweet::Data::Sphere2D::DataSpectral &i_prognostic_field,	//!< prognostic variables
			sweet::Data::Sphere2D::DataGrid &io_u,
			sweet::Data::Sphere2D::DataGrid &io_v,

			sweet::Data::Sphere2D::DataSpectral &o_prognostic_field,	//!< time updates

			double i_simulation_timestamp
	);

public:
	bool setup(
			sweet::Data::Sphere2D::Operators *io_ops
	);

	void runTimestep(
			std::vector<sweet::Data::Sphere2D::DataSpectral> &io_prognostic_fields,	//!< prognostic variables
			sweet::Data::Sphere2D::DataGrid &io_u,
			sweet::Data::Sphere2D::DataGrid &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp
	);



	virtual ~PDEAdvectionSphere2DTS_na_erk();
};

#endif
