/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_TIME_PDEADVECTIONSPHERE2DTS_NA_TRAJECTORIES_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_TIME_PDEADVECTIONSPHERE2DTS_NA_TRAJECTORIES_HPP


#include <programs/PDE_AdvectionSphere2D/time/PDEAdvectionSphere2DTS_BaseInterface.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Sphere2D/Operators_Sampler_Sphere2DDataGrid.hpp>
#include <sweet/SemiLagrangian/Sphere2D.hpp>
#include <sweet/Shacks/Dictionary.hpp>

#include "../PDEAdvectionSphere2DBenchmarksCombined.hpp"
#include <limits>


/*
 * This time integration method is based on given benchmark solutions
 *
 * See e.g. R. Nair, P. Lauritzen
 * "A class of deformational flow test cases for linear transport problems on the sphere2D"
 */
class PDEAdvectionSphere2DTS_na_trajectories	:
		public PDEAdvectionSphere2DTS_BaseInterface
{
	int timestepping_order;

	sweet::SemiLagrangian::Sphere2D semiLagrangian;
	sweet::Data::Sphere2D::Operators_Sampler_DataGrid sphere2DSampler;

	sweet::Data::Sphere2D::DataSpectral U_phi_prev;
	sweet::Data::Sphere2D::DataGrid U_u_prev, U_v_prev;

public:
	bool testImplementsTimesteppingMethod(const std::string &i_timestepping_method);

	std::string getStringId();

	void printImplementedTimesteppingMethods(
			std::ostream &o_ostream,
			const std::string &i_prefix
	);

public:
	bool setup(
			sweet::Data::Sphere2D::Operators *io_ops
	);

private:
	void runTimestep(
			std::vector<sweet::Data::Sphere2D::DataSpectral> &io_prognostic_fields,	//!< prognostic variables
			sweet::Data::Sphere2D::DataGrid &io_u,
			sweet::Data::Sphere2D::DataGrid &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp
	);

	void run_timestep_1(
			sweet::Data::Sphere2D::DataSpectral &io_prognostic_field,	//!< prognostic variables
			sweet::Data::Sphere2D::DataGrid &io_u,
			sweet::Data::Sphere2D::DataGrid &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp
	);



	virtual ~PDEAdvectionSphere2DTS_na_trajectories();
};

#endif
