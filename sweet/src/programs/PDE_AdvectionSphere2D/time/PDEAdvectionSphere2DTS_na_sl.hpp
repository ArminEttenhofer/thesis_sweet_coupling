/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_ADVECTIONSPHERE2D_TIME_PDEADVECTIONSPHERE2DTS_NA_SL_HPP
#define PROGRAMS_PDE_ADVECTIONSPHERE2D_TIME_PDEADVECTIONSPHERE2DTS_NA_SL_HPP


#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Sphere2D/Operators_Sampler_Sphere2DDataGrid.hpp>
#include <sweet/SemiLagrangian/Sphere2D.hpp>
#include <sweet/Shacks/Dictionary.hpp>

#include "../PDEAdvectionSphere2DBenchmarksCombined.hpp"
#include <limits>
#include "../time/PDEAdvectionSphere2DTS_BaseInterface.hpp"


class PDEAdvectionSphere2DTS_na_sl	:
		public PDEAdvectionSphere2DTS_BaseInterface
{
	int timestepping_order;

	sweet::SemiLagrangian::Sphere2D semiLagrangian;
	sweet::Data::Sphere2D::Operators_Sampler_DataGrid sphere2DSampler;

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
	void interpolate_departure_point_vec_3d(
			const sweet::Data::Sphere2D::DataSpectral &i_u,
			const sweet::Data::Sphere2D::DataSpectral &i_v,
			const sweet::Data::Sphere2D::DataSpectral &i_w,

			const sweet::Data::Vector::Vector<double> &i_pos_lon_D,
			const sweet::Data::Vector::Vector<double> &i_pos_lat_D,

			sweet::Data::Sphere2D::DataSpectral &o_u,
			sweet::Data::Sphere2D::DataSpectral &o_v,
			sweet::Data::Sphere2D::DataSpectral &o_w
	);

	void interpolate_departure_point_vec_uv(
			const sweet::Data::Sphere2D::DataGrid &i_u,
			const sweet::Data::Sphere2D::DataGrid &i_v,

			const sweet::Data::Vector::Vector<double> &i_pos_lon_D,
			const sweet::Data::Vector::Vector<double> &i_pos_lat_D,

			sweet::Data::Sphere2D::DataGrid &o_u,
			sweet::Data::Sphere2D::DataGrid &o_v
	);

	void runTimestep(
			std::vector<sweet::Data::Sphere2D::DataSpectral> &io_prog_fields,	//!< prognostic variables
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


	void run_timestep_2(
			std::vector<sweet::Data::Sphere2D::DataSpectral> &io_prog_fields,	//!< prognostic variables
			sweet::Data::Sphere2D::DataGrid &io_u,
			sweet::Data::Sphere2D::DataGrid &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp
	);



	void run_timestep_3(
			std::vector<sweet::Data::Sphere2D::DataSpectral> &io_prog_fields,	//!< prognostic variables
			sweet::Data::Sphere2D::DataGrid &io_u,
			sweet::Data::Sphere2D::DataGrid &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp
	);



	virtual ~PDEAdvectionSphere2DTS_na_sl();
};

#endif
