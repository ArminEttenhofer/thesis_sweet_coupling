/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */


#include "PDESWESphere2DTS_ln_erk_split_uv.hpp"




bool PDESWESphere2DTS_ln_erk_split_uv::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	/*
	 * l_na
	 */
	if (timestepping_method == "l_na_erk_split_uv")
		return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, true, true, true, false, false);

	if (timestepping_method == "l_na_erk_split_aa_uv")
		return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, true, true, true, false, true);

	/*
	 * l
	 */
	if (timestepping_method == "l_erk_split_uv")
		return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, true, true, false, false, false);

	if (timestepping_method == "l_erk_split_aa_uv")
		return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, true, true, false, false, true);

	/*
	 * ln
	 */
	if (timestepping_method == "ln_erk_split_uv")
		return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, true, true, true, true, false);

	if (timestepping_method == "ln_erk_split_aa_uv")
		return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order, true, true, true, true, true);

	SWEETErrorFatal("Should never happen");
	return false;
}


/*
 * Setup
 */
bool PDESWESphere2DTS_ln_erk_split_uv::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_order,	//!< order of RK time stepping method
		bool i_use_lg,
		bool i_use_lc,
		bool i_use_na,
		bool i_use_nr,

		bool i_antialiasing_for_each_term
)
{
	ops = io_ops;

	setupFG();

	timestepping_method = i_order;


	use_lg = i_use_lg;
	use_lc = i_use_lc;
	use_na = i_use_na;
	use_nr = i_use_nr;

	anti_aliasing_for_each_term = i_antialiasing_for_each_term;

	timestepping_order = i_order;
	return true;
}


/*
 * Main routine for method to be used in case of finite differences
 */
void PDESWESphere2DTS_ln_erk_split_uv::euler_timestep_update_lg(
		const sweet::Data::Sphere2D::DataSpectral &i_U_phi_pert,
		const sweet::Data::Sphere2D::DataSpectral &i_U_vrt,
		const sweet::Data::Sphere2D::DataSpectral &i_U_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_pert_t,
		sweet::Data::Sphere2D::DataSpectral &o_vrt_t,
		sweet::Data::Sphere2D::DataSpectral &o_div_t,

		double i_simulation_timestamp
)
{
	double gh0 = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0;

	const sweet::Data::Sphere2D::DataSpectral &U_phi_pert = i_U_phi_pert;
	const sweet::Data::Sphere2D::DataSpectral &U_div = i_U_div;


	sweet::Data::Sphere2D::DataGrid U_u_phys, U_v_phys;
	ops->vrtdiv_2_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

	o_phi_pert_t -= gh0*U_div;
	o_div_t -= ops->laplace(U_phi_pert);
}



/*
 * Main routine for method to be used in case of finite differences
 */
void PDESWESphere2DTS_ln_erk_split_uv::euler_timestep_update_lc(
		const sweet::Data::Sphere2D::DataSpectral &i_U_phi,
		const sweet::Data::Sphere2D::DataSpectral &i_U_vrt,
		const sweet::Data::Sphere2D::DataSpectral &i_U_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_pert_t,
		sweet::Data::Sphere2D::DataSpectral &o_vrt_t,
		sweet::Data::Sphere2D::DataSpectral &o_div_t,

		double i_simulation_timestamp
)
{
//	double gh0 = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0;


//	const sweet::Data::Sphere2D::DataSpectral &U_phi_pert = i_U_phi;
	//const sweet::Data::Sphere2D::DataSpectral &U_vrt = i_U_vrt;
//	const sweet::Data::Sphere2D::DataSpectral &U_div = i_U_div;


	sweet::Data::Sphere2D::DataGrid U_u_phys, U_v_phys;
	ops->vrtdiv_2_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);


	sweet::Data::Sphere2D::DataGrid fu_nl = fg*U_u_phys;
	sweet::Data::Sphere2D::DataGrid fv_nl = fg*U_v_phys;

	sweet::Data::Sphere2D::DataSpectral div, vrt;
	ops->uv_2_vrtdiv(fu_nl, fv_nl, vrt, div);

	o_vrt_t -= div;
	o_div_t += vrt;
}



/*
 * Main routine for method to be used in case of finite differences
 */
void PDESWESphere2DTS_ln_erk_split_uv::euler_timestep_update_na(
		const sweet::Data::Sphere2D::DataSpectral &i_U_phi,
		const sweet::Data::Sphere2D::DataSpectral &i_U_vrt,
		const sweet::Data::Sphere2D::DataSpectral &i_U_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_t,
		sweet::Data::Sphere2D::DataSpectral &o_vrt_t,
		sweet::Data::Sphere2D::DataSpectral &o_div_t,

		double i_simulation_timestamp
)
{
//	double gh0 = shackPDESWESphere2D->gravitation * shackPDESWESphere2D->h0;


	const sweet::Data::Sphere2D::DataSpectral &U_phi_pert = i_U_phi;
	//const sweet::Data::Sphere2D::DataSpectral &U_vrt = i_U_vrt;
	const sweet::Data::Sphere2D::DataSpectral &U_div = i_U_div;


	sweet::Data::Sphere2D::DataGrid U_u_phys, U_v_phys;
	ops->vrtdiv_2_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);



	sweet::Data::Sphere2D::DataGrid U_div_phys = U_div.toGrid();
	o_phi_t -= ops->V_dot_grad_scalar(U_u_phys, U_v_phys, U_div_phys, U_phi_pert.toGrid());

	/*
	 * Velocity
	 */
	sweet::Data::Sphere2D::DataGrid vrtg = i_U_vrt.toGrid();

	sweet::Data::Sphere2D::DataGrid u_nl = U_u_phys*vrtg;
	sweet::Data::Sphere2D::DataGrid v_nl = U_v_phys*vrtg;

	sweet::Data::Sphere2D::DataSpectral vrt, div;
	ops->uv_2_vrtdiv(u_nl, v_nl, vrt, div);
	o_vrt_t -= div;
	o_div_t += vrt;

	o_div_t -= ops->laplace(0.5*(U_u_phys*U_u_phys+U_v_phys*U_v_phys));
}



/*
 * Main routine for method to be used in case of finite differences
 */
void PDESWESphere2DTS_ln_erk_split_uv::euler_timestep_update_nr(
		const sweet::Data::Sphere2D::DataSpectral &i_U_phi,
		const sweet::Data::Sphere2D::DataSpectral &i_U_vrt,
		const sweet::Data::Sphere2D::DataSpectral &i_U_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_t,
		sweet::Data::Sphere2D::DataSpectral &o_vrt_t,
		sweet::Data::Sphere2D::DataSpectral &o_div_t,

		double i_simulation_timestamp
)
{
	o_phi_t -= sweet::Data::Sphere2D::DataSpectral(i_U_phi.toGrid()*i_U_div.toGrid());
}



/*
 * Main routine for method to be used in case of finite differences
 */
void PDESWESphere2DTS_ln_erk_split_uv::euler_timestep_set_tendencies(
		const sweet::Data::Sphere2D::DataSpectral &i_U_phi,
		const sweet::Data::Sphere2D::DataSpectral &i_U_vrt,
		const sweet::Data::Sphere2D::DataSpectral &i_U_div,

		sweet::Data::Sphere2D::DataSpectral &o_phi_t,
		sweet::Data::Sphere2D::DataSpectral &o_vrt_t,
		sweet::Data::Sphere2D::DataSpectral &o_div_t,

		double i_simulation_timestamp
)
{
	sweet::Data::Sphere2D::DataGrid U_u_phys, U_v_phys;
	ops->vrtdiv_2_uv(i_U_vrt, i_U_div, U_u_phys, U_v_phys);

	o_phi_t.spectral_setZero();
	o_vrt_t.spectral_setZero();
	o_div_t.spectral_setZero();


	if (anti_aliasing_for_each_term)
	{
		sweet::Data::Sphere2D::DataSpectral phi_tmp(i_U_phi.sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral vrt_tmp(i_U_vrt.sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral div_tmp(i_U_div.sphere2DDataConfig);


		/*
		 * See [SWEET]/doc/swe/swe_sphere2d_formulation/swe_on_sphere2d_formulation_in_sweet.pdf/lyx
		 */
		if (use_lg)
		{
			phi_tmp.spectral_setZero();
			vrt_tmp.spectral_setZero();
			div_tmp.spectral_setZero();

			euler_timestep_update_lg(
					i_U_phi, i_U_vrt, i_U_div,
					phi_tmp, vrt_tmp, div_tmp,
					i_simulation_timestamp);

			o_phi_t += phi_tmp.toGrid();
			o_vrt_t += vrt_tmp.toGrid();
			o_div_t += div_tmp.toGrid();
		}


		if (use_lc)
		{
			phi_tmp.spectral_setZero();
			vrt_tmp.spectral_setZero();
			div_tmp.spectral_setZero();

			euler_timestep_update_lc(
					i_U_phi, i_U_vrt, i_U_div,
					phi_tmp, vrt_tmp, div_tmp,
					i_simulation_timestamp);

			o_phi_t += phi_tmp.toGrid();
			o_vrt_t += vrt_tmp.toGrid();
			o_div_t += div_tmp.toGrid();
		}


		if (use_na)
		{
			phi_tmp.spectral_setZero();
			vrt_tmp.spectral_setZero();
			div_tmp.spectral_setZero();

			euler_timestep_update_na(
					i_U_phi, i_U_vrt, i_U_div,
					phi_tmp, vrt_tmp, div_tmp,
					i_simulation_timestamp);

			o_phi_t += phi_tmp.toGrid();
			o_vrt_t += vrt_tmp.toGrid();
			o_div_t += div_tmp.toGrid();
		}


		if (use_nr)
		{
			phi_tmp.spectral_setZero();
			vrt_tmp.spectral_setZero();
			div_tmp.spectral_setZero();

			euler_timestep_update_nr(
					i_U_phi, i_U_vrt, i_U_div,
					phi_tmp, vrt_tmp, div_tmp,
					i_simulation_timestamp);

			o_phi_t += phi_tmp.toGrid();
			o_vrt_t += vrt_tmp.toGrid();
			o_div_t += div_tmp.toGrid();
		}
	}
	else
	{

		/*
		 * See [SWEET]/doc/swe/swe_sphere2d_formulation/swe_on_sphere2d_formulation_in_sweet.pdf/lyx
		 */
		if (use_lg)
		{
			euler_timestep_update_lg(
					i_U_phi, i_U_vrt, i_U_div,
					o_phi_t, o_vrt_t, o_div_t,
					i_simulation_timestamp);
		}


		if (use_lc)
		{
			euler_timestep_update_lc(
					i_U_phi, i_U_vrt, i_U_div,
					o_phi_t, o_vrt_t, o_div_t,
					i_simulation_timestamp);
		}


		if (use_na)
		{
			euler_timestep_update_na(
					i_U_phi, i_U_vrt, i_U_div,
					o_phi_t, o_vrt_t, o_div_t,
					i_simulation_timestamp);
		}


		if (use_nr)
		{
			euler_timestep_update_nr(
					i_U_phi, i_U_vrt, i_U_div,
					o_phi_t, o_vrt_t, o_div_t,
					i_simulation_timestamp);
		}
	}
}




void PDESWESphere2DTS_ln_erk_split_uv::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		SWEETErrorFatal("Only constant time step size allowed");

	// standard time stepping
	timestepping_rk.runTimestep(
			this,
			&PDESWESphere2DTS_ln_erk_split_uv::euler_timestep_set_tendencies,	//!< pointer to function to compute euler time step updates
			io_phi, io_vrt, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}





PDESWESphere2DTS_ln_erk_split_uv::PDESWESphere2DTS_ln_erk_split_uv()
{
}



PDESWESphere2DTS_ln_erk_split_uv::~PDESWESphere2DTS_ln_erk_split_uv()
{
}

