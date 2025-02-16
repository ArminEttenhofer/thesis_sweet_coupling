/*
 * Author: Pedor Peixoto <ppeixoto@usp.br>
 * based on stuff from:
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *         
 */
#include "PDESWESphere2DTS_lg_0_lc_n_erk_bv.hpp"



/*
 * Barotropic vorticity equation implementation
 * Details from here: https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pjp/barotropic.pdf
 *
 *  The main prognostic variable will be vorticity (vrt)
 *  We are basically solving only the vorticity equation of the SWE equation written in vort-div formulation
 *  The flow is non-divergent, so div=0 everywhere, therefore, also, we don't need phi (fuild depth)
 *  We only need the streamfunction $\psi$,
 *  but not the the velocity potential (which is zero for this equation).
 * In summary:
 *  - Main prognostic: vrt
 *  - Zero variables: chi, div, phi_pert
 *  - Constant variables: phi_pert
 *  - Diagnostic variables: u, v (obtained from psi and chi)
 */



bool PDESWESphere2DTS_lg_0_lc_n_erk_bv::setup_auto(
		const std::string &i_timestepping_method,
		sweet::Data::Sphere2D::Operators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup_main(
			io_ops,
			shackPDESWETimeDisc->timestepping_order
		);
}

bool PDESWESphere2DTS_lg_0_lc_n_erk_bv::setup_main(
		const sweet::Data::Sphere2D::Operators *io_ops,
		int i_timestepping_order	//!< order of RK time stepping method
)
{
	ops = io_ops;

	timestepping_order = i_timestepping_order;
	timestep_size = shackTimestepControl->currentTimestepSize;

	setupFG();

	return true;
}


void PDESWESphere2DTS_lg_0_lc_n_erk_bv::runTimestep(
		sweet::Data::Sphere2D::DataSpectral &io_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &io_vrt,
		sweet::Data::Sphere2D::DataSpectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{

	/* Calculate velocities and stream function */
	//sweet::Data::Sphere2D::DataGrid ug(io_phi_pert.sphere2DDataConfig);
	//sweet::Data::Sphere2D::DataGrid vg(io_phi_pert.sphere2DDataConfig);

	//sweet::Data::Sphere2D::DataGrid vrtg = io_vrt.toGrid();
	//sweet::Data::Sphere2D::DataGrid divg = io_div.toGrid(); /* this should be zero! */

	//Sphere2DData_Spectral psi = ops->inv_laplace(io_vrt);
	//Sphere2DData_Spectral chi = ops->inv_laplace(io_div); /*this should be zero! */

	//ops->vrtdiv_2_uv(io_vrt, io_div, ug, vg);

	//ops->uv_2_vort(ug, vg);

	// standard time stepping RK
	timestepping_rk.runTimestep(
			this,
			&PDESWESphere2DTS_lg_0_lc_n_erk_bv::euler_timestep_update,	//!< pointer to function to compute euler time step updates
			io_phi_pert, io_vrt, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);

}

void PDESWESphere2DTS_lg_0_lc_n_erk_bv::euler_timestep_update(
		const sweet::Data::Sphere2D::DataSpectral &i_phi, //prog
		const sweet::Data::Sphere2D::DataSpectral &i_vrt, //prog
		const sweet::Data::Sphere2D::DataSpectral &i_div, //prog

		sweet::Data::Sphere2D::DataSpectral &o_phi_t, //updated with euler
		sweet::Data::Sphere2D::DataSpectral &o_vrt_t, //updated with euler
		sweet::Data::Sphere2D::DataSpectral &o_div_t, //updated with euler

		double i_simulation_timestamp
)
{
	//zero tendencies
	o_phi_t.spectral_setZero();
	o_vrt_t.spectral_setZero();
	o_div_t.spectral_setZero();


	// Calculate velocities in physical space
	sweet::Data::Sphere2D::DataGrid u_phys, v_phys;
	ops->vrtdiv_2_uv(i_vrt, i_div, u_phys, v_phys);

	/*
	 * Calculate absolute vorticity in physical space (vrt+f)
	 */
	sweet::Data::Sphere2D::DataGrid abs_vrtg = i_vrt.toGrid()+fg;
	//std::cout << "Vort" << std::endl;
	//i_vrt.spectral_print(6);

	// Nonlinear product (velocity * abs_vort)
	sweet::Data::Sphere2D::DataGrid u_nl = u_phys*abs_vrtg;
	sweet::Data::Sphere2D::DataGrid v_nl = v_phys*abs_vrtg;

	//nonlinear vort and divergence of (velocity * abs_vort)
	sweet::Data::Sphere2D::DataSpectral vrt, div; 
	ops->uv_2_vrtdiv(u_nl, v_nl, vrt, div);

	o_vrt_t -= div; //This is basically the tendency in the Barotropic Vorticity Eq.
	//std::cout << "Vort tendency" << std::endl;
	//o_vrt_t.spectral_print(6);
	
	//Keep div constant
	//o_div_t = o_div_t;

	//Phi stays constant
	//o_phi_t = o_phi_t; 

}


void PDESWESphere2DTS_lg_0_lc_n_erk_bv::printHelp()
{
	std::cout << "	Barotropic Vorticity Equation:" << std::endl;
	std::cout << "		+ lg_0_lc_n_erk_bv" << std::endl;
}

PDESWESphere2DTS_lg_0_lc_n_erk_bv::PDESWESphere2DTS_lg_0_lc_n_erk_bv()	:
		timestepping_order(-1)
{
}



PDESWESphere2DTS_lg_0_lc_n_erk_bv::~PDESWESphere2DTS_lg_0_lc_n_erk_bv()
{
}

