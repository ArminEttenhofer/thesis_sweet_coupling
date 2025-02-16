/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <programs/PDE_SWECart2D/TimeOld/SWE_Cart2D_TS_l_irk.hpp>
///#include <sweet/Data/Cart2D/Convert/Complex_DataSpectral_2_DataSpectral.hpp>
///#include <sweet/Data/Cart2D/Convert/DataSpectral_2_Complex_DataSpectral.hpp>
#include <sweet/Data/Cart2D/Convert/DataSpectral_2_Cart2DComplex_DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/Convert/DataSpectral_2_Cart2D_DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/DataSpectral.hpp>





/**
 * Solve SWE with implicit Euler time stepping
 *
 * U_t = L U(0)
 *
 * (U(tau) - U(0)) / tau = L U(tau)
 *
 * <=> U(tau) - U(0) = L U(tau) tau
 *
 * <=> U(tau) - L tau U(tau) = U(0)
 *
 * <=> (1 - L tau) U(tau) = U(0)
 *
 * <=> (1/tau - L) U(tau) = U(0)/tau
 */
void SWE_Cart2D_TS_l_irk::runTimestep(
		sweet::Data::Cart2D::DataSpectral &io_h,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_u,	//!< prognostic variables
		sweet::Data::Cart2D::DataSpectral &io_v,	//!< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETErrorFatal("SWE_Cart2D_TS_l_irk: only constant time step size allowed (Please set --dt)");

#if SWEET_USE_CART2D_SPECTRAL_SPACE

	sweet::Data::Cart2D::DataSpectral &eta0 = io_h;
	sweet::Data::Cart2D::DataSpectral &u0 = io_u;
	sweet::Data::Cart2D::DataSpectral &v0 = io_v;

	double alpha = 1.0/i_dt;

	eta0 *= alpha;
	u0 *= alpha;
	v0 *= alpha;

	// load kappa (k)
	double kappa = alpha*alpha + shackPDESWECart2D->cart2d_rotating_f0*shackPDESWECart2D->cart2d_rotating_f0;

	double eta_bar = shackPDESWECart2D->h0;
	double g = shackPDESWECart2D->gravitation;

	sweet::Data::Cart2D::DataSpectral rhs =
			(kappa/alpha) * eta0
			- eta_bar*(ops->diff_c_x(u0) + ops->diff_c_y(v0))
			- (shackPDESWECart2D->cart2d_rotating_f0*eta_bar/alpha) * (ops->diff_c_x(v0) - ops->diff_c_y(u0))
		;

	sweet::Data::Cart2D::DataSpectral lhs = (-g*eta_bar*(ops->diff2_c_x + ops->diff2_c_y)).spectral_addScalarAll(kappa);
	io_h = rhs.spectral_div_element_wise(lhs);

	sweet::Data::Cart2D::DataSpectral uh = u0 - g*ops->diff_c_x(io_h);
	sweet::Data::Cart2D::DataSpectral vh = v0 - g*ops->diff_c_y(io_h);

	io_u = alpha/kappa * uh     + shackPDESWECart2D->cart2d_rotating_f0/kappa * vh;
	io_v =    -shackPDESWECart2D->cart2d_rotating_f0/kappa * uh + alpha/kappa * vh;

#else

	sweet::Data::Cart2D::DataSpectralComplex eta0 = sweet::Data::Cart2D::Convert::Cart2DDataSpectral_2_Cart2DDataSpectralComplex::grid_convert(io_h);
	sweet::Data::Cart2D::DataSpectralComplex u0 = sweet::Data::Cart2D::Convert::Cart2DDataSpectral_2_Cart2DDataSpectralComplex::grid_convert(io_u);
	sweet::Data::Cart2D::DataSpectralComplex v0 = sweet::Data::Cart2D::Convert::Cart2DDataSpectral_2_Cart2DDataSpectralComplex::grid_convert(io_v);

	double alpha = 1.0/i_dt;

	eta0 *= alpha;
	u0 *= alpha;
	v0 *= alpha;

	// load kappa (k)
	double kappa = alpha*alpha + shackPDESWECart2D->cart2d_rotating_f0*shackPDESWECart2D->cart2d_rotating_f0;

	double eta_bar = shackPDESWECart2D->h0;
	double g = shackPDESWECart2D->gravitation;

	sweet::Data::Cart2D::DataSpectralComplex rhs =
			(kappa/alpha) * eta0
			- eta_bar*(opComplex.diff_c_x(u0) + opComplex.diff_c_y(v0))
			- (shackPDESWECart2D->cart2d_rotating_f0*eta_bar/alpha) * (opComplex.diff_c_x(v0) - opComplex.diff_c_y(u0))
		;

	sweet::Data::Cart2D::DataSpectralComplex lhs = (-g*eta_bar*(opComplex.diff2_c_x + opComplex.diff2_c_y)).spectral_addScalarAll(kappa);
	sweet::Data::Cart2D::DataSpectralComplex eta = rhs.spectral_div_element_wise(lhs);

	sweet::Data::Cart2D::DataSpectralComplex uh = u0 - g*opComplex.diff_c_x(eta);
	sweet::Data::Cart2D::DataSpectralComplex vh = v0 - g*opComplex.diff_c_y(eta);

	sweet::Data::Cart2D::DataSpectralComplex u1 = alpha/kappa * uh     + shackPDESWECart2D->cart2d_rotating_f0/kappa * vh;
	sweet::Data::Cart2D::DataSpectralComplex v1 =    -shackPDESWECart2D->cart2d_rotating_f0/kappa * uh + alpha/kappa * vh;

	io_h = sweet::Data::Cart2D::Convert::Cart2DDataSpectralComplex_2_Cart2DDataSpectral::grid_convert(eta);
	io_u = sweet::Data::Cart2D::Convert::Cart2DDataSpectralComplex_2_Cart2DDataSpectral::grid_convert(u1);
	io_v = sweet::Data::Cart2D::Convert::Cart2DDataSpectralComplex_2_Cart2DDataSpectral::grid_convert(v1);
#endif
}



/*
 * Setup
 */
bool SWE_Cart2D_TS_l_irk::setup(
	sweet::Data::Cart2D::Operators *io_ops,
	int i_order
)
{
	PDESWECart2DTS_BaseInterface::setup(io_ops);

	if (shackCart2DDataOps->space_grid_use_c_staggering)
		SWEETErrorFatal("Staggering not supported for l_irk");

	SWEET_ASSERT(i_order > 0);
	timestepping_order = i_order;

	if (timestepping_order != 1)
		SWEETErrorFatal("SWE_Cart2D_TS_l_irk: Only 1st order IRK is supported. Please set --timestepping-order 1.");

	return true;
}




/*
 * Setup
 */
bool SWE_Cart2D_TS_l_irk::setup(
	sweet::Data::Cart2D::Operators *io_ops
)
{
	return setup(io_ops, shackPDESWETimeDisc->timestepping_order);
}
