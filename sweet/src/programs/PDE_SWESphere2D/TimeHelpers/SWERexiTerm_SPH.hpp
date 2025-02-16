/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEHELPERS_SWEREXITERM_SPH_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEHELPERS_SWEREXITERM_SPH_HPP

#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Sphere2DComplex_DataSpectral.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2DComplex/Convert/DataSpectral_2_Sphere2D_DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/DataGrid.hpp>
#include <sweet/Data/Sphere2DComplex/DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/Operators.hpp>
#include <complex>

#include "SphBandedMatrix_GridComplex.hpp"



/**
 * REXI solver for SWE
 */
class SWERexiTerm_SPH
{
	//! SPH configuration
	const sweet::Data::Sphere2D::Config *sphere2DDataConfig;

	//! Solver for given alpha
	SphBandedMatrix_GridComplex sphSolverComplexDiv;

	sweet::Data::Sphere2DComplex::Operators opsComplex;

	//! REXI alpha
	std::complex<double> alpha;

	//! REXI beta
	std::complex<double> beta;


	//! (real) timestep size
	double timestep_size;

	//! earth radius
	double sphere2d_radius;

	//! inverse of earth radius
	//double ir;

	//! Coriolis omega
	double two_coriolis_omega;

	//! f0
	double f0;

	bool use_f_sphere2D;

	bool no_coriolis;

	//! Average geopotential
	double gh0;

	//! pseudo timestep size to use implicit solver for REXI
	std::complex<double> dt_implicit;

public:
	SWERexiTerm_SPH()	:
		sphere2DDataConfig(nullptr)
	{
	}




	/**
	 * REXI term formulated as backward Euler
	 */
	inline
	void solve_vectorinvariant_progphivortdiv(
			const sweet::Data::Sphere2D::DataSpectral &i_phi0,
			const sweet::Data::Sphere2D::DataSpectral &i_vrt0,
			const sweet::Data::Sphere2D::DataSpectral &i_div0,

			sweet::Data::Sphere2D::DataSpectral &o_phi,
			sweet::Data::Sphere2D::DataSpectral &o_vort,
			sweet::Data::Sphere2D::DataSpectral &o_div
	)
	{
		sweet::Data::Sphere2DComplex::DataSpectral phi1(sphere2DDataConfig);
		sweet::Data::Sphere2DComplex::DataSpectral vrt1(sphere2DDataConfig);
		sweet::Data::Sphere2DComplex::DataSpectral div1(sphere2DDataConfig);

		sweet::Data::Sphere2DComplex::DataSpectral phi0 = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(i_phi0);
		sweet::Data::Sphere2DComplex::DataSpectral vrt0 = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(i_vrt0);
		sweet::Data::Sphere2DComplex::DataSpectral div0 = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(i_div0);

		/*
		 * Preprocessing:
		 * U0* = U0 * beta/alpha
		 * dt_implicit = -timestep_size/alpha
		 */
		phi0 *= beta/alpha;
		div0 *= beta/alpha;
		vrt0 *= beta/alpha;

		if (no_coriolis)
		{
#if 0

			sweet::Data::Sphere2DComplex::DataSpectral rhs = div0 + opsComplex.implicit_L(phi0, dt_implicit);
			div1 = opsComplex.implicit_helmholtz(rhs, gh0*dt_implicit*dt_implicit, sphere2d_radius);

			phi1 = phi0 - dt_implicit*gh0*div1;
			vrt1 = vrt0;

#else

			std::complex<double> dt_two_omega = dt_implicit*0.0;

			sweet::Data::Sphere2DComplex::DataSpectral rhs = div0 + opsComplex.implicit_FJinv(vrt0, dt_two_omega) + opsComplex.implicit_L(phi0, dt_implicit);
			div1 = sphSolverComplexDiv.solve(rhs);
			phi1 = phi0 - dt_implicit*gh0*div1;
			vrt1 = opsComplex.implicit_Jinv(vrt0 - opsComplex.implicit_F(div1, dt_two_omega), dt_two_omega);

#endif
		}
		else if (use_f_sphere2D)
		{
			SWEETErrorFatal("TODO: This needs to be ported to the new formulation");
#if 0
			// TODO
			sweet::Data::Sphere2DComplex::DataSpectral rhs = dt_implicit*gh0*(div0 - f0*vrt0) + (1.0+f0*f0)*phi0;
			phi1 = rhs.spectral_solve_helmholtz_one(-gh0, sphere2d_radius);

			vrt1 = (vrt0 + f0*div1);
			div1 = -1.0/(dt_implicit*gh0)*(phi0 - phi1);
#endif
		}
		else
		{
			std::complex<double> dt_two_omega = dt_implicit*two_coriolis_omega;

			sweet::Data::Sphere2DComplex::DataSpectral rhs = div0 + opsComplex.implicit_FJinv(vrt0, dt_two_omega) + opsComplex.implicit_L(phi0, dt_implicit);
			div1 = sphSolverComplexDiv.solve(rhs);
			phi1 = phi0 - dt_implicit*gh0*div1;
			vrt1 = opsComplex.implicit_Jinv(vrt0 - opsComplex.implicit_F(div1, dt_two_omega), dt_two_omega);
		}

		o_phi = sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(phi1);
		o_vort = sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(vrt1);
		o_div = sweet::Data::Sphere2DComplex::Convert::DataSpectral_2_Sphere2D_DataSpectral::convert_real(div1);
	}


	/*
	 * Setup the SWE REXI solver with SPH
	 */
	bool setup_vectorinvariant_progphivortdiv(
			const sweet::Data::Sphere2D::Config *i_sphere2DDataConfigSolver,
			const sweet::Data::Sphere2D::Shack *i_shackSphere2DDataOps,

			const std::complex<double> &i_alpha,
			const std::complex<double> &i_beta,

			double i_radius,
			double i_coriolis_omega,
			double i_f0,
			double i_avg_geopotential,
			double i_timestepSize,

			bool i_use_f_sphere2D,
			bool i_no_coriolis
	)
	{
		sphere2DDataConfig = i_sphere2DDataConfigSolver;

		sphere2d_radius = i_radius;
		gh0 = i_avg_geopotential;
		no_coriolis = i_no_coriolis;
		timestep_size = i_timestepSize;
		use_f_sphere2D = i_use_f_sphere2D;

		if (use_f_sphere2D)
			f0 = i_f0;
		else
			two_coriolis_omega = 2.0*i_coriolis_omega;

		alpha = i_alpha;
		beta = i_beta;

		opsComplex.setup(sphere2DDataConfig, i_shackSphere2DDataOps);

		dt_implicit = -timestep_size/alpha;

		if (!no_coriolis)
		{
			if (!use_f_sphere2D)
			{
				std::complex<double> dt_two_omega = dt_implicit*two_coriolis_omega;

				sphSolverComplexDiv.setup(sphere2DDataConfig, 4);
				sphSolverComplexDiv.solver_addComponent_implicit_J(dt_two_omega);
				sphSolverComplexDiv.solver_addComponent_implicit_FJinvF(dt_two_omega);
				sphSolverComplexDiv.solver_addComponent_implicit_L(gh0*dt_implicit, dt_implicit, sphere2d_radius);
			}
		}
		else
		{
			std::complex<double> dt_two_omega = dt_implicit*0.0;

			sphSolverComplexDiv.setup(sphere2DDataConfig, 4);
			sphSolverComplexDiv.solver_addComponent_implicit_J(dt_two_omega);
			sphSolverComplexDiv.solver_addComponent_implicit_FJinvF(dt_two_omega);
			sphSolverComplexDiv.solver_addComponent_implicit_L(gh0*dt_implicit, dt_implicit, sphere2d_radius);

		}
		return true;
	}


};


#endif
