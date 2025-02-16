/*
 * SWE_Cart2D_Normal_Modes.hpp
 *
 *  Created on: 17 Nov 2019
 * Author: Pedro Peixoto <pedrosp@ime.usp.br>
 *
 *      based on previous implementation by Martin Schreiber in swe_cart2d.cpp
 *
 */

#ifndef PROGRAMS_PDE_SWECART2D_PDESWECART2D_NORMALMODES_HPP
#define PROGRAMS_PDE_SWECART2D_PDESWECART2D_NORMALMODES_HPP

#include <sweet/Data/Cart2D/DataGrid.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2D/Operators.hpp>
#include <sweet/Data/Cart2DComplex/DataSpectral.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/ExpIntegration/ExpFunction.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <functional>

#if SWEET_EIGEN
#include <Eigen/Eigenvalues>
#endif


#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/TimeTree/Shack.hpp>
#include <sweet/IO/Shack.hpp>
#include <programs/PDE_SWECart2D/Shack.hpp>

/**
 * SWE Cart2D normal mode
 */

namespace PDE_SWECart2D {
namespace NormalModes {

class NormalModes
{
public:

	sweet::Error::Base error;

	typedef double T;
	typedef std::complex<T> complex;

	sweet::Data::Cart2D::Shack *shackCart2DDataOps;
	sweet::TimeTree::Shack *shackTimestepControl;
	sweet::IO::Shack *shackIOData;

	PDE_SWECart2D::Shack *shackPDESWECart2D;

	NormalModes() :
		shackCart2DDataOps(nullptr),
		shackTimestepControl(nullptr),
		shackIOData(nullptr),
		shackPDESWECart2D(nullptr)
	{
	}

	NormalModes(
			const NormalModes &i_val
	)
	{
		shackCart2DDataOps = i_val.shackCart2DDataOps;
		shackTimestepControl = i_val.shackTimestepControl;
		shackIOData = i_val.shackIOData;
		shackPDESWECart2D = i_val.shackPDESWECart2D;
	}

	bool shackRegistration(sweet::Shacks::Dictionary &io_dict)
	{
		shackCart2DDataOps = io_dict.getAutoRegistration<sweet::Data::Cart2D::Shack>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(io_dict);

		shackTimestepControl = io_dict.getAutoRegistration<sweet::TimeTree::Shack>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(io_dict);

		shackIOData = io_dict.getAutoRegistration<sweet::IO::Shack>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(io_dict);

		shackPDESWECart2D = io_dict.getAutoRegistration<PDE_SWECart2D::Shack>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(io_dict);

		return true;
	}

	void add_normal_mode(
			std::size_t ik0,				//wavenumber in x
			std::size_t ik1,				// wavenumeber in y
			double geo_mode,    //Coeficient multiplying geostrophic mode
			double igwest_mode, //Coeficient multiplying west gravity mode
			double igeast_mode, //Coeficient multiplying east gravity mode
			sweet::Data::Cart2D::DataSpectral &io_h, // h: surface height (perturbation)
			sweet::Data::Cart2D::DataSpectral &io_u, // u: velocity in x-direction
			sweet::Data::Cart2D::DataSpectral &io_v // v: velocity in y-direction
	)
	{
		sweet::ExpIntegration::ExpFunction<T> rexiFunctions;

		const sweet::Data::Cart2D::Config *cart2DDataConfig = io_h.cart2DDataConfig;

		if (shackCart2DDataOps->space_grid_use_c_staggering)
			SWEETErrorFatal("Staggering not supported");
		
		//std::cout << "Adding mode to fields" <<std::endl;

		//Check if k0 is in correct sprectral area
		//std::cout << io_h.cart2DDataConfig->spectral_data_size[1] << std::endl;
		if ( ik1<0 || ik1 >= cart2DDataConfig->spectral_data_size[1]) 
			SWEETErrorFatal("Normal_mode: mode not within reach");

		if ( ik0<0 || ik0 >= cart2DDataConfig->spectral_data_size[0]) 
			SWEETErrorFatal("Normal_mode: mode not within reach");
	
		//Check for mirror effects
		T k0 = (T)ik0;
		T k1;
		if (ik1 < cart2DDataConfig->spectral_data_size[1]/2)
			k1 = (T)ik1;
		else
			k1 = (T)((int)ik1-(int)cart2DDataConfig->spectral_data_size[1]);
		
		complex v[3][3];
		complex lambda[3];

		sw_eigen_decomp(
				k0,				//wavenumber in x
				k1,				// wavenumeber in y
				false, // Direct EV matrix (not inverse)
				v , // EV matrix
				lambda // output eigen values */
		);

		complex U[3];
		// Set normal mode acording to desired wave type
		// These are weights for the modes
		U[0] = geo_mode;
		U[1] = igwest_mode;
		U[2] = igeast_mode;


		//Define normal mode as combination of eigen vectors
		complex UEV[3] = {0.0, 0.0, 0.0};
		for (int k = 0; k < 3; k++)
			for (int j = 0; j < 3; j++)
				UEV[k] += v[k][j] * U[j];

		//std::cout << "spectral before" << std::endl;
		//io_v.print_spectralIndex();

		complex h_add, u_add, v_add;
		h_add = io_h.spectral_get(ik1, ik0)+UEV[0];
		u_add = io_u.spectral_get(ik1, ik0)+UEV[1];
		v_add = io_v.spectral_get(ik1, ik0)+UEV[2];

		/* Add normal mode to data */
		io_h.spectral_set(ik1, ik0, h_add);
		io_u.spectral_set(ik1, ik0, u_add);
		io_v.spectral_set(ik1, ik0, v_add);

		io_h.spectral_zeroAliasingModes();
		io_u.spectral_zeroAliasingModes();
		io_v.spectral_zeroAliasingModes();

/*Debug output*/
#if 0

		std::cout << "EV matrix" <<std::endl;
		for (int j = 0; j < 3; j++)	{
			for (int i = 0; i < 3; i++)
				std::cout << v[j][i]<< " "; 
			std::cout <<std::endl;
		}
	
		std::cout << "Eigen values" <<std::endl;
		for (int j = 0; j < 3; j++)
			std::cout << lambda[j]<< " "; 
		std::cout <<std::endl;
		
		std::cout << "Adding normal mode" <<std::endl;
		for (int j = 0; j < 3; j++)
			std::cout << UEV[j]<< " "; 
		std::cout <<std::endl;
	
		std::cout <<ik0<< " " <<ik1<< " " << io_v.p_spectral_get(ik1, ik0) << UEV[2] << std::endl;

#endif
		return;
	}


	void convert_allspectralmodes_2_normalmodes(
			sweet::Data::Cart2D::DataSpectral &i_h, // h: surface height (perturbation)
			sweet::Data::Cart2D::DataSpectral &i_u, // u: velocity in x-direction
			sweet::Data::Cart2D::DataSpectral &i_v, // v: velocity in y-direction
			sweet::Data::Cart2D::DataSpectral &o_geo_mode,    //Output: Coeficients multiplying geostrophic mode
			sweet::Data::Cart2D::DataSpectral &o_igwest_mode, //Output: Coeficients multiplying west gravity mode
			sweet::Data::Cart2D::DataSpectral &o_igeast_mode //Output: Coeficients multiplying east gravity mode
	)
	{
		const sweet::Data::Cart2D::Config *cart2DDataConfig = i_h.cart2DDataConfig;

		o_geo_mode.spectral_setZero();
		o_igwest_mode.spectral_setZero();
		o_igeast_mode.spectral_setZero();

		complex geo_mode_c;
		complex igwest_mode_c;
		complex igeast_mode_c;
		for (std::size_t ik1 = 0; ik1 < cart2DDataConfig->spectral_data_size[1]; ik1++)
		{
			for (std::size_t ik0 = 0; ik0 < cart2DDataConfig->spectral_data_size[0]; ik0++)
			{
				
				convert_spectralmode_2_normalmode(
									ik0, ik1,
									i_h,
									i_u,
									i_v,
									geo_mode_c,
									igwest_mode_c,
									igeast_mode_c
							);
				o_geo_mode.spectral_set(ik1, ik0, geo_mode_c);
				o_igwest_mode.spectral_set(ik1, ik0, igwest_mode_c);
				o_igeast_mode.spectral_set(ik1, ik0, igeast_mode_c);

			}
		}	
		return;
	}

	void convert_spectralmode_2_normalmode(
			std::size_t ik0,				//wavenumber in x
			std::size_t ik1,				// wavenumeber in y
			sweet::Data::Cart2D::DataSpectral &i_h, // h: surface height (perturbation)
			sweet::Data::Cart2D::DataSpectral &i_u, // u: velocity in x-direction
			sweet::Data::Cart2D::DataSpectral &i_v, // v: velocity in y-direction
			complex &o_geo_mode,    //Output: Coeficient multiplying geostrophic mode
			complex &o_igwest_mode, //Output: Coeficient multiplying west gravity mode
			complex &o_igeast_mode //Output: Coeficient multiplying east gravity mode
	)
	{

		sweet::ExpIntegration::ExpFunction<T> rexiFunctions;

		const sweet::Data::Cart2D::Config *cart2DDataConfig = i_h.cart2DDataConfig;

		if (shackCart2DDataOps->space_grid_use_c_staggering)
			SWEETErrorFatal("Staggering not supported");
		
		//std::cout << "Adding mode to fields" <<std::endl;

		//Check if k0 is in correct sprectral area
		//std::cout << io_h.cart2DDataConfig->spectral_data_size[1] << std::endl;
		if ( ik1<0 || ik1 >= cart2DDataConfig->spectral_data_size[1]) 
			SWEETErrorFatal("Normal_mode: mode not within reach");

		if ( ik0<0 || ik0 >= cart2DDataConfig->spectral_data_size[0]) 
			SWEETErrorFatal("Normal_mode: mode not within reach");
	
		//Check for mirror effects
		T k0 = (T)ik0;
		T k1;
		if (ik1 < cart2DDataConfig->spectral_data_size[1]/2)
			k1 = (T)ik1;
		else
			k1 = (T)((int)ik1-(int)cart2DDataConfig->spectral_data_size[1]);
		
		complex v[3][3];
		complex lambda[3];

		sw_eigen_decomp(
				k0,				//wavenumber in x
				k1,				// wavenumeber in y
				true, // Inverse ev matrix
				v , // inverse EV matrix
				lambda // output eigen values */
		);

		complex U[3];
		// Set (h,u,v) spectral coeficients
		// These are weights for the modes
		U[0] = i_h.spectral_get(ik1, ik0);
		U[1] = i_u.spectral_get(ik1, ik0);
		U[2] = i_v.spectral_get(ik1, ik0);


		//Apply inverse EV matrix to obtain data in EV space
		complex UEV[3] = {0.0, 0.0, 0.0};
		for (int k = 0; k < 3; k++)
			for (int j = 0; j < 3; j++)
				UEV[k] += v[k][j] * U[j];

		//Return the modes
		o_geo_mode=UEV[0];
		o_igwest_mode=UEV[1];
		o_igeast_mode=UEV[2];
		//std::cout << " Geost: " << o_geo_mode<<std::endl;
		//std::cout << " IGWest: " << o_igwest_mode<<std::endl;
		//std::cout << " IGEast: " << o_igeast_mode<<std::endl;

		return;		
	}


	/* Get linear shallow water operator eigen decomposition */
	/* Default case considering the full linear operator */
	void sw_eigen_decomp(
			T k0,				//wavenumber in x
			T k1,				// wavenumeber in y
			bool i_inverse,		// Input true, returns inverse matriz, false: returns direct
			complex o_v[3][3],	// output eigen vector (direct or inverse)
			complex o_evalues[3]	// output eigen values (optional)
	)
	{
		sw_eigen_decomp(
					k0,
					k1,
					i_inverse,
					o_v,
					o_evalues,
					shackPDESWECart2D->cart2d_rotating_f0,
					shackPDESWECart2D->h0,
					shackPDESWECart2D->gravitation
				);
	}


	/* Get linear shallow water operator eigen decomposition */
	/* Passing f, h, g as arguments allow to compute eigendecomposition of L_g and L_c separately*/
	void sw_eigen_decomp(
			T k0,				//wavenumber in x
			T k1,				// wavenumeber in y
			bool i_inverse,		// Input true, returns inverse matriz, false: returns direct
			complex o_v[3][3],	// output eigen vector (direct or inverse)
			complex o_evalues[3],	// output eigen values (optional)
			T f,			// Coriolis
			T h,			// h0
			T g
	)
	{
		sweet::ExpIntegration::ExpFunction<T> rexiFunctions;
		bool i_evalues = false;
		if (o_evalues){
			i_evalues = true;
			//std::cout <<i_evalues<< " " << o_evalues[0]   << std::endl;
		}
		else{
			i_evalues = false;
		}
		//std::cout <<o_evalues[1]<<std::endl;
		//std::cout <<i_inverse<<std::endl;

		if (shackCart2DDataOps->space_grid_use_c_staggering)
			SWEETErrorFatal("Staggering not supported");
		
		complex I(0.0, 1.0);
		//std::cout << "Calculating EV for mode (" << k0 << ", " << k1 << ")" << std::endl;
		//std::cout << "hi" << std::endl;
		T s0 = shackCart2DDataOps->cart2d_domain_size[0];
		T s1 = shackCart2DDataOps->cart2d_domain_size[1];

		///T f = shackPDESWECart2D->cart2d_rotating_f0;
		///T h = shackPDESWECart2D->h0;
		///T g = shackPDESWECart2D->gravitation;

		T sqrt_h = rexiFunctions.l_sqrt(h);
		T sqrt_g = rexiFunctions.l_sqrt(g);

		complex b = -k0*I;	// d/dx exp(I*k0*x) = I*k0 exp(I*k0*x)
		complex c = -k1*I;

		b = b*rexiFunctions.pi2/s0;
		c = c*rexiFunctions.pi2/s1;

		/*
		 * Matrix with Eigenvectors (column-wise)
		 */
		complex v[3][3];
		complex v_inv[3][3];

		/*
		 * Eigenvalues
		 */
		complex lambda[3];

		///if (shackPDESWECart2D->cart2d_rotating_f0 == 0)
		if (f == 0)
		{
			/*
			 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c%7D,%7Bg*b,0,0%7D,%7Bg*c,0,0%7D%7D
			 */
			if (k0 == 0 && k1 == 0)
			{
				v[0][0] = 1;
				v[1][0] = 0;
				v[2][0] = 0;

				v[0][1] = 0;
				v[1][1] = 1;
				v[2][1] = 0;

				v[0][2] = 0;
				v[1][2] = 0;
				v[2][2] = 1;

				if (i_evalues){
					lambda[0] = 0;
					lambda[1] = 0;
					lambda[2] = 0;
				}
			}
			else if (k0 == 0)
			{
				v[0][0] = 0;
				v[1][0] = 1;
				v[2][0] = 0;

				v[0][1] = -sqrt_h/sqrt_g;
				v[1][1] = 0;
				v[2][1] = 1;

				v[0][2] = sqrt_h/sqrt_g;
				v[1][2] = 0;
				v[2][2] = 1;

				if (i_evalues){
					lambda[0] = 0;
					lambda[1] = -c*sqrt_g*sqrt_h;
					lambda[2] = c*sqrt_g*sqrt_h;;
				}
			}
			else if (k1 == 0)
			{
				/*
				 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c*0%7D,%7Bg*b,0,0%7D,%7Bg*c*0,0,0%7D%7D
				 */

				v[0][0] = 0;
				v[1][0] = 0;
				v[2][0] = 1;

				v[0][1] = -sqrt_h/sqrt_g;
				v[1][1] = 1;
				v[2][1] = 0;

				v[0][2] = sqrt_h/sqrt_g;
				v[1][2] = 1;
				v[2][2] = 0;

				if (i_evalues){
					lambda[0] = 0;
					lambda[1] = -b*sqrt_g*sqrt_h;
					lambda[2] = b*sqrt_g*sqrt_h;
				}
			}
			else
			{
				v[0][0] = 0;
				v[1][0] = -c/b;
				v[2][0] = 1.0;

				v[0][1] = -(sqrt_h*rexiFunctions.l_sqrtcplx(b*b + c*c))/(c*sqrt_g);
				v[1][1] = b/c;
				v[2][1] = 1.0;

				v[0][2] = (sqrt_h*rexiFunctions.l_sqrtcplx(b*b + c*c))/(c*sqrt_g);
				v[1][2] = b/c;
				v[2][2] = 1.0;

				if (i_evalues){
					lambda[0] = 0.0;
					lambda[1] = -rexiFunctions.l_sqrtcplx(b*b + c*c)*sqrt_h*sqrt_g;
					lambda[2] = rexiFunctions.l_sqrtcplx(b*b + c*c)*sqrt_h*sqrt_g;
				}
			}
		}
		else if (g == 0 && h == 0)
		{
			/*
			 * https://www.wolframalpha.com/input?i=eigenvector%7B%7B0%2C0%2C0%7D%2C%7B0%2C0%2Cf%7D%2C%7B0%2C-f%2C0%7D%7D
			 */
			v[0][0] = 0.;
			v[1][0] = -I;
			v[2][0] = 1.;

			v[0][1] = 0.;
			v[1][1] = I;
			v[2][1] = 1.;

			v[0][2] = 1.;
			v[1][2] = 0.;
			v[2][2] = 0.;

			if (i_evalues){
				lambda[0] = I * f;
				lambda[1] = -I * f;
				lambda[2] = 0.;
			}
		}
		else
		{
			if (k0 == 0 && k1 == 0)
			{
				/*
				 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,0,0%7D,%7B0,0,f%7D,%7B0,-f,0%7D%7D
				 */
				v[0][0] = 0;
				v[1][0] = -I;
				v[2][0] = 1;

				v[0][1] = 0;
				v[1][1] = I;
				v[2][1] = 1;

				v[0][2] = 1;
				v[1][2] = 0;
				v[2][2] = 0;

				if (i_evalues){
					lambda[0] = I*f;
					lambda[1] = -I*f;
					lambda[2] = 0;
				}
			}
			else if (k0 == 0)
			{
				/*
				 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b*0,h*c%7D,%7Bg*b*0,0,f%7D,%7Bg*c,-f,0%7D%7D
				 */
				v[0][0] = f/(c*g);
				v[1][0] = 1;
				v[2][0] = 0;

				v[0][1] = -(c*h)/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
				v[1][1] =  -f/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
				v[2][1] = 1;

				v[0][2] = (c*h)/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
				v[1][2] = f/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
				v[2][2] = 1;

				if (i_evalues){
					lambda[0] = 0;
					lambda[1] = -rexiFunctions.l_sqrtcplx(c*c*g*h-f*f);
					lambda[2] = rexiFunctions.l_sqrtcplx(c*c*g*h-f*f);
				}
			}
			else if (k1 == 0)
			{
					/*
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c*0%7D,%7Bg*b,0,f%7D,%7Bg*c*0,-f,0%7D%7D
					 */
				v[0][0] = -f/(b*g);
				v[1][0] = 0;
				v[2][0] = 1;

				v[0][1] = -(b*h)/f;
				v[1][1] = rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h)/f;
				v[2][1] = 1;

				v[0][2] = -(b*h)/f;
				v[1][2] = -rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h)/f;
				v[2][2] = 1;

				if (i_evalues){
					lambda[0] = 0;
					lambda[1] = -rexiFunctions.l_sqrtcplx(b*b*g*h-f*f);
					lambda[2] = rexiFunctions.l_sqrtcplx(b*b*g*h-f*f);
				}
			}
			else
				{
					/*
					 * Compute EV's of
					 * Linear operator
					 *
					 * [ 0  hb  hc ]
					 * [ gb  0   f ]
					 * [ gc -f   0 ]
					 *
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c%7D,%7Bg*b,0,f%7D,%7Bg*c,-f,0%7D%7D
					 */

					v[0][0] = -f/(b*g);
					v[1][0] = -c/b;
					v[2][0] = 1.0;

					v[0][1] = -(c*f*h + b*h*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h))/(b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[1][1] = -(f*f - b*b*g*h)/(b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[2][1] = 1.0;

					v[0][2] = -(-c*f*h + b*h*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h))/(-b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[1][2] =  -(-f*f + b*b*g*h)/(-b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[2][2] = 1.0;

					if (i_evalues){
						lambda[0] = 0.0;
						lambda[1] = -rexiFunctions.l_sqrtcplx(b*b*g*h + c*c*g*h - f*f);
						lambda[2] =  rexiFunctions.l_sqrtcplx(b*b*g*h + c*c*g*h - f*f);
					}
				}
		}

			/*
			 * Invert Eigenvalue matrix
			 */

		if (i_inverse){
			v_inv[0][0] =  (v[1][1]*v[2][2] - v[1][2]*v[2][1]);
			v_inv[0][1] = -(v[0][1]*v[2][2] - v[0][2]*v[2][1]);
			v_inv[0][2] =  (v[0][1]*v[1][2] - v[0][2]*v[1][1]);

			v_inv[1][0] = -(v[1][0]*v[2][2] - v[1][2]*v[2][0]);
			v_inv[1][1] =  (v[0][0]*v[2][2] - v[0][2]*v[2][0]);
			v_inv[1][2] = -(v[0][0]*v[1][2] - v[0][2]*v[1][0]);

			v_inv[2][0] =  (v[1][0]*v[2][1] - v[1][1]*v[2][0]);
			v_inv[2][1] = -(v[0][0]*v[2][1] - v[0][1]*v[2][0]);
			v_inv[2][2] =  (v[0][0]*v[1][1] - v[0][1]*v[1][0]);

			complex s = v[0][0]*v_inv[0][0] + v[0][1]*v_inv[1][0] + v[0][2]*v_inv[2][0];

			for (int j = 0; j < 3; j++)	{
				for (int i = 0; i < 3; i++)
					v_inv[j][i] /= s;
			}
			//Return inverse matrix
			for (int j = 0; j < 3; j++)	{
				for (int i = 0; i < 3; i++)
					o_v[j][i] = v_inv[j][i] ;
			}
		}	
		else{
			//Return direct matrix
			for (int j = 0; j < 3; j++)	{
				for (int i = 0; i < 3; i++)
					o_v[j][i] = v[j][i] ;
			}
		}
		if (i_evalues){
			for (int j = 0; j < 3; j++)	{
				o_evalues[j] = lambda[j] ;
			}
		}
		return;
}



	template <typename TCallbackClass>
	void normal_mode_analysis(
			sweet::Data::Cart2D::DataSpectral &io_prog_h_pert, // h: surface height (perturbation)
			sweet::Data::Cart2D::DataSpectral &io_prog_u, // u: velocity in x-direction
			sweet::Data::Cart2D::DataSpectral &io_prog_v, // v: velocity in y-direction
			int number_of_prognostic_variables,
			sweet::Shacks::Dictionary *shackDict, // Simulation variables
			TCallbackClass *i_class,
			bool(TCallbackClass::* const i_run_timestep_method)(void)
	)
	{

		const sweet::Data::Cart2D::Config *cart2DDataConfig = io_prog_h_pert.cart2DDataConfig;

		// dummy time step to get time step size
		if (shackTimestepControl->currentTimestepSize <= 0)
			SWEETErrorFatal("Normal mode analysis requires setting fixed time step size");

		/*
		 *
		 * Mode-wise normal mode analysis
		 *
		 *
		 */

		if (shackPDESWECart2D->normal_mode_analysis_generation == 4)
		{
#if SWEET_EIGEN
#if SWEET_USE_CART2D_SPECTRAL_DEALIASING
			SWEETErrorFatal("SWE_Cart2D_Normal_Modes: This test was build for linear or linearized models, so please compile without dealising --cart2d-spectral-dealiasing=disable.");
#endif


			/*
			 * Setup all output files
			 */
			const char* filename; //general filename
			char buffer_real[1024];

			if (shackIOData->outputFileName == "")
				filename = "output_%s_t%020.8f.csv";
			else
				filename = shackIOData->outputFileName.c_str();

			sprintf(buffer_real, filename, "normal_modes_cart2d", shackTimestepControl->currentTimestepSize*shackIOData->outputFormatTimeScale);
			std::ofstream file(buffer_real, std::ios_base::trunc);
			std::cout << "Writing normal mode analysis to files of the form '" << buffer_real << "'" << std::endl;

			//Positive inertia-gravity modes
			sprintf(buffer_real, filename, "normal_modes_cart2d_igpos", shackTimestepControl->currentTimestepSize*shackIOData->outputFormatTimeScale);
			std::ofstream file_igpos(buffer_real, std::ios_base::trunc);

			//Negative inertia-gravity modes
			sprintf(buffer_real, filename, "normal_modes_cart2d_igneg", shackTimestepControl->currentTimestepSize*shackIOData->outputFormatTimeScale);
			std::ofstream file_igneg(buffer_real, std::ios_base::trunc);

			//Geostrophic modes
			sprintf(buffer_real, filename, "normal_modes_cart2d_geo", shackTimestepControl->currentTimestepSize*shackIOData->outputFormatTimeScale);
			std::ofstream file_geo(buffer_real, std::ios_base::trunc);

			//std::cout << "WARNING: OUTPUT IS TRANSPOSED!" << std::endl;

			// use very high precision
			file << std::setprecision(20);
			file_igpos << std::setprecision(20);
			file_igneg << std::setprecision(20);
			file_geo << std::setprecision(20);

			file << "# dt " << shackTimestepControl->currentTimestepSize << std::endl;
			file << "# g " << shackPDESWECart2D->gravitation << std::endl;
			file << "# h " << shackPDESWECart2D->h0 << std::endl;
			file << "# f " << shackPDESWECart2D->cart2d_rotating_f0 << std::endl;

#if SWEET_USE_CART2D_SPECTRAL_SPACE
			int specmodes = cart2DDataConfig->getSpectralIterationRangeArea(0)+cart2DDataConfig->getSpectralIterationRangeArea(1);
			file << "# specnummodes " << specmodes << std::endl;
			file << "# specrealresx " << cart2DDataConfig->spectral_representable_modes[0] << std::endl;
			file << "# specrealresy " << cart2DDataConfig->spectral_representable_modes[1] << std::endl;
#endif

			file << "# physresx " << cart2DDataConfig->grid_res[0] << std::endl;
			file << "# physresy " << cart2DDataConfig->grid_res[1] << std::endl;
			file << "# normalmodegeneration " << shackPDESWECart2D->normal_mode_analysis_generation << std::endl;
			file << "# antialiasing ";
#if SWEET_USE_CART2D_SPECTRAL_DEALIASING
			file << 1;
#else
			file << 0;
#endif
			file << std::endl;

			sweet::Data::Cart2D::DataSpectral* prog[3] = {&io_prog_h_pert, &io_prog_u, &io_prog_v};

			int number_of_prognostic_variables = 3;
			//The basic state is with zero in all variables
			// The only non zero variable in the basic state is the total height
			//    for which the constant is added within runTimestep()
			io_prog_h_pert.spectral_setZero();
			io_prog_u.spectral_setZero();
			io_prog_v.spectral_setZero();

			//int num_timesteps = 1;

			// Timestep and perturbation
			double dt = shackTimestepControl->currentTimestepSize;
			double eps = dt;

			//Matrix representing discrete linear operator in spectral space
			Eigen::MatrixXcf A(3,3) ;
			//Eigen solver
			Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces;
			//Final eigenvalues
			std::complex<double> eval[3];

			//For each spectral mode
			//for (int r = 0; r < 2; r++) //only required to get the symmetric half of the spectrum
			//{
			int r = 0;

			for (std::size_t i = cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]; i < cart2DDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
			{
				std::cout << "." << std::flush;
				for (std::size_t j = cart2DDataConfig->spectral_data_iteration_ranges[r][1][0]; j < cart2DDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
				{
					//This is the mode to be analysed
					//std::cout << "Mode (i,j)= (" << i << " , " << j << ")" << std::endl;


					for (int outer_prog_id = 0; outer_prog_id < number_of_prognostic_variables; outer_prog_id++)
					{

						// reset time control
						shackTimestepControl->currentTimestepNr = 0;
						shackTimestepControl->currentSimulationTime = 0;

						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
							prog[inner_prog_id]->spectral_setZero();

						// activate mode via real coefficient
						prog[outer_prog_id]->spectral_set(j, i, 1.0);
						//Activate the symetric couterpart of the mode (only needed if j>0 )
						if (j > 0)
							prog[outer_prog_id]->spectral_set(cart2DDataConfig->spectral_data_size[1]-j, i, 1.0);

						/*
						 * RUN timestep
						 */
						////prog[outer_prog_id]->request_data_physical();
						(i_class->*i_run_timestep_method)();

						/*
						 * compute
						 * 1/dt * (U(t+1) - U(t))
						 */
						///////prog[outer_prog_id]->request_data_spectral();

						std::complex<double> val = prog[outer_prog_id]->spectral_get(j, i);
						val = val - 1.0; //subtract U(0) from mode
						prog[outer_prog_id]->spectral_set(j, i, val);

						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
							(*prog[inner_prog_id]) /= eps;

						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
						{
							A(inner_prog_id,outer_prog_id)=prog[inner_prog_id]->spectral_get(j, i);;
						}

					}

					//std::cout << "Lik matrix" << std::endl;
					//std::cout << A << std::endl;

					//std::cout << "Normal modes" << std::endl;
					ces.compute(A);
					for(int i=0; i<3; i++)
					{
						eval[i]=ces.eigenvalues()[i];
						//std::cout << "Eigenvalue " << i << " : " << eval[i].real() << " " <<eval[i].imag() << std::endl;

					}
					/* We will try to separate the modes in 3 types:
					 * -positive inertia-gravity (imag>f) - we will adopt coriolis f to test as if > zero, since the exact freq is sqrt(f^2+cK*K)
					 * -negative inertia-gravity (imag<-f)
					 * -negative inertia-gravity (imag aprox 0) - we will fit all other modes here
					 */
					int count_igpos=0;
					int count_igneg=0;
					int count_geo=0;
					for(int i=0; i<3; i++)
					{
						if (eval[i].imag() > 0.5 * shackPDESWECart2D->cart2d_rotating_f0)
						{
							//std::cout << "IG pos mode: " << eval[i].imag() << std::endl;
							//file_igpos << eval[i].imag();
							file_igpos << eval[i].real()<< "\t" << eval[i].imag();
							file_igpos << "\t";
							count_igpos++;
						}
						if (eval[i].imag() < - 0.5 * shackPDESWECart2D->cart2d_rotating_f0)
						{
							//std::cout << "IG neg mode: " << eval[i].imag() << std::endl;
							//file_igneg << eval[i].imag();
							file_igneg << eval[i].real()<< "\t" << eval[i].imag();
							file_igneg << "\t";
							count_igneg++;
						}
						if (eval[i].imag() >= - 0.5 * shackPDESWECart2D->cart2d_rotating_f0 && eval[i].imag() <=  0.5 * shackPDESWECart2D->cart2d_rotating_f0 )
						{
							//std::cout << "IG geo mode: " << eval[i].imag() << std::endl;
							//file_geo << eval[i].imag();
							file_geo << eval[i].real()<< "\t" << eval[i].imag();
							file_geo << "\t";
							count_geo++;
						}
					}
					//Check if we got the correct modes
					if ( count_igpos * count_igneg * count_geo > 0 )
					{
						count_igpos=0;
						count_igneg=0;
						count_geo=0;
					}
					else
					{
						SWEETErrorFatal("SWE_Cart2D_Normal_Modes: Could not separate modes!!");
					}

					//std::cout << "-------------------------" << std::endl;
				}
				file_igpos << std::endl;
				file_igneg << std::endl;
				file_geo << std::endl;
			}

			//}
			//std::cout << "-------------------------" << std::endl;
			//SWEETErrorFatal("still needs work...");
#else
			SWEETErrorFatal("SWE_Cart2D_Normal_Modes: Cannot test this without Eigen library. Please compile with --eigen=enable");
#endif
		}
		/*
		 * Do a normal mode analysis using perturbation, see
		 * Hillary Weller, John Thuburn, Collin J. Cotter,
		 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
		 */
		else
		{

			//runTimestep();
			const char* filename;
			char buffer_real[1024];

			if (shackIOData->outputFileName == "")
				filename = "output_%s_normalmodes.csv";
			else
				filename = shackIOData->outputFileName.c_str();


			sprintf(buffer_real, filename, "normal_modes_physical", shackTimestepControl->currentTimestepSize*shackIOData->outputFormatTimeScale);
			std::ofstream file(buffer_real, std::ios_base::trunc);
			std::cout << "Writing normal mode analysis to file '" << buffer_real << "'" << std::endl;

			std::cout << "WARNING: OUTPUT IS TRANSPOSED!" << std::endl;

			// use very high precision
			file << std::setprecision(20);

			sweet::Data::Cart2D::DataSpectral* prog[3] = {&io_prog_h_pert, &io_prog_u, &io_prog_v};

			/*
			 * Maximum number of prognostic variables
			 *
			 * Advection e.g. has only one
			 */
			if (number_of_prognostic_variables <= 0)
				SWEETErrorFatal("shackDict.pde.number_of_prognostic_variables must be set!");

			if (number_of_prognostic_variables == 3)
			{
				io_prog_h_pert.spectral_setZero();
				io_prog_u.spectral_setZero();
				io_prog_v.spectral_setZero();
			}
			else if (number_of_prognostic_variables == 1)
			{
				io_prog_h_pert.spectral_setZero();
			}
			else
			{
				SWEETErrorFatal("Not yet supported");
			}

#if 0
			if (i_shackDict.disc.timestepping_method == sweet::Shacks::Dictionary::Discretization::LEAPFROG_EXPLICIT)
			{
				SWEETErrorFatal("Not yet tested and supported");
				std::cout << "WARNING: Leapfrog time stepping doesn't make real sense since 1st step is based on RK-like method" << std::endl;
				std::cout << "We'll do two Leapfrog time steps here to take the LF errors into account!" << std::endl;
				std::cout << "Therefore, we also halve the time step size here" << std::endl;

				shackTimestepControl->currentTimestepSize = 0.5*i_shackDict.sim.CFL;
				i_shackDict.sim.CFL = -shackTimestepControl->currentTimestepSize;
			}
#endif

			int num_timesteps = 1;
			if (shackPDESWECart2D->normal_mode_analysis_generation >= 10)
			{
				if (shackTimestepControl->maxTimestepsNr > 0)
					num_timesteps = shackTimestepControl->maxTimestepsNr;
			}

			if (shackTimestepControl->maxSimulationTime > 0)
				file << "# t " << shackTimestepControl->maxSimulationTime << std::endl;
			else
				file << "# t " << (num_timesteps*(-shackTimestepControl->currentTimestepSize)) << std::endl;

			file << "# g " << shackPDESWECart2D->gravitation << std::endl;
			file << "# h " << shackPDESWECart2D->h0 << std::endl;
//			file << "# r " << shackSWECart2D->sphere2d_radius << std::endl;
			file << "# f " << shackPDESWECart2D->cart2d_rotating_f0 << std::endl;

#if SWEET_USE_CART2D_SPECTRAL_SPACE
			int specmodes = cart2DDataConfig->getSpectralIterationRangeArea(0)+cart2DDataConfig->getSpectralIterationRangeArea(1);
			file << "# specnummodes " << specmodes << std::endl;
			file << "# specrealresx " << cart2DDataConfig->spectral_representable_modes[0] << std::endl;
			file << "# specrealresy " << cart2DDataConfig->spectral_representable_modes[1] << std::endl;
#endif

			file << "# physresx " << cart2DDataConfig->grid_res[0] << std::endl;
			file << "# physresy " << cart2DDataConfig->grid_res[1] << std::endl;
			file << "# normalmodegeneration " << shackPDESWECart2D->normal_mode_analysis_generation << std::endl;
			file << "# antialiasing ";

#if SWEET_USE_CART2D_SPECTRAL_DEALIASING
			file << 1;
#else
			file << 0;
#endif

			file << std::endl;


			// iterate over all prognostic variables
			for (int outer_prog_id = 0; outer_prog_id < number_of_prognostic_variables; outer_prog_id++)
			{
				if (shackPDESWECart2D->normal_mode_analysis_generation == 1 || shackPDESWECart2D->normal_mode_analysis_generation == 11)
				{
					// iterate over physical space
					for (std::size_t outer_i = 0; outer_i < cart2DDataConfig->grid_number_elements; outer_i++)
					{
						// reset time control
						shackTimestepControl->currentTimestepNr = 0;
						shackTimestepControl->currentSimulationTime = 0;

						std::cout << "." << std::flush;

						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
							prog[inner_prog_id]->spectral_setZero();

						// activate mode
						///prog[outer_prog_id]->request_data_physical();
						///prog[outer_prog_id]->grid_space_data[outer_i] = 1;
						sweet::Data::Cart2D::DataGrid tmp = prog[outer_prog_id]->toGrid();
						tmp.grid_space_data[outer_i] = 1;
						prog[outer_prog_id]->loadCart2DDataGrid(tmp);

						/*
						 * RUN timestep
						 */

						(i_class->*i_run_timestep_method)();

						if (shackPDESWECart2D->normal_mode_analysis_generation == 1)
						{
							/*
							 * compute
							 * 1/dt * (U(t+1) - U(t))
							 */
							tmp = prog[outer_prog_id]->toGrid();
							tmp.grid_space_data[outer_i] -= 1.0;
							prog[outer_prog_id]->loadCart2DDataGrid(tmp);

							for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
								(*prog[inner_prog_id]) /= shackTimestepControl->currentTimestepSize;
						}

						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
						{
							tmp = prog[outer_prog_id]->toGrid();
							for (std::size_t k = 0; k < cart2DDataConfig->grid_number_elements; k++)
							{
								file << tmp.grid_space_data[k];
								if (inner_prog_id != number_of_prognostic_variables-1 || k != cart2DDataConfig->grid_number_elements-1)
									file << "\t";
								else
									file << std::endl;
							}
						}
					}
				}
#if 1
				else if (shackPDESWECart2D->normal_mode_analysis_generation == 3 || shackPDESWECart2D->normal_mode_analysis_generation == 13)
				{
#if !SWEET_USE_CART2D_SPECTRAL_SPACE
					SWEETErrorFatal("Only available with if cart2d spectral space is activated during compile time!");
#else

					// iterate over spectral space
					for (int r = 0; r < 2; r++)
					{

						for (std::size_t j = cart2DDataConfig->spectral_data_iteration_ranges[r][1][0]; j < cart2DDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
						{
							for (std::size_t i = cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]; i < cart2DDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
							{
								// reset time control
								shackTimestepControl->currentTimestepNr = 0;
								shackTimestepControl->currentSimulationTime = 0;

								std::cout << "." << std::flush;

								for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
									prog[inner_prog_id]->spectral_setZero();

								// activate mode via real coefficient
								prog[outer_prog_id]->spectral_set(j, i, 1.0);

								/*
								 * RUN timestep
								 */
								(i_class->*i_run_timestep_method)();


								if (shackPDESWECart2D->normal_mode_analysis_generation == 3)
								{
									/*
									 * compute
									 * 1/dt * (U(t+1) - U(t))
									 */
									///prog[outer_prog_id]->request_data_spectral();

									std::complex<double> val = prog[outer_prog_id]->spectral_get(j, i);
									val = val - 1.0;
									prog[outer_prog_id]->spectral_set(j, i, val);

									for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
										(*prog[inner_prog_id]) /= shackTimestepControl->currentTimestepSize;
								}


								for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
								{
									///prog[inner_prog_id]->request_data_spectral();

									/*
									 * REAL
									 */

									for (int r = 0; r < 2; r++)
									{
										for (std::size_t j = cart2DDataConfig->spectral_data_iteration_ranges[r][1][0]; j < cart2DDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
										{
											for (std::size_t i = cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]; i < cart2DDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
											{
												file << prog[inner_prog_id]->spectral_get(j, i).real();
												file << "\t";
											}
										}
									}
								}


								for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
								{
									/*
									 * IMAG
									 */
									int c = 0;
									for (int r = 0; r < 2; r++)
									{
										for (std::size_t j = cart2DDataConfig->spectral_data_iteration_ranges[r][1][0]; j < cart2DDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
										{
											for (std::size_t i = cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]; i < cart2DDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
											{
												file << prog[inner_prog_id]->spectral_get(j, i).imag();

												if (inner_prog_id != number_of_prognostic_variables-1 || c != specmodes-1)
													file << "\t";
												else
													file << std::endl;

												c++;
											}
										}
									}
								}
							}
						}
					}
#endif
				}
#else
				else if (shackPDESWECart2D->normal_mode_analysis_generation == 3 || shackPDESWECart2D->normal_mode_analysis_generation == 13)
				{
					sweet::Data::Cart2D::DataSpectralComplex t1(cart2DDataConfig);
					sweet::Data::Cart2D::DataSpectralComplex t2(cart2DDataConfig);
					sweet::Data::Cart2D::DataSpectralComplex t3(cart2DDataConfig);
					Cart2DDataComplex* prog_cplx[3] = {&t1, &t2, &t3};

					// iterate over spectral space
					for (std::size_t outer_i = 0; outer_i < cart2DDataConfig->spectral_complex_array_data_number_of_elements; outer_i++)
					{
						// reset time control
						shackTimestepControl->currentTimestepNr = 0;
						shackTimestepControl->currentSimulationTime = 0;

						std::cout << "." << std::flush;

						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
							prog_cplx[inner_prog_id]->spectral_setZero();

						// activate mode via real coefficient
						prog_cplx[outer_prog_id]->request_data_spectral();
						prog_cplx[outer_prog_id]->spectral_space_data[outer_i].real(1);

						// convert sweet::Data::Cart2D::DataSpectralComplex to Cart2DData
						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
						{
							*prog[inner_prog_id] = sweet::Data::Cart2D::Convert::Cart2DDataSpectralComplex_2_Cart2DDataSpectral::grid_convert(*prog_cplx[inner_prog_id]);
							prog[inner_prog_id]->spectral_zeroAliasingModes();
						}


						/*
						 * RUN timestep
						 */
						(i_class->*i_run_timestep_method)();

						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
						{
							prog[inner_prog_id]->spectral_zeroAliasingModes();
#warning "update this grid_convert maybe to spectral_convert"

							*prog_cplx[inner_prog_id] = Convert_Cart2DData_2_Cart2DDataComplex::grid_convert(*prog[inner_prog_id]);

							prog_cplx[inner_prog_id]->request_data_spectral();
						}

						if (shackPDESWECart2D->normal_mode_analysis_generation == 3)
						{
							/*
							 * compute
							 * 1/dt * (U(t+1) - U(t))
							 */
							prog_cplx[outer_prog_id]->request_data_spectral();
							prog_cplx[outer_prog_id]->spectral_space_data[outer_i] -= 1.0;

							for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
								prog_cplx[inner_prog_id]->operator*=(1.0/shackTimestepControl->currentTimestepSize);
						}


						// convert Cart2DData_SpectralComplex to Cart2DData
						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
						{
							prog_cplx[inner_prog_id]->request_data_spectral();

							/*
							 * REAL
							 */
							for (std::size_t k = 0; k < cart2DDataConfig->spectral_complex_array_data_number_of_elements; k++)
							{
								file << prog_cplx[inner_prog_id]->spectral_space_data[k].real();
								file << "\t";
							}

							/*
							 * IMAG
							 */
							for (std::size_t k = 0; k < cart2DDataConfig->spectral_complex_array_data_number_of_elements; k++)
							{
								file << prog_cplx[inner_prog_id]->spectral_space_data[k].imag();

								if (inner_prog_id != number_of_prognostic_variables-1 || k != cart2DDataConfig->spectral_complex_array_data_number_of_elements-1)
									file << "\t";
								else
									file << std::endl;
							}
						}
					}
				}
#endif
			}
		}
	}
};

}}

#endif
