/*
 * PDESWECart2DBench_NormalModes.hpp
 *
 *  Created on: 03 Nov 2019
 * Author: Pedro Peixoto <pedrosp@ime.usp.br>
 */
#ifndef PROGRAMS_PDE_SWECART2D_BENCHMARKS_PDESWECART2DBENCH_NORMALMODES_HPP
#define PROGRAMS_PDE_SWECART2D_BENCHMARKS_PDESWECART2DBENCH_NORMALMODES_HPP

#include <programs/PDE_SWECart2D/Benchmarks/BaseInterface.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2D/Operators.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

#include <cmath>
#include <sweet/ExpIntegration/ExpFunction.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <functional>

/**
 * SWE Cart2D normal mode benchmark
 */

		

namespace PDE_SWECart2D {
namespace Benchmarks {

class normal_modes	:
		public PDE_SWECart2D::Benchmarks::BaseInterface
{
public:
	std::string bcasename; //Benchmark case name
	std::size_t nwaves;   //number of waves to be added
	static const int maxwaves=10; //max number of waves
	std::size_t k0[maxwaves], k1[maxwaves];   // wavenumber to be used (-1,-1) refers to all wavenumbers
	double d0[maxwaves], dwest[maxwaves], deast[maxwaves]; //coefficients of normal modes eigen vectors

/**
 * SWE Cart2D normal mode benchmark helper functions
 */

#if SWEET_QUADMATH && 0
	typedef __float128 T;
#else
	typedef double T;
#endif
	typedef std::complex<T> complex;

public:
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

		////////////////io_h.request_data_spectral();
		////////////////io_u.request_data_spectral();
		////////////////io_v.request_data_spectral();

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
		//std::cout << "Eigen-values" <<std::endl;
		//std::cout <<lambda[0]<<std::endl;
		//std::cout <<lambda[1]<<std::endl;
		//std::cout <<lambda[2]<<std::endl;

		complex U[3];
		// Set normal mode acording to desired wave type
		// These are weights for the modes
		
		U[0] = geo_mode;
		U[1] = igwest_mode;
		U[2] = igeast_mode;
		//std::cout <<U[0]<<std::endl;
		//std::cout <<U[1]<<std::endl;
		//std::cout <<U[2]<<std::endl;
		//Define normal mode as combination of eigen vectors
		complex UEV[3] = {0.0, 0.0, 0.0};
		for (int k = 0; k < 3; k++)
			for (int j = 0; j < 3; j++)
				UEV[k] += v[k][j] * U[j];

		//std::cout << "spectral before" << std::endl;
		//io_v.print_spectralIndex();

		//Add normal mode to spectral space of variables
		// The scale factor is required to ensure the normal mode is dimensionally of order 1 (normalized)
		//     since in the backwards transform a scale factor is also applied.
		complex h_add, u_add, v_add;
		double scale_factor = ((double)(cart2DDataConfig->grid_data_size[0]*cart2DDataConfig->grid_data_size[1]));
		h_add = io_h.spectral_get(ik1, ik0)+UEV[0]*scale_factor;
		u_add = io_u.spectral_get(ik1, ik0)+UEV[1]*scale_factor;
		v_add = io_v.spectral_get(ik1, ik0)+UEV[2]*scale_factor;

		/* Add normal mode to data */
		io_h.spectral_set(ik1, ik0, h_add);
		io_u.spectral_set(ik1, ik0, u_add);
		io_v.spectral_set(ik1, ik0, v_add);

		io_h.spectral_zeroAliasingModes();
		io_u.spectral_zeroAliasingModes();
		io_v.spectral_zeroAliasingModes();

		//Request physical data, to ensure that mirroring is well well balanced (it may fill in the mirror mode)
		//////////////io_h.request_data_physical();
		//////////////io_u.request_data_physical();
		//////////////io_v.request_data_physical();

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
		
		std::cout << "Adding normal mode to wavenumber : (" <<ik0<< "," <<ik1<< ")" <<std::endl;
		std::cout << "h: " << UEV[0]<<std::endl; 
		std::cout << "u: " << UEV[1]<<std::endl; 
		std::cout << "v: " << UEV[2]<<std::endl; 
		

	#endif
		return;
	}


public:
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
		/////o_geo_mode.request_data_physical();	
		/////o_igwest_mode.request_data_physical();	
		/////o_igeast_mode.request_data_physical();	
		return;
	}


public:
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

		//////i_h.request_data_spectral();
		//////i_u.request_data_spectral();
		//////i_v.request_data_spectral();

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

public:
	void sw_eigen_decomp(
			T k0,				//wavenumber in x
			T k1,				// wavenumeber in y4
			bool i_inverse,		// Input true, returns inverse matriz, false: returns direct
			complex o_v[3][3],	// output eigen vector (direct or inverse)
			complex o_evalues[3]	// output eigen values (optional)
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

		T f = shackPDESWECart2D->cart2d_rotating_f0;
		T h = shackPDESWECart2D->h0;
		T g = shackPDESWECart2D->gravitation;

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

		
		if (shackPDESWECart2D->cart2d_rotating_f0 == 0)
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
		else
		{
			if (k0 == 0 && k1 == 0)
			{
				/*
					* http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,0,0%7D,%7B0,0,f%7D,%7B0,-f,0%7D%7D
					* Order need to be changed
					*/
				//west
				v[0][1] = 0;
				v[1][1] = -I;
				v[2][1] = 1;
				//east
				v[0][2] = 0;
				v[1][2] = I;
				v[2][2] = 1;
				//geost
				v[0][0] = 1;
				v[1][0] = 0;
				v[2][0] = 0;

				if (i_evalues){
					lambda[2] = I*f;
					lambda[1] = -I*f;
					lambda[0] = 0;
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
		std::cout << "EV matrix" <<std::endl;
		for (int j = 0; j < 3; j++)	{
			for (int i = 0; i < 3; i++)
				std::cout << v[j][i]<< " "; 
			std::cout <<std::endl;
		}
		*/

		//Normalize EV matrix (each EV has norm 1)
		complex evnorm[3];
		for (int j = 0; j < 3; j++)	{
			//Calculate EV norms
			evnorm[j]=0;
			for (int i = 0; i < 3; i++){
					//std::cout << j << ":" << evnorm[j] << std::endl;
					evnorm[j] += v[i][j]*std::conj(v[i][j]);
			}
		}
		for (int j = 0; j < 3; j++)	{
			for (int i = 0; i < 3; i++){
				//std::cout << k0 << " , " << k1 << " , " << i << " , " << j << " , " << v[i][j] <<std::endl;
				v[i][j] /= std::__complex_sqrt(evnorm[j]) ;
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

public:
	std::stringstream dump_normal_modes(
			sweet::Data::Cart2D::DataSpectral &i_mode    //Coeficients multiplying  mode
	)
	{
		const sweet::Data::Cart2D::Config *cart2DDataConfig = i_mode.cart2DDataConfig;

		std::stringstream buffer;
		buffer << std::setprecision(8);

		buffer << shackTimestepControl->currentTimestepNr;
		buffer << "\t" << shackTimestepControl->currentSimulationTime;
		const double zero=0.0;
		const double scale_factor =((double)(cart2DDataConfig->spectral_array_data_number_of_elements)); 
		//cart2DDataConfig->grid_data_size[0]*cart2DDataConfig->grid_data_size[1]));
		for (std::size_t ik1 = 0; ik1 < cart2DDataConfig->spectral_data_size[1]; ik1++)
		{
			for (std::size_t ik0 = 0; ik0 < cart2DDataConfig->spectral_data_size[0]; ik0++)
			{
				const std::complex<double> &value = i_mode.spectral_get(ik1, ik0);
				double norm = value.real()*value.real()+value.imag()*value.imag();
				norm=std::sqrt(norm/scale_factor);
				if (norm > 1.0e-15)
				{
					buffer << "\t" << norm;
				}
				else
				{
					buffer << "\t" << zero;
				}
			}
		}	
		return buffer;
	}

public:
	void dump_all_normal_modes(
			sweet::Data::Cart2D::DataSpectral &i_mode_geo,    //Coeficients multiplying  mode
			sweet::Data::Cart2D::DataSpectral &i_mode_igwest,    //Coeficients multiplying  mode
			sweet::Data::Cart2D::DataSpectral &i_mode_igeast    //Coeficients multiplying  mode
	)
	{
		const char* filename_geo = "output_nm_geo_evol.txt";
		const char* filename_igwest = "output_nm_igwest_evol.txt";
		const char* filename_igeast = "output_nm_igeast_evol.txt";

		const sweet::Data::Cart2D::Config *cart2DDataConfig = i_mode_geo.cart2DDataConfig;

		std::ofstream file0;	
		std::ofstream file1;	
		std::ofstream file2;	
		if (shackTimestepControl->currentTimestepNr == 0){
			file0.open(filename_geo, std::ofstream::out | std::ofstream::trunc);	
			file1.open(filename_igwest, std::ofstream::out | std::ofstream::trunc);	
			file2.open(filename_igeast, std::ofstream::out | std::ofstream::trunc);	
		}
		else
		{
			file0.open(filename_geo, std::ofstream::out | std::ofstream::app);	
			file1.open(filename_igwest, std::ofstream::out | std::ofstream::app);	
			file2.open(filename_igeast, std::ofstream::out | std::ofstream::app);	
		}

		std::stringstream buffer0, buffer1, buffer2;
		buffer0 << std::setprecision(8);
		buffer1 << std::setprecision(8);
		buffer2 << std::setprecision(8);

		//Headers
		if (shackTimestepControl->currentTimestepNr == 0){
			//header
			buffer0 << "n\t time";
			buffer1 << "n\t time";
			buffer2 << "n\t time";
			T k0, k1;
			for (std::size_t ik1 = 0; ik1 < cart2DDataConfig->spectral_data_size[1]; ik1++)
			{
				for (std::size_t ik0 = 0; ik0 < cart2DDataConfig->spectral_data_size[0]; ik0++)
				{
					k0 = (T) ik0;
					if (ik1 <= cart2DDataConfig->spectral_data_size[1]/2)
						k1 = (T)ik1;
					else
						k1 = (T)((int)ik1-(int)cart2DDataConfig->spectral_data_size[1]);

					//Geo modes
					buffer0 << "\t(" << k0 << "," <<k1<< ")";
					//West modes (adjust y_mode index)
					buffer1 << "\t(" << k0 << "," <<k1<< ")";
					//East mode (adjust x_mode index)
					buffer2 << "\t(" << (k0 == 0 ? std::abs(k0) : -k0)  << "," << (k1 == 0 ? std::abs(k1) : -k1) << ")";

				}
			}
			file0 << buffer0.str() << std::endl;
			file1 << buffer1.str() << std::endl;
			file2 << buffer2.str() << std::endl;
		}

		buffer0 = dump_normal_modes(i_mode_geo);
		file0 << buffer0.str() << std::endl;
		buffer0.str(std::string());

		buffer1 = dump_normal_modes(i_mode_igwest);
		file1 << buffer1.str() << std::endl;
		buffer1.str(std::string());

		buffer2 = dump_normal_modes(i_mode_igeast);
		file2 << buffer2.str() << std::endl;
		buffer2.str(std::string());

		file0.close();
		file1.close();
		file2.close();

		return ;
	}

private:
	void extract_bench_info(const std::string &bcase)
	{
		
		if (bcase==""){
			SWEETErrorFatal("PDESWECart2DBench_NormalModes: please choose the normal mode case with --benchmark-normal-mode-case [string] (see PDESWECart2DBench_NormalModes.hpp file)");
		};
		std::cout << bcase <<std::endl;

		//Convert parameter to words
		std::string str = bcase;
		std::replace( str.begin(), str.end(), '_', ' ');
		//std::cout << str<<std::endl;

		std::stringstream iss(str);
		iss >> bcasename;
		std::cout << "[MULE] benchmark_normal_modes.case:" << bcasename << std::endl;
		if (bcasename=="waves"){
			iss >> nwaves;
			if (nwaves>maxwaves){
				std::cout << "Waves:" <<nwaves<< " , maxwaves hardcoded:" <<maxwaves<<std::endl;
				SWEETErrorFatal("PDESWECart2DBench_NormalModes: Adjust maximun number of waves");
			}
			std::cout << "[MULE] simulation_benchmark_normal_modes.case: " <<bcasename<< std::endl;
			std::cout << "[MULE] simulation_benchmark_normal_modes.nwaves: " << nwaves << std::endl;
			
			//loop over waves
			for (int n = 0; n < (int)nwaves; n++){
				//get a single wave
				iss >> k0[n];
				iss >> k1[n];
				iss >> d0[n];
				iss >> dwest[n];
				iss >> deast[n];
				std::cout << "[MULE] simulation_benchmark_normal_modes.w" <<n<< ".k0: " << k0[n] << std::endl;
				std::cout << "[MULE] simulation_benchmark_normal_modes.w" <<n<< ".k1: " << k1[n] << std::endl;
				std::cout << "[MULE] simulation_benchmark_normal_modes.w" <<n<< ".d0: " << d0[n] << std::endl;
				std::cout << "[MULE] simulation_benchmark_normal_modes.w" <<n<< ".dwest: " << dwest[n] << std::endl;
				std::cout << "[MULE] simulation_benchmark_normal_modes.w" <<n<< ".deast: " << deast[n] << std::endl;
				//SWEETErrorFatal("PDESWECart2DBench_NormalModes: Adjust maximun number of waves");
			}
		}
		else{
			SWEETErrorFatal("PDESWECart2DBench_NormalModes: Please follow naming convention for nomal mode initialization: waves_N_k0_k1_d0_deast_dwest_k0_k1_d0_deast_dwest_k0_k1_d0_deast_dwest");
		}
		return;

	}

public:
	bool setupBenchmark(
			sweet::Data::Cart2D::DataSpectral &o_h,
			sweet::Data::Cart2D::DataSpectral &o_u,
			sweet::Data::Cart2D::DataSpectral &o_v
	)
	{
		std::cout << "Generating Normal Modes Initial Conditions: ";

		extract_bench_info(shackBenchmarks->benchmark_normal_modes_case);
		
		//Naming convention for normal_mode_cases (for arbitraty use):
		//  name_k0_k1_d0_deast_dwest
		// (k0,k1) are wave numbers (set to -1 for all wavenumbers)
		// d0, dwest, deast are numbers (floats) that are coefficients for different normal wave types

		// General convention (for N waves, repeat pattern N times)
		// waves_N_k0_k1_d0_deast_dwest_k0_k1_d0_deast_dwest_k0_k1_d0_deast_dwest

		//zero initial conditions
		const sweet::Data::Cart2D::Config *cart2DDataConfig = o_h.cart2DDataConfig;

		o_h.spectral_setZero();
		o_u.spectral_setZero();
		o_v.spectral_setZero();

		if (bcasename=="waves")
		{
			for (int n = 0; n<(int)nwaves; n++){
				//Set a single wavenumber with appropriate modes
				if (k0[n]<0 && k1[n] <0 ){
					std::cout << "Adding normal modes to all wavenumbers" <<std::endl;
					for (std::size_t ik1 = 0; ik1 < cart2DDataConfig->spectral_data_size[1]; ik1++)
					{
						for (std::size_t ik0 = 0; ik0 < cart2DDataConfig->spectral_data_size[0]; ik0++)
						{
							add_normal_mode(
										ik0, ik1,
										d0[n], dwest[n], deast[n],
										o_h, o_u, o_v
									);
						}
					}
				}
				else if (k0[n]>=0 && k1[n] >=0 ){
						add_normal_mode(
									k0[n], k1[n],
									d0[n], dwest[n], deast[n],
									o_h, o_u, o_v
							);
				}
				else
				{
					SWEETErrorFatal("PDESWECart2DBench_NormalModes: invalid wavenumber selection in --benchmark-normal-mode-case [string] (see PDESWECart2DBench_NormalModes.hpp file)");
				
				}
			}
		};
		std::cout << "[MULE] benchmark_normal_modes.h_max:" <<o_h.spectral_reduce_max_abs()<<std::endl;
		std::cout << "[MULE] benchmark_normal_modes.u_max:" <<o_u.spectral_reduce_max_abs()<<std::endl;
		std::cout << "[MULE] benchmark_normal_modes.v_max:" <<o_v.spectral_reduce_max_abs()<<std::endl;
		
		std::cout << "Normal mode initial conditions generated successfully! " << std::endl;

		return true;
	}

};


}}

#endif
