/*
 * Author: Pedro Peixoto <ppeixoto@usp.br>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_BAROTROPIC_VORT_MODES_HPP
#define PROGRAMS_PDE_SWESPHERE2D_BENCHMARKS_BAROTROPIC_VORT_MODES_HPP

#include <algorithm>
#include <sweet/Data/Sphere2D/Sphere2D.hpp>

#include "HelperGeostropicBalance.hpp"
#include "BaseInterface.hpp"
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>

namespace PDE_SWESphere2D {
namespace Benchmarks {

class barotropic_vort_modes	:
		public BaseInterface
{
	HelperGeostropicBalance helperGeostropicBalance;

	// Mode setup
	std::size_t maxmodes;   //number of waves to be added
	static const int maxmodeslimit=100; //max number of waves
	std::size_t nmode[maxmodeslimit], mmode[maxmodeslimit];   // spherical harmonic indexes
	double ampl[maxmodeslimit]; //coefficients of normal modes 

	std::string benchmark_name;


public:
	barotropic_vort_modes()
	{
	}

	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)
	{
		BaseInterface::shackRegistration(io_shackDict);
		helperGeostropicBalance.shackRegistration(io_shackDict);
		return true;
	}

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		bool found_bench = false;
		benchmark_name = i_benchmark_name;
		if (benchmark_name.find("barotropic_vort_modes") !=std::string::npos )
			found_bench = true;

		return
			benchmark_name == "barotropic_vort_modes"			||
			found_bench;
	}


	void setup_1_shackData()
	{
		/*
		 * barotropic vorticity initialization:
		 * uses normal modes, that are spherical harmonic modes
		 *
		 * WARNING: It is designed for the barotropic vorticity equations, but can be used for SWE
		 * considering a balanced solution for the full non-linear equations
		 *
		 */
		if (shackParallelization->isMPIRoot)
		{
			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Amplitudes of normal modes will be divided by Earth Radius !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;
		}

		shackPDESWESphere2D->sphere2d_rotating_coriolis_omega = 7.292e-5;
		shackPDESWESphere2D->gravitation = 9.80616;
		shackSphere2DDataOps->sphere2d_radius = 6.37122e6;
		shackPDESWESphere2D->h0 = 29400.0/shackPDESWESphere2D->gravitation;
	}

	void setup_2_withOps(
			sweet::Data::Sphere2D::Operators *io_ops
	)
	{
		ops = io_ops;

		helperGeostropicBalance.setup(ops);

	}


	void clear()
	{
		helperGeostropicBalance.clear();
	}


	std::string getHelp()
	{
		std::ostringstream stream;
		stream << "  BOROTROPIC VORTICITY MODES :" << std::endl;
		stream << "     'barotropic_vort_modes'" << std::endl;
		stream << "     'barotropic_vort_modes': standard test" << std::endl;
		stream << "     'barotropic_vort_modes_[N]_[n1]_[m1]_[v1]_...[nN]_[mN]_[vN]': mode choices" << std::endl;
		return stream.str();
	}

	void getInitialState(
		sweet::Data::Sphere2D::DataSpectral &o_phi_pert,
		sweet::Data::Sphere2D::DataSpectral &o_vrt,
		sweet::Data::Sphere2D::DataSpectral &o_div
	)
	{
		double a = shackSphere2DDataOps->sphere2d_radius;
		//double omega = shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;
		//double u0 = 2.0*M_PI*a/(12.0*24.0*60.0*60.0);
		//double alpha = 0;


		extract_bench_info(benchmark_name);
	

		o_vrt.spectral_setZero();
		o_div.spectral_setZero();
		
		//sweet::Data::Sphere2D::DataGrid ug, vg;

		// m=-M,...M are Fourier modes
		// n=|m|, ..., M are Legendre nodes
		// M max modes
		//                 n,m  (n>=m)
		// Set modes for the vorticity stream function
		//Sphere2DData_Spectral psi(ops->sphere2DDataConfig); // = inv_laplace(i_vrt);
	
		//psi.spectral_setZero();

		//Add mode values
		for (int n = 0; n < (int)maxmodes; n++){			
			if (nmode[n] < mmode[n]){
				if (shackParallelization->isMPIRoot)
				{
					std::cout << "Modes: n=" <<nmode[n]<< " , m=" <<mmode[n]<<std::endl;
				}
				SWEETErrorFatal("SWESphere2DBenchmark_barotropic_vort_modes: n cannot be smaller than m");	
			}

			// Only real part of amplitude is considered, so ampl[] is a vector of real values
			o_vrt.spectral_set(nmode[n],mmode[n],ampl[n]/a);
		}
		

		helperGeostropicBalance.computeGeostrophicBalance_nonlinear(
				o_vrt,
				o_div,
				o_phi_pert
		);
	}

private:
	void extract_bench_info(const std::string &bcase)
	{
		
		std::string basic_name = "barotropic_vort_modes";
		std::string bcase_code = bcase;

		bcase_code.replace(benchmark_name.find(basic_name),basic_name.length()+1,"");

		if (bcase_code==""){
			SWEETErrorFatal("SWESphere2DBenchmark_barotropic_vort_modes: please choose the normal mode case appending to benchmark name the code _[N]_[n1]_[m1]_[v1]_...[nN]_[mN]_[vN]");
		};

		std::cout << "[MULE] benchmark_barotropic_vort_modes.case:" << benchmark_name << std::endl;

		//Convert parameter to words
		std::string str = bcase_code;
		std::replace( str.begin(), str.end(), '_', ' ');

		std::stringstream iss(str);
		
		iss >> maxmodes;
		if (maxmodes>maxmodeslimit){
			if (shackParallelization->isMPIRoot)
			{
				std::cout << "Modes:" <<maxmodes<< " , maxmodes hardcoded:" <<maxmodeslimit<<std::endl;
			}
			SWEETErrorFatal("SWESphere2DBenchmark_barotropic_vort_modes: Adjust maximun number of waves");	
		}
		std::cout << "[MULE] benchmark_barotropic_vort_modes.code: " <<bcase_code<< std::endl;
		std::cout << "[MULE] benchmark_barotropic_vort_modes.maxmodes: " << maxmodes << std::endl;
		
		//loop over waves
		for (int n = 0; n < (int)maxmodes; n++){
			//get a single mode
			iss >> nmode[n];
			iss >> mmode[n];
			iss >> ampl[n];
			if (shackParallelization->isMPIRoot)
			{
				std::cout << "[MULE] benchmark_barotropic_vort_modes." <<n<< ".nmode: " << nmode[n] << std::endl;
				std::cout << "[MULE] benchmark_barotropic_vort_modes." <<n<< ".mmode: " << mmode[n] << std::endl;
				std::cout << "[MULE] benchmark_barotropic_vort_modes." <<n<< ".ampl: " << ampl[n] << std::endl;
			}
		}
		
		return;

	}

	void setup_topography() override
	{
		// initialize the topography
		shackPDESWEBenchmarks->topography.setup(ops->sphere2DDataConfig);
	}


};

}}

#endif
