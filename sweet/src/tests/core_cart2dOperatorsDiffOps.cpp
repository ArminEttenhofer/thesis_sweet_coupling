/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_SCONS_OPTIONS: --cart2d-spectral-space=enable
 */

#include <sweet/Tools/DefaultPrecompilerValues.hpp>

#if SWEET_GUI
#	error	"GUI not supported"
#endif


#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/Tools/ProgramArguments.hpp>

#include <ostream>
#include <cmath>


int main(
		int i_argc,
		char *i_argv[]
)
{
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	sweet::Data::Cart2D::Shack *shackCart2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Cart2D::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.printShackData();

	double prev_error_diff_x = 0;
	double prev_error_diff_y = 0;
	double prev_error_diff_z = 0;
	
	double prev_error_diff2_x = 0;
	double prev_error_diff2_y = 0;

	double freq_x = 2.0;
	double freq_y = 2.0;

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	std::size_t res_x = shackCart2DDataOps->space_res_physical[0];
	std::size_t res_y = shackCart2DDataOps->space_res_physical[1];

	std::size_t max_res = 2048;

	if (res_x > max_res || res_y > max_res)
		max_res = std::max(res_x, res_y);

	for (; res_x <= max_res && res_y <= max_res; res_x *= 2, res_y *= 2)
	{
		/*
		 * upper bound for roundoff errors.
		 *
		 * The upper bound of roundoff errors for a 1D FFT is given by sqrt(N) using
		 * an RMS norm, see "Roundoff Error Analysis of the Fast Fourier Transform" by George U. Ramos
		 *
		 * Let
		 * d: truncated roundoff errors
		 * T: matrix for Fourier transformation (not the FT)
		 *
		 * Then, according to the work mentioned above
		 * |T d|_rms <= sqrt(N) |d|rms
		 *
		 * Rearranging to Euclidian norm, this yields
		 * sqrt(1/N) |T d|_2 <= sqrt(N) * sqrt(1/N) |d|_2
		 *
		 * Cancelling the ugly sqrt terms, we get
		 * |T d|_2 <= sqrt(N) |d|_2
		 *
		 * So the max. error is increasing linearly to sqrt(N) with an increasing
		 * problem size N to be transformed via FT transformed.
		 *
		 * For a 2D problem, we assume the error thresholds to be additive correlated
		 * due to the 2D problem expressed by two successively 1D FFTs, hence the
		 * upper error thresholds should be related in an additive way.
		 *
		 * We can assume, that the error for a 2D problem is increasing linearly with sqrt(N):
		 *
		 * err ~ sqrt(Nx) + sqrt(Ny)
		 */
		double tolerance_increase = sqrt(res_x) + sqrt(res_y);

		/*
		 * error tolerance for machine accuracy
		 *
		 * We assume 1e-12 for double precision
		 */
		double eps = 1e-9*tolerance_increase;

		/*
		 * error tolerance for convergence
		 *
		 * Here, we are very patronizing due to flickering convergence for coarse solutions which
		 * are not really representable in the Fourier space where the discretization errors
		 * are dominating.
		 */
		double eps_convergence = 1e-3;


		std::cout << "*************************************************************" << std::endl;
		std::cout << "Testing operators with resolution " << res_x << " x " << res_y << std::endl;
		std::cout << "*************************************************************" << std::endl;
		std::size_t res[2] = {res_x, res_y};

		shackCart2DDataOps->space_res_physical[0] = res[0];
		shackCart2DDataOps->space_res_physical[1] = res[1];
		shackCart2DDataOps->space_res_spectral[0] = -1;
		shackCart2DDataOps->space_res_spectral[1] = -1;

		sweet::Data::Cart2D::Config cart2DDataConfig;
		cart2DDataConfig.setupAuto(shackCart2DDataOps);

		sweet::Data::Cart2D::Operators ops;
		ops.setup(cart2DDataConfig, shackCart2DDataOps);
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(ops);

		/*
		 * keep h in the outer regions to allocate it only once and avoid reinitialization of FFTW
		 */
		sweet::Data::Cart2D::DataGrid h(cart2DDataConfig);

		{
			std::cout << "**********************************************" << std::endl;
			std::cout << "> Resolution (" << res_x << "x" << res_y << ")" << std::endl;
			std::cout << "> Domain size (" << shackCart2DDataOps->cart2d_domain_size[0] << "x" << shackCart2DDataOps->cart2d_domain_size[1] << ")" << std::endl;
			std::cout << "**********************************************" << std::endl;
			std::cout << "error tol = " << eps << std::endl;
			std::cout << "**********************************************" << std::endl;
		}

		/**
		 * Tests for basic operators which are not amplifying the solution depending on the domain size
		 */
		{
			sweet::Data::Cart2D::DataGrid u(cart2DDataConfig);
			sweet::Data::Cart2D::DataGrid v(cart2DDataConfig);



#if 0
			std::cout << std::endl;
			std::cout << "op.diff_c_x" << std::endl;
			op.diff_c_x.print_spectralData_zeroNumZero();
			std::cout << std::endl;
			std::cout << "op.diff_c_y" << std::endl;
			op.diff_c_y.print_spectralData_zeroNumZero();

			std::cout << std::endl;
			std::cout << "op.diff2_c_x" << std::endl;
			op.diff2_c_x.print_spectralData_zeroNumZero();
			std::cout << std::endl;
			std::cout << "op.diff2_c_y" << std::endl;
			op.diff2_c_y.print_spectralData_zeroNumZero();
			exit(1);
#endif

			for (int j = 0; j < shackCart2DDataOps->space_res_physical[1]; j++)
			{
				for (int i = 0; i < shackCart2DDataOps->space_res_physical[0]; i++)
				{
					double x = ((double)i+0.5)/(double)shackCart2DDataOps->space_res_physical[0];
					double y = ((double)j+0.5)/(double)shackCart2DDataOps->space_res_physical[1];

#define FUN_ID	1

	#if FUN_ID==1
					u.grid_setValue(j, i, sin(freq_x*M_PI*x));
					v.grid_setValue(j, i, cos(freq_y*M_PI*y));
	#elif FUN_ID==2
					u.grid_setValue(j, i, sin(freq_x*M_PI*x));
					v.grid_setValue(j, i, 1.0/(cos(freq_y*M_PI*y)+2.0));
	#endif

					h.grid_setValue(
						j, i,
	#if FUN_ID==1
						sin(freq_x*M_PI*x)*cos(freq_y*M_PI*y)
	#elif FUN_ID==2
						sin(freq_x*M_PI*x)*sin(freq_x*M_PI*x)*cos(freq_y*M_PI*y)*cos(freq_y*M_PI*y)
	#elif FUN_ID==3
						sin(freq_x*M_PI*x)/(cos(freq_y*M_PI*y)+2.0)
	#endif
					);
				}
			}

			//sweet::Data::Cart2D::DataSpectral u_spec(u);
			//sweet::Data::Cart2D::DataSpectral v_spec(v);
			//double err_z = (u_spec*v_spec-h).toGrid().grid_reduce_rms_quad();
			double err_z = (u*v-h).grid_reduce_rms_quad();

#if FUN_ID == 1
			if (err_z > eps)
			{
				std::cerr << "SPEC: Error threshold exceeded for err_z!" << std::endl;
				exit(-1);
			}
#endif

			prev_error_diff_z = err_z;
		}

		std::cout << "TEST B: DONE" << std::endl;


		/**
		 * Tests for 1st order differential operator
		 */
		{
			sweet::Data::Cart2D::DataGrid u(cart2DDataConfig);
			sweet::Data::Cart2D::DataGrid v(cart2DDataConfig);
			sweet::Data::Cart2D::DataGrid h_diff_x(cart2DDataConfig);
			sweet::Data::Cart2D::DataGrid h_diff_y(cart2DDataConfig);

			for (int j = 0; j < shackCart2DDataOps->space_res_physical[1]; j++)
			{
				for (int i = 0; i < shackCart2DDataOps->space_res_physical[0]; i++)
				{
					double x = ((double)i+0.5)/(double)shackCart2DDataOps->space_res_physical[0];
					double y = ((double)j+0.5)/(double)shackCart2DDataOps->space_res_physical[1];

	#if FUN_ID==1
					u.grid_setValue(j, i, sin(freq_x*M_PI*x));
					v.grid_setValue(j, i, cos(freq_y*M_PI*y));
	#elif FUN_ID==2
					u.grid_setValue(j, i, sin(freq_x*M_PI*x));
					v.grid_setValue(j, i, 1.0/(cos(freq_y*M_PI*y)+2.0));
	#endif

					h.grid_setValue(
						j, i,
	#if FUN_ID==1
						sin(freq_x*M_PI*x)*cos(freq_y*M_PI*y)
	#elif FUN_ID==2
						sin(freq_x*M_PI*x)*sin(freq_x*M_PI*x)*cos(freq_y*M_PI*y)*cos(freq_y*M_PI*y)
	#elif FUN_ID==3
						sin(freq_x*M_PI*x)/(cos(freq_y*M_PI*y)+2.0)
	#endif
					);

					h_diff_x.grid_setValue(
						j, i,
	#if FUN_ID==1
						freq_x*M_PI*cos(freq_x*M_PI*x)*cos(freq_y*M_PI*y)/(double)shackCart2DDataOps->cart2d_domain_size[0]
	#elif FUN_ID==2
						2.0*sin(freq_x*M_PI*x)*std::pow(cos(freq_y*M_PI*y),2.0)*freq_x*M_PI*cos(freq_x*M_PI*x)/(double)shackCart2DDataOps->cart2d_domain_size[0]
	#elif FUN_ID==3
						freq_x*M_PI*cos(freq_x*M_PI*x)/(cos(freq_y*M_PI*y)+2.0)/(double)shackCart2DDataOps->cart2d_domain_size[0]
	#endif
					);

					h_diff_y.grid_setValue(
						j, i,
	#if FUN_ID==1
						-sin(freq_x*M_PI*x)*freq_y*M_PI*sin(freq_y*M_PI*y)/(double)shackCart2DDataOps->cart2d_domain_size[1]
	#elif FUN_ID==2
						-2.0*std::pow(std::sin(freq_x*M_PI*x),2.0)*std::cos(freq_y*M_PI*y)*freq_y*M_PI*std::sin(freq_y*M_PI*y)/(double)shackCart2DDataOps->cart2d_domain_size[1]
	#elif FUN_ID==3
						sin(freq_x*M_PI*x)*freq_y*M_PI*sin(freq_y*M_PI*y)/pow(cos(freq_y*M_PI*y)+2.0, 2.0)/(double)shackCart2DDataOps->cart2d_domain_size[1]
	#endif
					);
				}
			}


			double res_normalization = sqrt(1.0/(shackCart2DDataOps->space_res_physical[0]*shackCart2DDataOps->space_res_physical[1]));

			sweet::Data::Cart2D::DataSpectral h_spec(h);
//			h_spec.loadCart2DDataGrid(h);

			// normalization for diff = 2 pi / L
			double err_x = (ops.diff_c_x(h_spec)-h_diff_x).toGrid().grid_reduce_norm2()*res_normalization*shackCart2DDataOps->cart2d_domain_size[0]/(2.0*M_PI);
			double err_y = (ops.diff_c_y(h_spec)-h_diff_y).toGrid().grid_reduce_norm2()*res_normalization*shackCart2DDataOps->cart2d_domain_size[1]/(2.0*M_PI);
//			double err_z = (u*v-h).reduce_norm2()*res_normalization;

			{
				double conv_x = prev_error_diff_x/err_x;
				double conv_y = prev_error_diff_y/err_y;
//				double conv_z = prev_error_diff_z/err_z;
				std::cout << "error diff x = " << err_x << std::endl;
				std::cout << "error diff y = " << err_y << std::endl;
//				std::cout << "error diff z = " << err_z << std::endl;
				std::cout << "conv x = " << conv_x << std::endl;
				std::cout << "conv y = " << conv_y << std::endl;
//				std::cout << "conv z = " << conv_z << std::endl;

				if (conv_x != 0)
				if (abs(conv_x-4.0) > eps_convergence)
				{
					std::cerr << "Cart: Error threshold exceeded for conv_x, no convergence given!" << std::endl;
					exit(-1);
				}

				if (conv_y != 0)
				if (abs(conv_y-4.0) > eps_convergence)
				{
					std::cerr << "Cart: Error threshold exceeded for conv_y, no convergence given!" << std::endl;
					exit(-1);
				}
#if 0
				if (abs(err_z) > eps)
				{
					std::cerr << "Cart: Error threshold exceeded for err_z!" << std::endl;
					exit(-1);
				}
#endif
				prev_error_diff_x = err_x;
				prev_error_diff_y = err_y;
			}
		}

		std::cout << "TEST C: DONE" << std::endl;



		/**
		 * 2nd order differential operator
		 *
		 * note, that the function on which the 2nd diff operator is computed on has
		 * to be scaled up be a factor of domain_size^2, since e.g.
		 *
		 * diff(sin(2 pi x / size), x, x) = 4.0 pi^2 sin(2 pi x / size) / size^2
		 */
		{
			sweet::Data::Cart2D::DataGrid h_diff2_x(cart2DDataConfig);
			sweet::Data::Cart2D::DataGrid h_diff2_y(cart2DDataConfig);

			for (int j = 0; j < shackCart2DDataOps->space_res_physical[1]; j++)
			{
				for (int i = 0; i < shackCart2DDataOps->space_res_physical[0]; i++)
				{
					double x = ((double)i+0.5)/(double)shackCart2DDataOps->space_res_physical[0];
					double y = ((double)j+0.5)/(double)shackCart2DDataOps->space_res_physical[1];

					h.grid_setValue(
						j, i,
	#if FUN_ID==1
						sin(freq_x*M_PI*x)*cos(freq_y*M_PI*y)
	#elif FUN_ID==2
						sin(freq_x*M_PI*x)*sin(freq_x*M_PI*x)*cos(freq_y*M_PI*y)*cos(freq_y*M_PI*y)
	#elif FUN_ID==3
						sin(freq_x*M_PI*x)/(cos(freq_y*M_PI*y)+2.0)
	#endif
					);

					h_diff2_x.grid_setValue(
						j, i,
	#if FUN_ID==1
						freq_x*freq_x*M_PI*M_PI*(-1.0)*sin(freq_x*M_PI*x)*cos(freq_y*M_PI*y)/(shackCart2DDataOps->cart2d_domain_size[0]*shackCart2DDataOps->cart2d_domain_size[0])
	#elif FUN_ID==2
	//					2.0*sin(freq_x*M_PI*x)*std::pow(cos(freq_y*M_PI*y),2.0)*freq_x*M_PI*cos(freq_x*M_PI*x)/(double)parameters.sim.domain_size[0]
	#elif FUN_ID==3
	//					freq_x*M_PI*cos(freq_x*M_PI*x)/(cos(freq_y*M_PI*y)+2.0)/(double)parameters.sim.domain_size[0]
	#endif
					);

					h_diff2_y.grid_setValue(
						j, i,
	#if FUN_ID==1
						-sin(freq_x*M_PI*x)*freq_y*M_PI*freq_y*M_PI*cos(freq_y*M_PI*y)/(shackCart2DDataOps->cart2d_domain_size[1]*shackCart2DDataOps->cart2d_domain_size[1])
	#elif FUN_ID==2
	//					-2.0*std::pow(std::sin(freq_x*M_PI*x),2.0)*std::cos(freq_y*M_PI*y)*freq_y*M_PI*std::sin(freq_y*M_PI*y)/(double)parameters.sim.domain_size[1]
	#elif FUN_ID==3
	//					sin(freq_x*M_PI*x)*freq_y*M_PI*sin(freq_y*M_PI*y)/pow(cos(freq_y*M_PI*y)+2.0, 2.0)/(double)parameters.sim.domain_size[1]
	#endif
					);
				}
			}

			double normalization = sqrt(1.0/(shackCart2DDataOps->space_res_physical[0]*shackCart2DDataOps->space_res_physical[1]));

			sweet::Data::Cart2D::DataSpectral h_spec(h.cart2DDataConfig);
			h_spec.loadCart2DDataGrid(h);

			// diff2 normalization = 4.0 pi^2 / L^2
			double err2_x = (ops.diff2_c_x(h_spec) - h_diff2_x).toGrid().grid_reduce_norm2_quad()*normalization*(shackCart2DDataOps->cart2d_domain_size[0]*shackCart2DDataOps->cart2d_domain_size[0])/(4.0*M_PI*M_PI);
			double err2_y = (ops.diff2_c_y(h_spec) - h_diff2_y).toGrid().grid_reduce_norm2_quad()*normalization*(shackCart2DDataOps->cart2d_domain_size[1]*shackCart2DDataOps->cart2d_domain_size[1])/(4.0*M_PI*M_PI);

			{
				double conv2_x = prev_error_diff2_x/err2_x;
				double conv2_y = prev_error_diff2_y/err2_y;
				std::cout << "error diff2 x = " << err2_x << std::endl;
				std::cout << "error diff2 y = " << err2_y << std::endl;
				std::cout << "conv2 x = " << conv2_x << std::endl;
				std::cout << "conv2 y = " << conv2_y << std::endl;

				if (conv2_x != 0)
				if (abs(conv2_x-4.0) > eps_convergence)
				{
					std::cerr << "Cart: Error threshold exceeded for conv2_x, no convergence given!" << std::endl;
					exit(-1);
				}

				if (conv2_y != 0)
				if (abs(conv2_y-4.0) > eps_convergence)
				{
					std::cerr << "Cart: Error threshold exceeded for conv2_y, no convergence given!" << std::endl;
					exit(-1);
				}

				prev_error_diff2_x = err2_x;
				prev_error_diff2_y = err2_y;
			}
		}

		std::cout << "TEST D: DONE" << std::endl;
	}

	std::cout << "SUCCESSFULLY FINISHED" << std::endl;

	return 0;
}
