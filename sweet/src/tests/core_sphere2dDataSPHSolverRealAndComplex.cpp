/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * Include these files and directory for compilation
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWESphere2D/TimeOld/PDESWESphere2DTS_l_erk.cpp
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWESphere2D/TimeOld/PDESWESphere2DTS_lg_erk.cpp
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWESphere2D/TimeOld/PDESWESphere2DTS_l_irk.cpp
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWESphere2D/BenchmarksCombined.cpp
 *
 * MULE_SCONS_OPTIONS: --fortran-source=enable
 */

#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Sphere2DComplex_DataSpectral.hpp>
#include <sweet/Data/Sphere2D/DataSpectral.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/Data/Sphere2D/Shack.hpp>
#include <sweet/Data/Sphere2DComplex/Convert/DataSpectral_2_Sphere2D_DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/Operators.hpp>
#include <cmath>

#include "../programs/PDE_SWESphere2D/BenchmarksCombined.hpp"

#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>

#include "../programs/PDE_SWESphere2D/TimeOld/PDESWESphere2DTS_l_erk.hpp"
#include "../programs/PDE_SWESphere2D/TimeOld/PDESWESphere2DTS_l_irk.hpp"

#include "../programs/PDE_SWESphere2D/Shack.hpp"
#include "../programs/PDE_SWESphere2D/TimeHelpers/SphBandedMatrix_GridComplex.hpp"
#include "../programs/PDE_SWESphere2D/TimeHelpers/SphBandedMatrix_GridReal.hpp"


template<
	typename phys_value_type = double,
	typename sphere2d_data_spec_type = sweet::Data::Sphere2D::DataSpectral,
	typename sphere2d_data_phys_type = sweet::Data::Sphere2D::DataGrid,
	typename sphere2d_operators_type = sweet::Data::Sphere2D::Operators,
	typename sph_banded_solver_type = SphBandedMatrix_GridReal
>
class Test
{
	double double_precision_digits = -1.0;

	sweet::Data::Sphere2D::Operators *opsReal;
	sphere2d_operators_type *ops;

	sweet::Shacks::Dictionary *shackDict;
	PDE_SWESphere2D::Shack *shackPDESWESphere2D;
	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps;

	PDE_SWESphere2D::Benchmarks::BenchmarksCombined benchmarks;

public:
	sweet::Error::Base error;

	Test()	:
		ops(nullptr)
	{
	}

	~Test()
	{
		delete ops;
	}

public:
	bool shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;
		shackPDESWESphere2D = shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
		shackSphere2DDataOps = shackDict->getAutoRegistration<sweet::Data::Sphere2D::Shack>();

		benchmarks.setup_1_registerAllBenchmark();
		benchmarks.setup_2_shackRegistration(io_shackDict);

		return true;
	}

	bool setup(sweet::Data::Sphere2D::Operators *i_sphere2DOps)
	{
		opsReal = i_sphere2DOps;

		ops = new sphere2d_operators_type(i_sphere2DOps->sphere2DDataConfig, shackSphere2DDataOps);
		return true;
	}

	void check_error(double err, double rel_value)
	{
		if (err > double_precision_digits*rel_value || std::isnan(err))
		{
			std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
			std::cout << "+ err: " << err << std::endl;
			std::cout << "+ max_err_threshold: " << double_precision_digits*rel_value << std::endl;
			std::cout << "+ rel_value: " << rel_value << std::endl;
			SWEETErrorFatal("Error too high");
		}
	};


	void benchmark_setup_geostrophic_balance(
			sweet::Data::Sphere2D::DataSpectral &o_phi,
			sweet::Data::Sphere2D::DataSpectral &o_vrt,
			sweet::Data::Sphere2D::DataSpectral &o_div
	)
	{
		benchmarks.clear_3_benchmarkDetection();
		benchmarks.setup_3_benchmarkDetection("geostrophic_balance_linear_16");
		//benchmarks.setup_4_benchmarkSetup_1_withoutOps();
		benchmarks.setup_5_benchmarkSetup_2_withOps(opsReal);

		benchmarks.benchmark->getInitialState(o_phi, o_vrt, o_div);
	}

	void benchmark_setup_geostrophic_balance(
			sweet::Data::Sphere2DComplex::DataSpectral &o_phi,
			sweet::Data::Sphere2DComplex::DataSpectral &o_vrt,
			sweet::Data::Sphere2DComplex::DataSpectral &o_div
	)
	{

		sweet::Data::Sphere2D::DataSpectral phi(opsReal->sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral vrt(opsReal->sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral div(opsReal->sphere2DDataConfig);

		benchmark_setup_geostrophic_balance(phi, vrt, div);

		{
			// complex
			o_phi = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(phi);
			o_vrt = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(vrt);
			o_div = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(div);
		}

	}


	void benchmark_setup_pvd(
			sweet::Data::Sphere2D::DataSpectral &o_phi,
			sweet::Data::Sphere2D::DataSpectral &o_vrt,
			sweet::Data::Sphere2D::DataSpectral &o_div
	)
	{
		benchmarks.clear_3_benchmarkDetection();
		benchmarks.setup_3_benchmarkDetection("gaussian_bumps_pvd");
		//benchmarks.setup_4_benchmarkSetup_1_withoutOps();
		benchmarks.setup_5_benchmarkSetup_2_withOps(opsReal);


		benchmarks.benchmark->getInitialState(o_phi, o_vrt, o_div);
	}


	void benchmark_setup_pvd(
			sweet::Data::Sphere2DComplex::DataSpectral &o_phi,
			sweet::Data::Sphere2DComplex::DataSpectral &o_vrt,
			sweet::Data::Sphere2DComplex::DataSpectral &o_div
	)
	{

		sweet::Data::Sphere2D::DataSpectral phi(opsReal->sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral vrt(opsReal->sphere2DDataConfig);
		sweet::Data::Sphere2D::DataSpectral div(opsReal->sphere2DDataConfig);

		benchmark_setup_pvd(phi, vrt, div);

		{
			// complex
			o_phi = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(phi);
			o_vrt = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(vrt);
			o_div = sweet::Data::Sphere2D::Convert::DataSpectral_2_Sphere2DComplex_DataSpectral::convert(div);
		}
	}



public:
	void run_tests(
			sweet::Data::Sphere2D::Config *sphere2DDataConfig,
			phys_value_type &alpha
	)
	{
		shackDict->printProgramArguments();

		if (sizeof(phys_value_type) == sizeof(double))
		{
			double_precision_digits = 1e-10*2.0*sphere2DDataConfig->spectral_modes_m_max;
		}
		else
		{
			double_precision_digits = 1e-7*2.0*sphere2DDataConfig->spectral_modes_m_max;
		}

		/*
		 * Setup local simulation variables
		 */
		double gh0 = shackPDESWESphere2D->gravitation*shackPDESWESphere2D->h0;
		double sphere2d_radius = shackSphere2DDataOps->sphere2d_radius;

		double dt_implicit;

		if (gh0 == 0)
			dt_implicit = 1.0/sphere2d_radius;
		else
			dt_implicit = 1.0/std::sqrt(gh0)*sphere2d_radius;

		dt_implicit /= 4.0*2.0*sphere2DDataConfig->spectral_modes_n_max;

		double dt_two_omega = dt_implicit*2.0*shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;

		std::cout << "*********************************************************" << std::endl;
		std::cout << "* COMMON PARAMTERS" << std::endl;
		std::cout << "*********************************************************" << std::endl;
		std::cout << "dt_implicit: " << dt_implicit << std::endl;
		std::cout << "h0: " << shackPDESWESphere2D->h0 << std::endl;
		std::cout << "gh0: " << gh0 << std::endl;
		std::cout << "gravitation: " << shackPDESWESphere2D->gravitation << std::endl;
		std::cout << "sphere2d_rotating_coriolis_omega: " << shackPDESWESphere2D->sphere2d_rotating_coriolis_omega << std::endl;
		std::cout << "sphere2d_radius: " << shackSphere2DDataOps->sphere2d_radius << std::endl;


		/*
		 * Tests for stationary solution
		 */
		if (shackPDESWESphere2D->sphere2d_rotating_coriolis_omega == 0.0)
		{
			std::cout << std::endl;
			std::cout << "*********************************************************" << std::endl;
			std::cout << "* SKIPPING 'STATIONARY SOLUTION TESTS' (coriolis = 0 => only trivial stationary case exists)" << std::endl;
			std::cout << "*********************************************************" << std::endl;
		}
		else
		{
			std::cout << std::endl;
			std::cout << "*********************************************************" << std::endl;
			std::cout << "* STATIONARY SOLUTION TESTS" << std::endl;
			std::cout << "*********************************************************" << std::endl;


			sphere2d_data_spec_type phi(sphere2DDataConfig);
			sphere2d_data_spec_type vrt(sphere2DDataConfig);
			sphere2d_data_spec_type div(sphere2DDataConfig);

			benchmark_setup_geostrophic_balance(phi, vrt, div);

			sphere2d_data_spec_type foo, rhs;
			sph_banded_solver_type sphSolverTest;
			double dt = dt_implicit;
			double scalar = 1.234;
			double err = -1;

			double phi_max = phi.toGrid().grid_reduce_max_abs();
			double vrt_max = vrt.toGrid().grid_reduce_max_abs();
			double div_max = div.toGrid().grid_reduce_max_abs();

			if (phi_max == 0)
				phi_max = 1;

			if (vrt_max == 0)
				vrt_max = phi_max/(shackSphere2DDataOps->sphere2d_radius*shackSphere2DDataOps->sphere2d_radius);

			if (div_max == 0)
				div_max = phi_max/(shackSphere2DDataOps->sphere2d_radius*shackSphere2DDataOps->sphere2d_radius);

			std::cout << std::endl;
			std::cout << "phi_max: " << phi_max << std::endl;
			std::cout << "vrt_max: " << vrt_max << std::endl;
			std::cout << "div_max: " << div_max << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " J(x)-x" << std::endl;
			{
				foo = ops->implicit_J(vrt, dt_two_omega) - vrt;
				err = foo.toGrid().grid_reduce_max_abs();
				std::cout << " + vrt: 0 = " << err << std::endl;
				check_error(err, vrt_max);

				foo = ops->implicit_J(phi, dt_two_omega) -phi;
				err = foo.toGrid().grid_reduce_max_abs();
				std::cout << " + phi: 0 = " << err << std::endl;
				check_error(err, phi_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " x-Jinv(x)" << std::endl;
			{
				foo =vrt - ops->implicit_Jinv(vrt, dt_two_omega);
				err = foo.toGrid().grid_reduce_max_abs();
				std::cout << "vrt - Jinv(vrt) = 0 = " << err << std::endl;
				check_error(err, vrt_max);

				foo =phi - ops->implicit_Jinv(phi, dt_two_omega);
				err = foo.toGrid().grid_reduce_max_abs();
				std::cout << "phi - Jinv(phi) = 0 = " << err << std::endl;
				check_error(err, phi_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " FJinv(x)-F(Jinv(x))" << std::endl;
			{
				foo = ops->implicit_FJinv(vrt, dt_two_omega) - ops->implicit_F(ops->implicit_Jinv(vrt, dt_two_omega), dt_two_omega);
				err = foo.toGrid().grid_reduce_max_abs();
				std::cout << "FJinv(vrt) - F(Jinv(vrt)) = 0 = " << err << std::endl;
				check_error(err, vrt_max);

				foo = ops->implicit_FJinv(phi, dt_two_omega) - ops->implicit_F(ops->implicit_Jinv(phi, dt_two_omega), dt_two_omega);
				err = foo.toGrid().grid_reduce_max_abs();
				std::cout << "FJinv(phi) - F(Jinv(phi)) = 0 = " << err << std::endl;
				check_error(err, phi_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " -Linv(F(vrt)) - phi = 0" << std::endl;
			{
				foo = -ops->implicit_Linv(ops->implicit_F(vrt, dt_two_omega), dt) - phi;
				err = foo.toGrid().grid_reduce_max_abs();
				std::cout << "-Linv(F*vrt) - phi = 0 = " << err << std::endl;
				check_error(err, phi_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " -F(vrt) - L(phi) = 0" << std::endl;
			{
				foo = -ops->implicit_F(vrt, dt_two_omega) - ops->implicit_L(phi, dt);
				err = foo.toGrid().grid_reduce_max_abs();
				std::cout << "-F*vrt - L*phi = 0 = " << err << std::endl;
				check_error(err, div_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " invert(Jinv)" << std::endl;
			{
				sphSolverTest.setup(sphere2DDataConfig, 1);
				sphSolverTest.solver_addComponent_implicit_J(dt_two_omega);

				rhs = ops->implicit_J(vrt, dt_two_omega);
				foo = sphSolverTest.solve(rhs) -vrt;
				err = foo.toGrid().grid_reduce_max_abs();
				std::cout << " + vrt: 0 = " << err << std::endl;
				check_error(err, vrt_max);

				rhs = ops->implicit_J(phi, dt_two_omega);
				foo = sphSolverTest.solve(rhs) -phi;
				err = foo.toGrid().grid_reduce_max_abs();
				std::cout << " + phi: 0 = " << err << std::endl;
				check_error(err, phi_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " invert(F+I)" << std::endl;
			{
				sphSolverTest.setup(sphere2DDataConfig, 2);
				sphSolverTest.solver_addComponent_implicit_F(dt_two_omega);
				sphSolverTest.solver_addComponent_implicit_I(scalar);

				rhs = ops->implicit_F(vrt, dt_two_omega) + scalar*vrt;
				foo = sphSolverTest.solve(rhs) - vrt;
				err = foo.toGrid().grid_reduce_max_abs();
				std::cout << " + vrt: 0 = " << err << std::endl;
				check_error(err, vrt_max);

				rhs = ops->implicit_F(phi, dt_two_omega) + scalar*phi;
				foo = sphSolverTest.solve(rhs) -phi;
				err = foo.toGrid().grid_reduce_max_abs();
				std::cout << " + phi: 0 = " << err << std::endl;
				check_error(err, phi_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " b = F(a)+J(a)" << std::endl;
			std::cout << " c = FJinv(b)+J(b)" << std::endl;
			std::cout << " b' = invert(FJinv+J)*c" << std::endl;
			std::cout << " a' = invert(F+J)*b'" << std::endl;
			{
				sph_banded_solver_type sphSolverTestF_J;
				sphSolverTestF_J.setup(sphere2DDataConfig, 4);
				sphSolverTestF_J.solver_addComponent_implicit_F(dt_two_omega);
				sphSolverTestF_J.solver_addComponent_implicit_J(dt_two_omega);

				sph_banded_solver_type sphSolverTestFJinv_J;
				sphSolverTestFJinv_J.setup(sphere2DDataConfig, 4);
				sphSolverTestFJinv_J.solver_addComponent_implicit_FJinv(dt_two_omega);
				sphSolverTestFJinv_J.solver_addComponent_implicit_J(dt_two_omega);

				sphere2d_data_spec_type a = vrt;
				sphere2d_data_spec_type b = ops->implicit_F(a, dt_two_omega) + ops->implicit_J(a, dt_two_omega);
				sphere2d_data_spec_type c = ops->implicit_FJinv(b, dt_two_omega) + ops->implicit_J(b, dt_two_omega);
				sphere2d_data_spec_type b_ = sphSolverTestFJinv_J.solve(c);
				sphere2d_data_spec_type a_ = sphSolverTestF_J.solve(b_);

				err = (b-b_).toGrid().grid_reduce_max_abs();
				std::cout << " + b_: 0 = " << err << std::endl;

				err = (b-b_).toGrid().grid_reduce_max_abs();
				std::cout << " + a_: 0 = " << err << std::endl;
				check_error(err, a.toGrid().grid_reduce_max_abs());
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " invert(FJinvF+J)" << std::endl;
			{
				sphSolverTest.setup(sphere2DDataConfig, 4);
				sphSolverTest.solver_addComponent_implicit_FJinvF(dt_two_omega);
				sphSolverTest.solver_addComponent_implicit_J(dt_two_omega);

				rhs = ops->implicit_FJinv(ops->implicit_F(vrt, dt_two_omega), dt_two_omega) + ops->implicit_J(vrt, dt_two_omega);
				foo = sphSolverTest.solve(rhs) - vrt;

				err = foo.toGrid().grid_reduce_max_abs();
				std::cout << " + vrt: 0 = " << err << std::endl;
				check_error(err, vrt_max);

				rhs = ops->implicit_FJinv(ops->implicit_F(phi, dt_two_omega), dt_two_omega) + ops->implicit_J(phi, dt_two_omega);
				foo = sphSolverTest.solve(rhs) - phi;
				err = foo.toGrid().grid_reduce_max_abs();
				std::cout << " + phi: 0 = " << err << std::endl;
				check_error(err, phi_max);
			}
			std::cout << std::endl;

			std::cout << "============================================" << std::endl;
			std::cout << " FIN" << std::endl;
			std::cout << "============================================" << std::endl;
		}



		/*
		 * Tests with time-varying solution
		 */
		for (int i = 0; i < 8; i++)
		{
			std::cout << std::endl;
			std::cout << "*********************************************************" << std::endl;
			std::cout << "* TIME VARYING SOLUTION TESTS " << i << std::endl;
			std::cout << "*********************************************************" << std::endl;


			sphere2d_data_spec_type phi(sphere2DDataConfig);
			sphere2d_data_spec_type vrt(sphere2DDataConfig);
			sphere2d_data_spec_type div(sphere2DDataConfig);

			benchmark_setup_pvd(phi, vrt, div);

			double phi_order = phi.toGrid().grid_reduce_max_abs();
			double vrt_order = vrt.toGrid().grid_reduce_max_abs();
			double div_order = div.toGrid().grid_reduce_max_abs();

			if ((i & 1) == 0)
			{
				phi.spectral_setZero();
				phi_order = 1.0;
				std::cout << " + phi: set to zero" << std::endl;
			}
			else
			{
				std::cout << " + phi: " << phi_order << std::endl;
			}

			if ((i & 2) == 0)
			{
				vrt.spectral_setZero();
				vrt_order = phi_order/(shackSphere2DDataOps->sphere2d_radius*shackSphere2DDataOps->sphere2d_radius);
				std::cout << " + vrt: set to zero" << std::endl;
			}
			else
			{
				std::cout << " + vrt: " << vrt_order << std::endl;
			}

			if ((i & 4) == 0)
			{
				div.spectral_setZero();
				div_order = phi_order/(shackSphere2DDataOps->sphere2d_radius*shackSphere2DDataOps->sphere2d_radius);
				std::cout << " + div: set to zero" << std::endl;
			}
			else
			{
				std::cout << " + div: " << div_order << std::endl;
			}


			std::cout << std::endl;
			std::cout << "phi_order: " << phi_order << std::endl;
			std::cout << "vrt_order: " << vrt_order << std::endl;
			std::cout << "div_order: " << div_order << std::endl;


			double gh0 = shackPDESWESphere2D->gravitation*shackPDESWESphere2D->h0;
			double sphere2d_radius = shackSphere2DDataOps->sphere2d_radius;

			double dt_two_omega = dt_implicit*2.0*shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;

			sph_banded_solver_type sphSolverTest;

			sphere2d_data_spec_type& S1 = vrt;
			sphere2d_data_spec_type& S2 = div;
			sphere2d_data_spec_type& S3 = phi;

			sph_banded_solver_type sphSolverDiv;
			sphSolverDiv.setup(sphere2DDataConfig, 4);
			sphSolverDiv.solver_addComponent_implicit_J(dt_two_omega);
			sphSolverDiv.solver_addComponent_implicit_FJinvF(dt_two_omega);
			sphSolverDiv.solver_addComponent_implicit_L(gh0*dt_implicit, dt_implicit, sphere2d_radius);

			sphere2d_data_spec_type rhs = S2 + ops->implicit_FJinv(S1, dt_two_omega) + ops->implicit_L(S3, dt_implicit);
			sphere2d_data_spec_type div1 = sphSolverDiv.solve(rhs);

			sphere2d_data_spec_type phi1 = S3 - dt_implicit*gh0*div1;
			sphere2d_data_spec_type vrt1 = ops->implicit_Jinv(S1 - ops->implicit_F(div1, dt_two_omega), dt_two_omega);

			/*
			 * Testcase for correct solver
			 */
			{
				double gh0 = shackPDESWESphere2D->gravitation*shackPDESWESphere2D->h0;
				double dt_two_omega = dt_implicit*2.0*shackPDESWESphere2D->sphere2d_rotating_coriolis_omega;

				sphere2d_data_spec_type foo, rhs;
				sph_banded_solver_type sphSolverTest;
				double err = -1;

				std::cout << "============================================" << std::endl;
				std::cout << " CURL equation" << std::endl;
				{
					foo = ops->implicit_J(vrt1, dt_two_omega) + ops->implicit_F(div1, dt_two_omega) - S1;
					err = foo.toGrid().grid_reduce_max_abs();
					std::cout << " + err: 0 = " << err << std::endl;
					check_error(err, vrt_order);
				}

				std::cout << "============================================" << std::endl;
				std::cout << " DIV equation" << std::endl;
				{
					foo = ops->implicit_J(div1, dt_two_omega) - ops->implicit_F(vrt1, dt_two_omega) - ops->implicit_L(phi1, dt_implicit) - S2;
					err = foo.toGrid().grid_reduce_max_abs();
					std::cout << " + err: 0 = " << err << std::endl;
					check_error(err, div_order);
				}

				std::cout << "============================================" << std::endl;
				std::cout << " GEOPOT equation" << std::endl;
				{
					foo = phi1 + dt_implicit*gh0*div1 - S3;
					err = foo.toGrid().grid_reduce_max_abs();
					std::cout << " + err: 0 = " << err << std::endl;
					check_error(err, phi_order);
				}
				std::cout << std::endl;
			}
			std::cout << "Tests successful" << std::endl;

		}



		std::cout << "*********************************************************" << std::endl;
		std::cout << "* COMMON PARAMTERS" << std::endl;
		std::cout << "*********************************************************" << std::endl;
		std::cout << " + dt_implicit: " << dt_implicit << std::endl;
		std::cout << " + h0: " << shackPDESWESphere2D->h0 << std::endl;
		std::cout << " + gh0: " << gh0 << std::endl;
		std::cout << " + gravitation: " << shackPDESWESphere2D->gravitation << std::endl;
		std::cout << " + sphere2d_rotating_coriolis_omega: " << shackPDESWESphere2D->sphere2d_rotating_coriolis_omega << std::endl;
		std::cout << " + sphere2d_radius: " << shackSphere2DDataOps->sphere2d_radius << std::endl;
		std::cout << std::endl;


		std::cout << "*********************************************************" << std::endl;
		std::cout << "All tests successful" << std::endl;
		std::cout << "*********************************************************" << std::endl;
	}
};




template<
	typename phys_value_type = double,
	typename sphere2d_data_spec_type = sweet::Data::Sphere2D::DataSpectral,
	typename sphere2d_data_phys_type = sweet::Data::Sphere2D::DataGrid,
	typename sphere2d_operators_type = sweet::Data::Sphere2D::Operators,
	typename sph_banded_solver_type = double
>
class TestFB
{
	double double_precision_digits = -1.0;

	sweet::Shacks::Dictionary *shackDict;
	PDE_SWESphere2D::Shack *shackPDESWESphere2D;
	sweet::Data::Sphere2D::Shack *shackSphere2DDataOps;

	sweet::Data::Sphere2D::Operators *opsReal;
	sphere2d_operators_type *ops;

	PDESWESphere2DTS_l_erk l_erk;
	PDESWESphere2DTS_l_irk l_irk;

	PDE_SWESphere2D::Benchmarks::BenchmarksCombined benchmarks;

public:
	sweet::Error::Base error;

	TestFB()	:
		ops(nullptr)
	{
	}

	~TestFB()
	{
		delete ops;
	}

	bool shackRegistration(
		sweet::Shacks::Dictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackPDESWESphere2D = io_shackDict->getAutoRegistration<PDE_SWESphere2D::Shack>();
		shackSphere2DDataOps = io_shackDict->getAutoRegistration<sweet::Data::Sphere2D::Shack>();
		//shackTimestepControl = shackDict->getAutoRegistration<sweet::TimeTree::ShackTimestepControl>();
		//shackSphere2DDataOps = shackDict->getAutoRegistration<sweet::Data::Sphere2D::ShackSphere2DDataOps>();
		//shackPDESWETimeDisc = shackDict->getAutoRegistration<ShackTimeDiscretization>();
		//shackPDESWEBenchmark = shackDict->getAutoRegistration<PDE_SWESphere2D::Benchmarks::Shack>();

		benchmarks.setup_1_registerAllBenchmark();
		benchmarks.setup_2_shackRegistration(io_shackDict);

		l_erk.shackRegistration(shackDict);
		l_irk.shackRegistration(shackDict);
		return true;
	}


	bool setup(sweet::Data::Sphere2D::Operators *i_sphere2DOps)
	{
		opsReal = i_sphere2DOps;

#if 1
		benchmarks.clear_3_benchmarkDetection();
		benchmarks.setup_3_benchmarkDetection("gaussian_bumps_pvd");
		//benchmarks.setup_4_benchmarkSetup_1_withoutOps();
		benchmarks.setup_5_benchmarkSetup_2_withOps(opsReal);
#endif

		//l_irk.setup(opsReal, 1, dt_implicit);
		l_irk.setup(opsReal, 1, 1.0);
		l_erk.setup_main(opsReal, 1);

		ops = new sphere2d_operators_type(i_sphere2DOps->sphere2DDataConfig, shackSphere2DDataOps);

		return true;
	}

	void check_error(double err, double rel_value)
	{
		if (err > double_precision_digits*rel_value || std::isnan(err))
		{
			std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
			std::cout << "+ err: " << err << std::endl;
			std::cout << "+ max_err_threshold: " << double_precision_digits*rel_value << std::endl;
			std::cout << "+ rel_value: " << rel_value << std::endl;
			SWEETErrorFatal("Error too high");
		}
	};

public:
	void run_tests(
			sweet::Data::Sphere2D::Config *i_sphere2DDataConfig,
			phys_value_type &i_alpha
	)
	{
		shackDict->printProgramArguments();

		double_precision_digits = 1e-10*2.0*i_sphere2DDataConfig->spectral_modes_m_max;

		/*
		 * Setup local simulation variables
		 */
		double gh0 = shackPDESWESphere2D->gravitation*shackPDESWESphere2D->h0;
		double sphere2d_radius = shackSphere2DDataOps->sphere2d_radius;

		double dt_implicit;

		if (gh0 == 0)
			dt_implicit = 1.0/sphere2d_radius;
		else
			dt_implicit = 1.0/std::sqrt(gh0)*sphere2d_radius;

		dt_implicit /= 4.0*2.0*i_sphere2DDataConfig->spectral_modes_n_max;

		std::cout << "*********************************************************" << std::endl;
		std::cout << "* COMMON PARAMTERS" << std::endl;
		std::cout << "*********************************************************" << std::endl;
		std::cout << "dt_implicit: " << dt_implicit << std::endl;
		std::cout << "h0: " << shackPDESWESphere2D->h0 << std::endl;
		std::cout << "gh0: " << gh0 << std::endl;
		std::cout << "gravitation: " << shackPDESWESphere2D->gravitation << std::endl;
		std::cout << "sphere2d_rotating_coriolis_omega: " << shackPDESWESphere2D->sphere2d_rotating_coriolis_omega << std::endl;
		std::cout << "sphere2d_radius: " << shackSphere2DDataOps->sphere2d_radius << std::endl;

		/*
		 * Test forward/backward time stepping
		 *
		 * (I - dtexpl*L)^-1 (I + dtimpl*L) = I
		 */
		{
			std::cout << std::endl;
			std::cout << "*********************************************************" << std::endl;
			std::cout << "* FORWARD / BACKWARD TIME INTEGRATION" << std::endl;
			std::cout << "*********************************************************" << std::endl;

			sweet::Data::Sphere2D::DataSpectral phi(i_sphere2DDataConfig);
			sweet::Data::Sphere2D::DataSpectral vrt(i_sphere2DDataConfig);
			sweet::Data::Sphere2D::DataSpectral div(i_sphere2DDataConfig);

			benchmarks.benchmark->getInitialState(phi, vrt, div);

			//phi.spectral_setZero();
			vrt.spectral_setZero();
			div.spectral_setZero();

			double dt_implicit = 1.0/std::sqrt(gh0)*sphere2d_radius;
			dt_implicit = 0.5*i_sphere2DDataConfig->spectral_modes_n_max;
			std::cout << "dt_implicit: " << dt_implicit << std::endl;

			double phi_max = phi.toGrid().grid_reduce_max_abs();
			double vrt_max = vrt.toGrid().grid_reduce_max_abs();
			double div_max = div.toGrid().grid_reduce_max_abs();

			if (phi_max == 0)
				phi_max = 1;

			if (vrt_max == 0)
				vrt_max = phi_max/(shackSphere2DDataOps->sphere2d_radius*shackSphere2DDataOps->sphere2d_radius);

			if (div_max == 0)
				div_max = phi_max/(shackSphere2DDataOps->sphere2d_radius*shackSphere2DDataOps->sphere2d_radius);

			std::cout << std::endl;
			std::cout << "phi_max: " << phi_max << std::endl;
			std::cout << "vrt_max: " << vrt_max << std::endl;
			std::cout << "div_max: " << div_max << std::endl;
			std::cout << std::endl;

			{
				sweet::Data::Sphere2D::DataSpectral test_phi = phi;
				sweet::Data::Sphere2D::DataSpectral test_vrt = vrt;
				sweet::Data::Sphere2D::DataSpectral test_div = div;

				double dt_explicit = dt_implicit;
				double dt_implicit = -dt_explicit;		// Backward with flipped minus sigh

				l_erk.runTimestep(test_phi, test_vrt, test_div, dt_explicit);

				l_irk.runTimestep(test_phi, test_vrt, test_div, dt_implicit);

				double err_phi = (test_phi-phi).toGrid().grid_reduce_max_abs();
				double err_vrt = (test_vrt-vrt).toGrid().grid_reduce_max_abs();
				double err_div = (test_div-div).toGrid().grid_reduce_max_abs();

				std::cout << "err_phi: " << err_phi << std::endl;
				check_error(err_phi, phi_max);

				std::cout << "err_vrt: " << err_vrt << std::endl;
				check_error(err_vrt, vrt_max);

				std::cout << "err_div: " << err_div << std::endl;
				check_error(err_div, div_max);
			}
		}

		std::cout << "*********************************************************" << std::endl;
		std::cout << "* COMMON PARAMTERS" << std::endl;
		std::cout << "*********************************************************" << std::endl;
		std::cout << " + dt_implicit: " << dt_implicit << std::endl;
		std::cout << " + h0: " << shackPDESWESphere2D->h0 << std::endl;
		std::cout << " + gh0: " << gh0 << std::endl;
		std::cout << " + gravitation: " << shackPDESWESphere2D->gravitation << std::endl;
		std::cout << " + sphere2d_rotating_coriolis_omega: " << shackPDESWESphere2D->sphere2d_rotating_coriolis_omega << std::endl;
		std::cout << " + sphere2d_radius: " << shackSphere2DDataOps->sphere2d_radius << std::endl;
		std::cout << std::endl;


		std::cout << "*********************************************************" << std::endl;
		std::cout << "All tests successful" << std::endl;
		std::cout << "*********************************************************" << std::endl;
	}
};



int main(
		int i_argc,
		char *const i_argv[]
)
{
	for (int i = 0; i < 2; i++)
	{
		sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict(i_argc, i_argv);
		shackProgArgDict.setup();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

		sweet::Data::Sphere2D::Shack *shackSphere2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Sphere2D::Shack>();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

		shackProgArgDict.printShackData();

		sweet::Data::Sphere2D::Config sphere2DDataConfig;
		sphere2DDataConfig.setupAuto(shackSphere2DDataOps);
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(sphere2DDataConfig);

		sweet::Data::Sphere2D::Operators ops;
		ops.setup(&sphere2DDataConfig, shackSphere2DDataOps);
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(ops);

		// Stop overriding simulation variables for this test case
		//simVars.benchmark.benchmark_override_simvars = false;

		shackSphere2DDataOps->validateResolution();
		ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackSphere2DDataOps);


		if (i == 0)
		{
			/*
			 * Real-valued solver
			 *
			 * E.g. used for backward Euler
			 */
			double alpha_real = 3.0;

			std::cout << "***********************************************************" << std::endl;
			std::cout << "* REAL" << std::endl;
			std::cout << "***********************************************************" << std::endl;
			Test<
				double,
				sweet::Data::Sphere2D::DataSpectral,
				sweet::Data::Sphere2D::DataGrid,
				sweet::Data::Sphere2D::Operators,
				SphBandedMatrix_GridReal
			>t_real;
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(t_real);

			t_real.shackRegistration(&shackProgArgDict);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(t_real);

			t_real.setup(&ops);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(t_real);

			t_real.run_tests(&sphere2DDataConfig, alpha_real);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(t_real);

			std::cout << "***********************************************************" << std::endl;
			std::cout << "* REAL FB" << std::endl;
			std::cout << "***********************************************************" << std::endl;
			TestFB<
				double,
				sweet::Data::Sphere2D::DataSpectral,
				sweet::Data::Sphere2D::DataGrid,
				sweet::Data::Sphere2D::Operators,
				SphBandedMatrix_GridReal
			>t_real_fb;
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(t_real_fb);

			t_real_fb.shackRegistration(&shackProgArgDict);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(t_real_fb);

			t_real_fb.setup(&ops);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(t_real_fb);

			t_real_fb.run_tests(&sphere2DDataConfig, alpha_real);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(t_real_fb);
		}

		if (i == 1)
		{
			/*
			 * Complex-valued solver
			 *
			 * E.g. used for REXI
			 */
			std::cout << "***********************************************************" << std::endl;
			std::cout << "* COMPLEX" << std::endl;
			std::cout << "***********************************************************" << std::endl;

			std::complex<double> alpha_complex(1.0, 3.0);

			Test<
				std::complex<double>,
				sweet::Data::Sphere2DComplex::DataSpectral,
				sweet::Data::Sphere2DComplex::DataGrid,
				sweet::Data::Sphere2DComplex::Operators,
				SphBandedMatrix_GridComplex
			>t_complex;
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(t_complex);

			t_complex.shackRegistration(&shackProgArgDict);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(t_complex);

			t_complex.setup(&ops);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(t_complex);

			t_complex.run_tests(&sphere2DDataConfig, alpha_complex);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(t_complex);
		}
	}

}

