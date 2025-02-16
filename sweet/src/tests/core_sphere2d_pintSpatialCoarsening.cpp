/*
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWESphere2D/Benchmarks/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/PDE_SWESphere2D/BenchmarksCombined.cpp
 *
 * MULE_SCONS_OPTIONS: --parareal-sphere2d=enable
 * MULE_SCONS_OPTIONS: --sphere2d-spectral-space=enable
 *
 *  Created on: 26 Jul 2022
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */


#include <programs/PDE_AdvectionSphere2D/PDEAdvectionSphere2DBenchmarksCombined.hpp>
#include <programs/PDE_AdvectionSphere2D/PDEAdvectionSphere2DTimeSteppers.hpp>
#include <programs/PDE_AdvectionSphere2D/ProgramPDEAdvectionSphere2D.hpp>
#include <sweet/Data/Sphere2D/DataSpectral.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Sphere2DData_Spectral.hpp>
#include <sweet/Data/Sphere2D/Convert/DataGrid_2_Cart2D_DataGrid.hpp>
#include <sweet/Data/Sphere2D/Convert/DataSpectral_2_Cart2D_DataGrid.hpp>
#include <sweet/SemiLagrangian/Shack.hpp>
#include <sweet/SemiLagrangian/Sphere2D.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include "../programs/PDE_AdvectionSphere2D/time/ShackPDEAdvectionSphere2DTimeDiscretization.hpp"

#include "../programs/PDE_SWESphere2D/BenchmarksCombined.hpp"


void printTest(const std::string &test_name)
{
	std::cout << std::endl;
	std::cout << " ******************************************* " << std::endl;
	std::cout << " ****** " << test_name << std::endl;
	std::cout << " ******************************************* " << std::endl;
}

void printError(double val, const std::string &msg = "")
{
	std::cout << "    ----> ERROR " << msg << " : " << val << std::endl << std::endl;
}


class Data
{
public:
	sweet::Error::Base error;

	sweet::Data::Sphere2D::Config* sphere2DDataConfig;
	sweet::Data::Sphere2D::Operators* ops;
	PDE_SWESphere2D::Benchmarks::BenchmarksCombined sphere2DBenchmarks;
	sweet::DEPRECATED_pint::Parareal_GenericData* data = nullptr;

public:
	Data()
	{
	}

	~Data()
	{
		if (data)
		{
			delete data;
			data = nullptr;
		}
	}

	bool setup_1_shackRegister(sweet::Shacks::Dictionary *io_shackDict)
	{
		sphere2DBenchmarks.setup_1_registerAllBenchmark();
		sphere2DBenchmarks.setup_2_shackRegistration(io_shackDict);
		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(sphere2DBenchmarks);
		return true;
	}

	bool setup_2_benchmarkDetection()
	{
		sphere2DBenchmarks.setup_3_benchmarkDetection();
		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(sphere2DBenchmarks);
		return true;
	}

	bool setup_3_data(
			sweet::Data::Sphere2D::Config* i_sphere2DDataConfig,
			sweet::Data::Sphere2D::Operators* i_ops
	)
	{
		sphere2DDataConfig = i_sphere2DDataConfig;
		ops = i_ops;

		data = new sweet::DEPRECATED_pint::Parareal_GenericData_Sphere2DData_Spectral<3>;
		data->setup_data_config(sphere2DDataConfig);
		data->allocate_data();
		sweet::Data::Sphere2D::DataSpectral* phi_pert = data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[0];
		sweet::Data::Sphere2D::DataSpectral* vrt = data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[1];
		sweet::Data::Sphere2D::DataSpectral* div = data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[2];

		sphere2DBenchmarks.setup_4_benchmarkSetup_1_withoutOps();
		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(sphere2DBenchmarks);

		sphere2DBenchmarks.setup_5_benchmarkSetup_2_withOps(ops);
		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(sphere2DBenchmarks);

		sphere2DBenchmarks.benchmark->getInitialState(*phi_pert, *vrt, *div);
		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(sphere2DBenchmarks);

		return true;
	}

	void set_zero()
	{
		for (int i = 0; i < 3; i++)
			data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[i]->spectral_setZero();
	}

	void restrict(Data* i_data)
	{
		data->restrict(*i_data->data);
	}

	void restrict(Data& i_data)
	{
		data->restrict(*i_data.data);
	}

	void pad_zeros(Data* i_data)
	{
		data->pad_zeros(*i_data->data);
	}
	void pad_zeros(Data& i_data)
	{
		data->pad_zeros(*i_data.data);
	}
};




int main(int i_argc, char *i_argv[])
{

	double eps = 1e-15;


	int N_Hs[5] = {16, 32, 64, 128, 256};
	int N_Ls[5] = {8, 16, 32, 64, 128};

	for (int i_H = 0; i_H < 5; i_H++)
	{
		for (int i_L = 0; i_L < 5; i_L++)
		{
			sweet::Data::Sphere2D::Config sphere2DDataConfig_H;
			sweet::Data::Sphere2D::Config sphere2DDataConfig_L;


			int N_H = N_Hs[i_H];
			int N_L = N_Ls[i_L];

			if (N_L >= N_H)
				continue;

			std::cout << std::endl;
			std::cout << "-------------------------------------------" << std::endl;
			std::cout << "-------------------------------------------" << std::endl;
			std::cout << " TESTING FOR N_H = " << N_H << "; " << "N_L = " << N_L << std::endl;
			std::cout << "-------------------------------------------" << std::endl;
			std::cout << "-------------------------------------------" << std::endl;



			sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict(i_argc, i_argv);
			shackProgArgDict.setup();
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

			sweet::Data::Sphere2D::Shack *shackSphere2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Sphere2D::Shack>();
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);


			Data data_H;
			Data data_L;
			Data data_H_2_H;
			Data data_H_2_L;
			Data data_H_2_L_2_H;
			Data data_L_2_H;
			Data data_L_2_H_2_L;

			data_H.setup_1_shackRegister(&shackProgArgDict);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(data_H);
			data_L.setup_1_shackRegister(&shackProgArgDict);
			data_H_2_H.setup_1_shackRegister(&shackProgArgDict);
			data_H_2_L.setup_1_shackRegister(&shackProgArgDict);
			data_H_2_L_2_H.setup_1_shackRegister(&shackProgArgDict);
			data_L_2_H.setup_1_shackRegister(&shackProgArgDict);
			data_L_2_H_2_L.setup_1_shackRegister(&shackProgArgDict);


			shackProgArgDict.processProgramArguments();
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);


			data_H.setup_2_benchmarkDetection();
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(data_H);
			data_L.setup_2_benchmarkDetection();
			data_H_2_H.setup_2_benchmarkDetection();
			data_H_2_L.setup_2_benchmarkDetection();
			data_H_2_L_2_H.setup_2_benchmarkDetection();
			data_L_2_H.setup_2_benchmarkDetection();
			data_L_2_H_2_L.setup_2_benchmarkDetection();



			shackProgArgDict.printShackData();


			sweet::Data::Sphere2D::Shack shackSphere2DDataOps_H;
			shackSphere2DDataOps_H = *shackSphere2DDataOps;
			shackSphere2DDataOps_H.space_res_physical[0] = -1;
			shackSphere2DDataOps_H.space_res_physical[1] = -1;
			shackSphere2DDataOps_H.space_res_spectral[0] = N_H;
			shackSphere2DDataOps_H.space_res_spectral[1] = N_H;

			sweet::Data::Sphere2D::Shack shackSphere2DDataOps_L;
			shackSphere2DDataOps_L = *shackSphere2DDataOps;
			shackSphere2DDataOps_L.space_res_physical[0] = -1;
			shackSphere2DDataOps_L.space_res_physical[1] = -1;
			shackSphere2DDataOps_L.space_res_spectral[0] = N_L;
			shackSphere2DDataOps_L.space_res_spectral[1] = N_L;

			sphere2DDataConfig_H.setupAuto(shackSphere2DDataOps_H);
			sphere2DDataConfig_L.setupAuto(shackSphere2DDataOps_L);

			sweet::Data::Sphere2D::Operators ops_H(&sphere2DDataConfig_H, &shackSphere2DDataOps_H);
			sweet::Data::Sphere2D::Operators ops_L(&sphere2DDataConfig_L, &shackSphere2DDataOps_L);

			data_H.setup_3_data(&sphere2DDataConfig_H, &ops_H);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(data_H);
			data_L.setup_3_data(&sphere2DDataConfig_L, &ops_L);
			data_H_2_H.setup_3_data(&sphere2DDataConfig_H, &ops_H);
			data_H_2_L.setup_3_data(&sphere2DDataConfig_L, &ops_L);
			data_H_2_L_2_H.setup_3_data(&sphere2DDataConfig_H, &ops_H);
			data_L_2_H.setup_3_data(&sphere2DDataConfig_H, &ops_H);
			data_L_2_H_2_L.setup_3_data(&sphere2DDataConfig_L, &ops_L);

			// Error storage
			sweet::DEPRECATED_pint::Parareal_GenericData* error_H = new sweet::DEPRECATED_pint::Parareal_GenericData_Sphere2DData_Spectral<3>;
			sweet::DEPRECATED_pint::Parareal_GenericData* error_L = new sweet::DEPRECATED_pint::Parareal_GenericData_Sphere2DData_Spectral<3>;
			error_H->setup_data_config(&sphere2DDataConfig_H);
			error_L->setup_data_config(&sphere2DDataConfig_L);
			error_H->allocate_data();
			error_L->allocate_data();


			// Tests 1 and 2:
			data_H_2_H.set_zero();

			printTest("Test 1: high res -> high res (dummy restriction) ");
			data_H_2_H.restrict(data_H);
			*error_H = *data_H.data;
			*error_H -= *data_H_2_H.data;
			double error = error_H->spectral_reduce_maxAbs();
			printError(error);
			SWEET_ASSERT(error < eps);

			data_H_2_H.set_zero();
			printTest("Test 2: high res -> high res (dummy prolongation) ");
			data_H_2_H.pad_zeros(data_H);
			*error_H = *data_H.data;
			*error_H -= *data_H_2_H.data;
			error = error_H->spectral_reduce_maxAbs();
			printError(error);
			SWEET_ASSERT(error < eps);

			// Test 3:
			printTest("Test 3: high res -> low res (restriction) ");
			data_H_2_L.set_zero();
			data_H_2_L.restrict(data_H);
			*error_L = *data_L.data;
			*error_L -= *data_H_2_L.data;
			error = error_L->spectral_reduce_maxAbs();
			printError(error);
			SWEET_ASSERT(error < eps);

			// Test 4:
			printTest("Test 4: high res -> low res -> high res (restriction + prolongation) ");
			data_H_2_L_2_H.set_zero();
			data_H_2_L_2_H.pad_zeros(data_H_2_L);
			*error_H = *data_H.data;
			*error_H -= *data_H_2_L_2_H.data;
			error = error_H->spectral_reduce_maxAbs(N_L - 1);
			printError(error, "(Up to mode " + std::to_string(N_L - 1) + ")" );
			printError(error_H->spectral_reduce_maxAbs(N_L), "(Up to mode " + std::to_string(N_L) + ")" );
			SWEET_ASSERT(error < eps);

			// Test 5:
			printTest("Test 5: low res -> high res (prolongation) ");
			data_L_2_H.set_zero();
			data_L_2_H.restrict(data_L);
			*error_H = *data_H.data;
			*error_H -= *data_L_2_H.data;
			error = error_H->spectral_reduce_maxAbs(N_L - 1);
			printError(error, "(Up to mode " + std::to_string(N_L - 1) + ")" );
			printError(error_H->spectral_reduce_maxAbs(N_L), "(Up to mode " + std::to_string(N_L) + ")" );
			SWEET_ASSERT(error < eps);

			// Test 6:
			printTest("Test 6: low res -> high res -> low res (prolongation + restriction) ");
			data_L_2_H_2_L.set_zero();
			data_L_2_H_2_L.pad_zeros(data_L_2_H);
			*error_L = *data_L.data;
			*error_L -= *data_L_2_H_2_L.data;
			error = error_L->spectral_reduce_maxAbs();
			printError(error);
			SWEET_ASSERT(error < eps);

			delete error_H;
			delete error_L;

			std::cout << "-------------------------------------------" << std::endl;
			std::cout << "-------------------------------------------" << std::endl;
			std::cout << " TEST SUCCESSFUL" << std::endl;
			std::cout << "-------------------------------------------" << std::endl;
			std::cout << "-------------------------------------------" << std::endl;
			std::cout << std::endl;
		}
	}

	std::cout << std::endl;
	std::cout << " !!!!!!!!!!!!!!!!!! " << std::endl;
	std::cout << "  ALL TESTS PASSED " << std::endl;
	std::cout << " !!!!!!!!!!!!!!!!!! " << std::endl;
	std::cout << std::endl;
}
