/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <iostream>
#include <complex>

// Our shack directory to store different objects and get them back later on
#include <sweet/ExpIntegration/ExpFunction.hpp>
#include <sweet/ExpIntegration/REXI/REXI.hpp>
#include <sweet/ExpIntegration/REXI/REXICoefficients.hpp>
#include <sweet/ExpIntegration/REXI/REXICoefficientsSet.hpp>
#include <sweet/ExpIntegration/Shack.hpp>
#include <sweet/Shacks/Base.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/TimeTree/Shack.hpp>

typedef double T;
typedef std::complex<T> TComplex;

std::string function_name;

/**
 * Values and parameters to setup benchmarks simulations
 */
class ShackREXITest	:
		public sweet::Shacks::Base
{
public:
	std::string function_name;
	double lambda_real = std::numeric_limits<double>::infinity();
	double lambda_imag = std::numeric_limits<double>::infinity();
	std::complex<double> lambda = 0;
	int test_mode = 0;

	void printProgramArguments(const std::string& i_prefix = "") override
	{

		std::cout << "REXI APPROXIMATION PROGRAM:" << std::endl;
		std::cout << std::endl;
		std::cout << "	--function-name=..." << std::endl;
		std::cout << "	--lambda-real=..." << std::endl;
		std::cout << "	--lambda-imag=..." << std::endl;
		std::cout << "  --test-mode=..." << std::endl;
		std::cout << "        0: use standard time stepping" << std::endl;
		std::cout << "        1: always start from u(0) with increasing time step sizes" << std::endl;
		std::cout << std::endl;
	}

	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa) override
	{
		i_pa.getArgumentValueByKey("--function-name", function_name);
		i_pa.getArgumentValueByKey("--lambda-real", lambda_real);
		i_pa.getArgumentValueByKey("--lambda-imag", lambda_imag);
		i_pa.getArgumentValueByKey("--test-mode", test_mode);

		lambda = std::complex<double>(lambda_real, lambda_imag);

		return error.forwardWithPositiveReturn(i_pa.error);
	}

	void printShack(
		const std::string& i_prefix = ""
	) override
	{
		std::cout << std::endl;
		std::cout << "REXI APPROXIMATION PROGRAM:" << std::endl;
		std::cout << " + function_name: " << function_name << std::endl;
		std::cout << " + lambda_real: " << lambda_real << std::endl;
		std::cout << " + lambda_imag: " << lambda_imag << std::endl;
		std::cout << " + test_mode: " << test_mode << std::endl;
		std::cout << std::endl;
	}

	bool validateLambda()
	{
		if (std::abs(lambda) == 0)
			return error.set("Specify lambda coefficient of linear operators using --lambda=[float]");

		return true;
	}
};



int main(
		int i_argc,
		char *const i_argv[]
)
{
	sweet::Shacks::ProgramArgumentsDictionary shackDict(i_argc, i_argv);

	sweet::TimeTree::Shack *shackTimestepControl = shackDict.getAutoRegistration<sweet::TimeTree::Shack>();
	ShackREXITest *shackREXITest = shackDict.getAutoRegistration<ShackREXITest>();
	sweet::ExpIntegration::Shack *shackExpIntegration =  shackDict.getAutoRegistration<sweet::ExpIntegration::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackDict);

	shackDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackDict);

	shackTimestepControl->validateTimestepSize();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackTimestepControl);

	shackREXITest->validateLambda();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackREXITest);

	/*
	 * Load analytical function
	 */
	sweet::ExpIntegration::ExpFunction<double> expFunctions;
	expFunctions.setup(function_name);

	/*
	 * Load REXI coefficients from file
	 */
	sweet::ExpIntegration::REXI::REXICoefficientsSet<> rexiCoefficientsSet;

	if (shackExpIntegration->exp_method == "direct")
	{
		SWEETErrorFatal("Direct REXI mode not supported");
	}
	else if (shackExpIntegration->exp_method == "file")
	{
		rexiCoefficientsSet.setupFromFiles(shackExpIntegration->rexi_files);

		if (rexiCoefficientsSet.rexiCoefficientVector.size() == 0)
			SWEETErrorFatal("No REXI coefficient loaded");
	}
	else if (shackExpIntegration->exp_method == "terry" || shackExpIntegration->exp_method == "ci")
	{
		sweet::ExpIntegration::REXI::REXICoefficients<> rexiCoefficients;

		sweet::ExpIntegration::REXI::REXI<> rexi;
		rexi.load(shackExpIntegration, function_name, rexiCoefficients, 0);

		rexiCoefficientsSet.rexiCoefficientVector.push_back(rexiCoefficients);
	}
	else
	{
		SWEETErrorFatal("This REXI method is not supported");
	}

	std::cout << "+ test_mode: " << shackREXITest->test_mode << std::endl;

	shackREXITest->printShack();

	for (std::size_t i = 0; i < rexiCoefficientsSet.rexiCoefficientVector.size(); i++)
	{
		const std::string &function_name = rexiCoefficientsSet.rexiCoefficientVector[i].function_name;

		std::cout << "Running tests for function " << function_name << std::endl;

		// Load coefficients for function
		const sweet::ExpIntegration::REXI::REXICoefficients<> *rexiCoeffs = rexiCoefficientsSet.getByFunctionName(function_name);

		// Initial condition
		std::complex<double> U0(1.0, 0.0);

		// Current solution
		std::complex<double> U = U0;

		auto computeAndOutputError = [&](std::complex<double> &i_U) -> double
		{
			TComplex analU;
			expFunctions.eval(shackREXITest->lambda*shackTimestepControl->currentSimulationTime, analU);

			double error = std::abs(analU-U);
			std::cout <<
					"t=" << shackTimestepControl->currentSimulationTime << "\t"
					"error=" << error
					<< std::endl;

			return error;
		};

		for (	shackTimestepControl->currentSimulationTime = 0;
				shackTimestepControl->currentSimulationTime < shackTimestepControl->maxSimulationTime*(1.0-1e-12);
				shackTimestepControl->currentSimulationTime += shackTimestepControl->currentTimestepSize
		)
		{
			computeAndOutputError(U);

			// REXI time integration
			{
				std::complex<double> approx = rexiCoeffs->gamma*U;

				for (std::size_t i = 0; i < rexiCoeffs->alphas.size(); i++)
					approx += rexiCoeffs->betas[i]/(shackREXITest->lambda*shackTimestepControl->currentTimestepSize + rexiCoeffs->alphas[i])*U;

				U = approx;
			}
		}

		double error = computeAndOutputError(U);

		std::cout << "[MULE] error: " << error << std::endl;
	}

	return 0;
}
