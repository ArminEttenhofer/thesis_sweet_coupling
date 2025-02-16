/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *      
 * MULE_SCONS_OPTIONS: --cart2d-spectral-space=enable
 */

#include <sweet/Tools/DefaultPrecompilerValues.hpp>

#if SWEET_GUI
#	error	"GUI not supported"
#endif


#include <sweet/Error/Base.hpp>
#include <sweet/Data/Cart2D/Cart2D.hpp>
#include <sweet/Data/Cart2D/DataSampler.hpp>
#include <sweet/Shacks/ProgramArgumentsDictionary.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>
#include <sweet/Tools/ProgramArguments.hpp>


#include <sweet/Data/Cart2DComplex/DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/DataGrid.hpp>
#include <sweet/Data/Cart2D/Convert/DataSpectral_2_Cart2DComplex_DataSpectral.hpp>
#include <sweet/Data/Cart2DComplex/Convert/DataSpectral_2_Cart2D_DataSpectral.hpp>
#include <sweet/Data/Cart2D/Convert/DataGrid_2_Cart2DComlex_DataGrid.hpp>
#include <sweet/Data/Cart2DComplex/Convert/DataGrid_2_Cart2D_DataGrid.hpp>


class Core_cart2dSamplerInterpolation
{
public:
	sweet::Error::Base error;


	/*
	 * Just a class to store simulation data all together
	 */
	class Data
	{
	public:
		sweet::Error::Base error;

		sweet::Data::Cart2D::Config cart2DDataConfig;
		sweet::Data::Cart2D::Config cart2DDataConfigOversampling;
		sweet::Data::Cart2D::Operators ops;

		sweet::Data::Cart2D::DataSpectral prog_h;

		sweet::Data::Vector::Vector<double> posx_a, posy_a;

		double *gaussianCenter;
		double gaussianExpFac = 50.0;

		sweet::Data::Cart2D::DataSampler cart2DDataSampler;


		bool setup(
				sweet::Data::Cart2D::Shack *i_shackCart2DDataOps,
				sweet::Data::Cart2D::Shack *i_shackCart2DDataOpsOversampling,
				double i_gaussianCenter[]
		)
		{
			gaussianCenter = i_gaussianCenter;

			/*
			 * Setup Cart2D Data Config & Operators
			 */
			cart2DDataConfig.setupAuto(*i_shackCart2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2DDataConfig);

			ops.setup(cart2DDataConfig, *i_shackCart2DDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

			cart2DDataConfigOversampling.setupAuto(*i_shackCart2DDataOpsOversampling);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(cart2DDataConfigOversampling);

			prog_h.setup(cart2DDataConfig);

			sweet::Data::Cart2D::DataGrid prog_h_phys(cart2DDataConfig);

			prog_h_phys.grid_update_lambda_array_indices(
					[&](int i, int j, double &io_data)
				{
					double x = (double)i*(i_shackCart2DDataOps->cart2d_domain_size[0]/(double)i_shackCart2DDataOps->space_res_physical[0]);
					double y = (double)j*(i_shackCart2DDataOps->cart2d_domain_size[1]/(double)i_shackCart2DDataOps->space_res_physical[1]);

					io_data = gaussianValue(
							i_shackCart2DDataOps,
							gaussianCenter,
							x,
							y,
							gaussianExpFac
						);
				}
			);

			prog_h.loadCart2DDataGrid(prog_h_phys);

			posx_a.setup(cart2DDataConfigOversampling.grid_number_elements);
			posy_a.setup(cart2DDataConfigOversampling.grid_number_elements);

			// setup some test sampling points
			// we use 2 arrays - one for each sampling position

			posx_a.update_lambda_array_indices(
				[&](int idx, double &io_data)
				{
					int i = idx % cart2DDataConfigOversampling.grid_res[0];
					//int j = idx / cart2DDataConfig->grid_data_size[0];

					io_data = i_shackCart2DDataOps->cart2d_domain_size[0]*(double)i/(double)cart2DDataConfigOversampling.grid_res[0];
					SWEET_ASSERT(io_data >= 0);
					SWEET_ASSERT(io_data < i_shackCart2DDataOps->cart2d_domain_size[0]);
				}
			);
			posy_a.update_lambda_array_indices(
					[&](int idx, double &io_data)
				{
					//int i = idx % cart2DDataConfig->grid_data_size[0];
					int j = idx / cart2DDataConfigOversampling.grid_res[0];

					io_data = i_shackCart2DDataOps->cart2d_domain_size[1]*(double)j/(double)cart2DDataConfigOversampling.grid_res[1];

					SWEET_ASSERT(io_data >= -M_PI*0.5);
					SWEET_ASSERT(io_data < i_shackCart2DDataOps->cart2d_domain_size[1]);
				}
			);

			cart2DDataSampler.setup(i_shackCart2DDataOps->cart2d_domain_size, &cart2DDataConfig);

			return true;
		}

		void clear()
		{
			prog_h.clear();

			ops.clear();
			cart2DDataConfig.clear();
		}
	};

	// Simulation data
	Data data;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict;
	sweet::Data::Cart2D::Shack *shackCart2DDataOps;

	int interpolation_order = 3;

	double max_error;

public:
	Core_cart2dSamplerInterpolation(
			int i_argc,
			char *const * const i_argv
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackCart2DDataOps(nullptr)
	{
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN(shackProgArgDict);
	}

	bool setup(
			double i_gaussianCenter[],
			int i_interpolationOrder,
			int i_specModes[2],
			int i_specModesOversampled[2]
	)
	{
		/*
		 * SHACK: Register classes which we require
		 */
		shackCart2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Cart2D::Shack>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.setup();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

#if 0
		shackProgArgDict.printShackData();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);
#endif

		shackCart2DDataOps->space_res_spectral[0] = i_specModes[0];
		shackCart2DDataOps->space_res_spectral[1] = i_specModes[1];

		sweet::Data::Cart2D::Shack shackCart2DDataOpsOversampling = *shackCart2DDataOps;

		shackCart2DDataOpsOversampling.space_res_spectral[0] = i_specModesOversampled[0];
		shackCart2DDataOpsOversampling.space_res_spectral[1] = i_specModesOversampled[1];

#if 0
		std::cout << "Shack Data regular:" << std::endl;
		shackCart2DDataOps->printShack("  ");

		std::cout << "Shack Data oversampled:" << std::endl;
		shackCart2DDataOpsOversampling.printShack("  ");
#endif
		data.setup(
				shackCart2DDataOps,
				&shackCart2DDataOpsOversampling,
				i_gaussianCenter
		);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(data);

		return true;
	}


	void clear()
	{
		data.clear();

		shackCart2DDataOps = nullptr;
		shackProgArgDict.clear();
	}


	static
	double gaussianValue(
			sweet::Data::Cart2D::Shack *i_shackCart2DDataOps,
			double *i_center,
			double i_x, double i_y,
			double i_exp_fac
	)
	{
		double sx = i_shackCart2DDataOps->cart2d_domain_size[0];
		double sy = i_shackCart2DDataOps->cart2d_domain_size[1];

		// Gaussian
		double dx = i_x-i_center[0]*sx;
		double dy = i_y-i_center[1]*sy;

		if (dx > 0.5*i_shackCart2DDataOps->cart2d_domain_size[0])
			dx -= i_shackCart2DDataOps->cart2d_domain_size[0];
		else if (dx < -0.5*i_shackCart2DDataOps->cart2d_domain_size[0])
			dx += i_shackCart2DDataOps->cart2d_domain_size[0];

		if (dy > 0.5*i_shackCart2DDataOps->cart2d_domain_size[1])
			dy -= i_shackCart2DDataOps->cart2d_domain_size[1];
		else if (dy < -0.5*i_shackCart2DDataOps->cart2d_domain_size[1])
			dy += i_shackCart2DDataOps->cart2d_domain_size[1];

		dx /= sx;
		dy /= sy;

		return std::exp(-i_exp_fac*(dx*dx + dy*dy));
	}



	void run_tests()
	{
		sweet::Data::Vector::Vector<double> out_data;
		out_data.setup(data.posx_a.numberOfElements);

		if (interpolation_order == 2)
		{
			data.cart2DDataSampler.bilinear_scalar(
					data.prog_h.toGrid(),
					data.posx_a,
					data.posy_a,
					out_data
			);
		}
		else if (interpolation_order == 3)
		{
			data.cart2DDataSampler.bicubic_scalar(
					data.prog_h.toGrid(),
					data.posx_a,
					data.posy_a,
					out_data
			);
		}
		else
		{
			SWEETErrorFatal("Interpolation order not available");
		}

		/*
		 * Compute errors
		 */
		max_error = 0;
		for (std::size_t i = 0; i < data.posx_a.numberOfElements; i++)
		{
			double value = gaussianValue(
					shackCart2DDataOps,
					data.gaussianCenter,
					data.posx_a.data[i],
					data.posy_a.data[i],
					data.gaussianExpFac
				);

			max_error = std::max(max_error, std::abs(value - out_data.data[i]));
		}
	}


	bool should_quit()
	{
		return false;
	}
};




int main(int i_argc, char *i_argv[])
{
	/*
	 * SHACK: Register classes which we require
	 */

	sweet::Shacks::ProgramArgumentsDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	sweet::Data::Cart2D::Shack *shackCart2DDataOps = shackProgArgDict.getAutoRegistration<sweet::Data::Cart2D::Shack>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);


	int initialSpectralModes = shackCart2DDataOps->space_res_spectral[0];

	double gaussianCenterArray[6][2] = {
			{0.5, 0.5},
			{0.9, 0.4},
			{0.4, 0.9},
			{0.3, 0.8},
			{0.7, 0.3},
			{0.3, 0.7},
	};

	for (int i_gaussians = 2; i_gaussians < 6; i_gaussians++)
	{
		std::cout << "*********************************************************" << std::endl;
		std::cout << "* Running studies for Gaussian at " << gaussianCenterArray[i_gaussians][0] << ", " << gaussianCenterArray[i_gaussians][1] << std::endl;
		std::cout << "*********************************************************" << std::endl;

		for (int interpolation_order = 2; interpolation_order <= 3; interpolation_order++)
		{
			std::cout << std::endl;
			std::cout << "*********************************************************" << std::endl;
			std::cout << "* Running studies for interpolation of order " << interpolation_order << std::endl;
			std::cout << "*********************************************************" << std::endl;

			int oversamplingFactor = 13;
			std::cout << "Using oversampling of " << oversamplingFactor << std::endl;

			double prev_max_error = -1;
			for (int specModes = initialSpectralModes; specModes <= 256; specModes *= 2)
			{

				Core_cart2dSamplerInterpolation simulation(i_argc, i_argv);
				ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

				int specModes_[2] = {specModes, specModes};
				int specModesOversampled_[2] = {specModes*oversamplingFactor, specModes*oversamplingFactor};

				simulation.setup(
						gaussianCenterArray[i_gaussians],
						interpolation_order,
						specModes_,
						specModesOversampled_
					);
				ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

				simulation.run_tests();
				ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

				std::cout << "Error: " << simulation.max_error << std::endl;
				if (prev_max_error >= 0)
				{
					if (std::isnan(simulation.max_error))
						SWEETErrorFatal("NaN detected");

					//double conv = (prev_max_error - simulation.max_error) / simulation.max_error;
					double conv = prev_max_error / simulation.max_error;
					std::cout << "Convergence: " << conv << std::endl;

					if (conv*1.1 < std::pow(2.0, interpolation_order))
						SWEETErrorFatal("Convergence not given!");
				}
				prev_max_error = simulation.max_error;

				simulation.clear();
				ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(simulation);

			}
		}
	}


	std::cout << "FIN" << std::endl;
	return 0;
}
