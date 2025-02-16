/*
 * PInT_Common.hpp
 *
 *  Created on: 10 Jun 2022
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 */

#include <sweet/Error/Base.hpp>
#include <sweet/IO/Shack.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>

/*
 * Contains functions used by both Parareal and Xbraid (and possibly PFASST?) implementations in SWEET,
 * mostly error computation and file output functions.
 */

#ifndef INCLUDE_SWEET__DEPRECATED_PINT_PINT_COMMON_HPP
#define INCLUDE_SWEET__DEPRECATED_PINT_PINT_COMMON_HPP

/////#include <sweet/core/SimulationVariables.hpp>
#include <sweet/_DEPRECATED_pint/Parareal_GenericData.hpp>
#if SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR
	#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Scalar.hpp>
#elif SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
	#include <sweet/Data/Cart2D/Operators.hpp>
	#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Cart2DData_Spectral.hpp>
	#include <sweet/Data/Cart2D/Shack.hpp>
	#include <sweet/Data/Cart2D/GridMapping.hpp>
#elif SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D
	#include <sweet/_DEPRECATED_pint/Parareal_GenericData_Sphere2DData_Spectral.hpp>
	#include <sweet/Data/Sphere2D/Shack.hpp>
	#include <programs/PDE_SWESphere2D/Shack.hpp>
#endif

#include <map>

namespace sweet {
namespace DEPRECATED_pint {

class PInT_Common
{

protected:
	///SimulationVariables* simVars = nullptr;

/////////	sweet::Error::Base error;
/////////
/////////	sweet::IO::Shack* shackIOData;
/////////
/////////#if SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
/////////	// Grid Mapping (staggered grid)
/////////	sweet::Data::Cart2D::GridMapping gridMapping;
/////////#endif
/////////
/////////	// Operators and DataConfig
/////////#if SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
/////////	sweet::Data::Cart2D::Operators* base_op_cart2d;
/////////	sweet::Data::Cart2D::Config* base_cart2DDataConfig;
/////////	std::vector<sweet::Data::Cart2D::Operators*> op_cart2d;
/////////	std::vector<sweet::Data::Cart2D::Config*> cart2DDataConfig;
/////////	sweet::Data::Cart2D::Shack* shackCart2DDataOps;
/////////#elif SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D
/////////	sweet::Data::Sphere2D::Operators* base_op_sphere2D;
/////////	sweet::Data::Sphere2D::Operators* base_op_sphere2d_nodealiasing;
/////////	sweet::Data::Sphere2D::Config* base_sphere2DDataConfig;
/////////	std::vector<sweet::Data::Sphere2D::Operators*> op_sphere2D;
/////////	std::vector<sweet::Data::Sphere2D::Operators*> op_sphere2d_nodealiasing;
/////////	std::vector<sweet::Data::Sphere2D::Config*> sphere2DDataConfig;
/////////	sweet::Data::Sphere2D::Shack* shackSphere2DDataOps;
/////////	PDE_SWESphere2D::Shack* shackPDESWESphere2D;
/////////#endif

/////////	// list of SL schemes
/////////	std::vector<std::string> SL_tsm = {};


#if SWEET_PARAREAL_CART2D_BURGERS
	// required for computing analytical solution
	class BenchmarkErrors
	{
	public:
		// Max difference to initial conditions
		double benchmark_diff_u;
		double benchmark_diff_v;

		// Error measures L2 norm
		double benchmark_analytical_error_rms_u;
		double benchmark_analytical_error_rms_v;

		// Error measures max norm
		double benchmark_analytical_error_maxabs_u;
		double benchmark_analytical_error_maxabs_v;
	};
	Burgers_Cart2D_TimeSteppers* timeSteppersFineBurgers = nullptr;
	BenchmarkErrors benchmark;

	void set_tsm_burgers(
				Burgers_Cart2D_TimeSteppers* i_timeSteppersFineBurgers
	)
	{
		timeSteppersFineBurgers = i_timeSteppersFineBurgers;
	}
#endif

public:
	PInT_Common()
	{
	}

	~PInT_Common()
	{
	}


public:

	void setup(
			///SimulationVariables* i_simVars
			sweet::IO::Shack* i_shackIOData,
			sweet::Shacks::Base* i_shackDataOps
/////#if SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
/////			,
/////			sweet::Data::Cart2D::Shack* i_shackCart2DDataOps
/////#elif SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D
/////			,
/////			sweet::Data::Sphere2D::Shack* i_shackSphere2DDataOps
/////#endif
	)
	{
		////simVars = i_simVars;

/////////		shackIOData = i_shackIOData;
/////////#if SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
/////////		shackCart2DDataOps = i_shackCart2DDataOps;
/////////#elif SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D
/////////		shackSphere2DDataOps = i_shackSphere2DDataOps;
/////////#endif

/////////	#if SWEET_PARAREAL_CART2D_SWE || SWEET_XBRAID_CART2D_SWE
/////////		SL_tsm = { "l_cn_na_sl_nd_settls",
/////////				 "l_rexi_na_sl_nd_etdrk",
/////////				 "l_rexi_na_sl_nd_settls"
/////////				};
/////////	#elif SWEET_PARAREAL_CART2D_BURGERS || SWEET_XBRAID_CART2D_BURGERS
/////////		SL_tsm = { "l_cn_n_sl",
/////////				 "l_irk_n_sl",
/////////				 "l_irk_n_sl_forcing"
/////////				};
/////////	#elif SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D
/////////		SL_tsm = { "lg_exp_na_sl_lc_nr_etd_uv",
/////////				 "l_irk_na_sl_nr_settls_uv_only",
/////////				 "l_irk_na_sl_nr_settls_vd_only",
/////////				 "l_irk_na_sl_settls_uv_only",
/////////				 "l_irk_na_sl_settls_vd_only",
/////////				 "ln_sl_exp_settls_uv",
/////////				 "ln_sl_exp_settls_vd",
/////////				 "lg_exp_na_sl_lc_nr_etdrk_uv"
/////////				};
/////////	#endif


	}

	void output_residual_file(
			double res,
			int iteration_id
	)
	{

		char buffer[1024];

		const char* filename_template = "residual_iter%03d.csv";
		sprintf(buffer, filename_template, iteration_id);

		std::ofstream file(buffer, std::ios_base::trunc);
		file << std::setprecision(16);

		file << "#SWEET" << std::endl;
		file << "#FORMAT ASCII" << std::endl;
		file << "#PRIMITIVE SCALAR" << std::endl;

		file << res;

		file.close();

	}


	void output_data_file(
			sweet::DEPRECATED_pint::Parareal_GenericData* i_data,
			int iteration_id,
			int time_slice_id,
			double t
	)
	{
#if SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR
		double u_out;
		i_data->GenericData_Scalar_2_dataArrays(u_out);

		// Dump  data in csv, if output filename is not empty
		if (shackIOData->outputFileName.size() > 0)
		{
			std::string output_filenames = "";
			output_filenames = write_file_pint_scalar(u_out, "prog_u", iteration_id, t);
		}

#elif SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
	#if SWEET_PARAREAL_CART2D_BURGERS || SWEET_XBRAID_CART2D_BURGERS

		sweet::Data::Cart2D::DataSpectral dummy(cart2DDataConfig[0]);
		sweet::Data::Cart2D::DataSpectral u_out(cart2DDataConfig[0]);
		sweet::Data::Cart2D::DataSpectral v_out(cart2DDataConfig[0]);
		i_data->GenericData_Cart2DData_Spectral_2_dataArrays(u_out, v_out);

		sweet::Data::Cart2D::DataGrid u_out_phys = u_out.toGrid();
		sweet::Data::Cart2D::DataGrid v_out_phys = v_out.toGrid();

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
		// For output, variables need to be on unstaggered A-grid
		sweet::Data::Cart2D::DataGrid t_u(cart2DDataConfig[0]);
		sweet::Data::Cart2D::DataGrid t_v(cart2DDataConfig[0]);

		if (shackCart2DDataOps->space_grid_use_c_staggering) // Remap in case of C-grid
		{
			gridMapping.mapCtoA_u(u_out_phys, t_u);
			gridMapping.mapCtoA_v(v_out_phys, t_v);
		}
		else
		{
			t_u = u_out_phys;
			t_v = v_out_phys;
		}

		// Dump  data in csv, if output filename is not empty
		if (shackIOData->outputFileName.size() > 0)
		{
			std::string output_filenames = "";

			output_filenames = write_file_pint_cart2d(t_u, "prog_u", iteration_id, t);
			output_filenames += ";" + write_file_pint_cart2d(t_v, "prog_v", iteration_id, t);

			output_filenames += ";" + write_file_spec_pint_cart2d(u_out, "prog_u_spec", iteration_id, t);
			output_filenames += ";" + write_file_spec_pint_cart2d(v_out, "prog_v_spec", iteration_id, t);

		}

		write_file_spec_amp_phase_pint_cart2d(u_out, "prog_u", iteration_id, t);

		if (simVars->misc.compute_errors)
		{
			sweet::Data::Cart2D::DataSpectral ana = compute_errors2(u_out, v_out);

			write_file_pint_cart2d(ana.toGrid(),"analytical",iteration_id, t);
			write_file_spec_amp_phase_pint_cart2d(ana.toGrid(), "analytical", iteration_id, t);
		}

	#elif SWEET_PARAREAL_CART2D_SWE || SWEET_XBRAID_CART2D_SWE

		sweet::Data::Cart2D::DataSpectral h_out(cart2DDataConfig[0]);
		sweet::Data::Cart2D::DataSpectral u_out(cart2DDataConfig[0]);
		sweet::Data::Cart2D::DataSpectral v_out(cart2DDataConfig[0]);
		i_data->GenericData_Cart2DData_Spectral_2_dataArrays(h_out, u_out, v_out);

		sweet::Data::Cart2D::DataGrid h_out_phys = h_out.toGrid();
		sweet::Data::Cart2D::DataGrid u_out_phys = u_out.toGrid();
		sweet::Data::Cart2D::DataGrid v_out_phys = v_out.toGrid();

		// Save .vtk files for visualizing in paraview

#if 0
		// @ Joao: I (Martin) kicked this out for cleaning up the code

		std::ostringstream ss2;
		ss2 << "output_slice" << time_slice_id << "_iter" << iteration_id << ".vtk";
		std::string filename2 = ss2.str();
		h_out_phys.file_grid_saveData_vtk(filename2.c_str(), filename2.c_str());
#endif

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
		// For output, variables need to be on unstaggered A-grid
		sweet::Data::Cart2D::DataGrid t_h(cart2DDataConfig[0]);
		sweet::Data::Cart2D::DataGrid t_u(cart2DDataConfig[0]);
		sweet::Data::Cart2D::DataGrid t_v(cart2DDataConfig[0]);

		if (shackCart2DDataOps->space_grid_use_c_staggering) // Remap in case of C-grid
		{
			t_h = h_out_phys;
			gridMapping.mapCtoA_u(u_out_phys, t_u);
			gridMapping.mapCtoA_v(v_out_phys, t_v);
		}
		else
		{
			t_h = h_out_phys;
			t_u = u_out_phys;
			t_v = v_out_phys;
		}

		// Dump  data in csv, if output filename is not empty
		if (shackIOData->outputFileName.size() > 0)
		{
			std::string output_filenames = "";

			output_filenames = write_file_pint_cart2d(t_h, "prog_h_pert", iteration_id, t);
			output_filenames += ";" + write_file_pint_cart2d(t_u, "prog_u", iteration_id, t);
			output_filenames += ";" + write_file_pint_cart2d(t_v, "prog_v", iteration_id, t);

			output_filenames += ";" + write_file_pint_cart2d(op_cart2d[0]->ke(t_u,t_v).toGrid(),"diag_ke", iteration_id, t);

			output_filenames += ";" + write_file_spec_pint_cart2d(h_out, "prog_h_pert_spec", iteration_id, t);
			output_filenames += ";" + write_file_spec_pint_cart2d(u_out, "prog_u_spec", iteration_id, t);
			output_filenames += ";" + write_file_spec_pint_cart2d(v_out, "prog_v_spec", iteration_id, t);

			output_filenames += ";" + write_file_spec_pint_cart2d(op_cart2d[0]->ke(t_u,t_v).toGrid(),"diag_ke_spec", iteration_id, t);

			output_filenames += ";" + write_file_pint_cart2d(op_cart2d[0]->vort(t_u, t_v).toGrid(), "diag_vort", iteration_id, t);
			output_filenames += ";" + write_file_pint_cart2d(op_cart2d[0]->div(t_u, t_v).toGrid(), "diag_div", iteration_id, t);

			/////////if (compute_normal_modes){
			/////////	SWEETErrorFatal("TODO");
			/////////	///output_filenames += ";" + write_file_spec_parareal(normalmodes.geo, "nm_geo", iteration_id, output_initial_data);
			/////////	///output_filenames += ";" + write_file_spec_parareal(normalmodes.igwest, "nm_igwest", iteration_id, output_initial_data);
			/////////	///output_filenames += ";" + write_file_spec_parareal(normalmodes.igeast, "nm_igeast", iteration_id, output_initial_data);
			/////////}
			
		}
	#endif


#elif SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D

		sweet::Data::Sphere2D::DataSpectral phi_out(sphere2DDataConfig[0]);
		sweet::Data::Sphere2D::DataSpectral vrt_out(sphere2DDataConfig[0]);
		sweet::Data::Sphere2D::DataSpectral div_out(sphere2DDataConfig[0]);
		i_data->GenericData_Sphere2DData_Spectral_to_dataArrays(phi_out, vrt_out, div_out);

		sweet::Data::Sphere2D::DataGrid phi_out_phys = phi_out.toGrid();
		sweet::Data::Sphere2D::DataGrid vrt_out_phys = vrt_out.toGrid();
		sweet::Data::Sphere2D::DataGrid div_out_phys = div_out.toGrid();

		////////////////! Save .vtk files for visualizing in paraview
		///////////////std::ostringstream ss2;
		///////////////if (output_initial_data)
		///////////////	ss2 << "output_slice" << time_slice_id - 1 << "_iter" << iteration_id << ".vtk";
		///////////////else
		///////////////	ss2 << "output_slice" << time_slice_id << "_iter" << iteration_id << ".vtk";
		///////////////std::string filename2 = ss2.str();
		///////////////phi_out_phys.file_grid_saveData_vtk(filename2.c_str(), filename2.c_str());

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
		// Dump  data in csv, if output filename is not empty
		if (shackIOData->outputFileName.size() > 0)
		{
			if (shackIOData->outputFileMode == "csv")
			{
				std::string output_filename;

				sweet::Data::Sphere2D::DataSpectral h = phi_out_phys*(1.0/shackPDESWESphere2D->gravitation);
				h += shackPDESWESphere2D->h0;

				output_filename = write_file_csv_pint_sphere2D(h, t, "prog_h", iteration_id);
				//std::cout << " + " << output_filename << " (min: " << h.toGrid().grid_reduce_min() << ", max: " << h.toGrid().grid_reduce_max() << ")" << std::endl;

				sweet::Data::Sphere2D::DataGrid phi_phys = h.toGrid() * shackPDESWESphere2D->gravitation;
				sweet::Data::Sphere2D::DataSpectral phi(sphere2DDataConfig[0]);
				phi.loadSphere2DDataGrid(phi_phys);
				output_filename = write_file_csv_pint_sphere2D(phi, t, "prog_phi", iteration_id);

				output_filename = write_file_csv_pint_sphere2D(phi_out, t, "prog_phi_pert", iteration_id);
				//std::cout << " + " << output_filename << " (min: " << phi_out_phys.grid_reduce_min() << ", max: " << phi_out_phys.grid_reduce_max() << ")" << std::endl;
	
				sweet::Data::Sphere2D::DataGrid u(sphere2DDataConfig[0]);
				sweet::Data::Sphere2D::DataGrid v(sphere2DDataConfig[0]);
	
				op_sphere2D[0]->vrtdiv_2_uv(vrt_out_phys, div_out_phys, u, v);
	
				output_filename = write_file_csv_pint_sphere2D(u, t, "prog_u", iteration_id);
				//std::cout << " + " << output_filename << std::endl;
	
				output_filename = write_file_csv_pint_sphere2D(v, t, "prog_v", iteration_id);
				//std::cout << " + " << output_filename << std::endl;
	
				output_filename = write_file_csv_pint_sphere2D(vrt_out, t, "prog_vrt", iteration_id);
				//std::cout << " + " << output_filename << std::endl;
	
				output_filename = write_file_csv_pint_sphere2D(div_out, t, "prog_div", iteration_id);
				//std::cout << " + " << output_filename << std::endl;
	
				sweet::Data::Sphere2D::DataSpectral potvrt = (phi_out / shackPDESWESphere2D->gravitation)*vrt_out;
	
				output_filename = write_file_csv_pint_sphere2D(potvrt, t, "prog_potvrt", iteration_id);
				//std::cout << " + " << output_filename << std::endl;


				////output_filename = write_file_csv_pint_sphere2d_spec(phi_out, t, "prog_phi_pert", iteration_id);

			}
			else if (shackIOData->outputFileMode == "bin")
			{
				std::string output_filename;
	
				{
					output_filename = write_file_bin_pint_sphere2D(phi_out, t, "prog_phi_pert", iteration_id);
					sweet::Data::Sphere2D::DataGrid prog_phys = phi_out.toGrid();
	
					//std::cout << " + " << output_filename << " (min: " << prog_phys.grid_reduce_min() << ", max: " << prog_phys.grid_reduce_max() << ")" << std::endl;
				}
	
				{
					output_filename = write_file_bin_pint_sphere2D(vrt_out, t, "prog_vrt", iteration_id);
					sweet::Data::Sphere2D::DataGrid prog_phys = vrt_out.toGrid();
	
					//std::cout << " + " << output_filename << " (min: " << prog_phys.grid_reduce_min() << ", max: " << prog_phys.grid_reduce_max() << ")" << std::endl;
				}
	
				{
					output_filename = write_file_bin_pint_sphere2D(div_out, t, "prog_div", iteration_id);
					sweet::Data::Sphere2D::DataGrid prog_phys = div_out.toGrid();
	
					//std::cout << " + " << output_filename << " (min: " << prog_phys.grid_reduce_min() << ", max: " << prog_phys.grid_reduce_max() << ")" << std::endl;
				}
			}
			
		}


#endif
	};


#if SWEET_PARAREAL_CART2D_BURGERS || SWEET_XBRAID_CART2D_BURGERS
	// For Burgers
	sweet:Cart2DData_Spectral compute_errors2(
						const sweet::Data::Cart2D::DataSpectral &i_cart2DData_u,
						const sweet::Data::Cart2D::DataSpectral &i_cart2DData_v
	)
	{

		int analytic_solution;
		if (simVars->misc.compute_errors)
		{
			bool foundl = (simVars->disc.timestepping_method.find("l_")==0) || (simVars->disc.timestepping_method.find("_l_")!=std::string::npos);
			bool foundn = (simVars->disc.timestepping_method.find("n_")==0) || (simVars->disc.timestepping_method.find("_n_")!=std::string::npos);
			bool foundnl = (simVars->disc.timestepping_method.find("ln_")==0) || (foundl && foundn);
		
			if (foundnl)
				analytic_solution = 1;
			else if (foundl)
				analytic_solution = 2;
			else
				SWEETErrorFatal("Computing errors for this timestepping-method is not possible");
		}



		// Necessary to circumvent FFTW transformations on i_cart2DData_u and i_cart2DData_v, which would lead to errors
		sweet::Data::Cart2D::DataGrid u = i_cart2DData_u.toGrid();
		sweet::Data::Cart2D::DataGrid v = i_cart2DData_v.toGrid();

		////! Analytical solution at current time on original grid
		///Cart2DData_Spectral ts_u = t0_prog_u;
		///Cart2DData_Spectral ts_v = t0_prog_v;

		sweet::Data::Cart2D::DataSpectral ts_u(cart2DDataConfig[0]);
		sweet::Data::Cart2D::DataSpectral ts_v(cart2DDataConfig[0]);
		sweet::Data::Cart2D::DataGrid ts_u_phys(cart2DDataConfig[0]);
		sweet::Data::Cart2D::DataGrid ts_v_phys(cart2DDataConfig[0]);

		if (simVars->misc.compute_errors)
		{
			//if (simVars.setup.benchmark_id > 51 && simVars.setup.benchmark_id < 65)
			if (simVars->disc.timestepping_method.find("forcing")!=std::string::npos)
			{
				if (simVars->disc.space_grid_use_c_staggering)
				{
					ts_u_phys.grid_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							double x = (((double)i)/(double)simVars->disc.space_res_physical[0])*simVars->sim.cart2d_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars->disc.space_res_physical[1])*simVars->sim.cart2d_domain_size[1];
							io_data = BurgersValidationBenchmarks::return_u(*simVars, x, y);
						}
					);

					ts_v_phys.grid_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							io_data = 0.0;
#if 0
							double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.cart2d_domain_size[0];
							double y = (((double)j)/(double)simVars.disc.space_res_physical[1])*simVars.sim.cart2d_domain_size[1];
							io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
#endif
						}
					);
				}
				else
				{
					ts_u_phys.grid_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							double x = (((double)i+0.5)/(double)simVars->disc.space_res_physical[0])*simVars->sim.cart2d_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars->disc.space_res_physical[1])*simVars->sim.cart2d_domain_size[1];

							io_data = BurgersValidationBenchmarks::return_u(*simVars, x, y);
						}
					);

					ts_v_phys.grid_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							io_data = 0.0;
#if 0
							double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.cart2d_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.cart2d_domain_size[1];

							io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
#endif
						}
					);
				}
				ts_u.loadCart2DDataGrid(ts_u_phys);
				ts_v.loadCart2DDataGrid(ts_v_phys);
			}
			else //if (simVars.setup.benchmark_id == 70)
			{
				if (analytic_solution == 1)
				{
				   timeSteppersFineBurgers->ln_cole_hopf->runTimestep(
						 ts_u, ts_v,
						 //ts_u, ts_v,
						 simVars->timecontrol.current_simulation_time,
						 0
				   );
				}
				else if (analytic_solution == 2)
				{
				   timeSteppersFineBurgers->l_direct->runTimestep(
						 ts_u, ts_v,
						 //ts_u, ts_v,
						 simVars->timecontrol.current_simulation_time,
						 0
				   );
				}
			}
			benchmark.benchmark_analytical_error_rms_u = (ts_u-u).toGrid().grid_reduce_rms();
			benchmark.benchmark_analytical_error_rms_v = (ts_v-v).toGrid().grid_reduce_rms();

			benchmark.benchmark_analytical_error_maxabs_u = (ts_u-u).toGrid().grid_reduce_max_abs();
			benchmark.benchmark_analytical_error_maxabs_v = (ts_v-v).toGrid().grid_reduce_max_abs();

			return ts_u;
		}
		return nullptr;
	}

#endif


#if SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR
	/**
	 * Write file to data and return string of file name (parareal)
	 */
	std::string write_file_pint_scalar(
			const double &i_u,
			const char* i_name,	//!< name of output variable
			int iteration_id,
			double t
		)
	{
		char buffer[1024];


		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		sprintf(buffer, filename_template, i_name, t, iteration_id);

		std::ofstream file(buffer, std::ios_base::trunc);
		file << std::setprecision(16);

		file << "#SWEET" << std::endl;
		file << "#FORMAT ASCII" << std::endl;
		file << "#PRIMITIVE SCALAR" << std::endl;

		file << i_u;

		file.close();

		return buffer;
	}
#endif

#if SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D
	/**
	 * Write physical data to file and return string of file name (parareal)
	 */
	std::string write_file_csv_pint_sphere2D(
			const sweet::Data::Sphere2D::DataSpectral &i_sphere2DData,
			double t,
			const char* i_name,	//!< name of output variable
			int iteration_id,
			bool i_phi_shifted = false
		)
	{
		char buffer[1024];

		// create copy
		sweet::Data::Sphere2D::DataGrid sphere2DData = i_sphere2DData.toGrid();

		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		sprintf(buffer, filename_template, i_name, t * shackIOData->outputFormatTimeScale, iteration_id);

		if (i_phi_shifted)
			sphere2DData.grid_file_write_lon_pi_shifted(buffer, "vorticity, lon pi shifted");
		else
			sphere2DData.grid_file_write(buffer);

		return buffer;

	}

/////	/**
/////	 * Write spectral data to file and return string of file name (parareal)
/////	 */
/////	std::string write_file_csv_pint_sphere2d_spec(
/////			const sweet::Data::Sphere2D::DataSpectral &i_sphere2DData,
/////			double t,
/////			const char* i_name,	//!< name of output variable
/////			int iteration_id
/////		)
/////	{
/////		char buffer[1024];
/////
/////		const char* filename_template = "output_spec_%s_t%020.8f_iter%03d.csv";
/////		sprintf(buffer, filename_template, i_name, t * simVars->iodata.output_time_scale, iteration_id);
/////		i_sphere2DData.file_spectral_saveData_ascii(buffer);
/////		return buffer;
/////
/////	}


	/**
	 * Write spectral data to file and return string of file name
	 */
	std::string write_file_bin_pint_sphere2D(
			const sweet::Data::Sphere2D::DataSpectral &i_sphere2DData,
			double t,
			const char* i_name,
			int iteration_id
	)
	{
		char buffer[1024];

		sweet::Data::Sphere2D::DataSpectral sphere2DData(i_sphere2DData);
		//const char* filename_template = simVars.iodata.output_file_name.c_str();
		const char* filename_template = "output_%s_t%020.8f_iter%03d.sweet";
		sprintf(buffer, filename_template, i_name, t * shackIOData->outputFormatTimeScale, iteration_id);
		sphere2DData.file_write_binary_spectral(buffer);

		return buffer;
	}
#endif


#if SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
	/**
	 * Write file to data and return string of file name (parareal)
	 */
	std::string write_file_pint_cart2d(
			const sweet::Data::Cart2D::DataGrid &i_cart2DData,
			const char* i_name,	//!< name of output variable
			int iteration_id,
			double t
		)
	{
		char buffer[1024];

		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		sprintf(buffer, filename_template, i_name, t, iteration_id);
		i_cart2DData.file_grid_saveData_ascii(buffer);
		return buffer;
	}

	/**
	 * Write spectrum info to data and return string of file name (parareal)
	 */
	std::string write_file_spec_pint_cart2d(
			const sweet::Data::Cart2D::DataSpectral &i_cart2DData,
			const char* i_name,	//!< name of output variable
			int iteration_id,
			double t
		)
	{
		char buffer[1024];

		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		sprintf(buffer, filename_template, i_name, t, iteration_id);
		i_cart2DData.file_spectral_saveData_ascii(buffer);
		///i_cart2DData.file_spectral_abs_saveData_ascii(buffer);
		return buffer;
	}

	/**
	 * Write spectrum info to data and return string of file name (parareal)
	 */
	std::string write_file_spec_amp_phase_pint_cart2d(
			const sweet::Data::Cart2D::DataSpectral &i_cart2DData,
			const char* i_name,	//!< name of output variable
			int iteration_id,
			double t
		)
	{

		char buffer[1024];

		const char* filename_template = "output_%s_amp_phase_t%020.8f_iter%03d.csv";
		sprintf(buffer, filename_template, i_name, t, iteration_id);

		std::ofstream file(buffer, std::ios_base::trunc);
		file << std::setprecision(12);

		for (std::size_t x = 0; x < cart2DDataConfig[0]->spectral_data_size[0]; x++)
		{
			file << x << ", " << i_cart2DData.spectral_return_amplitude(0,x) << ", " << i_cart2DData.spectral_return_phase(0,x) << std::endl;
		}
		file.close();
		file.clear();

		return buffer;
	}

#endif

	/**
	 * Compute and store parareal errors during simulation
	 */
	void store_pint_error(
			sweet::DEPRECATED_pint::Parareal_GenericData* i_data,
			sweet::DEPRECATED_pint::Parareal_GenericData* pint_data_ref,
			int nvar,
			int iteration_id,
			int time_slice_id,
			double t,
			std::string path_ref,
			std::string base_solution,	// "ref" or "fine"
			std::string pint_type,
			int i_precision = 32
	)
	{

#if SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR
		if (iteration_id == 0)
		{
			// load ref file
			char buffer[1024];
			std::string i_name = "prog_u";
			const char* filename_template = shackIOData->outputFileName.c_str();
			///sprintf(buffer, filename_template, i_name.c_str(), timeframe_end);
			sprintf(buffer, filename_template, i_name.c_str(), t);
			std::string buffer2 = path_ref + "/" + std::string(buffer);

			double tmp;

			std::cout << path_ref << std::endl;
			std::cout << "loading DATA from " << buffer2 << std::endl;
			std::ifstream file(buffer2);
			for (int i = 0; i < 4; i++)
			{
				std::string line;
				std::getline(file, line);
				std::istringstream iss(line);
				std::vector<std::string> str_vector((std::istream_iterator<std::string>(iss)),
					std::istream_iterator<std::string>());

				if (i == 0)
				{
					SWEET_ASSERT(str_vector.size() == 1);
					SWEET_ASSERT(str_vector[0] == "#SWEET");
				}
				else if (i == 1)
				{
					SWEET_ASSERT(str_vector.size() == 2);
					SWEET_ASSERT(str_vector[0] == "#FORMAT");
					SWEET_ASSERT(str_vector[1] == "ASCII");
				}
				else if (i == 2)
				{
					SWEET_ASSERT(str_vector.size() == 2);
					SWEET_ASSERT(str_vector[0] == "#PRIMITIVE");
					SWEET_ASSERT(str_vector[1] == "SCALAR");
				}
				else if (i == 3)
				{
					SWEET_ASSERT(str_vector.size() == 1);
					tmp = stod(str_vector[0]);
				}
			}

			pint_data_ref->dataArrays_2_GenericData_Scalar(tmp);
		}

#elif SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
		sweet::Data::Cart2D::DataSpectral ref_data[] = { sweet::Data::Cart2D::DataSpectral(cart2DDataConfig[0]),
									 sweet::Data::Cart2D::DataSpectral(cart2DDataConfig[0]),
									 sweet::Data::Cart2D::DataSpectral(cart2DDataConfig[0])};

		for (int ivar = 0; ivar < nvar; ivar++)
		{
			std::string i_name;
			if (ivar == 0)
	#if SWEET_PARAREAL_CART2D_SWE || SWEET_XBRAID_CART2D_SWE
				i_name = "prog_h_pert";
	#elif SWEET_PARAREAL_CART2D_BURGERS || SWEET_XBRAID_CART2D_BURGERS
				i_name = "prog_u";
	#endif

			else if (ivar == 1)
	#if SWEET_PARAREAL_CART2D_SWE || SWEET_XBRAID_CART2D_SWE
				i_name = "prog_u";
	#elif SWEET_PARAREAL_CART2D_BURGERS || SWEET_XBRAID_CART2D_BURGERS
				i_name = "prog_v";
	#endif

			else if (ivar == 2)
				i_name = "prog_v";

			if (iteration_id == 0)
			{
				// load ref file
				char buffer[1024];
				const char* filename_template = shackIOData->outputFileName.c_str();
				sprintf(buffer, filename_template, i_name.c_str(), t);
				std::string buffer2 = path_ref + "/" + std::string(buffer);
				sweet::Data::Cart2D::DataGrid tmp(cart2DDataConfig[0]);
				tmp.file_grid_loadRefData_Parareal(buffer2.c_str());
				ref_data[ivar].loadCart2DDataGrid(tmp);

				// If necessary, interpolate to coarsest spatial grid
				if (	cart2DDataConfig[0]->grid_res[0] != ref_data[ivar].cart2DDataConfig->grid_res[0] ||
					cart2DDataConfig[0]->grid_res[1] != ref_data[ivar].cart2DDataConfig->grid_res[1]
				)
				{
					SWEETErrorFatal("TODO");
					//TODO
	///				/*
	///				 * setup sampler
	///				 */
	///				Cart2DDataSampler sampler2D;
	///				sampler2D.setup(simVars.sim.cart2d_domain_size, cart2DDataConfig);
	///		
	///		
	///					/*
	///					 * sample with BiLinear interpolation
	///					 */
	///					Cart2DData prog_h3_bilinear(cart2DDataConfig3);
	///		
	///					sampler2D.bilinear_scalar(
	///							prog_h_pert,	//!< input scalar field
	///							Convert_Cart2DData_2_ScalarDataArray::grid_convert(px),
	///							Convert_Cart2DData_2_ScalarDataArray::grid_convert(py),
	///							prog_h3_bilinear
	///					);
				}
				pint_data_ref->dataArrays_2_GenericData_Cart2DData_Spectral(
											ref_data[0],
											ref_data[1]
	#if SWEET_PARAREAL_CART2D_SWE || SWEET_XBRAID_CART2D_SWE
											, ref_data[2]
	#endif
										);
				}
			}

#elif SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D
		sweet::Data::Sphere2D::DataSpectral ref_data[] = { sweet::Data::Sphere2D::DataSpectral(sphere2DDataConfig[0]),
									  sweet::Data::Sphere2D::DataSpectral(sphere2DDataConfig[0]),
									  sweet::Data::Sphere2D::DataSpectral(sphere2DDataConfig[0])};

		for (int ivar = 0; ivar < nvar; ivar++)
		{
			std::string i_name;
			if (ivar == 0)
				i_name = "prog_phi_pert";
			else if (ivar == 1)
				i_name = "prog_vrt";
			else if (ivar == 2)
				i_name = "prog_div";

			if (shackIOData->outputFileMode == "csv")
			{
				if (iteration_id == 0)
				{
					// load ref file
					char buffer[1024];
					const char* filename_template = shackIOData->outputFileName.c_str();
					sprintf(buffer, filename_template, i_name.c_str(), t * shackIOData->outputFormatTimeScale);
					std::string buffer2 = path_ref + "/" + std::string(buffer);
					sweet::Data::Sphere2D::DataGrid tmp(sphere2DDataConfig[0]);
					tmp.file_grid_loadRefData_Parareal(buffer2.c_str());
					ref_data[ivar].loadSphere2DDataGrid(tmp);

					// If necessary, interpolate to coarsest spatial grid
					if (	sphere2DDataConfig[0]->grid_num_lat != ref_data[ivar].sphere2DDataConfig->grid_num_lat ||
						sphere2DDataConfig[0]->grid_num_lon != ref_data[ivar].sphere2DDataConfig->grid_num_lon
					)
					{
						SWEETErrorFatal("TODO");
						//TODO
	///					/*
	///					 * setup sampler
	///					 */
	///					Cart2DDataSampler sampler2D;
	///					sampler2D.setup(simVars.sim.cart2d_domain_size, cart2DDataConfig);
	///	
	///	
	///						/*
	///						 * sample with BiLinear interpolation
	///						 */
	///						Cart2DData prog_h3_bilinear(cart2DDataConfig3);
	///	
	///						sampler2D.bilinear_scalar(
	///								prog_h_pert,	//!< input scalar field
	///								Convert_Cart2DData_2_ScalarDataArray::grid_convert(px),
	///								Convert_Cart2DData_2_ScalarDataArray::grid_convert(py),
	///								prog_h3_bilinear
	///						);
					}
					pint_data_ref->dataArrays_2_GenericData_Sphere2DData_Spectral(ref_data[0], ref_data[1], ref_data[2]);
				}
			}
			else if (shackIOData->outputFileMode == "bin")
			{

				if (iteration_id == 0)
				{
					// load ref file
					char buffer[1024];
					const char* filename_template = shackIOData->outputFileName.c_str();
					sprintf(buffer, filename_template, i_name.c_str(), t * shackIOData->outputFormatTimeScale);
					std::string buffer2 = path_ref + "/" + std::string(buffer);
					ref_data[ivar].file_read_binary_spectral(buffer2);

					pint_data_ref->dataArrays_2_GenericData_Sphere2DData_Spectral(ref_data[0], ref_data[1], ref_data[2]);
				}

			}
			else
				SWEETErrorFatal("Invalid input data format.");
		}
#endif

		// COMPUTE AND STORE ERRORS
		for (int ivar = 0; ivar < nvar; ivar++)
		{
			int resx_data;
			int resy_data;

			double err_L1; // physical space
			double err_L2; // physical space
			double err_Linf; // physical space

			std::map<std::size_t, double> err_Linf_spectral;
			std::vector<std::size_t> rnorms;
			for (int ip = 0; ip <= 5; ip++)
			{
#if SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D
				int rnorm = cart2DDataConfig[0]->spectral_data_size[0] / std::pow(2, ip);
				if (rnorm >= 1)
					rnorms.push_back(rnorm);
#elif SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D
				int rnorm = sphere2DDataConfig[0]->spectral_modes_m_max / std::pow(2, ip);
				if (rnorm >= 8)
					rnorms.push_back(rnorm);
#endif
			}

			std::string i_name;

#if SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR
			i_name = "prog_u";
			double u_ref;
			pint_data_ref->GenericData_Scalar_2_dataArrays(u_ref);
			double err = std::abs(	i_data->get_pointer_to_data_Scalar()->simfields[ivar] -
						pint_data_ref->get_pointer_to_data_Scalar()->simfields[ivar]);
			err_L1 = err;
			err_L2 = err;
			err_Linf = err;


#elif SWEET_PARAREAL_CART2D || SWEET_XBRAID_CART2D

			if (ivar == 0)
	#if SWEET_PARAREAL_CART2D_SWE || SWEET_XBRAID_CART2D_SWE
				i_name = "prog_h_pert";
	#elif SWEET_PARAREAL_CART2D_BURGERS || SWEET_XBRAID_CART2D_BURGERS
				i_name = "prog_u";
	#endif

			else if (ivar == 1)
	#if SWEET_PARAREAL_CART2D_SWE || SWEET_XBRAID_CART2D_SWE
				i_name = "prog_u";
	#elif SWEET_PARAREAL_CART2D_BURGERS || SWEET_XBRAID_CART2D_BURGERS
				i_name = "prog_v";
	#endif

			else if (ivar == 2)
				i_name = "prog_v";

			resx_data = cart2DDataConfig[0]->grid_res[0];
			resy_data = cart2DDataConfig[0]->grid_res[1];

			sweet::Data::Cart2D::DataSpectral diff_spectral = *i_data->get_pointer_to_data_Cart2DData_Spectral()->simfields[ivar]
									- *pint_data_ref->get_pointer_to_data_Cart2DData_Spectral()->simfields[ivar];
			sweet::Data::Cart2D::DataGrid diff = i_data->get_pointer_to_data_Cart2DData_Spectral()->simfields[ivar]->toGrid() -
									pint_data_ref->get_pointer_to_data_Cart2DData_Spectral()->simfields[ivar]->toGrid();
			err_L1 = diff.grid_reduce_norm1() / (resx_data * resy_data);
			err_L2 = diff.grid_reduce_norm2() / std::sqrt(resx_data * resy_data);
			err_Linf = diff.grid_reduce_max_abs();


			// Spectral space
			double small = 1e-20;
			for (std::vector<std::size_t>::iterator it = rnorms.begin(); it != rnorms.end(); it++)
			{
				double norm_diff = std::sqrt(diff_spectral.spectral_reduce_max_abs(*it) );
				double norm_ref = std::sqrt(pint_data_ref->get_pointer_to_data_Cart2DData_Spectral()->simfields[ivar]->spectral_reduce_max_abs(*it) );
				if (norm_diff < small and norm_ref < small)
					err_Linf_spectral.emplace(std::make_pair(*it, 0.));
				else
					err_Linf_spectral.emplace(std::make_pair(*it, norm_diff / norm_ref ));
			}

#elif SWEET_PARAREAL_SPHERE2D || SWEET_XBRAID_SPHERE2D
			if (ivar == 0)
				i_name = "prog_phi_pert";
			else if (ivar == 1)
				i_name = "prog_vrt";
			else if (ivar == 2)
				i_name = "prog_div";

			resx_data = sphere2DDataConfig[0]->grid_num_lon;
			resy_data = sphere2DDataConfig[0]->grid_num_lat;

			sweet::Data::Sphere2D::DataSpectral diff_spectral = *i_data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[ivar]
									- *pint_data_ref->get_pointer_to_data_Sphere2DData_Spectral()->simfields[ivar];
			sweet::Data::Sphere2D::DataGrid diff = i_data->get_pointer_to_data_Sphere2DData_Spectral()->simfields[ivar]->toGrid() -
									pint_data_ref->get_pointer_to_data_Sphere2DData_Spectral()->simfields[ivar]->toGrid();
			err_L1 = diff.grid_reduce_norm1() / (resx_data * resy_data);
			err_L2 = diff.grid_reduce_norm2() / std::sqrt(resx_data * resy_data);
			err_Linf = diff.grid_reduce_max_abs();

			// Spectral space
			///double small = 1e-20;
			double small = 1e-16;
			for (std::vector<std::size_t>::iterator it = rnorms.begin(); it != rnorms.end(); it++)
			{
				double norm_diff = std::sqrt(diff_spectral.spectral_reduce_max_abs(*it));
				double norm_ref = std::sqrt(pint_data_ref->get_pointer_to_data_Sphere2DData_Spectral()->simfields[ivar]->spectral_reduce_max_abs(*it));
				if ( norm_diff < small && norm_ref < small )
					err_Linf_spectral.emplace(std::make_pair(*it, 0.));
				else
					err_Linf_spectral.emplace(std::make_pair(*it, norm_diff / norm_ref ));
			}

#endif

			// save physical errors in file
			char buffer_out[1024];

			///const char* filename_template_out = "parareal_error_%s_%s_t%020.8f_iter%03d.csv";
			std::string str = pint_type + "_error_%s_%s_t%020.8f_iter%03d.csv";
			const char* filename_template_out = str.c_str();
			sprintf(buffer_out, filename_template_out, base_solution.c_str(), i_name.c_str(), t * shackIOData->outputFormatTimeScale, iteration_id);

			std::ofstream file(buffer_out, std::ios_base::trunc);
			file << std::setprecision(i_precision);

			file << "#BASESOLUTION " << base_solution << " " << path_ref << std::endl;
			file << "#VAR " << i_name << std::endl;
			file << "#ITERATION " << iteration_id << std::endl;
			file << "#TIMESLICE " << time_slice_id << std::endl;
			file << "#TIMEFRAMEEND " << t  * shackIOData->outputFormatTimeScale << std::endl;
			file << "errL1 " << err_L1 << std::endl;
			file << "errL2 " << err_L2 << std::endl;
			file << "errLinf " << err_Linf << std::endl;

			file.close();


#if !(SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR)
			// save spectral errors in file
			char buffer_out_spec[1024];

			std::string str2 = pint_type + "_error_spec_%s_%s_t%020.8f_iter%03d.csv";
			const char* filename_template_out_spec = str2.c_str();
			sprintf(buffer_out_spec, filename_template_out_spec, base_solution.c_str(), i_name.c_str(), t * shackIOData->outputFormatTimeScale, iteration_id);

			std::ofstream file_spec(buffer_out_spec, std::ios_base::trunc);
			file_spec << std::setprecision(i_precision);

			file_spec << "#BASESOLUTION " << base_solution << " " << path_ref << std::endl;
			file_spec << "#VAR " << i_name << std::endl;
			file_spec << "#ITERATION " << iteration_id << std::endl;
			file_spec << "#TIMESLICE " << time_slice_id << std::endl;
			file_spec << "#TIMEFRAMEEND " << t  * shackIOData->outputFormatTimeScale << std::endl;
			for (std::vector<std::size_t>::iterator it = rnorms.begin(); it != rnorms.end(); it++)
				file_spec << "errLinf " << *it + 1 << " " << err_Linf_spectral.at(*it) << std::endl;

			file_spec.close();
#endif

		}
		pint_data_ref = nullptr;
		pint_data_ref = nullptr;
	}
};

}}

#endif
