/*
 * PInT_Common.hpp
 *
 *  Created on: 10 Jun 2022
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 */

#include <sweet/IO/ShackIOData.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Tools/DefaultPrecompilerValues.hpp>

/*
 * Contains functions used by both Parareal and Xbraid (and possibly PFASST?) implementations in SWEET,
 * mostly error computation and file output functions.
 */

#ifndef INCLUDE_SWEET_DEPRECATED_PINT_PINT_COMMON_HPP
#define INCLUDE_SWEET_DEPRECATED_PINT_PINT_COMMON_HPP

/////#include <sweet/core/SimulationVariables.hpp>
#include <sweet/DEPRECATED_pint/Parareal_GenericData.hpp>
#include <sweet/Data/GenericContainer/Base.hpp>
#if SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR
	#include <sweet/DEPRECATED_pint/Parareal_GenericData_Scalar.hpp>
#elif SWEET_PARAREAL_PLANE || SWEET_XBRAID_PLANE
	#include <sweet/Data/Plane2D/PlaneOperators.hpp>
	#include <sweet/DEPRECATED_pint/Parareal_GenericData_PlaneData_Spectral.hpp>
	#include <sweet/Data/Plane2D/ShackPlaneDataOps.hpp>
	#include <sweet/Data/Plane2D/PlaneDataGridMapping.hpp>
#elif SWEET_PARAREAL_SPHERE || SWEET_XBRAID_SPHERE
	#include <sweet/DEPRECATED_pint/Parareal_GenericData_SphereData_Spectral.hpp>
	#include <sweet/Data/Sphere2D/SphereOperators.hpp>
	#include <sweet/Data/Sphere2D/SphereOperatorsComplex.hpp>
	#include <sweet/Data/Sphere2D/ShackSphereDataOps.hpp>
	#include "src/programs/PDE_SWESphere/Shack.hpp"
	#include "src/programs/PDE_SWESphere/DataContainer.hpp"
#endif

#include <map>

namespace sweet {
namespace DEPRECATED_pint {

class PInT_Common
{

protected:
	///SimulationVariables* simVars = nullptr;

	sweet::Error::Base error;

	sweet::IO::ShackIOData* shackIOData;

#if SWEET_PARAREAL_PLANE || SWEET_XBRAID_PLANE
	// Grid Mapping (staggered grid)
	sweet::Data::Plane2D::PlaneDataGridMapping gridMapping;
#endif

	// Operators and DataConfig
#if SWEET_PARAREAL_PLANE || SWEET_XBRAID_PLANE
	sweet::Data::Plane2D::PlaneOperators* base_op_plane;
	sweet::Data::Plane2D::PlaneData_Config* base_planeDataConfig;
	std::vector<sweet::Data::Plane2D::PlaneOperators*> op_plane;
	std::vector<sweet::Data::Plane2D::PlaneData_Config*> planeDataConfig;
	sweet::Data::Plane2D::ShackPlaneDataOps* shackPlaneDataOps;
#elif SWEET_PARAREAL_SPHERE || SWEET_XBRAID_SPHERE
	sweet::Data::Sphere2D::SphereOperators* base_op_sphere;
	sweet::Data::Sphere2D::SphereOperatorsComplex* base_op_sphere_complex;
	sweet::Data::Sphere2D::SphereOperators* base_op_sphere_nodealiasing;
	sweet::Data::Sphere2D::SphereData_Config* base_sphereDataConfig;
	std::vector<sweet::Data::Sphere2D::SphereOperators*> op_sphere;
	std::vector<sweet::Data::Sphere2D::SphereOperatorsComplex*> op_sphere_complex;
	std::vector<sweet::Data::Sphere2D::SphereOperators*> op_sphere_nodealiasing;
	std::vector<sweet::Data::Sphere2D::SphereData_Config*> sphereDataConfig;
	sweet::Data::Sphere2D::ShackSphereDataOps* shackSphereDataOps;
	PDE_SWESphere::Shack* shackPDESWESphere;
#endif

	// list of SL schemes
	std::vector<std::string> SL_tsm = {};


#if SWEET_PARAREAL_PLANE_BURGERS
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
	Burgers_Plane_TimeSteppers* timeSteppersFineBurgers = nullptr;
	BenchmarkErrors benchmark;

	void set_tsm_burgers(
				Burgers_Plane_TimeSteppers* i_timeSteppersFineBurgers
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
			sweet::IO::ShackIOData* i_shackIOData
#if SWEET_PARAREAL_PLANE || SWEET_XBRAID_PLANE
			,
			sweet::Data::Plane2D::ShackPlaneDataOps* i_shackPlaneDataOps
#elif SWEET_PARAREAL_SPHERE || SWEET_XBRAID_SPHERE
			,
			sweet::Data::Sphere2D::ShackSphereDataOps* i_shackSphereDataOps
#endif
	)
	{
		////simVars = i_simVars;

		shackIOData = i_shackIOData;
#if SWEET_PARAREAL_PLANE || SWEET_XBRAID_PLANE
		shackPlaneDataOps = i_shackPlaneDataOps;
#elif SWEET_PARAREAL_SPHERE || SWEET_XBRAID_SPHERE
		shackSphereDataOps = i_shackSphereDataOps;
#endif

	#if SWEET_PARAREAL_PLANE_SWE || SWEET_XBRAID_PLANE_SWE
		SL_tsm = { "l_cn_na_sl_nd_settls",
				 "l_rexi_na_sl_nd_etdrk",
				 "l_rexi_na_sl_nd_settls"
				};
	#elif SWEET_PARAREAL_PLANE_BURGERS || SWEET_XBRAID_PLANE_BURGERS
		SL_tsm = { "l_cn_n_sl",
				 "l_irk_n_sl",
				 "l_irk_n_sl_forcing"
				};
	#elif SWEET_PARAREAL_SPHERE || SWEET_XBRAID_SPHERE
		SL_tsm = { "lg_exp_na_sl_lc_nr_etd_uv",
				 "l_irk_na_sl_nr_settls_uv_only",
				 "l_irk_na_sl_nr_settls_vd_only",
				 "l_irk_na_sl_settls_uv_only",
				 "l_irk_na_sl_settls_vd_only",
				 "ln_sl_exp_settls_uv",
				 "ln_sl_exp_settls_vd",
				 "lg_exp_na_sl_lc_nr_etdrk_uv"
				};
	#endif


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
			sweet::Data::GenericContainer::Base* i_data,
			int iteration_id,
			int time_slice_id,
			double t
	)
	{
		// Convert Container to PararealData
		sweet::DEPRECATED_pint::Parareal_GenericData* data;

#if SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR
		ODE_Scalar::DataContainer* data_container = (ODE_Scalar::DataContainer*) i_data;
		data = new sweet::DEPRECATED_pint::Parareal_GenericData_Scalar<1>;
		data->allocate_data();
		data->dataArrays_2_GenericData_Scalar(data_container->data[0]);
#elif SWEET_PARAREAL_PLANE || SWEET_XBRAID_PLANE
		PDE_SWEPlane::DataContainer* data_container = (PDE_SWEPlane::DataContainer*) i_data;
		data = new sweet::DEPRECATED_pint::Parareal_GenericData_PlaneData_Spectral<3>;
		data->setup_data_config(data_container->data[0].planeDataConfig);
		data->allocate_data();
		data->dataArrays_2_GenericData_PlaneData_Spectral(data_container->data[0], data_container->data[1], data_container->data[2]);
#elif SWEET_PARAREAL_SPHERE || SWEET_XBRAID_SPHERE
		PDE_SWESphere::DataContainer* data_container = (PDE_SWESphere::DataContainer*) i_data;
		data = new sweet::DEPRECATED_pint::Parareal_GenericData_SphereData_Spectral<3>;
		data->setup_data_config(data_container->data[0].sphereDataConfig);
		data->allocate_data();
		data->dataArrays_2_GenericData_SphereData_Spectral(data_container->data[0], data_container->data[1], data_container->data[2]);
#endif

		this->output_data_file(
						data,
						iteration_id,
						time_slice_id,
						t
					);
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
		if (shackIOData->output_file_name.size() > 0)
		{
			std::string output_filenames = "";
			output_filenames = write_file_pint_scalar(u_out, "prog_u", iteration_id, t);
		}

#elif SWEET_PARAREAL_PLANE || SWEET_XBRAID_PLANE
	#if SWEET_PARAREAL_PLANE_BURGERS || SWEET_XBRAID_PLANE_BURGERS

		sweet::Data::Plane2D::PlaneData_Spectral dummy(planeDataConfig[0]);
		sweet::Data::Plane2D::PlaneData_Spectral u_out(planeDataConfig[0]);
		sweet::Data::Plane2D::PlaneData_Spectral v_out(planeDataConfig[0]);
		i_data->GenericData_PlaneData_Spectral_2_dataArrays(u_out, v_out);

		sweet::Data::Plane2D::PlaneData_Physical u_out_phys = u_out.toPhys();
		sweet::Data::Plane2D::PlaneData_Physical v_out_phys = v_out.toPhys();

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
		// For output, variables need to be on unstaggered A-grid
		sweet::Data::Plane2D::PlaneData_Physical t_u(planeDataConfig[0]);
		sweet::Data::Plane2D::PlaneData_Physical t_v(planeDataConfig[0]);

		if (shackPlaneDataOps->space_grid_use_c_staggering) // Remap in case of C-grid
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
		if (shackIOData->output_file_name.size() > 0)
		{
			std::string output_filenames = "";

			output_filenames = write_file_pint_plane(t_u, "prog_u", iteration_id, t);
			output_filenames += ";" + write_file_pint_plane(t_v, "prog_v", iteration_id, t);

			output_filenames += ";" + write_file_spec_pint_plane(u_out, "prog_u_spec", iteration_id, t);
			output_filenames += ";" + write_file_spec_pint_plane(v_out, "prog_v_spec", iteration_id, t);

		}

		write_file_spec_amp_phase_pint_plane(u_out, "prog_u", iteration_id, t);

		if (simVars->misc.compute_errors)
		{
			sweet::Data::Plane2D::PlaneData_Spectral ana = compute_errors2(u_out, v_out);

			write_file_pint_plane(ana.toPhys(),"analytical",iteration_id, t);
			write_file_spec_amp_phase_pint_plane(ana.toPhys(), "analytical", iteration_id, t);
		}

	#elif SWEET_PARAREAL_PLANE_SWE || SWEET_XBRAID_PLANE_SWE

		sweet::Data::Plane2D::PlaneData_Spectral h_out(planeDataConfig[0]);
		sweet::Data::Plane2D::PlaneData_Spectral u_out(planeDataConfig[0]);
		sweet::Data::Plane2D::PlaneData_Spectral v_out(planeDataConfig[0]);
		i_data->GenericData_PlaneData_Spectral_2_dataArrays(h_out, u_out, v_out);

		sweet::Data::Plane2D::PlaneData_Physical h_out_phys = h_out.toPhys();
		sweet::Data::Plane2D::PlaneData_Physical u_out_phys = u_out.toPhys();
		sweet::Data::Plane2D::PlaneData_Physical v_out_phys = v_out.toPhys();

		// Save .vtk files for visualizing in paraview
		std::ostringstream ss2;
		ss2 << "output_slice" << time_slice_id << "_iter" << iteration_id << ".vtk";
		std::string filename2 = ss2.str();
		h_out_phys.file_physical_saveData_vtk(filename2.c_str(), filename2.c_str());

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
		// For output, variables need to be on unstaggered A-grid
		sweet::Data::Plane2D::PlaneData_Physical t_h(planeDataConfig[0]);
		sweet::Data::Plane2D::PlaneData_Physical t_u(planeDataConfig[0]);
		sweet::Data::Plane2D::PlaneData_Physical t_v(planeDataConfig[0]);

		if (shackPlaneDataOps->space_grid_use_c_staggering) // Remap in case of C-grid
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
		if (shackIOData->output_file_name.size() > 0)
		{
			std::string output_filenames = "";

			output_filenames = write_file_pint_plane(t_h, "prog_h_pert", iteration_id, t);
			output_filenames += ";" + write_file_pint_plane(t_u, "prog_u", iteration_id, t);
			output_filenames += ";" + write_file_pint_plane(t_v, "prog_v", iteration_id, t);

			output_filenames += ";" + write_file_pint_plane(op_plane[0]->ke(t_u,t_v).toPhys(),"diag_ke", iteration_id, t);

			output_filenames += ";" + write_file_spec_pint_plane(h_out, "prog_h_pert_spec", iteration_id, t);
			output_filenames += ";" + write_file_spec_pint_plane(u_out, "prog_u_spec", iteration_id, t);
			output_filenames += ";" + write_file_spec_pint_plane(v_out, "prog_v_spec", iteration_id, t);

			output_filenames += ";" + write_file_spec_pint_plane(op_plane[0]->ke(t_u,t_v).toPhys(),"diag_ke_spec", iteration_id, t);

			output_filenames += ";" + write_file_pint_plane(op_plane[0]->vort(t_u, t_v).toPhys(), "diag_vort", iteration_id, t);
			output_filenames += ";" + write_file_pint_plane(op_plane[0]->div(t_u, t_v).toPhys(), "diag_div", iteration_id, t);

			/////////if (compute_normal_modes){
			/////////	SWEETErrorFatal("TODO");
			/////////	///output_filenames += ";" + write_file_spec_parareal(normalmodes.geo, "nm_geo", iteration_id, output_initial_data);
			/////////	///output_filenames += ";" + write_file_spec_parareal(normalmodes.igwest, "nm_igwest", iteration_id, output_initial_data);
			/////////	///output_filenames += ";" + write_file_spec_parareal(normalmodes.igeast, "nm_igeast", iteration_id, output_initial_data);
			/////////}
			
		}
	#endif


#elif SWEET_PARAREAL_SPHERE || SWEET_XBRAID_SPHERE

		sweet::Data::Sphere2D::SphereData_Spectral phi_out(sphereDataConfig[0]);
		sweet::Data::Sphere2D::SphereData_Spectral vrt_out(sphereDataConfig[0]);
		sweet::Data::Sphere2D::SphereData_Spectral div_out(sphereDataConfig[0]);
		i_data->GenericData_SphereData_Spectral_to_dataArrays(phi_out, vrt_out, div_out);

		sweet::Data::Sphere2D::SphereData_Physical phi_out_phys = phi_out.toPhys();
		sweet::Data::Sphere2D::SphereData_Physical vrt_out_phys = vrt_out.toPhys();
		sweet::Data::Sphere2D::SphereData_Physical div_out_phys = div_out.toPhys();

		///////////////// Save .vtk files for visualizing in paraview
		///////////////std::ostringstream ss2;
		///////////////if (output_initial_data)
		///////////////	ss2 << "output_slice" << time_slice_id - 1 << "_iter" << iteration_id << ".vtk";
		///////////////else
		///////////////	ss2 << "output_slice" << time_slice_id << "_iter" << iteration_id << ".vtk";
		///////////////std::string filename2 = ss2.str();
		///////////////phi_out_phys.file_physical_saveData_vtk(filename2.c_str(), filename2.c_str());

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
		// Dump  data in csv, if output filename is not empty
		if (shackIOData->output_file_name.size() > 0)
		{
			if (shackIOData->output_file_mode == "csv")
			{
				std::string output_filename;

				sweet::Data::Sphere2D::SphereData_Spectral h = phi_out_phys*(1.0/shackPDESWESphere->gravitation);
				h += shackPDESWESphere->h0;

				output_filename = write_file_csv_pint_sphere(h, t, "prog_h", iteration_id);
				//std::cout << " + " << output_filename << " (min: " << h.toPhys().physical_reduce_min() << ", max: " << h.toPhys().physical_reduce_max() << ")" << std::endl;

				sweet::Data::Sphere2D::SphereData_Physical phi_phys = h.toPhys() * shackPDESWESphere->gravitation;
				sweet::Data::Sphere2D::SphereData_Spectral phi(sphereDataConfig[0]);
				phi.loadSphereDataPhysical(phi_phys);
				output_filename = write_file_csv_pint_sphere(phi, t, "prog_phi", iteration_id);

				output_filename = write_file_csv_pint_sphere(phi_out, t, "prog_phi_pert", iteration_id);
				//std::cout << " + " << output_filename << " (min: " << phi_out_phys.physical_reduce_min() << ", max: " << phi_out_phys.physical_reduce_max() << ")" << std::endl;
	
				sweet::Data::Sphere2D::SphereData_Physical u(sphereDataConfig[0]);
				sweet::Data::Sphere2D::SphereData_Physical v(sphereDataConfig[0]);
	
				op_sphere[0]->vrtdiv_2_uv(vrt_out_phys, div_out_phys, u, v);
	
				output_filename = write_file_csv_pint_sphere(u, t, "prog_u", iteration_id);
				//std::cout << " + " << output_filename << std::endl;
	
				output_filename = write_file_csv_pint_sphere(v, t, "prog_v", iteration_id);
				//std::cout << " + " << output_filename << std::endl;
	
				output_filename = write_file_csv_pint_sphere(vrt_out, t, "prog_vrt", iteration_id);
				//std::cout << " + " << output_filename << std::endl;
	
				output_filename = write_file_csv_pint_sphere(div_out, t, "prog_div", iteration_id);
				//std::cout << " + " << output_filename << std::endl;
	
				sweet::Data::Sphere2D::SphereData_Spectral potvrt = (phi_out / shackPDESWESphere->gravitation)*vrt_out;
	
				output_filename = write_file_csv_pint_sphere(potvrt, t, "prog_potvrt", iteration_id);
				//std::cout << " + " << output_filename << std::endl;


				////output_filename = write_file_csv_pint_sphere_spec(phi_out, t, "prog_phi_pert", iteration_id);

			}
			else if (shackIOData->output_file_mode == "bin")
			{
				std::string output_filename;
	
				{
					output_filename = write_file_bin_pint_sphere(phi_out, t, "prog_phi_pert", iteration_id);
					sweet::Data::Sphere2D::SphereData_Physical prog_phys = phi_out.toPhys();
	
					//std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
				}
	
				{
					output_filename = write_file_bin_pint_sphere(vrt_out, t, "prog_vrt", iteration_id);
					sweet::Data::Sphere2D::SphereData_Physical prog_phys = vrt_out.toPhys();
	
					//std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
				}
	
				{
					output_filename = write_file_bin_pint_sphere(div_out, t, "prog_div", iteration_id);
					sweet::Data::Sphere2D::SphereData_Physical prog_phys = div_out.toPhys();
	
					//std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
				}
			}
			
		}


#endif
	};


#if SWEET_PARAREAL_PLANE_BURGERS || SWEET_XBRAID_PLANE_BURGERS
	// For Burgers
	sweet:PlaneData_Spectral compute_errors2(
						const sweet::Data::Plane2D::PlaneData_Spectral &i_planeData_u,
						const sweet::Data::Plane2D::PlaneData_Spectral &i_planeData_v
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



		// Necessary to circumvent FFTW transformations on i_planeData_u and i_planeData_v, which would lead to errors
		sweet::Data::Plane2D::PlaneData_Physical u = i_planeData_u.toPhys();
		sweet::Data::Plane2D::PlaneData_Physical v = i_planeData_v.toPhys();

		///// Analytical solution at current time on original grid
		///PlaneData_Spectral ts_u = t0_prog_u;
		///PlaneData_Spectral ts_v = t0_prog_v;

		sweet::Data::Plane2D::PlaneData_Spectral ts_u(planeDataConfig[0]);
		sweet::Data::Plane2D::PlaneData_Spectral ts_v(planeDataConfig[0]);
		sweet::Data::Plane2D::PlaneData_Physical ts_u_phys(planeDataConfig[0]);
		sweet::Data::Plane2D::PlaneData_Physical ts_v_phys(planeDataConfig[0]);

		if (simVars->misc.compute_errors)
		{
			//if (simVars.setup.benchmark_id > 51 && simVars.setup.benchmark_id < 65)
			if (simVars->disc.timestepping_method.find("forcing")!=std::string::npos)
			{
				if (simVars->disc.space_grid_use_c_staggering)
				{
					ts_u_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							double x = (((double)i)/(double)simVars->disc.space_res_physical[0])*simVars->sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars->disc.space_res_physical[1])*simVars->sim.plane_domain_size[1];
							io_data = BurgersValidationBenchmarks::return_u(*simVars, x, y);
						}
					);

					ts_v_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							io_data = 0.0;
#if 0
							double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
							double y = (((double)j)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];
							io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
#endif
						}
					);
				}
				else
				{
					ts_u_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							double x = (((double)i+0.5)/(double)simVars->disc.space_res_physical[0])*simVars->sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars->disc.space_res_physical[1])*simVars->sim.plane_domain_size[1];

							io_data = BurgersValidationBenchmarks::return_u(*simVars, x, y);
						}
					);

					ts_v_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							io_data = 0.0;
#if 0
							double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];

							io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
#endif
						}
					);
				}
				ts_u.loadPlaneDataPhysical(ts_u_phys);
				ts_v.loadPlaneDataPhysical(ts_v_phys);
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
			benchmark.benchmark_analytical_error_rms_u = (ts_u-u).toPhys().physical_reduce_rms();
			benchmark.benchmark_analytical_error_rms_v = (ts_v-v).toPhys().physical_reduce_rms();

			benchmark.benchmark_analytical_error_maxabs_u = (ts_u-u).toPhys().physical_reduce_max_abs();
			benchmark.benchmark_analytical_error_maxabs_v = (ts_v-v).toPhys().physical_reduce_max_abs();

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
			const char* i_name,	///< name of output variable
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

#if SWEET_PARAREAL_SPHERE || SWEET_XBRAID_SPHERE
	/**
	 * Write physical data to file and return string of file name (parareal)
	 */
	std::string write_file_csv_pint_sphere(
			const sweet::Data::Sphere2D::SphereData_Spectral &i_sphereData,
			double t,
			const char* i_name,	///< name of output variable
			int iteration_id,
			bool i_phi_shifted = false
		)
	{
		char buffer[1024];

		// create copy
		sweet::Data::Sphere2D::SphereData_Physical sphereData = i_sphereData.toPhys();

		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		sprintf(buffer, filename_template, i_name, t * shackIOData->output_time_scale, iteration_id);

		if (i_phi_shifted)
			sphereData.physical_file_write_lon_pi_shifted(buffer, "vorticity, lon pi shifted");
		else
			sphereData.physical_file_write(buffer);

		return buffer;

	}

/////	/**
/////	 * Write spectral data to file and return string of file name (parareal)
/////	 */
/////	std::string write_file_csv_pint_sphere_spec(
/////			const sweet::Data::Sphere2D::SphereData_Spectral &i_sphereData,
/////			double t,
/////			const char* i_name,	///< name of output variable
/////			int iteration_id
/////		)
/////	{
/////		char buffer[1024];
/////
/////		const char* filename_template = "output_spec_%s_t%020.8f_iter%03d.csv";
/////		sprintf(buffer, filename_template, i_name, t * simVars->iodata.output_time_scale, iteration_id);
/////		i_sphereData.file_spectral_saveData_ascii(buffer);
/////		return buffer;
/////
/////	}


	/**
	 * Write spectral data to file and return string of file name
	 */
	std::string write_file_bin_pint_sphere(
			const sweet::Data::Sphere2D::SphereData_Spectral &i_sphereData,
			double t,
			const char* i_name,
			int iteration_id
	)
	{
		char buffer[1024];

		sweet::Data::Sphere2D::SphereData_Spectral sphereData(i_sphereData);
		//const char* filename_template = simVars.iodata.output_file_name.c_str();
		const char* filename_template = "output_%s_t%020.8f_iter%03d.sweet";
		sprintf(buffer, filename_template, i_name, t * shackIOData->output_time_scale, iteration_id);
		sphereData.file_write_binary_spectral(buffer);

		return buffer;
	}
#endif


#if SWEET_PARAREAL_PLANE || SWEET_XBRAID_PLANE
	/**
	 * Write file to data and return string of file name (parareal)
	 */
	std::string write_file_pint_plane(
			const sweet::Data::Plane2D::PlaneData_Physical &i_planeData,
			const char* i_name,	///< name of output variable
			int iteration_id,
			double t
		)
	{
		char buffer[1024];

		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		sprintf(buffer, filename_template, i_name, t, iteration_id);
		i_planeData.file_physical_saveData_ascii(buffer);
		return buffer;
	}

	/**
	 * Write spectrum info to data and return string of file name (parareal)
	 */
	std::string write_file_spec_pint_plane(
			const sweet::Data::Plane2D::PlaneData_Spectral &i_planeData,
			const char* i_name,	///< name of output variable
			int iteration_id,
			double t
		)
	{
		char buffer[1024];

		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		sprintf(buffer, filename_template, i_name, t, iteration_id);
		i_planeData.file_spectral_saveData_ascii(buffer);
		///i_planeData.file_spectral_abs_saveData_ascii(buffer);
		return buffer;
	}

	/**
	 * Write spectrum info to data and return string of file name (parareal)
	 */
	std::string write_file_spec_amp_phase_pint_plane(
			const sweet::Data::Plane2D::PlaneData_Spectral &i_planeData,
			const char* i_name,	///< name of output variable
			int iteration_id,
			double t
		)
	{

		char buffer[1024];

		const char* filename_template = "output_%s_amp_phase_t%020.8f_iter%03d.csv";
		sprintf(buffer, filename_template, i_name, t, iteration_id);

		std::ofstream file(buffer, std::ios_base::trunc);
		file << std::setprecision(12);

		for (std::size_t x = 0; x < planeDataConfig[0]->spectral_data_size[0]; x++)
		{
			file << x << ", " << i_planeData.spectral_return_amplitude(0,x) << ", " << i_planeData.spectral_return_phase(0,x) << std::endl;
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
			sweet::Data::GenericContainer::Base* i_data,
			sweet::Data::GenericContainer::Base* pint_data_ref,
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

		// Convert Container to PararealData
		sweet::DEPRECATED_pint::Parareal_GenericData* data;
		sweet::DEPRECATED_pint::Parareal_GenericData* data_ref;

#if SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR
		ODE_Scalar::DataContainer* data_container = (ODE_Scalar::DataContainer*) i_data;
		ODE_Scalar::DataContainer* data_container_ref = (ODE_Scalar::DataContainer*) pint_data_ref;

		data = new sweet::DEPRECATED_pint::Parareal_GenericData_ScalarData_Spectral<1>;
		data->allocate_data();

		data_ref = new sweet::DEPRECATED_pint::Parareal_GenericData_ScalarData_Spectral<1>;
		data_ref->allocate_data();

		data->dataArrays_2_GenericData_Scalar(data_container->data[0]);
		data_ref->dataArrays_2_GenericData_Scalar(data_container_ref->data[0]);
#elif SWEET_PARAREAL_PLANE || SWEET_XBRAID_PLANE
		PDE_SWEPlane::DataContainer* data_container = (PDE_SWEPlane::DataContainer*) i_data;
		PDE_SWEPlane::DataContainer* data_container_ref = (PDE_SWEPlane::DataContainer*) pint_data_ref;

		data = new sweet::DEPRECATED_pint::Parareal_GenericData_PlaneData_Spectral<3>;
		data->setup_data_config(data_container->data[0].planeDataConfig);
		data->allocate_data();

		data_ref = new sweet::DEPRECATED_pint::Parareal_GenericData_PlaneData_Spectral<3>;
		data_ref->setup_data_config(data_container_ref->data[0].planeDataConfig);
		data_ref->allocate_data();

		data->dataArrays_2_GenericData_PlaneData_Spectral(data_container->data[0], data_container->data[1], data_container->data[2]);
		data_ref->dataArrays_2_GenericData_PlaneData_Spectral(data_container_ref->data[0], data_container_ref->data[1], data_container_ref->data[2]);
#elif SWEET_PARAREAL_SPHERE || SWEET_XBRAID_SPHERE
		PDE_SWESphere::DataContainer* data_container = (PDE_SWESphere::DataContainer*) i_data;
		PDE_SWESphere::DataContainer* data_container_ref = (PDE_SWESphere::DataContainer*) pint_data_ref;

		data = new sweet::DEPRECATED_pint::Parareal_GenericData_SphereData_Spectral<3>;
		data->setup_data_config(data_container->data[0].sphereDataConfig);
		data->allocate_data();

		data_ref = new sweet::DEPRECATED_pint::Parareal_GenericData_SphereData_Spectral<3>;
		data_ref->setup_data_config(data_container_ref->data[0].sphereDataConfig);
		data_ref->allocate_data();

		data->dataArrays_2_GenericData_SphereData_Spectral(data_container->data[0], data_container->data[1], data_container->data[2]);
		data_ref->dataArrays_2_GenericData_SphereData_Spectral(data_container_ref->data[0], data_container_ref->data[1], data_container_ref->data[2]);
#endif

		store_pint_error(
				data,
				data_ref,
				nvar,
				iteration_id,
				time_slice_id,
				t,
				path_ref,
				base_solution,	// "ref" or "fine"
				pint_type,
				i_precision
		);

	}


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
			const char* filename_template = shackIOData->output_file_name.c_str();
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
					assert(str_vector.size() == 1);
					assert(str_vector[0] == "#SWEET");
				}
				else if (i == 1)
				{
					assert(str_vector.size() == 2);
					assert(str_vector[0] == "#FORMAT");
					assert(str_vector[1] == "ASCII");
				}
				else if (i == 2)
				{
					assert(str_vector.size() == 2);
					assert(str_vector[0] == "#PRIMITIVE");
					assert(str_vector[1] == "SCALAR");
				}
				else if (i == 3)
				{
					assert(str_vector.size() == 1);
					tmp = stod(str_vector[0]);
				}
			}

			pint_data_ref->dataArrays_2_GenericData_Scalar(tmp);
		}

#elif SWEET_PARAREAL_PLANE || SWEET_XBRAID_PLANE
		sweet::Data::Plane2D::PlaneData_Spectral ref_data[] = { sweet::Data::Plane2D::PlaneData_Spectral(planeDataConfig[0]),
									 sweet::Data::Plane2D::PlaneData_Spectral(planeDataConfig[0]),
									 sweet::Data::Plane2D::PlaneData_Spectral(planeDataConfig[0])};

		for (int ivar = 0; ivar < nvar; ivar++)
		{
			std::string i_name;
			if (ivar == 0)
	#if SWEET_PARAREAL_PLANE_SWE || SWEET_XBRAID_PLANE_SWE
				i_name = "prog_h_pert";
	#elif SWEET_PARAREAL_PLANE_BURGERS || SWEET_XBRAID_PLANE_BURGERS
				i_name = "prog_u";
	#endif

			else if (ivar == 1)
	#if SWEET_PARAREAL_PLANE_SWE || SWEET_XBRAID_PLANE_SWE
				i_name = "prog_u";
	#elif SWEET_PARAREAL_PLANE_BURGERS || SWEET_XBRAID_PLANE_BURGERS
				i_name = "prog_v";
	#endif

			else if (ivar == 2)
				i_name = "prog_v";

			if (iteration_id == 0)
			{
				// load ref file
				char buffer[1024];
				const char* filename_template = shackIOData->output_file_name.c_str();
				sprintf(buffer, filename_template, i_name.c_str(), t);
				std::string buffer2 = path_ref + "/" + std::string(buffer);
				sweet::Data::Plane2D::PlaneData_Physical tmp(planeDataConfig[0]);
				tmp.file_physical_loadRefData_Parareal(buffer2.c_str());
				ref_data[ivar].loadPlaneDataPhysical(tmp);

				// If necessary, interpolate to coarsest spatial grid
				if (	planeDataConfig[0]->physical_res[0] != ref_data[ivar].planeDataConfig->physical_res[0] ||
					planeDataConfig[0]->physical_res[1] != ref_data[ivar].planeDataConfig->physical_res[1]
				)
				{
					SWEETErrorFatal("TODO");
					//TODO
	///				/*
	///				 * setup sampler
	///				 */
	///				PlaneDataSampler sampler2D;
	///				sampler2D.setup(simVars.sim.plane_domain_size, planeDataConfig);
	///		
	///		
	///					/*
	///					 * sample with BiLinear interpolation
	///					 */
	///					PlaneData prog_h3_bilinear(planeDataConfig3);
	///		
	///					sampler2D.bilinear_scalar(
	///							prog_h_pert,	///< input scalar field
	///							Convert_PlaneData_2_ScalarDataArray::physical_convert(px),
	///							Convert_PlaneData_2_ScalarDataArray::physical_convert(py),
	///							prog_h3_bilinear
	///					);
				}
				pint_data_ref->dataArrays_2_GenericData_PlaneData_Spectral(
											ref_data[0],
											ref_data[1]
	#if SWEET_PARAREAL_PLANE_SWE || SWEET_XBRAID_PLANE_SWE
											, ref_data[2]
	#endif
										);
				}
			}

#elif SWEET_PARAREAL_SPHERE || SWEET_XBRAID_SPHERE
		sweet::Data::Sphere2D::SphereData_Spectral ref_data[] = { sweet::Data::Sphere2D::SphereData_Spectral(sphereDataConfig[0]),
									  sweet::Data::Sphere2D::SphereData_Spectral(sphereDataConfig[0]),
									  sweet::Data::Sphere2D::SphereData_Spectral(sphereDataConfig[0])};

		for (int ivar = 0; ivar < nvar; ivar++)
		{
			std::string i_name;
			if (ivar == 0)
				i_name = "prog_phi_pert";
			else if (ivar == 1)
				i_name = "prog_vrt";
			else if (ivar == 2)
				i_name = "prog_div";

			if (shackIOData->output_file_mode == "csv")
			{
				if (iteration_id == 0)
				{
					// load ref file
					char buffer[1024];
					const char* filename_template = shackIOData->output_file_name.c_str();
					sprintf(buffer, filename_template, i_name.c_str(), t * shackIOData->output_time_scale);
					std::string buffer2 = path_ref + "/" + std::string(buffer);
					sweet::Data::Sphere2D::SphereData_Physical tmp(sphereDataConfig[0]);
					tmp.file_physical_loadRefData_Parareal(buffer2.c_str());
					ref_data[ivar].loadSphereDataPhysical(tmp);

					// If necessary, interpolate to coarsest spatial grid
					if (	sphereDataConfig[0]->physical_num_lat != ref_data[ivar].sphereDataConfig->physical_num_lat ||
						sphereDataConfig[0]->physical_num_lon != ref_data[ivar].sphereDataConfig->physical_num_lon
					)
					{
						SWEETErrorFatal("TODO");
						//TODO
	///					/*
	///					 * setup sampler
	///					 */
	///					PlaneDataSampler sampler2D;
	///					sampler2D.setup(simVars.sim.plane_domain_size, planeDataConfig);
	///	
	///	
	///						/*
	///						 * sample with BiLinear interpolation
	///						 */
	///						PlaneData prog_h3_bilinear(planeDataConfig3);
	///	
	///						sampler2D.bilinear_scalar(
	///								prog_h_pert,	///< input scalar field
	///								Convert_PlaneData_2_ScalarDataArray::physical_convert(px),
	///								Convert_PlaneData_2_ScalarDataArray::physical_convert(py),
	///								prog_h3_bilinear
	///						);
					}
					pint_data_ref->dataArrays_2_GenericData_SphereData_Spectral(ref_data[0], ref_data[1], ref_data[2]);
				}
			}
			else if (shackIOData->output_file_mode == "bin")
			{

				if (iteration_id == 0)
				{
					// load ref file
					char buffer[1024];
					const char* filename_template = shackIOData->output_file_name.c_str();
					sprintf(buffer, filename_template, i_name.c_str(), t * shackIOData->output_time_scale);
					std::string buffer2 = path_ref + "/" + std::string(buffer);
					ref_data[ivar].file_read_binary_spectral(buffer2);

					pint_data_ref->dataArrays_2_GenericData_SphereData_Spectral(ref_data[0], ref_data[1], ref_data[2]);
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
#if SWEET_PARAREAL_PLANE || SWEET_XBRAID_PLANE
				int rnorm = planeDataConfig[0]->spectral_data_size[0] / std::pow(2, ip);
				if (rnorm >= 1)
					rnorms.push_back(rnorm);
#elif SWEET_PARAREAL_SPHERE || SWEET_XBRAID_SPHERE
				int rnorm = sphereDataConfig[0]->spectral_modes_m_max / std::pow(2, ip);
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


#elif SWEET_PARAREAL_PLANE || SWEET_XBRAID_PLANE

			if (ivar == 0)
	#if SWEET_PARAREAL_PLANE_SWE || SWEET_XBRAID_PLANE_SWE
				i_name = "prog_h_pert";
	#elif SWEET_PARAREAL_PLANE_BURGERS || SWEET_XBRAID_PLANE_BURGERS
				i_name = "prog_u";
	#endif

			else if (ivar == 1)
	#if SWEET_PARAREAL_PLANE_SWE || SWEET_XBRAID_PLANE_SWE
				i_name = "prog_u";
	#elif SWEET_PARAREAL_PLANE_BURGERS || SWEET_XBRAID_PLANE_BURGERS
				i_name = "prog_v";
	#endif

			else if (ivar == 2)
				i_name = "prog_v";

			resx_data = planeDataConfig[0]->physical_res[0];
			resy_data = planeDataConfig[0]->physical_res[1];

			sweet::Data::Plane2D::PlaneData_Spectral diff_spectral = *i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[ivar]
									- *pint_data_ref->get_pointer_to_data_PlaneData_Spectral()->simfields[ivar];
			sweet::Data::Plane2D::PlaneData_Physical diff = i_data->get_pointer_to_data_PlaneData_Spectral()->simfields[ivar]->toPhys() -
									pint_data_ref->get_pointer_to_data_PlaneData_Spectral()->simfields[ivar]->toPhys();
			err_L1 = diff.physical_reduce_norm1() / (resx_data * resy_data);
			err_L2 = diff.physical_reduce_norm2() / std::sqrt(resx_data * resy_data);
			err_Linf = diff.physical_reduce_max_abs();


			// Spectral space
			double small = 1e-20;
			for (std::vector<std::size_t>::iterator it = rnorms.begin(); it != rnorms.end(); it++)
			{
				double norm_diff = std::sqrt(diff_spectral.spectral_reduce_max_abs(*it) );
				double norm_ref = std::sqrt(pint_data_ref->get_pointer_to_data_PlaneData_Spectral()->simfields[ivar]->spectral_reduce_max_abs(*it) );
				if (norm_diff < small and norm_ref < small)
					err_Linf_spectral.emplace(std::make_pair(*it, 0.));
				else
					err_Linf_spectral.emplace(std::make_pair(*it, norm_diff / norm_ref ));
			}

#elif SWEET_PARAREAL_SPHERE || SWEET_XBRAID_SPHERE
			if (ivar == 0)
				i_name = "prog_phi_pert";
			else if (ivar == 1)
				i_name = "prog_vrt";
			else if (ivar == 2)
				i_name = "prog_div";

			resx_data = sphereDataConfig[0]->physical_num_lon;
			resy_data = sphereDataConfig[0]->physical_num_lat;

			sweet::Data::Sphere2D::SphereData_Spectral diff_spectral = *i_data->get_pointer_to_data_SphereData_Spectral()->simfields[ivar]
									- *pint_data_ref->get_pointer_to_data_SphereData_Spectral()->simfields[ivar];
			sweet::Data::Sphere2D::SphereData_Physical diff = i_data->get_pointer_to_data_SphereData_Spectral()->simfields[ivar]->toPhys() -
									pint_data_ref->get_pointer_to_data_SphereData_Spectral()->simfields[ivar]->toPhys();
			err_L1 = diff.physical_reduce_norm1() / (resx_data * resy_data);
			err_L2 = diff.physical_reduce_norm2() / std::sqrt(resx_data * resy_data);
			err_Linf = diff.physical_reduce_max_abs();

			// Spectral space
			///double small = 1e-20;
			double small = 1e-16;
			for (std::vector<std::size_t>::iterator it = rnorms.begin(); it != rnorms.end(); it++)
			{
				double norm_diff = std::sqrt(diff_spectral.spectral_reduce_max_abs(*it));
				double norm_ref = std::sqrt(pint_data_ref->get_pointer_to_data_SphereData_Spectral()->simfields[ivar]->spectral_reduce_max_abs(*it));
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
			sprintf(buffer_out, filename_template_out, base_solution.c_str(), i_name.c_str(), t * shackIOData->output_time_scale, iteration_id);

			std::ofstream file(buffer_out, std::ios_base::trunc);
			file << std::setprecision(i_precision);

			file << "#BASESOLUTION " << base_solution << " " << path_ref << std::endl;
			file << "#VAR " << i_name << std::endl;
			file << "#ITERATION " << iteration_id << std::endl;
			file << "#TIMESLICE " << time_slice_id << std::endl;
			file << "#TIMEFRAMEEND " << t  * shackIOData->output_time_scale << std::endl;
			file << "errL1 " << err_L1 << std::endl;
			file << "errL2 " << err_L2 << std::endl;
			file << "errLinf " << err_Linf << std::endl;

			file.close();


#if !(SWEET_PARAREAL_SCALAR || SWEET_XBRAID_SCALAR)
			// save spectral errors in file
			char buffer_out_spec[1024];

			std::string str2 = pint_type + "_error_spec_%s_%s_t%020.8f_iter%03d.csv";
			const char* filename_template_out_spec = str2.c_str();
			sprintf(buffer_out_spec, filename_template_out_spec, base_solution.c_str(), i_name.c_str(), t * shackIOData->output_time_scale, iteration_id);

			std::ofstream file_spec(buffer_out_spec, std::ios_base::trunc);
			file_spec << std::setprecision(i_precision);

			file_spec << "#BASESOLUTION " << base_solution << " " << path_ref << std::endl;
			file_spec << "#VAR " << i_name << std::endl;
			file_spec << "#ITERATION " << iteration_id << std::endl;
			file_spec << "#TIMESLICE " << time_slice_id << std::endl;
			file_spec << "#TIMEFRAMEEND " << t  * shackIOData->output_time_scale << std::endl;
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
