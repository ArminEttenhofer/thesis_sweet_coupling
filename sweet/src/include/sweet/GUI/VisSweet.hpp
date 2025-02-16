/*
 * VisSweet.hpp
 *
 *  Created on: 30 Jun 2015
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef INCLUDE_SWEET_GUI_VISSWEET_HPP
#define INCLUDE_SWEET_GUI_VISSWEET_HPP

#include <limits>
#include <string>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <libgl/draw/GlDrawCube.hpp>
#include <libgl/shaders/shader_blinn/CShaderBlinn.hpp>
#include <libgl/VisualizationEngine.hpp>
#include <libgl/core/GlTexture.hpp>
#include <libgl/draw/GlDrawQuad.hpp>
#include <sweet/GUI/VisSweetHUD.hpp>

#ifndef SWEET_USE_SPHERE2D_SPECTRAL_SPACE
#define SWEET_USE_SPHERE2D_SPECTRAL_SPACE 1
#endif

#if SWEET_USE_SPHERE2D_SPECTRAL_SPACE
	#include <libgl/draw/GlDrawSphereSph.hpp>
#endif
#include <libgl/hud/GlFreeType.hpp>
#include <libgl/hud/GlRenderOStream.hpp>

#include <sweet/Data/Cart2D/DataGrid.hpp>

namespace sweet {
namespace GUI {

class SimulationGUICallbacks
{
public:
	/**
	 * postprocessing of frame: do time stepping
	 */
	virtual
	void vis_post_frame_processing(int i_num_iterations) = 0;

	virtual
	void vis_getDataArray(
			const sweet::Data::Cart2D::DataGrid **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data,
			double *o_vis_min,
			double *o_vis_max,
			bool *vis_reset
	) = 0;

	virtual
	const std::string vis_getStatusString(bool &o_replace_commas_with_newline) = 0;


	virtual
	void vis_pause() = 0;


	virtual
	void vis_keypress(int i_key) = 0;

	virtual
	bool should_quit() = 0;

	virtual
	bool runTimestep() = 0;

	virtual
	bool reset() = 0;

	~SimulationGUICallbacks() {}
};



/**
 * A visualization class specifically designed for SWEET applications
 */
class VisSweet	:
		public libgl::VisualizationEngine::ProgramCallbacks
{
	/**
	 * Simulation class from the SWEET-using application
	 *
	 * Certain interfaces have to be implemented, see other programs for these interfaces
	 */
	SimulationGUICallbacks *simCallbacks;

	libgl::GlDrawQuad *glDrawQuad;

#if SWEET_USE_SPHERE2D_SPECTRAL_SPACE
	libgl::GlDrawSphereSph *glDrawSphereSph;
#endif

	libgl::GlTexture *glTexture = nullptr;
	unsigned char *texture_data = nullptr;

	bool hud_visible = true;
	VisSweetHUD *visSweetHUD;

	int sim_runs_per_frame = 1;

	libgl::VisualizationEngine *visualizationEngine;

	double vis_min = 0;
	double vis_max = 0;

	int viewport_width = -1, viewport_height = -1;

	double font_size = -1;

	int screenshot_number = 0;

	void vis_setup(libgl::VisualizationEngine *i_visualizationEngine)
	{
		visualizationEngine = i_visualizationEngine;

		glDrawQuad = new libgl::GlDrawQuad;
#if SWEET_USE_SPHERE2D_SPECTRAL_SPACE
		glDrawSphereSph = nullptr;
#endif

		update_font();

		visSweetHUD = new VisSweetHUD;
		visSweetHUD->setup();
		visSweetHUD->assembleGui();
		visSweetHUD->setHudVisibility(hud_visible);

		screenshot_number = 0;
	}



	bool vis_shouldQuit()
	{
		return simCallbacks->should_quit();
	}



	void vis_render()
	{
		const sweet::Data::Cart2D::DataGrid *ro_visCart2DData;
		double aspect_ratio = 0;
		int render_primitive = 0;
		void *bogus_data;

		vis_min = std::numeric_limits<double>::infinity();
		vis_max = std::numeric_limits<double>::infinity();
		bool reset = false;


		simCallbacks->vis_getDataArray(
				&ro_visCart2DData,
				&aspect_ratio,
				&render_primitive,
				&bogus_data,
				&vis_min,
				&vis_max,
				&reset
		);

		sweet::Data::Cart2D::DataGrid &visData = (sweet::Data::Cart2D::DataGrid&)*ro_visCart2DData;

		if (std::isinf(vis_min))
		{
			vis_min = visData.grid_reduce_min();
			vis_max = visData.grid_reduce_max();

			vis_max = std::max(vis_max, vis_min+1e-20);	//< avoid numerical issues if min == max
		}


		if (glTexture == nullptr || reset)
		{
			delete glTexture;
			glTexture = new libgl::GlTexture(GL_TEXTURE_2D, GL_RED, GL_RED, GL_UNSIGNED_BYTE);
			glTexture->bind();
			glTexture->resize(visData.cart2DDataConfig->grid_data_size[0], visData.cart2DDataConfig->grid_data_size[1]);
			glTexture->unbind();

			texture_data = new unsigned char[visData.cart2DDataConfig->grid_number_elements];
		}

		double real_delta = vis_max-vis_min;
		double inv_delta = 1.0/real_delta;

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (	std::size_t i = 0;
				i < visData.cart2DDataConfig->grid_number_elements;
				i++
		)
		{
			double value = (visData.grid_space_data[i]-vis_min)*inv_delta;
			value *= 255.0;

			texture_data[i] = (unsigned char)std::min(255.0, std::max(0.0, value));//			texture_data[i] = 128;
		}

		glTexture->bind();
		glTexture->setData(texture_data);

			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

			visualizationEngine->engineState->commonShaderPrograms.shaderTexturizeRainbowmap.use();

			if (render_primitive == 0)
			{
				double scale_x = 1.0;
				double scale_y = aspect_ratio;

				if (aspect_ratio > 1.0)
				{
					scale_x /= aspect_ratio;
					scale_y = 1.0;
				}

				visualizationEngine->engineState->commonShaderPrograms.shaderTexturizeRainbowmap.pvm_matrix_uniform.set(
						visualizationEngine->engineState->matrices.pvm*libgl::scale((float)scale_x, (float)scale_y, (float)1.0)
				);

				glDrawQuad->renderWithoutProgram();
			}
			else
			{
#if SWEET_USE_SPHERE2D_SPECTRAL_SPACE
				libgl::VisualizationEngine::EngineState::Matrices &m = visualizationEngine->engineState->matrices;
				visualizationEngine->engineState->commonShaderPrograms.shaderTexturizeRainbowmap.pvm_matrix_uniform.set(
						m.pvm
				);

				/*
				 * Initialize on the first time we do this
				 */
				if (glDrawSphereSph == nullptr)
				{
					sweet::Data::Sphere2D::Config *sphere2DDataConfig = (sweet::Data::Sphere2D::Config*)bogus_data;

					glDrawSphereSph = new libgl::GlDrawSphereSph;
					glDrawSphereSph->initSphere(sphere2DDataConfig);
				}

				glDrawSphereSph->renderWithoutProgram();
#endif
			}

			visualizationEngine->engineState->commonShaderPrograms.shaderTexturizeRainbowmap.disable();

		glTexture->unbind();

		if (hud_visible)
		{
			//glFreeType->viewportChanged(visualizationEngine->renderWindow->window_width, visualizationEngine->renderWindow->window_height);

			bool o_replace_comma_with_linebreaks = true;
			std::string status_string = simCallbacks->vis_getStatusString(o_replace_comma_with_linebreaks);
			if (o_replace_comma_with_linebreaks)
				std::replace(status_string.begin(), status_string.end(), ',', '\n');

			//glFreeType->setPosition(10, visualizationEngine->renderWindow->window_height-font_size-10);
			//glFreeType->renderString(status_string.c_str());

			visSweetHUD->render();
		}

		// execute simulation time step
		simCallbacks->vis_post_frame_processing(sim_runs_per_frame);

		if (visSweetHUD->take_screenshot_series && visSweetHUD->run_simulation_timesteps)
		{
			std::ostringstream oss;
			oss << "output_screenshot_" << std::setw(8) << std::setfill('0') << screenshot_number << ".bmp";
			visualizationEngine->renderWindow->saveScreenshot(oss.str());
			screenshot_number++;
		}
	}


	const std::string vis_getStatusString(bool &o_replace_commas_with_newline)
	{
		std::ostringstream ss;
		ss <<  simCallbacks->vis_getStatusString(o_replace_commas_with_newline);
		ss.precision(12);
		ss << ", vis min/max: " << vis_min << ", " << vis_max;
		return ss.str();
	}

	void update_font()
	{
		// Approx. desired font size on standard display and standard viewport height
		font_size = 14.0;

		// rescale with viewport size
		font_size *= std::max(viewport_width, viewport_height)/(2.0*800.0);

		// Limit font size
		font_size = std::min(50.0, font_size);
		font_size = std::max(14.0, font_size);

		// Incorporate HiDPI information
		float ddpi;
		SDL_GetDisplayDPI(0, &ddpi, nullptr, nullptr);
		font_size *= ddpi/100.0;
	}


	void vis_viewportChanged(int i_width, int i_height)
	{
		viewport_width = i_width;
		viewport_height = i_height;

		update_font();
	}


	void vis_keypress(char i_key)
	{
		if (i_key >= '1' && i_key <= '9')
		{
			sim_runs_per_frame = std::pow(2, i_key-'1');
			return;
		}

		switch(i_key)
		{
		case 'r':
			simCallbacks->reset();
			screenshot_number = 0;
			break;

		case ' ':
			simCallbacks->vis_pause();
			visSweetHUD->run_simulation_timesteps = !visSweetHUD->run_simulation_timesteps;
			break;

		case SDLK_BACKSPACE:
			hud_visible = !hud_visible;
			visSweetHUD->setHudVisibility(hud_visible);
			break;

		case 'j':
			simCallbacks->runTimestep();
			break;

		default:
			simCallbacks->vis_keypress(i_key);
			break;
		}
	}


	bool vis_mouse_motion(
			int i_x,
			int i_y
	)
	{
		visSweetHUD->mouse_motion(i_x, viewport_height-i_y);
		return false;
	}

	bool vis_mouse_button_up(
			int i_button
	)
	{
		visSweetHUD->mouse_button_up(i_button);
		return false;
	}

	bool vis_mouse_button_down(
			int i_button
	)
	{
		visSweetHUD->mouse_button_down(i_button);
		return false;
	}

	bool vis_mouse_wheel(
			int i_x,
			int i_y
	)
	{
		visSweetHUD->mouse_wheel(i_x, viewport_height-i_y);
		return false;
	}

	void vis_shutdown()
	{
		delete visSweetHUD;

		delete [] texture_data;
		delete glTexture;

#if SWEET_USE_SPHERE2D_SPECTRAL_SPACE
		delete glDrawSphereSph;
#endif
		delete glDrawQuad;
	}



public:
	VisSweet(SimulationGUICallbacks &i_simCallbacks)	:
		glDrawQuad(nullptr),
#if SWEET_USE_SPHERE2D_SPECTRAL_SPACE
		glDrawSphereSph(nullptr),
#endif

		visSweetHUD(nullptr)
	{
		simCallbacks = &i_simCallbacks;

		libgl::VisualizationEngine(this, "SWEET");

	}

};

}}

#endif
