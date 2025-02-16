/*
 * Cart2DOperators.hpp
 *
 *  Created on: 30 Jun 2015
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef INCLUDE_SWEET_DATA_CART2D_OPERATORS_HPP
#define INCLUDE_SWEET_DATA_CART2D_OPERATORS_HPP


#if SWEET_USE_CART2D_SPECTRAL_SPACE
//	#include <sweet/Data/Cart2D/Cart2DDataComplex.hpp>
#endif

#include <sweet/Error/Base.hpp>
#include <sweet/Error/Fatal.hpp>
#include <sweet/Data/Cart2D/DataSpectral.hpp>
#include <sweet/Data/Cart2D/Config.hpp>
#include <sweet/Data/Cart2D/Shack.hpp>

namespace sweet {
namespace Data {
namespace Cart2D {


/*!
 * \brief Operators for a Cart2D grid such as computing the n-th derivative
 */
class Operators
{
public:
	Error::Base error;

	Config *cart2DDataConfig;

	///! differential operators (central / forward / backward)
	DataSpectral diff_c_x, diff_c_y;
	DataSpectral diff_f_x, diff_f_y;
	DataSpectral diff_b_x, diff_b_y;

	DataSpectral diff2_c_x, diff2_c_y;

	DataGrid avg_f_x, avg_f_y;
	DataGrid avg_b_x, avg_b_y;

	DataGrid shift_left;
	DataGrid shift_right;
	DataGrid shift_up;
	DataGrid shift_down;



	Operators()	:
		cart2DDataConfig(nullptr),

		diff_c_x(1),
		diff_c_y(1),

		diff_f_x(1),
		diff_f_y(1),
		diff_b_x(1),
		diff_b_y(1),

		diff2_c_x(1),
		diff2_c_y(1),

		avg_f_x(1),
		avg_f_y(1),
		avg_b_x(1),
		avg_b_y(1),

		shift_left(1),
		shift_right(1),
		shift_up(1),
		shift_down(1)
	{

	}


	Operators(
		Config *i_cart2DDataConfig,		//!< data config setup for spectral transformations
		const double i_domain_size[2],			//!< domain size
		bool i_use_spectral_basis_diffs = true	//!< use spectral differentiation (d/dx e^ix)
	)	 __attribute__ ((deprecated)):
		cart2DDataConfig(i_cart2DDataConfig),

		diff_c_x(i_cart2DDataConfig),
		diff_c_y(i_cart2DDataConfig),

		diff_f_x(i_cart2DDataConfig),
		diff_f_y(i_cart2DDataConfig),
		diff_b_x(i_cart2DDataConfig),
		diff_b_y(i_cart2DDataConfig),

		diff2_c_x(i_cart2DDataConfig),
		diff2_c_y(i_cart2DDataConfig),

		avg_f_x(i_cart2DDataConfig),
		avg_f_y(i_cart2DDataConfig),
		avg_b_x(i_cart2DDataConfig),
		avg_b_y(i_cart2DDataConfig),

		shift_left(i_cart2DDataConfig),
		shift_right(i_cart2DDataConfig),
		shift_up(i_cart2DDataConfig),
		shift_down(i_cart2DDataConfig)
	{
		_setup(i_domain_size, i_use_spectral_basis_diffs);
	}

	Operators(
		Config *i_cart2DDataConfig,		//!< data config setup for spectral transformations
		Shack *i_cart2DDataOps
	)	:
		cart2DDataConfig(i_cart2DDataConfig),

		diff_c_x(i_cart2DDataConfig),
		diff_c_y(i_cart2DDataConfig),

		diff_f_x(i_cart2DDataConfig),
		diff_f_y(i_cart2DDataConfig),
		diff_b_x(i_cart2DDataConfig),
		diff_b_y(i_cart2DDataConfig),

		diff2_c_x(i_cart2DDataConfig),
		diff2_c_y(i_cart2DDataConfig),

		avg_f_x(i_cart2DDataConfig),
		avg_f_y(i_cart2DDataConfig),
		avg_b_x(i_cart2DDataConfig),
		avg_b_y(i_cart2DDataConfig),

		shift_left(i_cart2DDataConfig),
		shift_right(i_cart2DDataConfig),
		shift_up(i_cart2DDataConfig),
		shift_down(i_cart2DDataConfig)
	{
		_setup(i_cart2DDataOps->cart2d_domain_size, i_cart2DDataOps->space_use_spectral_basis_diffs);
	}

	/*!
	 * D2, e.g. for viscosity
	 */
	DataSpectral diff2(
			const DataSpectral &i_dataArray
	)
	{
		return diff2_c_x(i_dataArray) + diff2_c_y(i_dataArray);
	}



	/*!
	 *        __2
	 * apply  \/  operator (aka Laplace)
	 */
	DataSpectral laplace(
			const DataSpectral &i_a
	)
	{
		return diff2_c_x(i_a)+diff2_c_y(i_a);
	}



	/*!
	 * Vorticity
	 *
	 * vort(a,b) = db/dx - da/dy
	 */
	DataSpectral vort(
			const DataSpectral &a,
			const DataSpectral &b
	)
	{
		return diff_c_x(b) - diff_c_y(a);
	}



	/*!
	 * Divergence
	 *
	 * div(a,b) = da/dx + db/dy
	 */
	DataSpectral div(
			const DataSpectral &a,
			const DataSpectral &b
	)
	{
		return diff_c_x(a) + diff_c_y(b);
	}

	/*!
	 * kinetic energy
	 *
	 * ke(a,b) = 0.5*(a^2+b^2)
	 */
	DataSpectral ke(
			const DataSpectral &a,
			const DataSpectral &b
	)
	{
		return 0.5*(a*a+b*b);
	}



	/*
	 * Compute Arakawa Jacobian
	 * See A. Arakawa, V. R. Lamb, "A potential enstrophy and energy conserving scheme for the shallow water equations"
	 *
	 * J(a,b) = da/dx db/dy - da/dy db/dx
	 */
	DataSpectral J(
			const DataSpectral &a,
			const DataSpectral &b
	)
	{
		return diff_c_x(a)*diff_c_y(b) - diff_c_y(a)*diff_c_x(b);
	}


	/*
	 * Compute time derivative of Arakawa Jacobian
	 *
	 * J(a,b)_t = (da/dx db/dy - da/dy db/dx)_t
	 */
	DataSpectral J_t(
			const DataSpectral &a,
			const DataSpectral &b,
			const DataSpectral &a_t,
			const DataSpectral &b_t
	)
	{
		return	  diff_c_x(a_t)*diff_c_y(b)
				+ diff_c_x(a)*diff_c_y(b_t)
				- diff_c_y(a_t)*diff_c_x(b)
				- diff_c_y(a)*diff_c_x(b_t);
	}



	/*!
	 *        __
	 * apply  \/ .  operator
	 */
	DataSpectral diff_dot(
			const DataSpectral &i_a
	)
	{
		return diff_c_x(i_a)+diff_c_y(i_a);
	}



	/*!
	 * Diff N operator for hyperviscosity, see
	 * "Numerical Techniques for Global Atmospheric Models", page 500
	 */
	DataSpectral diffN_x(
			const DataSpectral &io_u,
			int i_order
	)
	{
		if (i_order == 0)
			return io_u;

		DataSpectral tu = io_u;

		for (int i = 0; i < i_order/2; i++)
			tu = diff2_c_x(tu);

		if (i_order & 1)
			tu = diff_c_x(tu);

		return tu;
	}


	/*!
	 * Diff N operator for hyperviscosity, see
	 * "Numerical Techniques for Global Atmospheric Models", page 500
	 */
	DataSpectral diffN_y(
			const DataSpectral &io_v,
			int i_order
	)
	{
		if (i_order == 0)
			return io_v;

		DataSpectral tv = io_v;

		for (int i = 0; i < i_order/2; i++)
			tv = diff2_c_y(tv);

		if (i_order & 1)
			tv = diff_c_y(tv);

		return tv;
	}


#if SWEET_USE_CART2D_SPECTRAL_SPACE
	/*!
	 * Diffusion or hyperviscosity coefficients
	 * Simply calculates the spectral coefficients
	 *  see "Numerical Techniques for Global Atmospheric Models", page 500
	 *
	 * i_order (q) needs to be even!!! (second or forth order usually)
	 *
	 * Returns operator D^q
	 *
	 */
	DataSpectral diffusion_coefficient(
			int i_order
	)
	{
		//Check if even
		SWEET_ASSERT( i_order % 2 == 0);
		SWEET_ASSERT( i_order > 0);
		DataSpectral out = diff2_c_x+diff2_c_y;

		/*
		 * Always use negative sign for hyperdiffusion to allow using always positive viscosity
		 */
		for (int i = 1; i < i_order/2; i++)
			out = -(diff2_c_x(out)+diff2_c_y(out));

		return out;
	}

	/*!
	 * Calculates implicit diffusion (applies 1/(1-mu*dt*D^q) to spectrum)
	 *  see "Numerical Techniques for Global Atmospheric Models", page 500
	 *
	 * i_order (q) needs to be even!!! (second or forth order usually)
	 * i_coef is mu*dt
	 *
	 * Only works in spectral space
	 *
	 */
	DataSpectral implicit_diffusion(
			const DataSpectral &i_data,
			double i_coef,
			int i_order
	)
	{
		DataSpectral out=i_data;

		// Get diffusion coefficients (these are the -mu*dt*D^q, where q is the order
		DataSpectral diff = -i_coef*diffusion_coefficient(i_order);

		// Add 1 to get denominator
		diff = diff.spectral_addScalarAll(1.0);

		// Invert
		diff = diff.spectral_invert();
		// apply to data
		out=diff(out);
		return out;
	}
#endif


	bool setup(

		Config &i_cart2DDataConfig,		//!< data config setup for spectral transformations
		const double i_domain_size[2],			//!< domain size
		bool i_use_spectral_basis_diffs = true	//!< use spectral differentiation (d/dx e^ix)
	)
	{
		return setup(
				&i_cart2DDataConfig,
				i_domain_size,
				i_use_spectral_basis_diffs
			);
	}

	void clear()
	{
		diff_c_x.clear();
		diff_c_y.clear();

		diff_f_x.clear();
		diff_f_y.clear();
		diff_b_x.clear();
		diff_b_y.clear();

		diff2_c_x.clear();
		diff2_c_y.clear();

		avg_f_x.clear();
		avg_f_y.clear();
		avg_b_x.clear();
		avg_b_y.clear();

		shift_left.clear();
		shift_right.clear();
		shift_up.clear();
		shift_down.clear();

		cart2DDataConfig = nullptr;
	}


	bool setup(
		Config *i_cart2DDataConfig,		//!< data config setup for spectral transformations
		const double i_domain_size[2],			//!< domain size
		bool i_use_spectral_basis_diffs = true	//!< use spectral differentiation (d/dx e^ix)
	)
	{
		SWEET_ASSERT(cart2DDataConfig == nullptr);
		cart2DDataConfig = i_cart2DDataConfig;

		diff_c_x.setup(cart2DDataConfig);
		diff_c_y.setup(cart2DDataConfig);

		diff_f_x.setup(cart2DDataConfig);
		diff_f_y.setup(cart2DDataConfig);
		diff_b_x.setup(cart2DDataConfig);
		diff_b_y.setup(cart2DDataConfig);

		diff2_c_x.setup(cart2DDataConfig);
		diff2_c_y.setup(cart2DDataConfig);

		avg_f_x.setup(cart2DDataConfig);
		avg_f_y.setup(cart2DDataConfig);
		avg_b_x.setup(cart2DDataConfig);
		avg_b_y.setup(cart2DDataConfig);

		shift_left.setup(cart2DDataConfig);
		shift_right.setup(cart2DDataConfig);
		shift_up.setup(cart2DDataConfig);
		shift_down.setup(cart2DDataConfig);

		return _setup(i_domain_size, i_use_spectral_basis_diffs);
	}

	bool setup(
			Config *i_cart2DDataConfig,		//!< data config setup for spectral transformations
			Shack *i_shackCart2DDataOps
	)
	{
		return setup(
				i_cart2DDataConfig,
				i_shackCart2DDataOps->cart2d_domain_size,
				i_shackCart2DDataOps->space_use_spectral_basis_diffs
		);
	}


	bool setup(
			Config &i_cart2DDataConfig,		//!< data config setup for spectral transformations
			Shack &i_shackCart2DDataOps
	)
	{
		return setup(
				&i_cart2DDataConfig,
				i_shackCart2DDataOps.cart2d_domain_size,
				i_shackCart2DDataOps.space_use_spectral_basis_diffs
		);
	}

	bool setup(
			Config &i_cart2DDataConfig,		//!< data config setup for spectral transformations
			Shack *i_shackCart2DDataOps
	)
	{
		return setup(
				&i_cart2DDataConfig,
				i_shackCart2DDataOps->cart2d_domain_size,
				i_shackCart2DDataOps->space_use_spectral_basis_diffs
		);
	}

	bool _setup(
			const double i_domain_size[2],
			bool i_use_spectral_basis_diffs
	)
	{

		double h[2] = {
				(double)i_domain_size[0] / (double)cart2DDataConfig->grid_res[0],
				(double)i_domain_size[1] / (double)cart2DDataConfig->grid_res[1]
		};

		double avg_f_x_kernel[3][3] = {
				{0,0,0},
				{0,1,1},
				{0,0,0},
		};
		avg_f_x.kernel_stencil_setup(avg_f_x_kernel, 0.5);

		double avg_f_y_kernel[3][3] = {
				{0,1,0},
				{0,1,0},
				{0,0,0},
		};
		avg_f_y.kernel_stencil_setup(avg_f_y_kernel, 0.5);

		double avg_b_x_kernel[3][3] = {
				{0,0,0},
				{1,1,0},
				{0,0,0},
		};
		avg_b_x.kernel_stencil_setup(avg_b_x_kernel, 0.5);

		double avg_b_y_kernel[3][3] = {
				{0,0,0},
				{0,1,0},
				{0,1,0},
		};
		avg_b_y.kernel_stencil_setup(avg_b_y_kernel, 0.5);

/////////////////////////////////////////////////////////////////////

		double shift_left_kernel[3][3] = {
				{0,0,0},
				{0,0,1},
				{0,0,0},
		};
		shift_left.kernel_stencil_setup(shift_left_kernel);

		double shift_right_kernel[3][3] = {
				{0,0,0},
				{1,0,0},
				{0,0,0},
		};
		shift_right.kernel_stencil_setup(shift_right_kernel);

		double shift_up_kernel[3][3] = {
				{0,0,0},
				{0,0,0},
				{0,1,0},
		};
		shift_up.kernel_stencil_setup(shift_up_kernel);

		double shift_down_kernel[3][3] = {
				{0,1,0},
				{0,0,0},
				{0,0,0},
		};
		shift_down.kernel_stencil_setup(shift_down_kernel);

/////////////////////////////////////////////////////////////////////

		if (i_use_spectral_basis_diffs)
		{
			/*
			 * setup spectral differential operators
			 * 		diff(e(ix), x)
			 */
			// Assume, that errors are linearly depending on the resolution
			// see test_spectral_ops.cpp

#if !SWEET_USE_CART2D_SPECTRAL_SPACE
			std::cerr << "Activate spectral space during compile time to use spectral diffs. Otherwise, the convolution would be freakingly expensive" << std::endl;
			SWEET_ASSERT(false);
			exit(-1);
#else

			/*
			 * Note, that there's a last column which is set to 0 (Nyquist freq, noise in signal)
			 * PXT: removed this setting to zero (changed < to <=), because of 2nd and higher order differentiation
			 * MaS: changed <= to < for the x-axis because of invalid memory access...
			 */
			diff_c_x.spectral_setZero();

			for (int j = cart2DDataConfig->spectral_data_iteration_ranges[0][1][0]; j < (int)cart2DDataConfig->spectral_data_iteration_ranges[0][1][1]; j++)
			{
				for (int i = cart2DDataConfig->spectral_data_iteration_ranges[0][0][0]; i < (int)cart2DDataConfig->spectral_data_iteration_ranges[0][0][1]; i++)
				{
					std::complex<double> data(0.0, ((double)i*2.0*M_PI/(double)i_domain_size[0]));
					diff_c_x.spectral_set(j, i, data);
				}
			}

			for (int j = cart2DDataConfig->spectral_data_iteration_ranges[1][1][0]; j < (int)cart2DDataConfig->spectral_data_iteration_ranges[1][1][1]; j++)
			{
				for (int i = cart2DDataConfig->spectral_data_iteration_ranges[1][0][0]; i < (int)cart2DDataConfig->spectral_data_iteration_ranges[1][0][1]; i++)
				{
					std::complex<double> data(0.0, ((double)i*2.0*M_PI/(double)i_domain_size[0]));
					diff_c_x.spectral_set(j, i, data);
				}
			}

			/*
			 * DIFF operator in y axis
			 */
			diff_c_y.spectral_setZero();

			for (int j = cart2DDataConfig->spectral_data_iteration_ranges[0][1][0]; j < (int)cart2DDataConfig->spectral_data_iteration_ranges[0][1][1]; j++)
			{
				for (int i = cart2DDataConfig->spectral_data_iteration_ranges[0][0][0]; i < (int)cart2DDataConfig->spectral_data_iteration_ranges[0][0][1]; i++)
				{
					std::complex<double> data(0, (double)((double)j*2.0*M_PI/(double)i_domain_size[1]));
					diff_c_y.spectral_set(j, i, data);
				}
			}


			for (int j = cart2DDataConfig->spectral_data_iteration_ranges[1][1][0]; j < (int)cart2DDataConfig->spectral_data_iteration_ranges[1][1][1]; j++)
			{
				for (int i = cart2DDataConfig->spectral_data_iteration_ranges[1][0][0]; i < (int)cart2DDataConfig->spectral_data_iteration_ranges[1][0][1]; i++)
				{
					std::complex<double> data(0, -(double)((double)(cart2DDataConfig->spectral_data_size[1]-j)*2.0*M_PI/(double)i_domain_size[1]));
					diff_c_y.spectral_set(j, i, data);
				}
			}


			/*
			 * TODO: WARNING! These operators are setup in Cartesian space,
			 * hence they are not as accurate as spectral operators
			 */
			double d_f_x_kernel[3][3] = {
					{0,0,0},
					{0,-1,1},
					{0,0,0}
			};
			diff_f_x.toGrid().kernel_stencil_setup(d_f_x_kernel, 1.0/h[0]);

			double d_f_y_kernel[3][3] = {
					{0,1,0},
					{0,-1,0},
					{0,0,0},
			};
			diff_f_y.toGrid().kernel_stencil_setup(d_f_y_kernel, 1.0/h[1]);


			double d_b_x_kernel[3][3] = {
					{0,0,0},
					{-1,1,0},
					{0,0,0}
			};
			diff_b_x.toGrid().kernel_stencil_setup(d_b_x_kernel, 1.0/h[0]);

			double d_b_y_kernel[3][3] = {
					{0,0,0},
					{0,1,0},
					{0,-1,0},
			};
			diff_b_y.toGrid().kernel_stencil_setup(d_b_y_kernel, 1.0/h[1]);


			/*
			 * 2nd order differential operators
			 */
			diff2_c_x.spectral_setZero();
			for (int r = 0; r < 2; r++)
			{
				for (int j = cart2DDataConfig->spectral_data_iteration_ranges[r][1][0]; j < (int)cart2DDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
				{
					for (int i = cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]; i < (int)cart2DDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
					{
						std::complex<double> data = diff_c_x.spectral_get(j, i);
						diff2_c_x.spectral_set(j, i, data*data);
					}
				}
			}

			diff2_c_y.spectral_setZero();
			for (int r = 0; r < 2; r++)
			{
				for (int j = cart2DDataConfig->spectral_data_iteration_ranges[r][1][0]; j < (int)cart2DDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
				{
					for (int i = cart2DDataConfig->spectral_data_iteration_ranges[r][0][0]; i < (int)cart2DDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
					{
						std::complex<double> data = diff_c_y.spectral_get(j, i);
						diff2_c_y.spectral_set(j, i, data*data);
					}
				}
			}

			diff_c_x.spectral_debugCheckForZeroAliasingModes();
			diff_c_y.spectral_debugCheckForZeroAliasingModes();
			diff2_c_x.spectral_debugCheckForZeroAliasingModes();
			diff2_c_y.spectral_debugCheckForZeroAliasingModes();
#endif
		}
		else
		{
			DataGrid tmp(cart2DDataConfig);

			double diff1_x_kernel[3][3] = {
					{0,0,0},
					{-1.0,0,1.0},
					{0,0,0}
			};
			tmp.kernel_stencil_setup(diff1_x_kernel, 1.0/(2.0*h[0]));
			diff_c_x.loadCart2DDataGrid(tmp);

			double diff1_y_kernel[3][3] = {
					{0,1.0,0},	// higher y coordinate
					{0,0,0},
					{0,-1.0,0},	// lower y coordinate
			};
			tmp.kernel_stencil_setup(diff1_y_kernel, 1.0/(2.0*h[1]));
			diff_c_y.loadCart2DDataGrid(tmp);

			double d_f_x_kernel[3][3] = {
					{0,0,0},
					{0,-1,1},
					{0,0,0}
			};
			tmp.kernel_stencil_setup(d_f_x_kernel, 1.0/h[0]);
			diff_f_x.loadCart2DDataGrid(tmp);

			double d_f_y_kernel[3][3] = {
					{0,1,0},
					{0,-1,0},
					{0,0,0},
			};
			tmp.kernel_stencil_setup(d_f_y_kernel, 1.0/h[1]);
			diff_f_y.loadCart2DDataGrid(tmp);


			double d_b_x_kernel[3][3] = {
					{0,0,0},
					{-1,1,0},
					{0,0,0}
			};
			tmp.kernel_stencil_setup(d_b_x_kernel, 1.0/h[0]);
			diff_b_x.loadCart2DDataGrid(tmp);

			double d_b_y_kernel[3][3] = {
					{0,0,0},
					{0,1,0},
					{0,-1,0},
			};
			tmp.kernel_stencil_setup(d_b_y_kernel, 1.0/h[1]);
			diff_b_y.loadCart2DDataGrid(tmp);


			double diff2_x_kernel[3][3] = {
					{0,0,0},
					{1.0,-2.0,1.0},
					{0,0,0}
				};
			tmp.kernel_stencil_setup(diff2_x_kernel, 1.0/(h[0]*h[0]));
			diff2_c_x.loadCart2DDataGrid(tmp);

			double diff2_y_kernel[3][3] = {
					{0,1.0,0},
					{0,-2.0,0},
					{0,1.0,0}
			};
			tmp.kernel_stencil_setup(diff2_y_kernel, 1.0/(h[1]*h[1]));
			diff2_c_y.loadCart2DDataGrid(tmp);
		}

		return true;
	}

};

}}}

#endif
