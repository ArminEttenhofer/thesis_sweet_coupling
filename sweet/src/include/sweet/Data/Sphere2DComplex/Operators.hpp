/*
 *  Created on: 31 Aug 2016
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_DATA_SPHERE2DCOMPLEX_OPERATORS_HPP
#define INCLUDE_SWEET_DATA_SPHERE2DCOMPLEX_OPERATORS_HPP

#include <sweet/Error/Base.hpp>
#include <sweet/Data/Sphere2D/Config.hpp>
#include <sweet/Data/Sphere2DComplex/DataSpectral.hpp>
#include <sweet/Data/Sphere2DComplex/DataGrid.hpp>
#include <sweet/Data/Sphere2D/Helpers/Helpers_SPHIdentities.hpp>
#include <sweet/Data/Sphere2D/Shack.hpp>
#include <sweet/Data/Sphere2DComplex/DataSpectral.hpp>
#include <sweet/LibMath/shtns_inc.hpp>

#if SWEET_THREADING
#include <omp.h>
#endif

namespace sweet {
namespace Data {
namespace Sphere2DComplex {


class Operators	:
		public Sphere2D::Helpers_SPHIdentities
{
	friend Sphere2D::Config;

public:
	Error::Base error;

	const Sphere2D::Config *sphere2DDataConfig;

private:
	double r;			///! radius
	double ir;			///! 1/radius


public:
	Operators(
		const Sphere2D::Config *i_sphere2DDataConfig,
		const Sphere2D::Shack *i_shackSphere2DDataOps
	);


public:
	Operators();


public:
	void setup(
			const Sphere2D::Config *i_sphere2DDataConfig,
			const Sphere2D::Shack *i_shackSphere2DDataOps
	);


	/*!
	 * Solve a Helmholtz problem given by
	 *
	 * (I + b D^2) x = rhs
	 */
public:
	DataSpectral implicit_helmholtz(
			const Sphere2DComplex::DataSpectral &i_sphere2d_data,
			const std::complex<double> &i_b,
			double i_radius
	)	const;

	/*!
	 * Compute multiplication with "J" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
public:
	DataSpectral implicit_J(
			const Sphere2DComplex::DataSpectral &i_sphere2d_data,
			const std::complex<double>& i_dt_two_omega
	)	const;


	/*!
	 * Compute multiplication with "J^-1" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
public:
	DataSpectral implicit_Jinv(
			const Sphere2DComplex::DataSpectral &i_sphere2d_data,
			const std::complex<double>& i_dt_two_omega
	)	const;


	/*!
	 * Compute multiplication with "F" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
public:
	DataSpectral implicit_F(
			const Sphere2DComplex::DataSpectral &i_sphere2d_data,
			const std::complex<double>& i_dt_two_omega
	)	const;



	/*!
	 * Compute multiplication with "F" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
public:
	DataSpectral implicit_FJinv(
			const Sphere2DComplex::DataSpectral &i_sphere2d_data,
			const std::complex<double>& i_dt_two_omega
	)	const;


	/*!
	 * Compute multiplication with "L" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
public:
	DataSpectral implicit_L(
			const Sphere2DComplex::DataSpectral &i_sphere2d_data,
			const std::complex<double>& i_dt
	)	const;


	/*!
	 * Compute multiplication with "Linv" linear operator used for implicit time integration
	 * see Temperton "Coriolis Terms in SL spectral models"
	 */
public:
	DataSpectral implicit_Linv(
			const Sphere2DComplex::DataSpectral &i_sphere2d_data,
			const std::complex<double>& i_dt
	)	const;

	/*!
	 * Compute differential along longitude
	 *
	 * d/d lambda f(lambda,mu)
	 */
public:
	DataSpectral diff_lon(
			const Sphere2DComplex::DataSpectral &i_sph_data
	)	const;

	/*!
	 * Compute
	 * mu*F(\lambda,\mu)
	 */
public:
	DataSpectral mu(
			const Sphere2DComplex::DataSpectral &i_sph_data
	);

	/*!
	 * Compute
	 * mu*F(\lambda,\mu)
	 */
public:
	DataSpectral mu2(
			const Sphere2DComplex::DataSpectral &i_sph_data
	)	const;


	/*!
	 * Laplace operator
	 */
public:
	DataSpectral laplace(
			const Sphere2DComplex::DataSpectral &i_sph_data
	)	const;


	/*!
	 * Convert vorticity/divergence field to u,v velocity field
	 */
public:
	void vrtdiv_2_uv(
			const Sphere2DComplex::DataSpectral &i_vrt,
			const Sphere2DComplex::DataSpectral &i_div,
			Sphere2DComplex::DataGrid &o_u,
			Sphere2DComplex::DataGrid &o_v

	)	const;


	/*!
	 * Convert u,v velocity field to vorticity/divergence field
	 */
public:
	void uv_2_vrtdiv(
			const Sphere2DComplex::DataGrid &i_u,
			const Sphere2DComplex::DataGrid &i_v,
			Sphere2DComplex::DataSpectral &o_vrt,
			Sphere2DComplex::DataSpectral &o_div

	)	const;


	/*!
	 * Laplace operator
	 */
public:
	DataSpectral inv_laplace(
			const Sphere2DComplex::DataSpectral &i_sph_data,
			double i_radius = -1
	)	const;
};

}}}

#endif
