/*
 * SPHSolver.hpp
 *
 *  Created on: 24 Aug 2016
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMEHELPERS_SPHBANDEDMATRIX_GRIDREAL_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMEHELPERS_SPHBANDEDMATRIX_GRIDREAL_HPP

#include <sweet/Data/Sphere2D/Helpers/Helpers_SPHIdentities.hpp>
#include <sweet/LibMath/BandedMatrixGridReal.hpp>
#include <sweet/LibMath/BandedMatrixSolverComplex.hpp>



/*!
 * phi(lambda,mu) denotes the solution
 *
 * WARNING: This is only for real-valued physical space!
 */
class SphBandedMatrix_GridReal	:
		sweet::Data::Sphere2D::Helpers_SPHIdentities
{
	typedef std::complex<double> T;

public:
	/*!
	 * Matrix on left-hand side
	 */
	sweet::LibMath::BandedMatrixGridReal<T> lhs;


	/*!
	 * Buffer
	 */
	sweet::LibMath::BandedMatrixGridReal<T> lhsBuffer;

	/*!
	 * SPH configuration
	 */
	const sweet::Data::Sphere2D::Config *sphere2DDataConfig;
	
	/*!
	 * True if setup was called before
	 */
	bool isSetup;


	/*!
	 * Setup the SPH solver
	 */
public:
	bool setup(
			const sweet::Data::Sphere2D::Config *i_sphere2DConfig,		//!< Handler to sphere2DDataConfig
			int i_halosize_offdiagonal	//!< Size of the halo around. A value of 2 allocates data for 5 diagonals.
	)
	{
		sphere2DDataConfig = i_sphere2DConfig;

		lhs.setup(sphere2DDataConfig, i_halosize_offdiagonal);
		lhsBuffer.setup(sphere2DDataConfig, i_halosize_offdiagonal);

		isSetup = true;

		return true;
	}


	SphBandedMatrix_GridReal()	:
		sphere2DDataConfig(nullptr)
	{
	}



	~SphBandedMatrix_GridReal()
	{
		clear();
	}


	void clear()
	{
		isSetup = false;
	}

	void solver_setZero()
	{
		lhs.zeroAll();
	}

	/*!
	 * Solver for
	 * 	a*phi(lambda,mu)
	 */
	void solver_addComponent_scalar_phi(
			const std::complex<double> &i_value
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, 0, i_value);
			}
		}
	}



	/*!
	 * Solver for
	 * 	mu*phi(lambda,mu)
	 */
	void solver_addComponent_mu_phi(
			const std::complex<double> &i_scalar = 1.0
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, -1, R(n-1,m)*i_scalar);
				lhs.rowElement_add(row, n, m, +1, S(n+1,m)*i_scalar);
			}
		}
	}




	/*!
	 * Solver for
	 * 	(1-mu*mu)*d/dmu phi(lambda,mu)
	 */
	void solver_addComponent_one_minus_mu_mu_diff_mu_phi()
	{
		SWEET_ASSERT(isSetup);

		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, -1, (-(double)n+1.0)*R(n-1,m));
				lhs.rowElement_add(row, n, m, +1, ((double)n+2.0)*S(n+1,m));
			}
		}
	}



	/*!
	 * Solver for
	 * 	scalar*phi(lambda,mu)
	 */
	void solver_addComponent_rexi_z1(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		SWEET_ASSERT(isSetup);

		solver_addComponent_scalar_phi(i_scalar);
	}



	/*!
	 * Solver for
	 * 	mu^2*phi(lambda,mu)
	 */
	void solver_addComponent_rexi_z2(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		SWEET_ASSERT(isSetup);

		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				lhs.rowElement_add(row, n, m, -2, i_scalar*A(n-2,m));
				lhs.rowElement_add(row, n, m,  0, i_scalar*B(n,m));
				lhs.rowElement_add(row, n, m, +2, i_scalar*C(n+2,m));
			}
		}
	}



	/*!
	 * Solver for
	 * 	mu^4*phi(lambda,mu)
	 */
	void solver_addComponent_rexi_z3(
			const std::complex<double> &i_scalar,
			double i_r_not_required
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				lhs.rowElement_add(row, n, m, -4, i_scalar*(A(n-2,m)*A(n-4.0,m)));
				lhs.rowElement_add(row, n, m, -2, i_scalar*(A(n-2,m)*B(n-2,m) + B(n,m)*A(n-2,m)));
				lhs.rowElement_add(row, n, m,  0, i_scalar*(A(n-2,m)*C(n,m) + B(n,m)*B(n,m) + C(n+2,m)*A(n,m)));
				lhs.rowElement_add(row, n, m, +2, i_scalar*(B(n,m)*C(n+2,m) + C(n+2,m)*B(n+2,m)));
				lhs.rowElement_add(row, n, m, +4, i_scalar*(C(n+2,m)*C(n+4,m)));
			}
		}
	}


	/*!
	 * Solver for
	 * Z4 := grad_j(mu) * grad_i(phi)
	 */
	void solver_addComponent_rexi_z4(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m,  0, 1.0/(i_r*i_r)*i_scalar*T(0, m));
			}
		}
	}



	/*!
	 * Solver for
	 * Z4robert := grad_j(mu) * grad_i(phi)
	 *           = im/(r*r) * F_n^m(Phi)
	 */
	void solver_addComponent_rexi_z4robert(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			std::complex<double> fac = 1.0/(i_r*i_r)*i_scalar*T(0, m);
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m,  0, fac);
			}
		}
	}



	/*!
	 * Solver for
	 * Z5 := grad_j(mu) * mu^2 * grad_i(phi)
	 */
	void solver_addComponent_rexi_z5(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			T fac = std::complex<double>(1.0/(i_r*i_r))*i_scalar*T(0, m);

			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m,  -2, fac*A(n-2,m)	);
				lhs.rowElement_add(row, n, m,   0, fac*B(n+0,m)	);
				lhs.rowElement_add(row, n, m,  +2, fac*C(n+2,m)	);
			}
		}
	}


	/*!
	 * Solver for
	 *
	 * Z5robert := grad_j(mu) * mu^2 * grad_i(phi)
	 *           = i*m/(r*r) * F_n^m(mu^2 \phi)
	 */
	void solver_addComponent_rexi_z5robert(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			std::complex<double> fac = (i_scalar/(i_r*i_r))*std::complex<double>(0, m);

			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, -2, fac*A(n-2,m));
				lhs.rowElement_add(row, n, m,  0, fac*B(n,m));
				lhs.rowElement_add(row, n, m, +2, fac*C(n+2,m));
			}
		}
	}



	/*!
	 * Solver for
	 * Z6 := grad_j(mu) * mu * grad_j(phi)
	 */
	void solver_addComponent_rexi_z6(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m,  -2, -1.0/(i_r*i_r)*i_scalar*(D(n-1,m)*R(n-2,m) + (E(n,m)-3.0)*A(n-2,m)));
				lhs.rowElement_add(row, n, m,   0, -1.0/(i_r*i_r)*i_scalar*(1.0+D(n-1,m)*S(n,m)+(E(n,m)-3.0)*B(n,m)));
				lhs.rowElement_add(row, n, m,  +2, -1.0/(i_r*i_r)*i_scalar*(E(n,m)-3.0)*C(n+2,m));
			}
		}
	}



	/*!
	 * Solver for
	 * Z6robert := grad_j(mu) * mu * grad_j(phi)
	 *           =
	 */
	void solver_addComponent_rexi_z6robert(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		SWEET_ASSERT(isSetup);

		/*
		 * First part
		 */
		// phi
		solver_addComponent_rexi_z1(-i_scalar/(i_r*i_r), i_r);

		// mu^2*phi
		solver_addComponent_rexi_z2(3.0*i_scalar/(i_r*i_r), i_r);


		/*
		 * Second part
		 */
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			std::complex<double> fac = -1.0/(i_r*i_r)*i_scalar;
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				lhs.rowElement_add(row, n, m,  -2, fac*(G(n-1,m)*R(n-2,m))	);
				lhs.rowElement_add(row, n, m,   0, fac*(G(n-1,m)*S(n,m) + H(n+1,m)*R(n,m))	);
				lhs.rowElement_add(row, n, m,  +2, fac*(H(n+1,m)*S(n+2,m))	);
			}
		}
	}



	/*!
	 * Solver for
	 * Z7 := laplace(phi)
	 */
	void solver_addComponent_rexi_z7(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add(row, n, m, 0, -1.0/(i_r*i_r)*i_scalar*(double)n*((double)n+1.0));
			}
		}
	}


	/*!
	 * Solver for
	 * Z8 := mu*mu*laplace(phi)
	 */
	void solver_addComponent_rexi_z8(
			const std::complex<double> &i_scalar,
			double i_r
	)
	{
		SWEET_ASSERT(isSetup);

		std::complex<double> fac = (1.0/(i_r*i_r))*i_scalar;

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				// eq (D) in REXI SH document section 8.1.11 ver. 15
				std::complex<double> a = fac*(-(double)n*((double)n+1.0));
				lhs.rowElement_add(row, n, m, -2, a*A(n-2,m));
				lhs.rowElement_add(row, n, m,  0, a*B(n,m));
				lhs.rowElement_add(row, n, m, +2, a*C(n+2,m));

				// eq (Ba+Aa) in REXI SH document section 8.1.11 ver. 15
				lhs.rowElement_add(row, n, m, -2, 4.0*fac*( R(n-2,m)*G(n-1,m)) );
				lhs.rowElement_add(row, n, m,  0, 4.0*fac*( S(n,m)*G(n-1,m) + R(n,m)*H(n+1,m)) );
				lhs.rowElement_add(row, n, m, +2, 4.0*fac*( S(n+2,m)*H(n+1,m)) );

				// eq (Ab) in REXI SH document section 8.1.11 ver. 15
				lhs.rowElement_add(row, n, m, 0, 2.0*fac);

				lhs.rowElement_add(row, n, m, -2, -6.0*fac*A(n-2,m));
				lhs.rowElement_add(row, n, m,  0, -6.0*fac*B(n,m));
				lhs.rowElement_add(row, n, m, +2, -6.0*fac*C(n+2,m));
			}
		}
	}


	/*!
	 * Solver for implicit time integration
	 */
	void solver_addComponent_implicit_FJinvF(
			const double &i_dt_two_omega
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				// out of boundary check for P(n-2, m)
				if (n-2 >= m)
				{
					Tcomplex a =
							i_dt_two_omega * i_dt_two_omega
							* implicit_f_minus(n, m)
							/ implicit_J_scalar(n-1, m, i_dt_two_omega)
							* implicit_f_minus(n-1, m);

					lhs.rowElement_add_NEW(row, n, m, -2, a);
				}


				{
					// out of boundary check for P(n, m) (at least something like that)
					if (n > m)
					{
						Tcomplex a =
								i_dt_two_omega * i_dt_two_omega
								* implicit_f_minus(n, m)
								/ implicit_J_scalar(n-1, m, i_dt_two_omega)
								* implicit_f_plus(n-1, m);

						lhs.rowElement_add_NEW(row, n, m, 0, a);
					}

					// out of boundary check for P(n, m) (at least something like that)
					if (n < sphere2DDataConfig->spectral_modes_n_max)
					{
						Tcomplex a =
								i_dt_two_omega * i_dt_two_omega
								* implicit_f_plus(n, m)
								/ implicit_J_scalar(n+1, m, i_dt_two_omega)
								* implicit_f_minus(n+1, m);

						lhs.rowElement_add_NEW(row, n, m, 0, a);
					}

				}

				// out of boundary check for P(n+2, m)
				if (n+2 <= sphere2DDataConfig->spectral_modes_n_max)
				{
					Tcomplex a =
							i_dt_two_omega * i_dt_two_omega
							* implicit_f_plus(n, m)
							/ implicit_J_scalar(n+1, m, i_dt_two_omega)
							* implicit_f_plus(n+1, m);

					lhs.rowElement_add_NEW(row, n, m, 2, a);
				}
			}
		}
	}



	/*!
	 * Solver for implicit time integration
	 */
	void solver_addComponent_implicit_F(
			const double &i_dt_two_imega
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				// out of boundary check for P(n-1, m)
				if (n-1 >= m)
				{
					Tcomplex a = implicit_f_minus(n, m) * i_dt_two_imega;
					lhs.rowElement_add_NEW(row, n, m, -1, a);
				}

				// out of boundary check for P(n+1, m)
				if (n+1 <= sphere2DDataConfig->spectral_modes_n_max)
				{
					Tcomplex a = implicit_f_plus(n, m) * i_dt_two_imega;
					lhs.rowElement_add_NEW(row, n, m, 1, a);
				}
			}
		}
	}



	/*!
	 * Solver for implicit time integration
	 */
	void solver_addComponent_implicit_FJinv(
			const double &i_dt_two_imega
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{

			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				// out of boundary check for P(n-1, m)
				if (n-1 >= m)
				{
					Tcomplex a =	i_dt_two_imega
									* implicit_f_minus(n, m)
									/ implicit_J_scalar(n-1, m, i_dt_two_imega);

					lhs.rowElement_add_NEW(row, n, m, -1, a);
				}

				// out of boundary check for P(n+1, m)
				if (n+1 <= sphere2DDataConfig->spectral_modes_n_max)
				{
					// sort out div/0 from implicit_J_scalar
					//if (n > 0)
					{
						Tcomplex a =	i_dt_two_imega
										* implicit_f_plus(n, m)
										/ implicit_J_scalar(n+1, m, i_dt_two_imega);

						lhs.rowElement_add_NEW(row, n, m, 1, a);
					}
				}
			}
		}
	}



	/*!
	 * Solver for implicit time integration
	 */
	void solver_addComponent_implicit_J(
			const double &i_dt_two_imega
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				Tcomplex a = implicit_J_scalar(n, m, i_dt_two_imega);
				lhs.rowElement_add_NEW(row, n, m, 0, a);
			}
		}
	}


	/*!
	 * Add the identity matrix times a scalar
	 */
	void solver_addComponent_implicit_I(
			const double &i_scalar
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);
				lhs.rowElement_add_NEW(row, n, m, 0, i_scalar);
			}
		}
	}



	/*!
	 * Solver for implicit time integration
	 */
	void solver_addComponent_implicit_L(
			const double i_scalar,
			const double i_dt,		// note, that the dt is also in L. That's not a bug
			const double i_radius
	)
	{
		SWEET_ASSERT(isSetup);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				T *row = lhs.getMatrixRow(n, m);

				Tcomplex a = i_scalar*i_dt*(n*(n+1.0))/(i_radius*i_radius);
				lhs.rowElement_add_NEW(row, n, m, 0, a);
			}
		}
	}


	/*!
	 * Output matrix
	 */
	void print()
	{
		std::cout << "SphBandedMatrix_GridReal:" << std::endl;
		std::cout << " + lhs:" << std::endl;
		lhs.print();
	}


	/*!
	 * Apply the solver matrix.
	 *
	 * This function is intended to be used for debugging.
	 * This only multiplies the i_x values with the matrix.
	 * Use solve(...) to solve for the matrix
	 */
	sweet::Data::Sphere2D::DataSpectral apply(
			const sweet::Data::Sphere2D::DataSpectral &i_x	//!< solution to be searched
	)
	{
		SWEET_ASSERT(isSetup);

		sweet::Data::Sphere2D::DataSpectral out(sphere2DDataConfig);

		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			std::size_t idx = sphere2DDataConfig->getArrayIndexByModes(m, m);
			for (int n = m; n <= sphere2DDataConfig->spectral_modes_n_max; n++)
			{
				out.spectral_space_data[idx] = 0;

				std::complex<double> *row = lhs.getMatrixRow(n, m);
				for (int i = 0; i < lhs.num_diagonals; i++)
				{
					int delta = i-lhs.halosize_off_diagonal;

					int n_rel = n+delta;

					if (n_rel < 0 ||  m < 0)
						continue;

					if (n_rel > sphere2DDataConfig->spectral_modes_n_max)
						continue;

					if (m > sphere2DDataConfig->spectral_modes_m_max)
						continue;

					if (m > n_rel)
						continue;

					out.spectral_space_data[idx] += lhs.rowElement_getRef(row, n, m, delta)*i_x.spectral_get_(n_rel, m);
				}

				idx++;
			}
		}

		return out;
	}



	sweet::Data::Sphere2D::DataSpectral solve(
			const sweet::Data::Sphere2D::DataSpectral &i_rhs
	)	const
	{
		SWEET_ASSERT(isSetup);

		sweet::Data::Sphere2D::DataSpectral out(sphere2DDataConfig);

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (int m = 0; m <= sphere2DDataConfig->spectral_modes_m_max; m++)
		{
			int idx = sphere2DDataConfig->getArrayIndexByModes(m,m);

			int num_rows = sphere2DDataConfig->spectral_modes_n_max+1-m;

			std::memcpy(
					(void*)&lhsBuffer.data[idx*lhs.num_diagonals],
					(void*)&lhs.data[idx*lhs.num_diagonals],
					sizeof(std::complex<double>)*num_rows*lhs.num_diagonals
				);

			std::memcpy(
					(void*)&out.spectral_space_data[idx],
					(void*)&i_rhs.spectral_space_data[idx],
					sizeof(std::complex<double>)*num_rows
				);

			sweet::LibMath::BandedMatrixSolverComplex::solve(
							&lhsBuffer.data[idx*lhs.num_diagonals],
							num_rows,
							lhs.halosize_off_diagonal,
							&out.spectral_space_data[idx]
					);
		}

		return out;
	}
};


#endif
