/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_EXPINTEGRATION_EXPFUNCTION_HPP
#define INCLUDE_SWEET_EXPINTEGRATION_EXPFUNCTION_HPP

#include <iostream>
#include <complex>
#include <typeinfo>
#include <string>
#include <sweet/Parallelization/openmp_helper.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Error/Fatal.hpp>


/*!
 * When to switch to special handling of evaluation
 */
#define EXP_FUNCTIONS_SPECIAL_THRESHOLD	1

/*!
 * Use recursive phi definition if 1, explicit formula if 0
 */
#define EXP_FUNCTIONS_USE_REC_PHI 0

/*!
 * 0: Series
 * 1: Cauchy contour integration
 */
#define EXP_FUNCTIONS_SPECIAL_EVAL_TYPE	1

/*!
 * Max number of iterations for series-based evaluation
 */
#define EXP_FUNCTIONS_SPECIAL_EVAL_SERIES_MAX_ITERS_DEFAULT 20

/*!
 * Cauchy contour integral parameters
 */
#define EXP_FUNCTIONS_SPECIAL_EVAL_CAUCHY_N			32
#define EXP_FUNCTIONS_SPECIAL_EVAL_CAUCHY_RADIUS	2



namespace sweet {
namespace ExpIntegration {

/*!
 * This class implements various EXP functions.
 *
 * See Cox-Matthews paper for the provided EXP functions
 */

template <typename T = double>
class ExpFunction
{
public:
	Error::Base error;

	std::string functionName;

private:
	typedef std::complex<T> CT;

	enum fun_id_enum
	{
		INVALID = 0,

		PHI0 = 1,
		PHI1,
		PHI2,
		PHI3,
		PHI4,
		PHI5,
		PHI6,
		PHI7,
		PHI8,
		PHI9,
		PHI10,
		PHI11,
		PHI12,
		PHI13,
		PHI14,
		PHI15,
		PHI16,
		PHI17,
		PHI18,
		PHI19,
		PHI20,

		UPS1,
		UPS2,
		UPS3,

		PSI1,
		PSI2,
		PSI3
	};

public:
	fun_id_enum functionId;


public:
	T eps_phi;
	T eps_ups;
	T pi2;

	void setup_constvars()
	{
#if SWEET_QUADMATH
		if (typeid(T) == typeid(__float128))
		{
			eps_phi = 1e-10;
			eps_ups= 1e-10;

			static char *sp;

			pi2 = (__float128)2.0*strtoflt128("3.1415926535897932384626433832795029", &sp);
		}
		else
#endif
		if (typeid(T) == typeid(double))
		{
			eps_phi = 1e-10;
			eps_ups= 1e-10;
			pi2 = (double)2.0*(double)M_PI;
		}
		else
		{
			SWEETErrorFatal("Type not supported");
		}
	}



public:
	ExpFunction()	:
		functionId(INVALID)
	{
		setup_constvars();
	}


public:
	bool setup(
			const std::string &i_functionName
	)
	{
		// get phi functions
		if (i_functionName[1] == 'h')
		{
			functionId = static_cast<fun_id_enum>(std::stoi(i_functionName.substr(3)) + 1);			
		}

		else if (i_functionName == "ups1")
			functionId = UPS1;
		else if (i_functionName == "ups2")
			functionId = UPS2;
		else if (i_functionName == "ups3")
			functionId = UPS3;

		/*
		 * Semi-Lag phi functions (phi0 factored out) - see sl-rexi paper
		 * This is not the psi from some of Martin's presentation slides!
		 */
		else if (i_functionName == "psi1")
			functionId = PSI1;
		else if (i_functionName == "psi2")
			functionId = PSI2;
		else if (i_functionName == "psi3")
			functionId = PSI3;

		else
		{
			std::ostringstream ss;
			ss << "The function '" << i_functionName << "' is not supported!";
			error.set(ss.str());
			return false;
		}

		functionName = i_functionName;
		return true;
	}

	bool isSetup()
	{
		return functionId != INVALID;
	}

#if SWEET_QUADMATH
	/**************************************************************
	 * __float128 TYPES
	 **************************************************************/
	static std::complex<__float128> l_expcplx(const std::complex<__float128> &i_value)
	{
		__complex128 z;
		__real__ z = i_value.real();
		__imag__ z = i_value.imag();

		__complex128 val = cexpq(z);

		return std::complex<__float128>(crealq(val), cimagq(val));
	}

	inline
	static T l_sqrt(const __float128 &i_value)
	{
		return sqrtq(i_value);
	}

	static const std::complex<__float128> l_sqrtcplx(const std::complex<__float128> &i_value)
	{
		__complex128 z;
		__real__ z = i_value.real();
		__imag__ z = i_value.imag();

		__complex128 val = csqrtq(z);

		return std::complex<double>(crealq(val), cimagq(val));
	}
#endif


	/**************************************************************
	 * double TYPES
	 * Might suffer of numerical double precision limited effects
	 **************************************************************/

	std::complex<double> l_expcplx(const std::complex<double> &i_value)
	{
		return std::exp(i_value);
	};

	T l_sqrt(const double &i_value)
	{
		return std::sqrt(i_value);
	};

	std::complex<double> l_sqrtcplx(const std::complex<double> &i_value)
	{
		return std::sqrt(i_value);
	};


	int factorial(
			int N
	)
	{
		int retval = 1;

		for (; N > 0; N--)
			retval = retval * N;

		return retval;
	}



	/*
	 * \phi_N: Recursive computations
	 *
	 * ATTENTION: There's a singularity close to 0!
	 */
	std::complex<T> phiNRec(
		int n,
		const std::complex<T> &z
	)
	{
		if (n == 0)
			return l_expcplx(z);

		return (phiNRec(n-1, z) - (T)1.0/(T)factorial(n-1))/z;
	}




	/*
	 * \phi_N: Series-based computation
	 *
	 * This avoids the singularity close to 0
	 */
	std::complex<T> phiNSeries(
		int n,
		const std::complex<T> &K,
		int max_iters = EXP_FUNCTIONS_SPECIAL_EVAL_SERIES_MAX_ITERS_DEFAULT
	)
	{
		std::complex<T> pow_K = 1;
		T facn = factorial(n);
		std::complex<T> retval = pow_K/facn;

		for (int i = 1; i < max_iters; i++)
		{
			pow_K *= K;
			facn *= (n+i);

			retval += pow_K/facn;
		}

		return retval;
	}

	std::complex<T> phiNDirect(
		int i_phi_N,
		const std::complex<T> &i_K
	)
	{
		switch (i_phi_N)
		{
			case 0:
				return std::exp(i_K);
			case 1:
				return (std::exp(i_K) - 1.)/i_K;
			case 2:
				return (std::exp(i_K) - 1. - i_K)/(i_K*i_K);
			case 3:
				return (2.*std::exp(i_K) - 2. - 2. * i_K - i_K * i_K)/(2. * i_K * i_K * i_K);
			case 4:
				return (6.*std::exp(i_K) - 6. - 6. * i_K - 3. * i_K * i_K - i_K * i_K * i_K)/(6. * i_K * i_K * i_K * i_K);
			case 5:
				return (24.*std::exp(i_K) - 24. - 24.*i_K - 12. * i_K * i_K - 4.*i_K*i_K*i_K - i_K*i_K*i_K*i_K)/(24.*i_K*i_K*i_K*i_K*i_K);
			case 6:
				return (120.*std::exp(i_K) - 120. - 120.*i_K - 60.*i_K*i_K - 20.*i_K*i_K*i_K - 5.*i_K*i_K*i_K*i_K - i_K*i_K*i_K*i_K*i_K)/(120.*i_K*i_K*i_K*i_K*i_K*i_K);
			case 7:
				return (720.*std::exp(i_K) - 720. - 720.*i_K - 360.*i_K*i_K - 120.*i_K*i_K*i_K - 30.*i_K*i_K*i_K*i_K - 6.*i_K*i_K*i_K*i_K*i_K - i_K*i_K*i_K*i_K*i_K*i_K)/(720.*i_K*i_K*i_K*i_K*i_K*i_K*i_K);
			case 8:
				return (5040.*std::exp(i_K) - 5040. - 5040.*i_K - 2520.*i_K*i_K - 840.*i_K*i_K*i_K - 210.*i_K*i_K*i_K*i_K - 42.*i_K*i_K*i_K*i_K*i_K - 7.*i_K*i_K*i_K*i_K*i_K*i_K - i_K*i_K*i_K*i_K*i_K*i_K*i_K)/(5040.*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K);
			case 9:
				return (40320.*std::exp(i_K) - 40320. - 40320.*i_K - 20160.*i_K*i_K - 6720.*i_K*i_K*i_K - 1680.*i_K*i_K*i_K*i_K - 336.*i_K*i_K*i_K*i_K*i_K - 56.*i_K*i_K*i_K*i_K*i_K*i_K - 8.*i_K*i_K*i_K*i_K*i_K*i_K*i_K - i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K)/(40320.*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K);
			case 10:
				return (362880.*std::exp(i_K) - 362880. - 362880.0*i_K - 181440.0*i_K*i_K - 60480.0*i_K*i_K*i_K - 15120.0*i_K*i_K*i_K*i_K - 3024.0*i_K*i_K*i_K*i_K*i_K - 504.0*i_K*i_K*i_K*i_K*i_K*i_K - 72.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 9.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K)/(362880.*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K);
			case 11:
				return (3628800.*std::exp(i_K) - 3628800. - 3628800.0*i_K - 1814400.0*i_K*i_K - 604800.0*i_K*i_K*i_K - 151200.0*i_K*i_K*i_K*i_K - 30240.0*i_K*i_K*i_K*i_K*i_K - 5040.0*i_K*i_K*i_K*i_K*i_K*i_K - 720.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 90.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 10.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K)/(3628800.*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K);
			case 12:
				return (39916800.*std::exp(i_K) - 39916800. - 39916800.0*i_K - 19958400.0*i_K*i_K - 6652800.0*i_K*i_K*i_K - 1663200.0*i_K*i_K*i_K*i_K - 332640.0*i_K*i_K*i_K*i_K*i_K - 55440.0*i_K*i_K*i_K*i_K*i_K*i_K - 7920.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 990.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 110.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 11.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K)/(39916800.*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K);
			case 13 :
				return (479001600.0*std::exp(i_K) - 479001600.0 - 479001600.0*i_K - 239500800.0*i_K*i_K - 79833600.0*i_K*i_K*i_K - 19958400.0*i_K*i_K*i_K*i_K - 3991680.0*i_K*i_K*i_K*i_K*i_K - 665280.0*i_K*i_K*i_K*i_K*i_K*i_K - 95040.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 11880.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1320.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 132.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 12.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K)/(479001600.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K);
			case 14 :
				return (6227020800.0*std::exp(i_K) - 6227020800.0 - 6227020800.0*i_K - 3113510400.0*i_K*i_K - 1037836800.0*i_K*i_K*i_K - 259459200.0*i_K*i_K*i_K*i_K - 51891840.0*i_K*i_K*i_K*i_K*i_K - 8648640.0*i_K*i_K*i_K*i_K*i_K*i_K - 1235520.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 154440.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 17160.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1716.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 156.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 13.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K)/(6227020800.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K);
			case 15 :
				return (87178291200.0*std::exp(i_K) - 87178291200.0 - 87178291200.0*i_K - 43589145600.0*i_K*i_K - 14529715200.0*i_K*i_K*i_K - 3632428800.0*i_K*i_K*i_K*i_K - 726485760.0*i_K*i_K*i_K*i_K*i_K - 121080960.0*i_K*i_K*i_K*i_K*i_K*i_K - 17297280.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 2162160.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 240240.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 24024.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 2184.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 182.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 14.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K)/(87178291200.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K);
			case 16 :
				return (1307674368000.0*std::exp(i_K) - 1307674368000.0 - 1307674368000.0*i_K - 653837184000.0*i_K*i_K - 217945728000.0*i_K*i_K*i_K - 54486432000.0*i_K*i_K*i_K*i_K - 10897286400.0*i_K*i_K*i_K*i_K*i_K - 1816214400.0*i_K*i_K*i_K*i_K*i_K*i_K - 259459200.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 32432400.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 3603600.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 360360.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 32760.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 2730.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 210.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 15.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K)/(1307674368000.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K);
			case 17 :
				return (20922789888000.0*std::exp(i_K) - 20922789888000.0 - 20922789888000.0*i_K - 10461394944000.0*i_K*i_K - 3487131648000.0*i_K*i_K*i_K - 871782912000.0*i_K*i_K*i_K*i_K - 174356582400.0*i_K*i_K*i_K*i_K*i_K - 29059430400.0*i_K*i_K*i_K*i_K*i_K*i_K - 4151347200.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 518918400.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 57657600.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 5765760.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 524160.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 43680.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 3360.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 240.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 16.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K)/(20922789888000.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K);
			case 18 :
				return (355687428096000.0*std::exp(i_K) - 355687428096000.0 - 355687428096000.0*i_K - 177843714048000.0*i_K*i_K - 59281238016000.0*i_K*i_K*i_K - 14820309504000.0*i_K*i_K*i_K*i_K - 2964061900800.0*i_K*i_K*i_K*i_K*i_K - 494010316800.0*i_K*i_K*i_K*i_K*i_K*i_K - 70572902400.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 8821612800.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 980179200.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 98017920.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 8910720.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 742560.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 57120.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 4080.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 272.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 17.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K)/(355687428096000.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K);
			case 19 :
				return (6402373705728000.0*std::exp(i_K) - 6402373705728000.0 - 6402373705728000.0*i_K - 3201186852864000.0*i_K*i_K - 1067062284288000.0*i_K*i_K*i_K - 266765571072000.0*i_K*i_K*i_K*i_K - 53353114214400.0*i_K*i_K*i_K*i_K*i_K - 8892185702400.0*i_K*i_K*i_K*i_K*i_K*i_K - 1270312243200.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 158789030400.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 17643225600.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1764322560.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 160392960.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 13366080.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1028160.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 73440.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 4896.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 306.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 18.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K)/(6402373705728000.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K);
			case 20 :
				return (1.21645100408832e+17*std::exp(i_K) - 1.21645100408832e+17 - 1.21645100408832e+17*i_K - 6.0822550204416e+16*i_K*i_K - 2.0274183401472e+16*i_K*i_K*i_K - 5068545850368000.0*i_K*i_K*i_K*i_K - 1013709170073600.0*i_K*i_K*i_K*i_K*i_K - 168951528345600.0*i_K*i_K*i_K*i_K*i_K*i_K - 24135932620800.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 3016991577600.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 335221286400.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 33522128640.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 3047466240.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 253955520.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 19535040.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1395360.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 93024.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 5814.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 342.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 19.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K - 1.0*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K)/(1.21645100408832e+17*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K*i_K);
			default:
				SWEETErrorFatal("Phi_n explicit formula not implemented yet for such high order...");
				exit(-1);
		}
	}

	/*
	 * \phi_N: Buvoli et al Cauchy-based computation
	 */
	std::complex<T> phiNmyCauchy(
		int i_phi_N,
		const std::complex<T> &i_K
	)
	{
		if (i_phi_N == 0)
			return std::exp(i_K);

		int P = EXP_FUNCTIONS_SPECIAL_EVAL_CAUCHY_N;
		T R = EXP_FUNCTIONS_SPECIAL_EVAL_CAUCHY_RADIUS;

		std::complex<T> retval = 0;
		for (int k = 0; k < P; k++)
		{
			T theta_j = (pi2 * k / (T)P);
			retval += phiNDirect(i_phi_N, R * std::exp(std::complex<double>(0., theta_j)) + i_K);
		}

		return retval / (T)P;
	}

	std::complex<T> phiNCauchy(
		int i_phi_N,
		const std::complex<T> &i_K
	)
	{
		if (i_phi_N == 0)
			return std::exp(i_K);

		int N = EXP_FUNCTIONS_SPECIAL_EVAL_CAUCHY_N;
		T r = EXP_FUNCTIONS_SPECIAL_EVAL_CAUCHY_RADIUS;

		std::complex<T> retval = 0;

#if EXP_FUNCTIONS_USE_REC_PHI
		for (int j = 0; j < N; j++)
		{
			/*
			 * TODO:
			 * Figure out why the half-shifted version is a must-do
			 * for the upsN functions (forgot the reason for it)
			 */

			T theta_j = (T)(pi2*(j + 0.5)/N);
				
			// else
			// {
			// 	// allow points directly on axes
			// 	theta_j = (T)pi2*((T)j)/(T)N;
			// }

			// sampling position of support point
			std::complex<T> pos = r*std::exp(std::complex<double>(0, theta_j));

			// position
			std::complex<T> alpha = pos;
			std::complex<T> beta = -phiNRec(i_phi_N, alpha)*pos/(T)N;

			retval += beta/(i_K - alpha);
		}
#else
		for (int k = 0; k < N; k++)
		{
			T theta_j = (pi2 * k / (T)N);
			retval += phiNDirect(i_phi_N, r * std::exp(std::complex<double>(0., theta_j)) + i_K);
		}
		retval /= (T)N;
#endif
		return retval;
	}


	/*
	 * \phi_N: default caller switching between recursive formulation and series
	 */
	std::complex<T> phiN(
			int i_phi_N,
			const std::complex<T> &z
	)
	{
		if (std::abs(z) < EXP_FUNCTIONS_SPECIAL_THRESHOLD)
		{
#if EXP_FUNCTIONS_SPECIAL_EVAL_TYPE == 0
			return phiNSeries(i_phi_N, z);
#else
			return phiNmyCauchy(i_phi_N, z);
#endif
		}

		//return phiNRec(i_phi_N, z);
		return phiNDirect(i_phi_N, z);
	}



	std::complex<T> upsNDirect(
			int n,
			const std::complex<T> &z
	)
	{
		switch(n)
		{
		case 1:
			// http://www.wolframalpha.com/input/?i=(-4-K%2Bexp(K)*(4-3K%2BK*K))%2F(K*K*K)
			return (-4.0-z+l_expcplx(z)*(4.0-3.0*z+z*z)) / (z*z*z);

		case 2:
			// http://www.wolframalpha.com/input/?i=(2%2BK%2Bexp(K)*(-2%2BK))%2F(K*K*K)
			return (2.0+z+l_expcplx(z)*(-2.0+z)) / (z*z*z);

		case 3:
			// http://www.wolframalpha.com/input/?i=(-4-3*K-K*K%2Bexp(K)*(4-K))%2F(K*K*K)
			return (-4.0-3.0*z-z*z+l_expcplx(z)*(4.0-z)) / (z*z*z);
		}

		SWEETErrorFatal("ups number not supported!");
		return 0;
	}

	std::complex<T> ups1Series(
			const std::complex<T> &K,
			int max_iters = EXP_FUNCTIONS_SPECIAL_EVAL_SERIES_MAX_ITERS_DEFAULT
	)
	{
		std::complex<T> pow_K = 1;
		T fac_denom = 6;	//factorial(3);
		std::complex<T> retval = pow_K/fac_denom;

		for (int l = 1; l < max_iters; l++)
		{
			pow_K *= K;
			fac_denom *= (l+3);
			retval += pow_K*(l+1.0)*(l+1.0)/fac_denom;
		}
		return retval;
	}

	std::complex<T> ups2Series(
			const std::complex<T> &K,
			int max_iters = EXP_FUNCTIONS_SPECIAL_EVAL_SERIES_MAX_ITERS_DEFAULT
	)
	{
		std::complex<T> pow_K = 1;
		T facn = 6;//factorial(3);
		std::complex<T> retval = 1.0/2.0;

		retval += (K-2.0)*pow_K/facn;

		for (int l = 1; l < max_iters; l++)
		{
			pow_K *= K;
			facn *= (l+3);
			retval += (K-2.0)*pow_K/facn;
		}
		return retval;
	}

	std::complex<T> ups3Series(
			const std::complex<T> &K,
			int max_iters = EXP_FUNCTIONS_SPECIAL_EVAL_SERIES_MAX_ITERS_DEFAULT
	)
	{
		std::complex<T> pow_K = 1;
		T facn = 6;//factorial(3);
		std::complex<T> retval = -1.0/2.0;

		retval += (4.0-K)*pow_K/facn;

		for (int l = 1; l < max_iters; l++)
		{
			pow_K *= K;
			facn *= (l+3);
			retval += (4.0-K)*pow_K/facn;
		}
		return retval;
	}

	std::complex<T> upsNSeries(
			int n,
			const std::complex<T> &z,
			int max_iters = EXP_FUNCTIONS_SPECIAL_EVAL_SERIES_MAX_ITERS_DEFAULT
	)
	{
		switch(n)
		{
		case 1:
			return ups1Series(z, max_iters);

		case 2:
			return ups2Series(z, max_iters);

		case 3:
			return ups3Series(z, max_iters);
		}

		SWEETErrorFatal("This Upsilon function is not supported!");
		return 0;
	}


	std::complex<T> upsNCauchy(
			int i_ups_N,
			const std::complex<T> &i_K,
			int max_iters = EXP_FUNCTIONS_SPECIAL_EVAL_SERIES_MAX_ITERS_DEFAULT
	)
	{
		int N = EXP_FUNCTIONS_SPECIAL_EVAL_CAUCHY_N;
		T r = EXP_FUNCTIONS_SPECIAL_EVAL_CAUCHY_RADIUS;

		std::complex<T> retval = 0;
		for (int j = 0; j < N; j++)
		{
			/*
			 * TODO:
			 * Figure out why the half-shifted version is a must-do
			 * for the upsN functions (forgot the reason for it)
			 */
			T theta_j;
			if (1)
			{
				// avoid points directly on axes for upsN functions
				theta_j = (T)pi2*((T)j+(T)0.5)/(T)N;
			}
			else
			{
				// allow points directly on axes
				theta_j = (T)pi2*((T)j)/(T)N;
			}

			// sampling position of support point
			std::complex<T> pos = r*std::exp(std::complex<double>(0, theta_j));

			// position
			std::complex<T> alpha = pos;
			std::complex<T> beta = -upsNDirect(i_ups_N, alpha)*pos/(T)N;

			retval += beta/(i_K - alpha);
		}

		return retval;
	}


	/*
	 * \phi_N: default caller switching between recursive formulation and series
	 */
	std::complex<T> upsN(
			int N,
			const std::complex<T> &z
	)
	{
		T linf = z.real()*z.real() + z.imag()*z.imag();
		if (linf < EXP_FUNCTIONS_SPECIAL_THRESHOLD)
		{
#if EXP_FUNCTIONS_SPECIAL_EVAL_TYPE == 0
			return upsNSeries(N, z);
#else
			return upsNCauchy(N, z);
#endif
		}

		return upsNDirect(N, z);
	}




	/*
	 * Semi-Lagrangian psi functions
	 */
	std::complex<T> psiN(
			int n,
			const std::complex<T> &i_K
	)
	{
		switch(n)
		{
		case 1:	// SL psi1
			// psi1(z) = phi1(-z)
			return phiN(1, -i_K);
			break;


		case 2:	// SL psi2
			// psi2(z)=-phi2(-z)+phi1(-z)
			return -phiN(2, -i_K) + phiN(1, -i_K);
			break;


		case 3:	// SL psi3
#if 1
			SWEETErrorFatal("TODO: Redo this with e.g. series treatment");
#else
			if (lamdt < eps_phi)
				//					if (lamdt*lamdt*lamdt < expFunctions.eps_phi)
			{
				return 1.0/(2.0*3.0);
			}
			else
			{
				return (l_expcplx(i_K) - std::complex<T>(1.0) - i_K - i_K*i_K)/(i_K*i_K*i_K)
						- (l_expcplx(i_K) - std::complex<T>(1.0) - i_K)/(i_K*i_K)
						+ std::complex<T>(0.5)*(l_expcplx(i_K) - std::complex<T>(1.0))/i_K;
			}
#endif
			break;


		default:
			SWEETErrorFatal("This psi number is not yet supported");
		}

		return -1;
	}


	/*
	 * Workaround for error in gcc version 8 and 7.
	 * ... defined in discarded section ...
	 */
#ifndef __INTEL_COMPILER
#if __GNUC__ != 8 && __GNUC__ != 7
	SWEET_OMP_DECLARE_SIMD
#endif
#endif
	void eval(
		const std::complex<T> &i_K,
		std::complex<T> &o_retval
	)
	{
		// shortcut for phies
		if (functionId < 21 && functionId > 0) {
			o_retval = phiN(functionId - 1, i_K);
			return;
		}

		switch(functionId)
		{

		case UPS1:
			o_retval = upsN(1, i_K);
			return;

		case UPS2:
			o_retval = upsN(2, i_K);
			return;

		case UPS3:
			o_retval = upsN(3, i_K);
			return;



		case PSI1:
			o_retval = psiN(1, i_K);
			return;

		case PSI2:
			o_retval = psiN(2, i_K);
			return;

		case PSI3:
			o_retval = psiN(3, i_K);
			return;

		default:
			;
		}

		if (functionId == INVALID)
			SWEETErrorFatal("Exponential function likely not initialized (call setup())!");

		std::ostringstream ss;
		ss << "The phi function with id '" << functionId << "' is not supported!";
		error.set(ss.str());
		SWEETErrorFatal(ss.str().c_str());
	}

};

}}

#endif
