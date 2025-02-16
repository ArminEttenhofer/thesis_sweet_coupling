#ifndef PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_SHACK_HPP
#define PROGRAMS_ODE_GENERIC_DE_DAHLQUIST_SHACK_HPP

#include <complex>
#include <sweet/Tools/ProgramArguments.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/Base.hpp>


namespace ODE_Generic {
namespace DE_Dahlquist {

class Shack	:
		public sweet::Shacks::Base
{
	typedef std::complex<double> T;
public:

	/**
	 * Coefficients for
	 *
	 * d/dt u(t) = \lambda_1 * u(t) + \lambda_2 * u(t) + \lambda_3 * u(t) + mu * exp(phi*t)
	 */
	T lambda1;
	T lambda2;
	T lambda3;

	T mu;
	T phi;

	Shack()
	{
		lambda1.real(0.0);
		lambda1.imag(10.0);

		lambda2.real(0.0);
		lambda2.imag(0.0);

		lambda3.real(0.0);
		lambda3.imag(0.0);

		mu.real(0.0);
		mu.imag(0.0);

		phi.real(0.0);
		phi.imag(0.0);
	}



	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << "ODEGeneric_DE_Dahlquist parameters:" << std::endl;
		std::cout << i_prefix << "	--ode-dahlquist-lambda1 [float],[float]	set value of lambda1, 2nd optional float is imaginary part" << std::endl;
		std::cout << i_prefix << "	--ode-dahlquist-lambda2 [float],[float]	set value of lambda2, 2nd optional float is imaginary part" << std::endl;
		std::cout << i_prefix << "	--ode-dahlquist-lambda3 [float],[float]	set value of lambda3, 2nd optional float is imaginary part" << std::endl;
		std::cout << i_prefix << "	--ode-dahlquist-mu [float],[float]	set value of mu, 2nd optional float is imaginary part" << std::endl;
		std::cout << i_prefix << "	--ode-dahlquist-phi [float],[float]	set value of phi, 2nd optional float is imaginary part" << std::endl;
	}


	void printShack(const std::string& i_prefix = "")
	{
		std::cout << i_prefix << "ODEGeneric_DE_Dahlquist parameters:" << std::endl;
		std::cout << i_prefix << " + lambda1: " << lambda1 << std::endl;
		std::cout << i_prefix << " + lambda2: " << lambda2 << std::endl;
		std::cout << i_prefix << " + lambda3: " << lambda3 << std::endl;
		std::cout << i_prefix << " + mu: " << mu << std::endl;
		std::cout << i_prefix << " + phi: " << phi << std::endl;
	}



	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--ode-dahlquist-lambda1", lambda1);
		i_pa.getArgumentValueByKey("--ode-dahlquist-lambda2", lambda2);
		i_pa.getArgumentValueByKey("--ode-dahlquist-lambda3", lambda3);
		i_pa.getArgumentValueByKey("--ode-dahlquist-mu", mu);
		i_pa.getArgumentValueByKey("--ode-dahlquist-phi", phi);
		
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_pa);
		return true;
	}
};


}}

#endif
