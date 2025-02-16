/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_EXPINTEGRATION_REXI_REXICOEFFICIENTS_HPP
#define INCLUDE_SWEET_EXPINTEGRATION_REXI_REXICOEFFICIENTS_HPP

#include <sweet/Error/Base.hpp>
#include <sweet/Dict/Dict.hpp>
#include <sweet/LibMath/DQStuff.hpp>
#include <vector>
#include <complex>
#include <fstream>

namespace sweet {
namespace ExpIntegration {
namespace REXI {

template <typename T = double>
class REXICoefficients
{
public:
	typedef std::complex<T> TComplex;

public:
	Error::Base error;

	std::vector<TComplex> alphas;
	std::vector<TComplex> betas;
	TComplex gamma;

	std::string filename;
	std::string function_name;


	/*1
	 * Constructor
	 */
	REXICoefficients()	:
		gamma(0)
	{
	}

	/*!
	 * Load REXI coefficients from given file
	 */
	bool load_from_file(
			const std::string &i_filename
	)
	{
		Dict::Dict dict;
		dict.fileLoad(i_filename);

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dict);

		filename = i_filename;
		dict.get("function_name", function_name);
		dict.get("gamma", gamma);

		Dict::TypesArrayND<1, std::complex<double>> alphas_;
		dict.get("alphas", alphas_);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dict);
		alphas = alphas_.dataVector();

		Dict::TypesArrayND<1, std::complex<double>> betas_;
		dict.get("betas", betas_);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dict);
		betas = betas_.dataVector();

		return true;
	}
};

}}}

#endif
