/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef INCLUDE_SWEET_EXPINTEGRATION_REXI_REXICOEFFICIENTSSET_HPP
#define INCLUDE_SWEET_EXPINTEGRATION_REXI_REXICOEFFICIENTSSET_HPP

#include <sweet/ExpIntegration/REXI/REXICoefficients.hpp>
#include <vector>
#include <complex>
#include <sweet/Tools/StringSplit.hpp>

namespace sweet {
namespace ExpIntegration {
namespace REXI {

template <typename T = double>
class REXICoefficientsSet
{
public:
	typedef std::complex<T> TComplex;

	std::vector< REXICoefficients<> > rexiCoefficientVector;


	/**
	 * Load REXI coefficients from filenames
	 */
	void setupFromFiles(
			const std::string &i_rexi_filenames
	)
	{
		rexiCoefficientVector.clear();

		// do simply nothing if there's no file name
		if (i_rexi_filenames == "")
			return;

		std::vector<std::string> rexi_filenames = sweet::Tools::StringSplit::split(i_rexi_filenames, ",");

		for (auto iter = rexi_filenames.begin(); iter != rexi_filenames.end(); iter++)
		{
			std::string& rexi_filename = *iter;

			std::vector<std::string> split2 = sweet::Tools::StringSplit::split(rexi_filename, ":");

			if (split2.size() == 0)
			{
				SWEETErrorFatal("Strange things....");
			}
			else if (split2.size() == 1)
			{
				REXICoefficients<T> rexiCoefficients;
				rexiCoefficients.load_from_file(split2[0]);
			}
			else if (split2.size() == 2)
			{
				REXICoefficients<T> rexiCoefficients;
				rexiCoefficients.load_from_file(split2[1]);

				if (rexiCoefficients.filename != "")
					if (rexiCoefficients.function_name != split2[0])
						SWEETErrorFatal("Function name mismatch!");

				rexiCoefficientVector.push_back(rexiCoefficients);
			}
			else
			{
				SWEETErrorFatal("Too many split variable names");
			}
		}
	}


	const REXICoefficients<> *getByFunctionName(
			const std::string &i_function_name
	)	const
	{
		for (auto iter = rexiCoefficientVector.begin(); iter != rexiCoefficientVector.end(); iter++)
		{
			if (iter->function_name == i_function_name)
				return &(*iter);
		}

		SWEETErrorFatal("Not found");
		return nullptr;
	}
};

}}}

#endif
