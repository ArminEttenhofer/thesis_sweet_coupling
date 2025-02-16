/*
 * ShackSDC.hpp
 */

#ifndef INCLUDE_SWEET_SDC_SHACK_HPP
#define INCLUDE_SWEET_SDC_SHACK_HPP

#include <string>
#include <iostream>
#include <sweet/Dict/Dict.hpp>
#include <sweet/Shacks/Base.hpp>
#include <sweet/Tools/ProgramArguments.hpp>

namespace sweet {
namespace SDC {

/*!
 * \brief Shack for SDC
 */
class Shack	:
		public sweet::Shacks::Base
{
public:
	std::string fileName;

#if 0
	/// Nodes values (between 0 and 1)
	sweet::Dict::TypesArrayND<1, double> tauNodes;

	/// Quadrature weights
	sweet::Dict::TypesArrayND<1, double> weights;

	/// Collocation matrix
	sweet::Dict::TypesArrayND<2, double> qMat;

	/// QDelta matrix for implicit sweep
	sweet::Dict::TypesArrayND<2, double> qMatDeltaI;

	/// QDelta matrix for explicit sweep
	sweet::Dict::TypesArrayND<2, double> qMatDeltaE;

	/// QDelta matrix for initial (implicit) sweep
	sweet::Dict::TypesArrayND<2, double> qMatDelta0;

	/// Number of iterations (sweeps)
	sweet::Dict::Dict::int64 nIter;

	/// Type of initial sweep to use
	std::string preSweep;

	/// Whether or not use the diagonal implementation
	sweet::Dict::Dict::int64 diagonal;

	/// Whether or not use collocation update for end point
	sweet::Dict::Dict::int64 useEndUpdate;

	/// Unique string ID
	std::string idString;
#endif

	/// Whether or not activate time parallelization
	bool runParallel_DiagonalQDelta0Matrix;

	/// Default parameters for SDC shack
	Shack(){
#if 0
		int nNodes = 3;
		nIter = 3;
		diagonal = false;
		preSweep = "COPY";
		useEndUpdate = false;
		idString = "M3_RADAU-RIGHT_K3_BE_FE_COPY";

		// RADAU-RIGHT nodes, weights quadrature matrix
		tauNodes.resize(nNodes);
		const double _nodes[] = {
			0.15505102572168, 0.64494897427832, 1.
		};
		tauNodes = _nodes;

		weights.resize(nNodes);
		const double _weights[] = {
			0.3764030627004656, 0.51248582618842650, 0.1111111111111111
		};
		weights = _weights;

		qMat.resize(nNodes, nNodes);
		const double _qMatrix[] = {
			0.1968154772236567, -0.06553542585019642,  0.02377097434821968,
			0.394424314739085,   0.2920734116652353,  -0.04154875212600038,
			0.3764030627004656,  0.5124858261884265,   0.1111111111111079
		};
		qMat = _qMatrix;

		// BE for implicit sweep
		qMatDeltaI.resize(nNodes, nNodes);
		const double _qDeltaI[] = {
			0.15505102572168, 0.,         0.,
 			0.15505102572168, 0.48989794855664, 0.,
			0.15505102572168, 0.48989794855664, 0.35505102572168        
		};
		qMatDeltaI = _qDeltaI;

		// FE for explicit sweep
		qMatDeltaE.resize(nNodes, nNodes);
		const double _qDeltaE[] = {
			0.,         0.,         0.,
 			0.48989794855664, 0.,         0.,
			0.48989794855664, 0.35505102572168, 0.        
		};
		qMatDeltaE = _qDeltaE;

		// BEpar for initial sweep
		qMatDelta0.resize(nNodes, nNodes);
		const double _qDelta0[] = {
			0.15505102572168, 0.			  , 0.,
 			0.				, 0.64494897427832, 0.,
			0.				, 0.			  , 1.        
		};
		qMatDelta0 = _qDelta0;
#endif
		runParallel_DiagonalQDelta0Matrix = false;
	}

	void printProgramArguments(const std::string& i_prefix = "") override
	{
		std::cout << std::endl;
		std::cout << "SDC option:" << std::endl;
		std::cout << "	--sdc-file 	   [path]  SDC parameters in Dict format" << std::endl;
		std::cout << "	--sdc-parallel [bool]  wether or not activate time parallelization" << std::endl;
		std::cout << std::endl;
	}

#if 0
	bool loadFromFile(const std::string &i_fileName)
	{
		sweet::Dict::Dict params(i_fileName);
		params.get("nodes", tauNodes);
		params.get("weights", weights);
		params.get("qMatrix", qMat);
		params.get("qDeltaI", qMatDeltaI);
		params.get("qDeltaE", qMatDeltaE);
		params.get("qDelta0", qMatDelta0);
		params.get("nIter", nIter);
		params.get("diagonal", diagonal);
		params.get("preSweep", preSweep);
		params.get("useEndUpdate", useEndUpdate);
		params.get("idString", idString);

		ERROR_CHECK_COND_RETURN_BOOLEAN(params);

		return true;
	}
#endif

	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa) override
	{
		i_pa.getArgumentValueByKey("--sdc-parallel", runParallel_DiagonalQDelta0Matrix);

		i_pa.getArgumentValueByKey("--sdc-file", fileName);
#if 0
		if (i_pa.getArgumentValueByKey("--sdc-file", fileName))
		{
			loadFromFile(fileName);
		}
#endif

		return error.forwardWithPositiveReturn(i_pa.error);
	}

	void printShack(
		const std::string& i_prefix = ""
	) override
	{
		std::cout << std::endl;
		std::cout << "SDC:" << std::endl;
#if 0
		std::cout << " + M (number of nodes): " << tauNodes.size() << std::endl;
		std::cout << " + nodes: " << tauNodes << std::endl;
		std::cout << " + weights: " << weights << std::endl;
		std::cout << " + qMatrix: " << qMat << std::endl;
		std::cout << " + qDeltaI: " << qMatDeltaI << std::endl;
		std::cout << " + qDeltaE: " << qMatDeltaE << std::endl;
		std::cout << " + qDelta0: " << qMatDelta0 << std::endl;
		std::cout << " + nIter: " << nIter << std::endl;
		std::cout << " + diagonal: " << diagonal << std::endl;
		std::cout << " + preSweep: " << preSweep << std::endl;
		std::cout << " + useEndUpdate: " << useEndUpdate << std::endl;
#endif
		std::cout << " + parallel: " << runParallel_DiagonalQDelta0Matrix << std::endl;
		std::cout << std::endl;
	}
};


}}


#endif
