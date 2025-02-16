#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_SDC_BASE_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_SDC_BASE_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/Dict/TypesArrayND.hpp>
#include <sweet/SDC/Shack.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>
#include <vector>
#include <string>


namespace sweet {
namespace TimeTree {
namespace InteriorNodes {


/*!
 * A base class for all SDC methods to provide the same common ground
 */
class SDC_Coefficients
{
public:
	/*
	 * Utility attributes
 	 */

	//! Convenience typedef for vectors
	typedef sweet::Dict::TypesArrayND<1, double> Vec;
	//! Convenience typedef for matrix
	typedef sweet::Dict::TypesArrayND<2, double> Mat;
	//! Convenience typedef for 3D array
	typedef sweet::Dict::TypesArrayND<3, double> Arr3D;

	//! Filename where coefficients have been loaded from
	std::string fileName;
	//! Parameter shack 
	sweet::SDC::Shack* shackSDC;


	/**
	 * Collocation parameters
	 */

	//! Subtime stamps \in [0;1]
	Vec tauNodes;
	//! Subtime stamps \in [0;1], including 0 and 1
	Vec tau;
	//! Number of quadrature nodes
	int M;
	//! Distance between nodes including relative distances to 0 at the beginning and 1 at the end
	Vec deltaTau;
	//! Weights if final end update with quadrature is used (if useEndUpdate)
	Vec weights;
	//! Quadrature matrix (zero-to-node formulation)
	Mat qMat;
	//! Quadrature matrix (node-to-node formulation)
	Mat sMat;
	//! Interpolation coefficients (end-point update)
	Vec hCoeffs;


	/**
	 * SDC parameters
	 */

	//! Number of SDC iterations (excluding initial sweep)
	int nIters;
	//! Implicit QDelta matrices
	Arr3D qMatDeltaI;
	//! Explicit QDelta matrix
	Mat qMatDeltaE;
	//! dtau coefficients for explicit QDelta
	Vec dtauE;
	//! QDelta for initial sweep (implicit term)
	Mat qMatDelta0;
	//! dtau coefficients for initial (implicit) QDelta
	Vec dtau0;
	//! Enumeration for pre sweep types
	enum EnumPreSweep {
		PRE_SWEEP_INVALID 	= 0,
		PRE_SWEEP_QDELTA 	= 1,
		PRE_SWEEP_COPY 		= 2,
		PRE_SWEEP_ZEROS 	= 3
	};
	//! Type of pre sweep
	EnumPreSweep preSweep;
	//! Enumeration for post sweep types
	enum EnumPostSweep {
		POST_SWEEP_INVALID 			= 0,
		POST_SWEEP_LASTNODE 		= 1,
		POST_SWEEP_QUADRATURE 		= 2,
		POST_SWEEP_INTERPOLATION 	= 3
	};
	//! Type of post sweep
	EnumPostSweep postSweep;


	/**
	 * SDC description variables
	 */

	//! Whether or not qDelta matrices are parallel
	bool diagonal;
	//! Theoretical order of the implicit term
	int orderI;
	//! Theoretical order of the explicit term
	int orderE;
	//! Wether the left node is 0
	bool leftIsNode;
	//! Wether the right node is 0
	bool rightIsNode;
	//! Whether or not activate parallel computation for diagonal sweeps
	bool parallel;
	//! Unique ID String
	std::string idString;


public:
	SDC_Coefficients()	:
		// Utility attributes
		fileName(""),
		shackSDC(nullptr),
		// Collocation parameters
		tauNodes(),
		tau(),
		M(-1),
		deltaTau(),
		weights(),
		qMat(),
		sMat(),
		hCoeffs(),
	
		// SDC parameters
		nIters(-1),
		qMatDeltaI(),
		qMatDeltaE(),
		dtauE(),
		qMatDelta0(),
		dtau0(),
		preSweep(PRE_SWEEP_INVALID),
		postSweep(POST_SWEEP_INVALID),
		// SDC description variables
		diagonal(false),
		orderI(-1),
		orderE(-1),
		leftIsNode(false),
		rightIsNode(false),
		parallel(false),
		idString("UNDEF")
	{}

	~SDC_Coefficients()
	{}

	/*!
	 * Copy constructor (using default copying each attributes)
	 */
	SDC_Coefficients(const SDC_Coefficients &i_src) = default;

	bool loadSDCCoefficientsFromFile(
			const std::string &i_fileName = ""
	)
	{
		fileName = i_fileName;

		sweet::Dict::Dict params(fileName);
		ERROR_CHECK_COND_RETURN_BOOLEAN(params);

		// Collocation parameters
		params.get("tauNodes", tauNodes);
		params.get("tau", tau);
		tau.setOffset(1);
		M = tauNodes.size();
		params.get("deltaTau", deltaTau);
		params.get("weights", weights);
		params.get("qMatrix", qMat);
		params.get("sMatrix", sMat);
		params.get("hCoeffs", hCoeffs);
		

		// SDC parameters
		params.get("nIter", nIters);
		params.get("qDeltaI", qMatDeltaI);
		params.get("qDeltaE", qMatDeltaE);
		params.get("dtauE", dtauE);
		params.get("qDelta0", qMatDelta0);
		params.get("dtau0", dtau0);
		
		// Pre sweep
		std::string preSweepString;
		params.get("preSweep", preSweepString);
		if (preSweepString == "QDELTA")
		{
			preSweep = PRE_SWEEP_QDELTA;	
		}
		else if (preSweepString == "COPY")
		{
			preSweep = PRE_SWEEP_COPY;
		}
		else if (preSweepString == "ZEROS")
		{
			preSweep = PRE_SWEEP_ZEROS;
		}
		else
		{
			preSweep = PRE_SWEEP_INVALID;
			SWEETErrorFatal("Invalid pre sweep type '"+preSweepString+"'");
		}

		// Post sweep
		std::string postSweepString;
		params.get("postSweep", postSweepString);
		if (postSweepString == "LASTNODE")
		{
			postSweep = POST_SWEEP_LASTNODE;	
		}
		else if (postSweepString == "QUADRATURE")
		{
			postSweep = POST_SWEEP_QUADRATURE;
		}
		else if (postSweepString == "INTERPOLATION")
		{
			postSweep = POST_SWEEP_INTERPOLATION;
		}
		else
		{
			postSweep = POST_SWEEP_INVALID;
			SWEETErrorFatal("Invalid post sweep type '"+postSweepString+"'");
		}
		
		// SDC description variables
		params.get("diagonal", diagonal);
		params.get("orderI", orderI);
		params.get("orderE", orderE);
		leftIsNode = (tauNodes(0) == 0.0);
		rightIsNode = (tauNodes(M-1) == 1.0);
		params.get("idString", idString);

		ERROR_CHECK_COND_RETURN_BOOLEAN(params);

		return true;
	}


public:
	void printSDCInformation(const std::string &i_Prefix = "")
	{
		std::cout << i_Prefix << "fileName: " << fileName << std::endl;
		std::cout << i_Prefix << "idString: " << idString << std::endl;
		std::cout << i_Prefix << "nNodes: " << M << std::endl;
		std::cout << i_Prefix << "weights: " << weights << std::endl;
		std::cout << i_Prefix << "tau: " << tau << std::endl;
		std::cout << i_Prefix << "deltaTau: " << deltaTau << std::endl;
		std::cout << i_Prefix << "qMat: " << qMat << std::endl;
		std::cout << i_Prefix << "sMat: " << sMat << std::endl;
		std::cout << i_Prefix << "hCoeffs: " << hCoeffs << std::endl;
		std::cout << i_Prefix << "nIter: " << nIters << std::endl;
		std::cout << i_Prefix << "qMatDeltaI: " << qMatDeltaI << std::endl;
		std::cout << i_Prefix << "qMatDeltaE: " << qMatDeltaE << std::endl;
		std::cout << i_Prefix << "dtauE: " << dtauE << std::endl;
		std::cout << i_Prefix << "qMatDelta0: " << qMatDelta0 << std::endl;
		std::cout << i_Prefix << "dtau0: " << dtau0 << std::endl;
		std::cout << i_Prefix << "preSweep: " << preSweep << std::endl;
		std::cout << i_Prefix << "postSweep: " << postSweep << std::endl;
		std::cout << i_Prefix << "diagonal: " << diagonal << std::endl;
		std::cout << i_Prefix << "orderI: " << orderI << std::endl;
		std::cout << i_Prefix << "orderE: " << orderE << std::endl;
		std::cout << i_Prefix << "parallel: " << parallel << std::endl;
		std::cout << std::endl;
	}
};

}}}

#endif
