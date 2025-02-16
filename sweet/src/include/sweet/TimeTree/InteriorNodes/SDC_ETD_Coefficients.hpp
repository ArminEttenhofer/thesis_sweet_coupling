#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_SDC_ETD_COEFFICIENTS_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_SDC_ETD_COEFFICIENTS_HPP

#include <vector>
#include <string>
#include <sweet/Dict/TypesArrayND.hpp>
#include <sweet/SDC/Shack.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>

namespace sweet {
namespace TimeTree {
namespace InteriorNodes {

/*!
 * Special class for ETDSDC coefficients (see SDC_Coefficients for a generic SDC case)
 */
class SDC_ETD_Coefficients
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
	// Always used with lastNode -> no quadrature matrices

	/**
	 * SDC parameters
	 */
	//! Number of SDC iterations (excluding initial sweep)
	int nIters;
	//! Aijl matrix of finite differences
	Arr3D A;

	//! Enumeration for pre sweep types
	enum EnumPreSweep {
		PRE_SWEEP_INVALID 	= 0,
		PRE_SWEEP_ZEROS		= 1
	};
	//! Type of pre sweep
	EnumPreSweep preSweep;
	//! Enumeration for post sweep types
	enum EnumPostSweep {
		POST_SWEEP_INVALID 			= 0,
		POST_SWEEP_LASTNODE 		= 1
	};
	//! Type of post sweep
	EnumPostSweep postSweep;

	/**
	 * SDC description variables
	 */
	//! Unique ID String
	std::string idString;

public:
	SDC_ETD_Coefficients()	:
		// Utility attributes
		fileName(""),
		shackSDC(nullptr),
		// Collocation parameters
		tauNodes(),
		tau(),
		M(-1),
		deltaTau(),
		// SDC parameters
		nIters(-1),
		preSweep(PRE_SWEEP_ZEROS),
		postSweep(POST_SWEEP_LASTNODE),
		// SDC description variables
		idString("UNDEF")
	{}

	~SDC_ETD_Coefficients()
	{}

	/*!
	 * Copy constructor (using default copying each attributes)
	 */
	SDC_ETD_Coefficients(const SDC_ETD_Coefficients &i_src) = default;

	bool loadSDCCoefficientsFromFile(
			const std::string &i_fileName = ""
	)
	{
		fileName = i_fileName;

		sweet::Dict::Dict params(fileName);
		ERROR_CHECK_COND_RETURN_BOOLEAN(params);

		// Check if file is SDC_ETD_Coeffs, not SDC_Coeffs
		// would throw an error "key A not found"
		params.get("A", A);

		// Collocation parameters
		params.get("tauNodes", tauNodes);
		params.get("tau", tau);
		tau.setOffset(1);
		M = tauNodes.size();
		params.get("deltaTau", deltaTau);

		// SDC parameters
		params.get("nIter", nIters);

		// Pre and post sweep always stubbed
		preSweep = PRE_SWEEP_ZEROS;
		postSweep = POST_SWEEP_LASTNODE;	
		
		// SDC description variables
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
		std::cout << i_Prefix << "tau: " << tau << std::endl;
		std::cout << i_Prefix << "deltaTau: " << deltaTau << std::endl;
		std::cout << i_Prefix << "nIter: " << nIters << std::endl;
		std::cout << i_Prefix << "preSweep: " << preSweep << std::endl;
		std::cout << i_Prefix << "postSweep: " << postSweep << std::endl;
		std::cout << std::endl;
	}
};

}}}

#endif
