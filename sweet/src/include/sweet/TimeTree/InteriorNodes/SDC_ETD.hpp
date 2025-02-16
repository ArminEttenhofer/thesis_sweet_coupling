#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_SDC_ETD_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_SDC_ETD_HPP

#include <vector>
#include <string>
#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>
#include <sweet/TimeTree/InteriorNodes/SDC_ETD_Coefficients.hpp>

#include <sweet/ExpIntegration/ExpFunction.hpp>

namespace sweet {
namespace TimeTree {
namespace InteriorNodes {

/*!
 * Convinience class for linear DE term L(*)
 */
class LinearTerm {
	public:
		bool active;

		///! temporary node which can be used during setup
		std::shared_ptr<TimeTree_Node_Base> tmp_TimeTreeNode;

		// time tree indices
		int evalStartIdx_PreSweep_phi0;
		int evalStartIdx_PreSweep_phi1;
		int evalStartIdx_MainSweep_phi0;
		int evalStartIdx_MainSweep_phi1;
		int evalStartIdx_MainSweep_W;

		LinearTerm()	:
			active(false),
			evalStartIdx_PreSweep_phi0(-1),
			evalStartIdx_PreSweep_phi1(-1),
			evalStartIdx_MainSweep_phi0(-1),
			evalStartIdx_MainSweep_phi1(-1),
			evalStartIdx_MainSweep_W(-1)
		{}
	};

/*!
 * Convinience class for nonlinear DE term N(*)
 */
class NonlinearTerm {
	public:
		bool active;

		///! temporary node which can be used during setup
		std::shared_ptr<TimeTree_Node_Base> tmp_TimeTreeNode;

		// time tree indices
		int evalStartIdx_PreSweep;
		int evalStartIdx_MainSweep;

		//! Storage for nonlinear tendencies
		sweet::Data::GenericContainer::Base** tendencies_k0;
		sweet::Data::GenericContainer::Base** tendencies_k1;

		//! Storage for a correction sum Wi;i+1
		sweet::Data::GenericContainer::Base* accW;

		NonlinearTerm()	:
			active(false),
			evalStartIdx_PreSweep(-1),
			evalStartIdx_MainSweep(-1),
			tendencies_k0(nullptr),
			tendencies_k1(nullptr),
			accW(nullptr)
		{}
	};

/*!
  \brief Exponential Time Differencing SDC (From Buvoli 15' paper)
 */
class SDC_ETD	:
	public TimeTree_Node_InteriorHelper<SDC_ETD>,
	public SDC_ETD_Coefficients
{
private:
	//! Number of phi functions to use (M + 1)
	int _nPhi;

	//! Filename of SWEETDict with SDC coefficients
	std::string fileNameWithSDCCoefficients;

	// DE_Terms
	LinearTerm linearTerm;
	NonlinearTerm nonlinearTerm;

	//! Intermediate states during SDC
	sweet::Data::GenericContainer::Base* tmpU;
	sweet::Data::GenericContainer::Base* tmp;
	sweet::Data::GenericContainer::Base* tmp2;

	//! Provisional solution storage and intermediate states
	sweet::Data::GenericContainer::Base** provisional_U;

public:
	SDC_ETD()	:
		_nPhi(-1),
		tmpU(nullptr),
		tmp(nullptr),
		tmp2(nullptr)
	{
		setEvalAvailable(EVAL_INTEGRATION);
	}

	SDC_ETD(
			const SDC_ETD &i_src
	)	:
		TimeTree_Node_InteriorHelper<SDC_ETD>(i_src),
		SDC_ETD_Coefficients(i_src),
		tmpU(nullptr),
		tmp(nullptr),
		tmp2(nullptr)
	{
		_nPhi = i_src._nPhi;
		setEvalAvailable(EVAL_INTEGRATION);
	}

	~SDC_ETD()
	{
		clear();
	}

	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)	override
	{
		shackSDC = io_shackDict->getAutoRegistration<sweet::SDC::Shack>();

		return TimeTree_Node_InteriorHelper::shackRegistration(io_shackDict);
	}

	const std::vector<std::string> getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("SDC_ETD");
		retval.push_back("SDCETD");
		retval.push_back("ETD_SDC");
		retval.push_back("ETDSDC");
		return retval;
	}

	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override
	{
		o_ostream << i_prefix << "InteriorNode: 'SDCETD':" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Usage: ETDSDC(DeTermLinear,DeTermNonlinear)" << std::endl;
		o_ostream << i_prefix << "           Compute ETDSDC time integration (From Buvoli 15')." << std::endl;
		o_ostream << i_prefix << "           DeTermLinear has to support exponential, DeTermNonLinear - tendencies" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Parameters:" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "    - file=[str]:" << std::endl;
		o_ostream << i_prefix << "        Filename with ETD SDC coefficients" << std::endl;
		o_ostream << i_prefix << "    - r=[double]:" << std::endl;
		o_ostream << i_prefix << "        Residual-based stopping criteria based on lmax" << std::endl;
		o_ostream << i_prefix << "    - linear=[str]:" << std::endl;
		o_ostream << i_prefix << "    - l=[str]:" << std::endl;
		o_ostream << i_prefix << "        DETerm to treat linearly (supporting exponential)" << std::endl;
		o_ostream << i_prefix << "    - nonlinear=[str]:" << std::endl;
		o_ostream << i_prefix << "    - n=[str]:" << std::endl;
		o_ostream << i_prefix << "        DETerm to treat nonlinearly (supporting tendencies)" << std::endl;
		o_ostream << i_prefix << std::endl;

		return true;
	}
	
	bool setupArgumentInternals()
	{
		if (fileNameWithSDCCoefficients == "")
			fileNameWithSDCCoefficients = shackSDC->fileName;

		if (fileNameWithSDCCoefficients == "")
			fileNameWithSDCCoefficients = "params_SDC.sweet";

		if (!linearTerm.active && !nonlinearTerm.active)
			return error.set("No term set for SDC, stopping here!");

		return true;
	}

	bool setupTreeNodeByFunction(
			std::shared_ptr<TimeTreeIR::Function> &i_function,
			TimeTreeIR_2_TimeTreeNodes &i_tsAssemblation
	)	override
	{
		for (auto iter = i_function->arguments.begin(); iter != i_function->arguments.end(); iter++)
		{
			TimeTreeIR::Argument *a = iter->get();

			switch(a->argType)
			{
			case TimeTreeIR::Argument::ARG_TYPE_KEY_VALUE:
				if (a->key == "file" || a->key == "filename")
				{
					a->getValue(fileNameWithSDCCoefficients);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
					break;
				}

			case TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION:
				if (a->key == "l" || a->key == "linear")
				{
					if (linearTerm.active)
						return error.set("Linear term already set"+a->getNewLineDebugMessage());

					if (a->argType == TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION)
					{
						i_tsAssemblation.assembleTimeTreeNodeByFunction(
								a->function,
								linearTerm.tmp_TimeTreeNode
							);
					}
					else
					{
						i_tsAssemblation.assembleTimeTreeNodeByName(
								a->value,
								linearTerm.tmp_TimeTreeNode
							);
					}
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);

					linearTerm.active = true;
					break;
				}

				if (a->key == "n" || a->key == "nonlinear")
				{
					if (nonlinearTerm.active)
						return error.set("Nonlinear term already set"+a->getNewLineDebugMessage());

					if (a->argType == TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION)
					{
						i_tsAssemblation.assembleTimeTreeNodeByFunction(
								a->function,
								nonlinearTerm.tmp_TimeTreeNode
							);
					}
					else
					{
						i_tsAssemblation.assembleTimeTreeNodeByName(
								a->value,
								nonlinearTerm.tmp_TimeTreeNode
							);
					}

					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);

					nonlinearTerm.active = true;
					break;
				}

				return error.set("key_function not supported!"+a->getNewLineDebugMessage());


			case TimeTreeIR::Argument::ARG_TYPE_FUNCTION:
			case TimeTreeIR::Argument::ARG_TYPE_VALUE:
				if (linearTerm.active)
				{
					if (nonlinearTerm.active)
					{
						return error.set("Linear and nonlinear terms already set"+a->getNewLineDebugMessage());
					}
					else
					{
						nonlinearTerm.tmp_TimeTreeNode = std::shared_ptr<TimeTree_Node_Base>();

						if (a->argType == TimeTreeIR::Argument::ARG_TYPE_FUNCTION)
						{
							i_tsAssemblation.assembleTimeTreeNodeByFunction(
									a->function,
									nonlinearTerm.tmp_TimeTreeNode
								);
						}
						else
						{
							i_tsAssemblation.assembleTimeTreeNodeByName(
									a->value,
									nonlinearTerm.tmp_TimeTreeNode
								);
						}

						ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);

						nonlinearTerm.active = true;
						continue;
					}
				}
				else
				{
					linearTerm.tmp_TimeTreeNode = std::shared_ptr<TimeTree_Node_Base>();

					if (a->argType == TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION)
					{
						i_tsAssemblation.assembleTimeTreeNodeByFunction(
								a->function,
								linearTerm.tmp_TimeTreeNode
							);
					}
					else
					{
						i_tsAssemblation.assembleTimeTreeNodeByName(
								a->value,
								linearTerm.tmp_TimeTreeNode
							);
					}

					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);

					linearTerm.active = true;
					continue;
				}
				break;

			default:
				SWEETErrorFatal("Internal error");
				return error.set("Internal error");
			}
		}
		// provide debug message in case that something goes wrong with the arguments
		setDebugMessage(i_function->getDebugMessage());

		return setupArgumentInternals();
	}


	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper
	) override
	{
		// load SDC coefficients (sMat=aijl; nNodes=M; nIter=N) 
		loadSDCCoefficientsFromFile(fileNameWithSDCCoefficients);
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);
		SWEET_ASSERT(M >= 1);

		// Print out config in MULE format (for further processing later on)
		printMule();

		// Count phies: [0; M] used overall, [1; M] for W
		_nPhi = M + 1;

		TimeTree_Node_Base::registerTimeStepperEval(
				i_evalType,
				o_timeStepper
			);
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

		/*********************************************************************
		 * Setup time tree nodes
		 *********************************************************************/

		// Count time tree nodes
		int counterTimeTreeNodes = 0;

		if (linearTerm.active)
		{
			// phi0, phi1 pre-sweep and main-sweep, excluding tau = 1.
			counterTimeTreeNodes += 2 * (M - 1) + 2 * (M - 1);
		}

		if (nonlinearTerm.active)
		{
			// N(*) pre-sweep and main-sweep
			counterTimeTreeNodes += M + M;
		}

		if (linearTerm.active)
		{
			// phi1:M(N(*)) main-sweep
			counterTimeTreeNodes += (_nPhi - 1) * M  * (M - 1);
		}

		// Time stepper leaf node setup
		_timeTreeNodes.resize(counterTimeTreeNodes);
		_evalFuns.resize(counterTimeTreeNodes);

		// Setup evaluations
		int counterTimeTreeNodesCheck = 0;

		/*
		 * Pre-sweep
		 */
		if (nonlinearTerm.active)
		{
			// Get nonlinear tendencies N(*)
			nonlinearTerm.evalStartIdx_PreSweep = counterTimeTreeNodesCheck;
			for (int i = 0; i < M; i++)
			{
				int idx = nonlinearTerm.evalStartIdx_PreSweep + i;
				_timeTreeNodes[idx] = nonlinearTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck += M;
		}

		if (linearTerm.active)
		{
			// Get provisional solution: phi0 and phi1
			linearTerm.evalStartIdx_PreSweep_phi0 = counterTimeTreeNodesCheck;
			for (int i = 0; i < M - 1; i++)
			{
				int idx = linearTerm.evalStartIdx_PreSweep_phi0 + i;
				_timeTreeNodes[idx] = linearTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupByKeyValue("ExpIntegrationFunction", "phi0");
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EXPONENTIAL, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck += (M - 1);

			linearTerm.evalStartIdx_PreSweep_phi1 = counterTimeTreeNodesCheck;
			for (int i = 0; i < M - 1; i++)
			{
				int idx = linearTerm.evalStartIdx_PreSweep_phi1 + i;
				_timeTreeNodes[idx] = linearTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupByKeyValue("ExpIntegrationFunction", "phi1");
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EXPONENTIAL, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck += (M - 1);
		}

		/*
		 * Main-sweep
		 */
		if (nonlinearTerm.active)
		{
			// Get nonlinear tendencies N(*)
			nonlinearTerm.evalStartIdx_MainSweep = counterTimeTreeNodesCheck;

			for (int i = 0; i < M; i++)
			{
				int idx = nonlinearTerm.evalStartIdx_MainSweep + i;
				_timeTreeNodes[idx] = nonlinearTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck += M;
		}

		if (linearTerm.active)
		{
			// Get phi0 and phi1
			linearTerm.evalStartIdx_MainSweep_phi0 = counterTimeTreeNodesCheck;
			for (int i = 0; i < M - 1; i++)
			{
				int idx = linearTerm.evalStartIdx_MainSweep_phi0 + i;
				_timeTreeNodes[idx] = linearTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupByKeyValue("ExpIntegrationFunction", "phi0");
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EXPONENTIAL, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck += (M - 1);

			linearTerm.evalStartIdx_MainSweep_phi1 = counterTimeTreeNodesCheck;
			for (int i = 0; i < M - 1; i++)
			{
				int idx = linearTerm.evalStartIdx_MainSweep_phi1 + i;
				_timeTreeNodes[idx] = linearTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupByKeyValue("ExpIntegrationFunction", "phi1");
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EXPONENTIAL, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck+= (M - 1);
		}

		if (linearTerm.active)
		{
			// Calculate Wi;i+1: get phi[1:M](N(*))
			linearTerm.evalStartIdx_MainSweep_W = counterTimeTreeNodesCheck;

			for (int i = 0; i < M - 1; i++)
			{	
				for (int l = 0; l < M; l++)
				{
					for (int j = 0; j < _nPhi - 1; j++)
					{
						// Add to the sum phi_{j+1}
						int idx = linearTerm.evalStartIdx_MainSweep_W + i*M*(_nPhi-1) + l*(_nPhi - 1) + j;
						_timeTreeNodes[idx] = linearTerm.tmp_TimeTreeNode->getInstanceCopy();
					
						_timeTreeNodes[idx]->setupByKeyValue("ExpIntegrationFunction", "phi" + std::to_string(j+1));
						_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EXPONENTIAL, &_evalFuns[idx]);
						ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
					}
				}
			}
			counterTimeTreeNodesCheck += M * (_nPhi - 1) * (M - 1);
		}

		SWEET_ASSERT(counterTimeTreeNodes == counterTimeTreeNodesCheck);

		/*********************************************************************
		 * Setup data containers
		 *********************************************************************/

		// Count time tree nodes
		int counterDataContainer = 0;

		if (nonlinearTerm.active)
		{
			// Containers for tendencies k0 and k1, and sum Wi;i+1
			counterDataContainer += M;
			counterDataContainer += M;
			counterDataContainer += 1;
		}

		// Temporary buffers: tmpU, tmp1 and tmp2
		counterDataContainer += 3;
		
		// Container for provisional_U
		counterDataContainer += M;

		// Allocate containers
		_tmpDataContainer.resize(counterDataContainer);
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();

		// Setup convenience access arrays
		int counterDataContainerCheck = 0;

		if (nonlinearTerm.active)
		{
			nonlinearTerm.tendencies_k0 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			nonlinearTerm.tendencies_k1 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			nonlinearTerm.accW = _tmpDataContainer[counterDataContainerCheck];
        	counterDataContainerCheck += 1;
		}

		tmpU = _tmpDataContainer[counterDataContainerCheck];
		counterDataContainerCheck += 1;

		tmp = _tmpDataContainer[counterDataContainerCheck];
		counterDataContainerCheck += 1;

		tmp2 = _tmpDataContainer[counterDataContainerCheck];
		counterDataContainerCheck += 1;

		provisional_U = &_tmpDataContainer[counterDataContainerCheck];
        counterDataContainerCheck += M;

		SWEET_ASSERT(counterDataContainer == counterDataContainerCheck);
		return true;
	}

	void clear() override
	{
		TimeTree_Node_InteriorHelper::clear();
	}

	std::shared_ptr<TimeTree_Node_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<TimeTree_Node_Base>(new SDC_ETD(*this));
	}

	bool setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		/*
		 * Pre-sweep
		 */
		if (nonlinearTerm.active)
		{
			for (int i = 0; i < M - 1; i++)
			{
				int idx = nonlinearTerm.evalStartIdx_PreSweep + i;
				double dt = _dt*deltaTau(i);
				
				_timeTreeNodes[idx]->setTimeStepSize(dt);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
		}

		if (linearTerm.active)
		{
			for (int i = 0; i < M - 1; i++)
			{
				int idx = linearTerm.evalStartIdx_PreSweep_phi1 + i;
				double dt = _dt*deltaTau(i+1);
				
				_timeTreeNodes[idx]->setTimeStepSize(dt);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}

			for (int i = 0; i < M - 1; i++)
            {
                int idx = linearTerm.evalStartIdx_PreSweep_phi0 + i;
                double dt = _dt*deltaTau(i+1);

                _timeTreeNodes[idx]->setTimeStepSize(dt);
                ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
            }

		}

		/*
		 * Main-sweep
		 */
		if (nonlinearTerm.active)
		{
			for (int i = 0; i < M; i++)
			{
				int idx = nonlinearTerm.evalStartIdx_MainSweep + i;
				double dt = _dt*deltaTau(i);
				
				_timeTreeNodes[idx]->setTimeStepSize(dt);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
		}

		if (linearTerm.active)
		{
			for (int i = 0; i < M - 1; i++)
			{
				int idx = linearTerm.evalStartIdx_MainSweep_phi0 + i;
				double dt = _dt*deltaTau(i+1);
				
				_timeTreeNodes[idx]->setTimeStepSize(dt);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}

			for (int i = 0; i < M - 1; i++)
			{
				int idx = linearTerm.evalStartIdx_MainSweep_phi1 + i;
				double dt = _dt*deltaTau(i+1);
				
				_timeTreeNodes[idx]->setTimeStepSize(dt);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
		}

		if (linearTerm.active)
		{
			for (int i = 0; i < M - 1; i++)
			{
				double dt = _dt*deltaTau(i+1);
				for (int l = 0; l < M; l++)
				{
					for (int j = 0; j < _nPhi - 1; j++)
					{
						int idx = linearTerm.evalStartIdx_MainSweep_W + i*M*(_nPhi-1) + l*(_nPhi - 1) + j;
					
						_timeTreeNodes[idx]->setTimeStepSize(dt);
						ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
					}
				}
			}
		}
		return true;
	}

	//! Swap k+1 and k values for next iteration (swap pointers)
	void _swapK()
	{
		if (nonlinearTerm.active)
		{
			std::swap(nonlinearTerm.tendencies_k0, nonlinearTerm.tendencies_k1);
		}
	}

private:
	//! Provisional solution
	void _sweep_Node2Node_Pre(
			const sweet::Data::GenericContainer::Base& i_U,
			double i_simulationTime
	)
	{
		// Start with U_0
		tmpU->op_setVector(i_U);

		// Propagate provisional solution: U = phi0(hL) * U + hi phi1(hL) * N(U)
		for (size_t m = 0; m < M; m++)
		{
			double t = i_simulationTime + _dt*tau(m-1);
			double dTau = _dt*deltaTau(m);

			// Save N(U) for U_0
			if (nonlinearTerm.active && m == 0)
				evalTimeStepper(nonlinearTerm.evalStartIdx_PreSweep, *tmpU, *nonlinearTerm.tendencies_k0[0], t);

			if (m != 0)
			{
				if (linearTerm.active)
				{
					// U = phi0(hL) * U + hi phi1(hL) * N(U)
					//                       ========
					if (nonlinearTerm.active)
						evalTimeStepper(linearTerm.evalStartIdx_PreSweep_phi1 + m - 1, *nonlinearTerm.tendencies_k0[m - 1], *tmp, t);

					// U = phi0(hL) * U + hi phi1(hL) * N(U)
					//     ========
					evalTimeStepper(linearTerm.evalStartIdx_PreSweep_phi0 + m - 1, *tmpU, *tmp2, t);

					// U = phi0(hL) * U + hi phi1(hL) * N(U)
					//                    ==================
					if (nonlinearTerm.active)
						tmp2->op_addScalarMulVector(dTau, *tmp);	
					
					tmpU->swap(*tmp2);
				} else {
					// if linearTerm is inactive (l=0) phi0 = phi1 = 1.

					// U = U + hi * N(U)
					//     =============
					if (nonlinearTerm.active)
						tmpU->op_addScalarMulVector(dTau, *nonlinearTerm.tendencies_k0[m - 1]);
				}
			}
			
			// U = phi0(hL) * U + hi phi1(hL) * N(U)
			//                                  ====
			if (nonlinearTerm.active && m != 0)
                evalTimeStepper(nonlinearTerm.evalStartIdx_PreSweep + m, *tmpU, *nonlinearTerm.tendencies_k0[m], t);

			// Store provisional solution
			provisional_U[m]->op_setVector(*tmpU);
		}
	}

	//! Main SDC sweep
	void _sweep_Node2Node_Main(
			const sweet::Data::GenericContainer::Base& i_U,
			double i_simulationTime
	) {
		// Start with provisional_U[0], which is U_0 too
		tmpU->op_setVector(*provisional_U[0]);

		// Iterate over nodes with node-to-node scheme: U = phi0(hL) * U + hi phi1(hL) * [N^{k+1}(U) - N^k(U)] + sum_l sum_j h a_jl^i phi_j+1(hL) N^k(U_l)
		for (size_t m = 0; m < M; m++)
		{
			double t = i_simulationTime + _dt*tau(m-1);
			double dTau = _dt*deltaTau(m);

			if (linearTerm.active && m != 0)
			{
				//  U = phi0(hL) * U + hi phi1(hL) * [N^{k+1}(U) - N^k(U)] + sum_l sum_j h a_jl^i phi_j+1(hL) N^k(U_l)
				//      ============
				evalTimeStepper(linearTerm.evalStartIdx_MainSweep_phi0 + m - 1, *provisional_U[m-1], *tmp, t);
				tmpU->op_setVector(*tmp);
			}

			if (nonlinearTerm.active && m != 0)
			{
				//  U = phi0(hL) * U + hi phi1(hL) * [N^{k+1}(U) - N^k(U)] + sum_l sum_j h a_jl^i phi_j+1(hL) N^k(U_l)
				//                                    ====================
				tmp->op_setVectorPlusScalarMulVector(*nonlinearTerm.tendencies_k1[m-1], -1, *nonlinearTerm.tendencies_k0[m-1]);
				
				//  U = phi0(hL) * U + hi phi1(hL) * [N^{k+1}(U) - N^k(U)] + sum_l sum_j h a_jl^i phi_j+1(hL) N^k(U_l)
				//                        ===============================
				if (linearTerm.active)
				{
					evalTimeStepper(linearTerm.evalStartIdx_MainSweep_phi1 + m - 1, *tmp, *tmp2, t);
					tmp2->op_mulScalar(dTau);
				} else {
					// if linearTerm is inactive (l=0) phi0 = phi1 = 1.
					tmp2->op_setVector(*tmp);
					tmp2->op_mulScalar(dTau);
				}

				tmpU->op_addVector(*tmp2);

				// Accumulate Wi;i+1
				nonlinearTerm.accW->op_setZero();

				if (linearTerm.active)
				{
					// Iterate through timesteps N^k(U_l)
					for (size_t l = 0; l < M; l++)
					{
						// Iterate through phies_{j+1}
						for (size_t j = 0; j < _nPhi - 1; j++)
						{
							//  U = phi0(hL) * U + hi phi1(hL) * [N^{k+1}(U) - N^k(U)] + sum_l sum_j h a_jl^i phi_j+1(hL) N^k(U_l)
				            //                                                                       ========
							double aij = A(m-1, l, j) * dTau;

							//  U = phi0(hL) * U + hi phi1(hL) * [N^{k+1}(U) - N^k(U)] + sum_l sum_j h a_jl^i phi_j+1(hL) N^k(U_l)
				            //                                                                                ====================
							int idx = linearTerm.evalStartIdx_MainSweep_W + (m-1)*M*(_nPhi-1) + l*(_nPhi-1) + j;
							evalTimeStepper(idx, *nonlinearTerm.tendencies_k0[l], *tmp2, t);

							//  U = phi0(hL) * U + hi phi1(hL) * [N^{k+1}(U) - N^k(U)] + sum_l sum_j h a_jl^i phi_j+1(hL) N^k(U_l)
				            //                                                                       ============================
							tmp2->op_mulScalar(aij);
							nonlinearTerm.accW->op_addVector(*tmp2);
						}
					}
				} else {
					// if linearTerm is inactive phi_{j+1} = 1/(j+1)!

					// Iterate through timesteps N^k(U_l)
					for (size_t l = 0; l < M; l++)
					{
						double phi_j1 = 1.0;

						// Iterate through phies_{j+1}
						for (size_t j = 0; j < _nPhi - 1; j++)
						{
							//  U = phi0(hL) * U + hi phi1(hL) * [N^{k+1}(U) - N^k(U)] + sum_l sum_j h a_jl^i phi_j+1(hL) N^k(U_l)
				            //                                                                       ========
							double aij = A(m-1, l, j) * dTau;

							//  Update factorial
							phi_j1 /= (j + 1);

							//  U = phi0(hL) * U + hi phi1(hL) * [N^{k+1}(U) - N^k(U)] + sum_l sum_j h a_jl^i phi_j+1(hL) N^k(U_l)
				            //                                                                                ====================
							tmp2->op_setVector(*nonlinearTerm.tendencies_k0[l]);
							tmp2->op_mulScalar(aij * phi_j1);

							//  U = phi0(hL) * U + hi phi1(hL) * [N^{k+1}(U) - N^k(U)] + sum_l sum_j h a_jl^i phi_j+1(hL) N^k(U_l)
				            //                                                                       ============================
							nonlinearTerm.accW->op_addVector(*tmp2);
						}
					}
				}
				
				tmpU->op_addVector(*nonlinearTerm.accW);
			}

			// Store sweep state U_m
			provisional_U[m]->op_setVector(*tmpU);	

			// Save tendencies of updated state N_{k+1}(U_m)
			if (nonlinearTerm.active)
				evalTimeStepper(nonlinearTerm.evalStartIdx_MainSweep + m, *provisional_U[m], *nonlinearTerm.tendencies_k1[m], t);
		}
		if (nonlinearTerm.active)
		{
			_swapK();
		}
	}

	void _sweep_Post(
			const sweet::Data::GenericContainer::Base& i_U,
			sweet::Data::GenericContainer::Base& o_U,
			double i_timeStamp
	) {
		// Copy last step
		o_U.op_setVector(*tmpU);
		return;
	}

	bool _eval_integration(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)	override
	{
		// initialize nodes values and state
		_sweep_Node2Node_Pre(i_U, i_simulationTime);

		// perform sweeps
		for (size_t k = 0; k < nIters; k++) {
			_sweep_Node2Node_Main(i_U, i_simulationTime);
		}

		_sweep_Post(i_U, o_U, i_simulationTime);

		return true;
	}

	void printMule(const std::string &i_prefixStr = "")
	{
		std::cout << i_prefixStr << "SDC_ETD:" << std::endl;
		printSDCInformation("[MULE] sdc.");
	}

	void print(const std::string &i_prefixStr = "")
	{
		std::string newPrefix = i_prefixStr + "  ";
		std::cout << newPrefix << "SDC_ETD(" << std::endl;
		std::cout << newPrefix << "  file: " << fileNameWithSDCCoefficients << std::endl;
		printSDCInformation(newPrefix);
	}
};

}}}

#endif
