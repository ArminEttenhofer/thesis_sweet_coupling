#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_SDC_CLASSIC_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_SDC_CLASSIC_HPP

#include <vector>
#include <string>
#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/InteriorNodes/SDC_Coefficients.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>


namespace sweet {
namespace TimeTree {
namespace InteriorNodes {


class SDC_Classic	:
	public TimeTree_Node_InteriorHelper<SDC_Classic>,
	public SDC_Coefficients
{
private:
	//! File name of SWEETDict with SDC coefficients
	std::string fileNameWithSDCCoefficients;


	/*!
	 * Time integrator term
	 */
	class TimeIntegratorTerm {
	public:
		bool active;

		//! temporary node which can be used during setup
		std::shared_ptr<TimeTree_Node_Base> tmp_TimeTreeNode;
		std::shared_ptr<TimeTree_Node_Base> tmp_TimeTreeNodeTendencies;

		int evalStartIdx_Tendencies;
		int evalStartIdx_TimeIntegrator_PreSweep;
		int evalStartIdx_TimeIntegrator_MainSweep;

		//! Storage space for implicit tendencies
		sweet::Data::GenericContainer::Base** dataContainer_Tendencies_k0;
		sweet::Data::GenericContainer::Base** dataContainer_Tendencies_k1;

		//! Storage space for implicit tendencies
		sweet::Data::GenericContainer::Base** dataContainer_Integral_U0;
		sweet::Data::GenericContainer::Base** dataContainer_Integral_U1;


		TimeIntegratorTerm()	:
			active(false),
			evalStartIdx_Tendencies(-1),
			evalStartIdx_TimeIntegrator_PreSweep(-1),
			evalStartIdx_TimeIntegrator_MainSweep(-1),
			dataContainer_Tendencies_k0(nullptr),
			dataContainer_Tendencies_k1(nullptr),
			dataContainer_Integral_U0(nullptr),
			dataContainer_Integral_U1(nullptr)
		{}
	};
	TimeIntegratorTerm timeIntegratorTerm;


	/*!
	 * Explicit term
	 */
	class ExplicitTerm {
	public:
		bool active;

		///! temporary node which can be used during setup
		std::shared_ptr<TimeTree_Node_Base> tmp_TimeTreeNode;

		int evalStartIdx_Tendencies;

		// Storage space for explicit tendencies
		sweet::Data::GenericContainer::Base** dataContainer_Tendencies_k0;
		sweet::Data::GenericContainer::Base** dataContainer_Tendencies_k1;

		ExplicitTerm()	:
			active(false),
			evalStartIdx_Tendencies(-1),
			dataContainer_Tendencies_k0(nullptr),
			dataContainer_Tendencies_k1(nullptr)
		{}
	};
	ExplicitTerm explicitTerm;


	/*!
	 * Implicit term
	 */
	class ImplicitTerm {
	public:
		bool active;

		//! temporary node which can be used during setup
		std::shared_ptr<TimeTree_Node_Base> tmp_TimeTreeNode;

		int evalStartIdx_Tendencies;
		int evalStartIdx_EulerBackward_PreSweep;
		int evalStartIdx_EulerBackward_MainSweep;

		//! Storage space for implicit tendencies
		sweet::Data::GenericContainer::Base** dataContainer_Tendencies_k0;
		sweet::Data::GenericContainer::Base** dataContainer_Tendencies_k1;


		ImplicitTerm()	:
			active(false),
			evalStartIdx_Tendencies(-1),
			evalStartIdx_EulerBackward_PreSweep(-1),
			evalStartIdx_EulerBackward_MainSweep(-1),
			dataContainer_Tendencies_k0(nullptr),
			dataContainer_Tendencies_k1(nullptr)
		{}
	};
	ImplicitTerm implicitTerm;




	//! Intermediate states during SDC
	sweet::Data::GenericContainer::Base* dataContainer_intermediateState;

	//! Temporary state
	sweet::Data::GenericContainer::Base *tmp;
	sweet::Data::GenericContainer::Base *tmp2;

	//! Special state which is used to access the last updated state in case of no endpoint update
	sweet::Data::GenericContainer::Base* dataContainer_lastUpdateState;

public:
	SDC_Classic()	:
		dataContainer_intermediateState(nullptr),
		tmp(nullptr),
		tmp2(nullptr),
		dataContainer_lastUpdateState(nullptr)
	{
		setEvalAvailable(EVAL_INTEGRATION);
	}

	SDC_Classic(
			const SDC_Classic &i_src
	)	:
		TimeTree_Node_InteriorHelper<SDC_Classic>(i_src),
		SDC_Coefficients(i_src),
		dataContainer_intermediateState(nullptr),
		tmp(nullptr),
		tmp2(nullptr),
		dataContainer_lastUpdateState(nullptr)
	{
		setEvalAvailable(EVAL_INTEGRATION);
	}

	~SDC_Classic()
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


	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("SDC_CLASSIC");
		retval.push_back("SDCCLASSIC");
		retval.push_back("SDC_Classic");
		retval.push_back("SDCClassic");
		return retval;
	}

	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override
	{
		o_ostream << i_prefix << "InteriorNode: 'SDCClassic':" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Usage: SDCClassic(DeTermImplicit,DeTermExplicit,[parameter1],...)" << std::endl;
		o_ostream << i_prefix << "           Compute Classical SDC IMEX time integration (From Dutt et al.'s paper, using NodesToNodes)." << std::endl;
		o_ostream << i_prefix << "           DeTermStiff has to support 'tendencies' and 'backward_euler'" << std::endl;
		o_ostream << i_prefix << "           DeTermNonStiff has to support 'tendencies'" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Parameters:" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "    - file=[str]:" << std::endl;
		o_ostream << i_prefix << "        Filename with SDC coefficients" << std::endl;
		o_ostream << i_prefix << "    - r=[double]:" << std::endl;
		o_ostream << i_prefix << "        Residual-based stopping criteria based on lmax" << std::endl;
		o_ostream << i_prefix << "    - explicit=[str]:" << std::endl;
		o_ostream << i_prefix << "    - e=[str]:" << std::endl;
		o_ostream << i_prefix << "        DETerm to treat explicitly" << std::endl;
		o_ostream << i_prefix << "    - implicit=[str]:" << std::endl;
		o_ostream << i_prefix << "    - i=[str]:" << std::endl;
		o_ostream << i_prefix << "        DETerm to treat implicitly" << std::endl;
		o_ostream << i_prefix << "    - timeIntegrator=[str]:" << std::endl;
		o_ostream << i_prefix << "    - ti=[str]:" << std::endl;
		o_ostream << i_prefix << "        DETerm to treat as a time integrator" << std::endl;
		o_ostream << i_prefix << std::endl;

		return true;
	}

	bool _setupArgumentInternals()
	{
		if (fileNameWithSDCCoefficients == "")
			fileNameWithSDCCoefficients = shackSDC->fileName;

		if (fileNameWithSDCCoefficients == "")
			fileNameWithSDCCoefficients = "params_SDC.sweet";

		if (!implicitTerm.active && !explicitTerm.active && !timeIntegratorTerm.active)
			return error.set("No term set for SDC, stopping here!");

		if (timeIntegratorTerm.active)
		{
			if (!timeIntegratorTerm.tmp_TimeTreeNodeTendencies)
				return error.set("Using a time integrator also requires setting timeIntegratorTendencies=... or tit=...");
		}

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

				/*
				 * NO BREAK MISSING HERE!
				 * We continue processing things
				 */

			case TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION:
				if (a->key == "i" || a->key == "implicit")
				{
					if (implicitTerm.active)
						return error.set("Implicit term already set"+a->getNewLineDebugMessage());

					if (a->argType == TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION)
					{
						i_tsAssemblation.assembleTimeTreeNodeByFunction(
								a->function,
								implicitTerm.tmp_TimeTreeNode
							);
					}
					else
					{
						i_tsAssemblation.assembleTimeTreeNodeByName(
								a->value,
								implicitTerm.tmp_TimeTreeNode
							);
					}
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);

					implicitTerm.active = true;
					break;
				}

				if (a->key == "e" || a->key == "explicit")
				{
					if (explicitTerm.active)
						return error.set("Explicit term already set"+a->getNewLineDebugMessage());

					if (a->argType == TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION)
					{
						i_tsAssemblation.assembleTimeTreeNodeByFunction(
								a->function,
								explicitTerm.tmp_TimeTreeNode
							);
					}
					else
					{
						i_tsAssemblation.assembleTimeTreeNodeByName(
								a->value,
								explicitTerm.tmp_TimeTreeNode
							);
					}

					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);

					explicitTerm.active = true;
					break;
				}


				if (a->key == "ti" || a->key == "timeIntegrator")
				{
					if (timeIntegratorTerm.active)
						return error.set("Explicit term already set"+a->getNewLineDebugMessage());

					if (a->argType == TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION)
					{
						i_tsAssemblation.assembleTimeTreeNodeByFunction(
								a->function,
								timeIntegratorTerm.tmp_TimeTreeNode
							);
					}
					else
					{
						i_tsAssemblation.assembleTimeTreeNodeByName(
								a->value,
								timeIntegratorTerm.tmp_TimeTreeNode
							);
					}

					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);

					timeIntegratorTerm.active = true;
					break;
				}

				if (a->key == "tit" || a->key == "timeIntegratorTendencies")
				{
					if (!timeIntegratorTerm.active)
						return error.set("You first need to set timeIntegrator=... before setting the tendencies"+a->getNewLineDebugMessage());

					if (a->argType == TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION)
					{
						i_tsAssemblation.assembleTimeTreeNodeByFunction(
								a->function,
								timeIntegratorTerm.tmp_TimeTreeNodeTendencies
							);
					}
					else
					{
						i_tsAssemblation.assembleTimeTreeNodeByName(
								a->value,
								timeIntegratorTerm.tmp_TimeTreeNodeTendencies
							);
					}

					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
					break;
				}

				return error.set("key_function not supported!"+a->getNewLineDebugMessage());


			case TimeTreeIR::Argument::ARG_TYPE_FUNCTION:
			case TimeTreeIR::Argument::ARG_TYPE_VALUE:
				if (implicitTerm.active)
				{
					if (explicitTerm.active)
					{
						return error.set("linear and explicit terms already set"+a->getNewLineDebugMessage());
					}
					else
					{
						explicitTerm.tmp_TimeTreeNode = std::shared_ptr<TimeTree_Node_Base>();

						if (a->argType == TimeTreeIR::Argument::ARG_TYPE_FUNCTION)
						{
							i_tsAssemblation.assembleTimeTreeNodeByFunction(
									a->function,
									explicitTerm.tmp_TimeTreeNode
								);
						}
						else
						{
							i_tsAssemblation.assembleTimeTreeNodeByName(
									a->value,
									explicitTerm.tmp_TimeTreeNode
								);
						}

						ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);

						explicitTerm.active = true;
						continue;
					}
				}
				else
				{
					implicitTerm.tmp_TimeTreeNode = std::shared_ptr<TimeTree_Node_Base>();

					if (a->argType == TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION)
					{
						i_tsAssemblation.assembleTimeTreeNodeByFunction(
								a->function,
								implicitTerm.tmp_TimeTreeNode
							);
					}
					else
					{
						i_tsAssemblation.assembleTimeTreeNodeByName(
								a->value,
								implicitTerm.tmp_TimeTreeNode
							);
					}

					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);

					implicitTerm.active = true;
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

		return _setupArgumentInternals();
	}


	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper
	) override
	{
		/*********************************************************************
		 * We first load all SDC coefficients
		 *********************************************************************/
		loadSDCCoefficientsFromFile(fileNameWithSDCCoefficients);
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);
		SWEET_ASSERT(M >= 1);

		// Print out config in MULE format (for further processing later on)
		printMule();

		// Sanity check
		if (!rightIsNode && postSweep == POST_SWEEP_LASTNODE)
			return error.set("Error! Cannot use post sweep LASTNODE if not rightIsNode");

		/*********************************************************************
		 * First we care about the eval request of the caller
		 *********************************************************************/
		TimeTree_Node_Base::registerTimeStepperEval(
				i_evalType,
				o_timeStepper
			);
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);


		/*********************************************************************
		 * Setup time tree nodes
		 *********************************************************************/
		int counterTimeTreeNodes = 0;

		if (timeIntegratorTerm.active)
		{
			counterTimeTreeNodes += 1;

			// For first estimate and also sweeps
			counterTimeTreeNodes += 2*M;
		}

		if (explicitTerm.active)
		{
			counterTimeTreeNodes += 1;
		}

		if (implicitTerm.active)
		{
			// One for explicit evaluation
			counterTimeTreeNodes += 1;

			// For first estimate and also sweeps
			counterTimeTreeNodes += 2*M;
		}

		/*
		 * Time stepper leaf node setup
		 */
		_timeTreeNodes.resize(counterTimeTreeNodes);
		_evalFuns.resize(counterTimeTreeNodes);



		/*
		 * Setup evaluations
		 */
		int counterTimeTreeNodesCheck = 0;



		if (timeIntegratorTerm.active)
		{
			timeIntegratorTerm.evalStartIdx_Tendencies = counterTimeTreeNodesCheck;

			/*
			 * Tendency evaluations
			 */
			int idx = timeIntegratorTerm.evalStartIdx_Tendencies;
			_timeTreeNodes[idx] = timeIntegratorTerm.tmp_TimeTreeNodeTendencies->getInstanceCopy();
			_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);

			counterTimeTreeNodesCheck += 1;


			/*
			 * Integration for pre sweep
			 */
			timeIntegratorTerm.evalStartIdx_TimeIntegrator_PreSweep = counterTimeTreeNodesCheck;
			for (int i = 0; i < M; i++)
			{
				int idx = timeIntegratorTerm.evalStartIdx_TimeIntegrator_PreSweep + i;
				_timeTreeNodes[idx] = timeIntegratorTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_INTEGRATION, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck += M;


			/*
			 * Integration for main sweeps
			 */
			timeIntegratorTerm.evalStartIdx_TimeIntegrator_MainSweep = counterTimeTreeNodesCheck;
			for (int i = 0; i < M; i++)
			{
				int idx = timeIntegratorTerm.evalStartIdx_TimeIntegrator_MainSweep + i;
				_timeTreeNodes[idx] = timeIntegratorTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_INTEGRATION, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck += M;

			timeIntegratorTerm.tmp_TimeTreeNode.reset();
		}


		if (explicitTerm.active)
		{
			explicitTerm.evalStartIdx_Tendencies = counterTimeTreeNodesCheck;

			int idx = explicitTerm.evalStartIdx_Tendencies;
			_timeTreeNodes[idx] = explicitTerm.tmp_TimeTreeNode->getInstanceCopy();
			_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
			counterTimeTreeNodesCheck += 1;

			explicitTerm.tmp_TimeTreeNode.reset();
		}


		if (implicitTerm.active)
		{
			implicitTerm.evalStartIdx_Tendencies = counterTimeTreeNodesCheck;

			/*
			 * Tendency evaluations
			 */
			int idx = implicitTerm.evalStartIdx_Tendencies;
			_timeTreeNodes[idx] = implicitTerm.tmp_TimeTreeNode->getInstanceCopy();
			_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);

			counterTimeTreeNodesCheck += 1;


			/*
			 * Backward Euler time steps
			 */
			implicitTerm.evalStartIdx_EulerBackward_PreSweep = counterTimeTreeNodesCheck;
			for (int i = 0; i < M; i++)
			{
				int idx = implicitTerm.evalStartIdx_EulerBackward_PreSweep + i;
				_timeTreeNodes[idx] = implicitTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EULER_BACKWARD, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck += M;


			/*
			 * Time integrator
			 */
			implicitTerm.evalStartIdx_EulerBackward_MainSweep = counterTimeTreeNodesCheck;
			for (int i = 0; i < M; i++)
			{
				int idx = implicitTerm.evalStartIdx_EulerBackward_MainSweep + i;
				_timeTreeNodes[idx] = implicitTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EULER_BACKWARD, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck += M;

			implicitTerm.tmp_TimeTreeNode.reset();
		}


		/*********************************************************************
		 * Setup data containers
		 *********************************************************************/

		/*
		 * We need data containers for
		 *
		 *  - 2 iterations (k0/k1, nNodes each)
		 *  - solving & storing the quadrature results (nNodes)
		 *  - implicit and explicit terms (nNodes each)
		 */
		SWEET_ASSERT(M >= 1);

		int counterDataContainer = 0;

		if (explicitTerm.active)
			counterDataContainer += 2*M;

		if (implicitTerm.active)
			counterDataContainer += 2*M;

		if (timeIntegratorTerm.active)
		{
			counterDataContainer += 2*M;
			counterDataContainer += 2*M;
		}

		// current state
		counterDataContainer += 1;

		// temporary buffer 1 and 2
		counterDataContainer += 2;


		/*
		 * Setup data containers
		 */
		_tmpDataContainer.resize(counterDataContainer);
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();

		/*
		 * Setup convenience access arrays
		 */
		int counterDataContainerCheck = 0;

		if (timeIntegratorTerm.active)
		{
			timeIntegratorTerm.dataContainer_Tendencies_k0 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			timeIntegratorTerm.dataContainer_Tendencies_k1 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			timeIntegratorTerm.dataContainer_Integral_U0 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			timeIntegratorTerm.dataContainer_Integral_U1 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;
		}

		if (explicitTerm.active)
		{
			explicitTerm.dataContainer_Tendencies_k0 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			explicitTerm.dataContainer_Tendencies_k1 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;
		}

		if (implicitTerm.active)
		{
			implicitTerm.dataContainer_Tendencies_k0 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			implicitTerm.dataContainer_Tendencies_k1 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;
		}

		/*
		 * Without parallel evaluation we just need one container
		 */
		dataContainer_intermediateState = _tmpDataContainer[counterDataContainerCheck];
		counterDataContainerCheck += 1;

		tmp = _tmpDataContainer[counterDataContainerCheck];
		counterDataContainerCheck += 1;

		tmp2 = _tmpDataContainer[counterDataContainerCheck];
		counterDataContainerCheck += 1;

		SWEET_ASSERT(counterDataContainer == counterDataContainerCheck);

		return true;

	}


	void clear() override
	{
		TimeTree_Node_InteriorHelper::clear();
	}


	std::shared_ptr<TimeTree_Node_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<TimeTree_Node_Base>(new SDC_Classic(*this));
	}


	bool setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		/*
		 * Set time step size for time integrator terms
		 */
		if (timeIntegratorTerm.active)
		{
			for (int i = 0; i < M; i++)
			{
				double idx = timeIntegratorTerm.evalStartIdx_TimeIntegrator_PreSweep+i;
				double dt = _dt*deltaTau(i);

				_timeTreeNodes[idx]->setTimeStepSize(dt);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}

			for (int i = 0; i < M; i++)
			{
				double idx = timeIntegratorTerm.evalStartIdx_TimeIntegrator_MainSweep+i;
				double dt = _dt*deltaTau(i);

				_timeTreeNodes[idx]->setTimeStepSize(dt);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
		}


		/*
		 * Set time step size for explicit terms
		 */
		if (explicitTerm.active)
		{
			// Nothing to do since it's just about tendencies
			//_timeTreeNodes[explicitTerm.evalStartIdx_Tendencies]->setTimeStepSize(i_dt);
		}


		/*
		 * Set time step size for implicit terms
		 */
		if (implicitTerm.active)
		{
			for (int i = 0; i < M; i++)
			{
				double idx = implicitTerm.evalStartIdx_EulerBackward_PreSweep+i;
				double dt = _dt*deltaTau(i);

				_timeTreeNodes[idx]->setTimeStepSize(dt);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}

			for (int i = 0; i < M; i++)
			{
				double idx = implicitTerm.evalStartIdx_EulerBackward_MainSweep+i;
				double dt = _dt*deltaTau(i);

				_timeTreeNodes[idx]->setTimeStepSize(dt);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
		}


		return true;
	}


	/*!
	 * Compute time tendencies at quadrature point "m"
	 */
	void _timeTendencies(
			sweet::Data::GenericContainer::Base *io_state,
			int m,	// Tau node id
			double i_simulationTime
	)
	{
		double t = i_simulationTime + _dt*tau(m);

		if (timeIntegratorTerm.active)
		{
			int idx = timeIntegratorTerm.evalStartIdx_Tendencies;
			evalTimeStepper(idx, *io_state, *timeIntegratorTerm.dataContainer_Tendencies_k1[m], t);
		}

		if (explicitTerm.active)
		{
			int idx = explicitTerm.evalStartIdx_Tendencies;
			evalTimeStepper(idx, *io_state, *explicitTerm.dataContainer_Tendencies_k1[m], t);
		}

		if (implicitTerm.active)
		{
			int idx = implicitTerm.evalStartIdx_Tendencies;
			evalTimeStepper(idx, *io_state, *implicitTerm.dataContainer_Tendencies_k1[m], t);
		}
	}


	void _swapK()
	{
		/*
		 * Swap k+1 and k values for next iteration (or end-point update)
		 *
		 * (Just swap the pointers into the vector structure)
		 */

		if (timeIntegratorTerm.active)
		{
			std::swap(timeIntegratorTerm.dataContainer_Tendencies_k0, timeIntegratorTerm.dataContainer_Tendencies_k1);
			std::swap(timeIntegratorTerm.dataContainer_Integral_U0, timeIntegratorTerm.dataContainer_Integral_U1);
		}

		if (implicitTerm.active)
			std::swap(implicitTerm.dataContainer_Tendencies_k0, implicitTerm.dataContainer_Tendencies_k1);

		if (explicitTerm.active)
			std::swap(explicitTerm.dataContainer_Tendencies_k0, explicitTerm.dataContainer_Tendencies_k1);
	}


private:
	void _sweep_Node2Node_Pre(
			const sweet::Data::GenericContainer::Base& i_U,
			double i_simulationTime
	)
	{
		sweet::Data::GenericContainer::Base *U = dataContainer_intermediateState;

		// Start with U_0
		U->op_setVector(i_U);

		/*
		 * Integrate from \tau_{m-1} to \tau{m}
		 */
		for (size_t m = 0; m < M; m++)
		{
			/*
			 * Integrate from \tau_{m-1} to \tau_{m}
			 */

			/*
			 * We make this match the main sweeper updates.
			 *
			 * For the explicit time integrator, we need the time tendencies for the first node
			 */
			if (explicitTerm.active)
			{
				double t = i_simulationTime + _dt*tau(m-1);
				if (m == 0)
					evalTimeStepper(explicitTerm.evalStartIdx_Tendencies, *U, *tmp2, t);
			}


			/*
			 * Process time integrator first, because it's based on the current state U
			 */
			if (timeIntegratorTerm.active)
			{
				double t = i_simulationTime + _dt*tau(m-1);

				/*
				 * If timeIntegratorTerm is active, we do it Strang-splitting like
				 * First the time integration, then based on the solution, do the explicit solver
				 */
				int idx = timeIntegratorTerm.evalStartIdx_TimeIntegrator_PreSweep+m;
				evalTimeStepper(idx, *U, *tmp, t);

				/*
				 * Now, we need to store the update of the time integrator "U1-U0"
				 */
				timeIntegratorTerm.dataContainer_Integral_U1[m]->op_setVectorPlusScalarMulVector(*tmp, -1, *U);

				U->swap(*tmp);
			}

			/*
			 * Now we do the 2nd part of the time integrator
			 */
			if (explicitTerm.active)
			{
				if (m == 0)
					U->op_addScalarMulVector(_dt*deltaTau(m), *tmp2);
				else
					U->op_addScalarMulVector(_dt*deltaTau(m), *explicitTerm.dataContainer_Tendencies_k1[m-1]);
			}

			if (implicitTerm.active)
			{
				double t = i_simulationTime + _dt*tau(m);
				evalTimeStepper(implicitTerm.evalStartIdx_EulerBackward_PreSweep+m, *U, *tmp, t);

				U->swap(*tmp);
			}



			/*
			 * Tendencies are required if
			 * - end update (quadrature is used)
			 * - there are some sweeps (hence quadrature of residual)
			 */
			bool allTendenciesRequired = (postSweep == POST_SWEEP_QUADRATURE) || (nIters > 0);
			allTendenciesRequired = true;


			double t = i_simulationTime + _dt*tau(m);

			/*
			 * We also need tendencies for the explicit evaluations
			 */
			if (explicitTerm.active)
				if ((m < M-1) || allTendenciesRequired)
					evalTimeStepper(explicitTerm.evalStartIdx_Tendencies, *U, *explicitTerm.dataContainer_Tendencies_k1[m], t);


			/*
			 * Update tendency evaluations at \tau_m (for residual computation and next time step)
			 */
			if (allTendenciesRequired)
			{
				if (implicitTerm.active)
					evalTimeStepper(implicitTerm.evalStartIdx_Tendencies, *U, *implicitTerm.dataContainer_Tendencies_k1[m], t);

				if (timeIntegratorTerm.active)
					evalTimeStepper(timeIntegratorTerm.evalStartIdx_Tendencies, *U, *timeIntegratorTerm.dataContainer_Tendencies_k1[m], t);
			}

			if (m == M-1)
				dataContainer_lastUpdateState = U;
		}

		_swapK();
	}


	/*!
	 * Residual reconstruction (Q matrix update)
	 */
private:
	void _sweepMain_Node2Node_AddResidual(
			sweet::Data::GenericContainer::Base *ts_state,
			int m		//! Quadrature node
	) {
		// Estimate residual by adding quadrature terms
		for (size_t j = 0; j < M; j++)
		{
			double a = _dt*sMat(m, j);

			if (timeIntegratorTerm.active)
				ts_state->op_addScalarMulVector(a, *timeIntegratorTerm.dataContainer_Tendencies_k0[j]);

			if (explicitTerm.active)
				ts_state->op_addScalarMulVector(a, *explicitTerm.dataContainer_Tendencies_k0[j]);

			if (implicitTerm.active)
				ts_state->op_addScalarMulVector(a, *implicitTerm.dataContainer_Tendencies_k0[j]);
		}
	}


	/*!
	 * Integrator corrections
	 *
	 * From \tau_{m-1} to {m}
	 */
private:
	void _sweep_Node2Node_Main(
			sweet::Data::GenericContainer::Base *ts_state,	//! State to update
			double i_simulationTime,	//! Simulation time
			int m						//! Quadrature node
	) {
		// Get buffer for state vector where the residual quadrature has been stored to before

		double t = i_simulationTime + _dt*tau(m-1);
		double dTau = _dt*deltaTau(m);

		if (timeIntegratorTerm.active)
		{
			if (m != 0)
			{
				evalTimeStepper(timeIntegratorTerm.evalStartIdx_TimeIntegrator_MainSweep+m, *ts_state, *tmp, t);

				/*
				 * Now, we need to store the update of the time integrator "U1-U0"
				 */
				timeIntegratorTerm.dataContainer_Integral_U1[m]->op_setVectorPlusScalarMulVector(*tmp, -1, *ts_state);
				ts_state->swap(*tmp);

				ts_state->op_subVector(*timeIntegratorTerm.dataContainer_Integral_U0[m]);
			}
		}

		/*
		 * Add residual
		 */
		_sweepMain_Node2Node_AddResidual(ts_state, m);

		if (explicitTerm.active)
		{
			/*
			 * explicit update is f(u^{k+1}_0) - f(u^{k}_0), hence cancels out for m=0
			 */
			if (m != 0)
			{
				// Last implicit term from iteration k
				ts_state->op_addScalarMulVector(dTau, *explicitTerm.dataContainer_Tendencies_k1[m-1]);
				ts_state->op_addScalarMulVector(-dTau, *explicitTerm.dataContainer_Tendencies_k0[m-1]);
			}
		}

		if (implicitTerm.active)
		{
			// Tendencies from "m"
			ts_state->op_addScalarMulVector(-dTau, *implicitTerm.dataContainer_Tendencies_k0[m]);

			// Backward Euler, evaluated at the next point in time!
			double t = i_simulationTime + _dt*tau(m);
			evalTimeStepper(implicitTerm.evalStartIdx_EulerBackward_MainSweep+m, *ts_state, *tmp, t);
			ts_state->swap(*tmp);
		}
	}


private:
	void _sweepMain_Tendencies(
			sweet::Data::GenericContainer::Base *ts_state,
			double i_simulationTime,
			int m			//! Quadrature node
	) {
		double t = i_simulationTime + _dt*tau(m);

		if (implicitTerm.active)
			evalTimeStepper(implicitTerm.evalStartIdx_Tendencies, *ts_state, *implicitTerm.dataContainer_Tendencies_k1[m], t);

		if (explicitTerm.active)
			evalTimeStepper(explicitTerm.evalStartIdx_Tendencies, *ts_state, *explicitTerm.dataContainer_Tendencies_k1[m], t);

		if (timeIntegratorTerm.active)
			evalTimeStepper(timeIntegratorTerm.evalStartIdx_Tendencies, *ts_state, *timeIntegratorTerm.dataContainer_Tendencies_k1[m], t);
	}


	/*!
	 * Main SDC sweeps
	 */
private:
	void _sweep_Node2Node_Main(
			const sweet::Data::GenericContainer::Base& i_U,
			double i_simulationTime,
			int k
	) {
		// Initialize state with U0
		dataContainer_intermediateState->op_setVector(i_U);

		/*
		 * Iterate over nodes with node-to-node scheme
		 */
		for (size_t m = 0; m < M; m++)
		{
			_sweep_Node2Node_Main(dataContainer_intermediateState, i_simulationTime, m);

			//if (useEndUpdate && (k < nIters-1))
			_sweepMain_Tendencies(dataContainer_intermediateState, i_simulationTime, m);
		}

		_swapK();
	}


private:
	void _sweep_Post(
			const sweet::Data::GenericContainer::Base& i_U,
			sweet::Data::GenericContainer::Base& o_U,
			double i_timeStamp
	) {
		if (postSweep == POST_SWEEP_LASTNODE)
		{
			o_U.op_setVector(*dataContainer_lastUpdateState);
			return;
		}

		/*
		 * Use quadrature for end update
		 */

		/*
		 * Use quadrature
		 */
		const Vec& w = weights;

		// Get buffer for state vector
		sweet::Data::GenericContainer::Base *ts_state = dataContainer_intermediateState;


		// Compute collocation update

		// Initialize with U0
		ts_state->op_setVector(i_U);

		for (size_t j = 0; j < M; j++) {
			double a = _dt*w(j);

			if (implicitTerm.active)
				ts_state->op_addScalarMulVector(a, *implicitTerm.dataContainer_Tendencies_k0[j]);

			if (explicitTerm.active)
				ts_state->op_addScalarMulVector(a, *explicitTerm.dataContainer_Tendencies_k0[j]);

			if (timeIntegratorTerm.active)
				ts_state->op_addScalarMulVector(a, *timeIntegratorTerm.dataContainer_Tendencies_k0[j]);
		}

		o_U.op_setVector(*ts_state);
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
			_sweep_Node2Node_Main(i_U, i_simulationTime, k);
		}

		_sweep_Post(i_U, o_U, i_simulationTime);

		return true;
	}


	void printMule(const std::string &i_prefixStr = "")
	{
		std::cout << i_prefixStr << "SDC_Classic:" << std::endl;
		printSDCInformation("[MULE] sdc.");
	}

	void print(const std::string &i_prefixStr = "")
	{
		std::string newPrefix = i_prefixStr + "  ";
		std::cout << newPrefix << "SDC_Classic(" << std::endl;
		std::cout << newPrefix << "  file: " << fileNameWithSDCCoefficients << std::endl;
		printSDCInformation(newPrefix);
	}
};

}}}

#endif
