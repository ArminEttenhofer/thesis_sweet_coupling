/*
 * SDC_IMEX.cpp
 *
 *      Authors: Thibaut LUNET <thibaut.lunet@tuhh.de>, Martin SCHREIBER <schreiberx@gmail.com>
 *
 * Changelog:
 * 	2023-07-21: Martin converted Thibaut's code to TimeTree
 */

#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_SDC_MARTIN_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_SDC_MARTIN_HPP

#include <vector>
#include <string>
#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/InteriorNodes/SDC_Coefficients.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>


namespace sweet {
namespace TimeTree {
namespace InteriorNodes {


class SDC_Martin	:
	public TimeTree_Node_InteriorHelper<SDC_Martin>,
	public SDC_Coefficients
{
private:
	//! File name of SWEETDict with SDC coefficients
	std::string fileNameWithSDCCoefficients;

	/*!
	 * Implicit term
	 */
	class ImplicitTerm {
	public:
		bool active;

		//! temporary node which can be used during setup
		std::shared_ptr<TimeTree_Node_Base> tmp_TimeTreeNode;

		int evalStartIdx_Tendencies;

		int evalStartIdx_EulerBackward_QMatrix0;
		int evalStartIdx_EulerBackward_QMatrixI;

		//! Storage space for linear tendencies
		sweet::Data::GenericContainer::Base** dataContainer_Tendencies_k0;
		sweet::Data::GenericContainer::Base** dataContainer_Tendencies_k1;

		ImplicitTerm()	:
			active(false),
			evalStartIdx_Tendencies(-1),
			evalStartIdx_EulerBackward_QMatrix0(-1),
			evalStartIdx_EulerBackward_QMatrixI(-1),
			dataContainer_Tendencies_k0(nullptr),
			dataContainer_Tendencies_k1(nullptr)
		{}
	};
	ImplicitTerm implicitTerm;


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
	 * Time integrator term
	 */
	class TimeIntegratorTerm {
	public:
		bool active;

		//! temporary node which can be used during setup
		std::shared_ptr<TimeTree_Node_Base> tmp_TimeTreeNode;
		std::shared_ptr<TimeTree_Node_Base> tmp_TimeTreeNodeTendencies;

		int evalStartIdx_Tendencies;

		int evalStartIdx_TimeIntegrator_QMatrix0;
		int evalStartIdx_TimeIntegrator_QMatrixI;

		//! Storage space for linear tendencies
		sweet::Data::GenericContainer::Base** dataContainer_Tendencies_k0;
		sweet::Data::GenericContainer::Base** dataContainer_Tendencies_k1;

		//! Storage space for linear tendencies
		sweet::Data::GenericContainer::Base** dataContainer_Integration_U0;
		sweet::Data::GenericContainer::Base** dataContainer_Integration_U1;

		TimeIntegratorTerm()	:
			active(false),
			evalStartIdx_Tendencies(-1),
			evalStartIdx_TimeIntegrator_QMatrix0(-1),
			evalStartIdx_TimeIntegrator_QMatrixI(-1),
			dataContainer_Tendencies_k0(nullptr),
			dataContainer_Tendencies_k1(nullptr)
		{}
	};
	TimeIntegratorTerm timeIntegratorTerm;


	//! Intermediate states during SDC
	sweet::Data::GenericContainer::Base** dataContainer_intermediateStates;

public:
	SDC_Martin()	:
		dataContainer_intermediateStates(nullptr)
	{
		setEvalAvailable(EVAL_INTEGRATION);
	}

	SDC_Martin(
			const SDC_Martin &i_src
	)	:
		TimeTree_Node_InteriorHelper<SDC_Martin>(i_src),
		SDC_Coefficients(i_src),
		dataContainer_intermediateStates(nullptr)
	{
		setEvalAvailable(EVAL_INTEGRATION);
	}

	~SDC_Martin()
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
		retval.push_back("SDC_MARTIN");
		retval.push_back("SDCMARTIN");
		retval.push_back("SDC_Martin");
		retval.push_back("SDCMartin");
		return retval;
	}

	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override
	{
		o_ostream << i_prefix << "InteriorNode: 'SDCExplicit':" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Usage: SDCExplicit(DeTermImplicit,DeTermExplicit,[parameter1],...)" << std::endl;
		o_ostream << i_prefix << "           Compute IMEX SDC time integration." << std::endl;
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

		runParallelDiagonalQDelta0Matrix = shackSDC->runParallel_DiagonalQDelta0Matrix;

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

				if (a->key == "residualThreshold" || a->key == "r")
				{
					a->getValue(residualStopThreshold);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);

					residualStopActive = true;
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
		if (runParallelDiagonalQDelta0Matrix && !diagonal)
			return error.set("Error! Diagonal implementation not activated, but parallelization requested!");

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

		if (implicitTerm.active)
		{
			// One for explicit evaluation
			if (runParallelDiagonalQDelta0Matrix)
				counterTimeTreeNodes += M;
			else
				counterTimeTreeNodes += 1;

			/*
			 * TODO: avoid this if COPY is used
			 */
			// For first estimate and also sweeps
			counterTimeTreeNodes += 2*M;
		}

		if (explicitTerm.active)
		{
			if (runParallelDiagonalQDelta0Matrix)
				counterTimeTreeNodes += M;
			else
				counterTimeTreeNodes += 1;
		}

		if (timeIntegratorTerm.active)
		{
			// One for explicit evaluation
			if (runParallelDiagonalQDelta0Matrix)
				counterTimeTreeNodes += M;
			else
				counterTimeTreeNodes += 1;

			/*
			 * TODO: avoid this if COPY is used
			 */
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


		if (implicitTerm.active)
		{
			implicitTerm.evalStartIdx_Tendencies = counterTimeTreeNodesCheck;

			/*
			 * Tendency evaluations
			 */
			if (runParallelDiagonalQDelta0Matrix)
			{
				for (int i = 0; i < M; i++)
				{
					int idx = implicitTerm.evalStartIdx_Tendencies + i;
					_timeTreeNodes[idx] = implicitTerm.tmp_TimeTreeNode->getInstanceCopy();
					_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
				}
				counterTimeTreeNodesCheck += M;
			}
			else
			{
				int idx = implicitTerm.evalStartIdx_Tendencies;
				_timeTreeNodes[idx] = implicitTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);

				counterTimeTreeNodesCheck += 1;
			}


			/*
			 * Backward Euler time steps for Q0
			 */
			implicitTerm.evalStartIdx_EulerBackward_QMatrix0 = counterTimeTreeNodesCheck;
			for (int i = 0; i < M; i++)
			{
				int idx = implicitTerm.evalStartIdx_EulerBackward_QMatrix0 + i;
				_timeTreeNodes[idx] = implicitTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EULER_BACKWARD, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck += M;

			/*
			 * Time integrator
			 */
			implicitTerm.evalStartIdx_EulerBackward_QMatrixI = counterTimeTreeNodesCheck;
			for (int i = 0; i < M; i++)
			{
				int idx = implicitTerm.evalStartIdx_EulerBackward_QMatrixI + i;
				_timeTreeNodes[idx] = implicitTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EULER_BACKWARD, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck += M;

			implicitTerm.tmp_TimeTreeNode.reset();
		}

		if (explicitTerm.active)
		{
			explicitTerm.evalStartIdx_Tendencies = counterTimeTreeNodesCheck;

			if (runParallelDiagonalQDelta0Matrix)
			{
				for (int i = 0; i < M; i++)
				{
					int idx = explicitTerm.evalStartIdx_Tendencies + i;
					_timeTreeNodes[idx] = explicitTerm.tmp_TimeTreeNode->getInstanceCopy();
					_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
				}
				counterTimeTreeNodesCheck += M;
			}
			else
			{
				int idx = explicitTerm.evalStartIdx_Tendencies;
				_timeTreeNodes[idx] = explicitTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
				counterTimeTreeNodesCheck += 1;
			}

			explicitTerm.tmp_TimeTreeNode.reset();
		}


		if (timeIntegratorTerm.active)
		{
			timeIntegratorTerm.evalStartIdx_Tendencies = counterTimeTreeNodesCheck;

			/*
			 * Tendency evaluations
			 */
			if (runParallelDiagonalQDelta0Matrix)
			{
				for (int i = 0; i < M; i++)
				{
					int idx = timeIntegratorTerm.evalStartIdx_Tendencies + i;
					_timeTreeNodes[idx] = timeIntegratorTerm.tmp_TimeTreeNodeTendencies->getInstanceCopy();
					_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
				}
				counterTimeTreeNodesCheck += M;
			}
			else
			{
				int idx = timeIntegratorTerm.evalStartIdx_Tendencies;
				_timeTreeNodes[idx] = timeIntegratorTerm.tmp_TimeTreeNodeTendencies->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);

				counterTimeTreeNodesCheck += 1;
			}


			/*
			 * Integration for setup
			 */
			timeIntegratorTerm.evalStartIdx_TimeIntegrator_QMatrix0 = counterTimeTreeNodesCheck;
			for (int i = 0; i < M; i++)
			{
				int idx = timeIntegratorTerm.evalStartIdx_TimeIntegrator_QMatrix0 + i;
				_timeTreeNodes[idx] = timeIntegratorTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_INTEGRATION, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck += M;

			/*
			 * Integration for sweeps
			 */
			timeIntegratorTerm.evalStartIdx_TimeIntegrator_QMatrixI = counterTimeTreeNodesCheck;
			for (int i = 0; i < M; i++)
			{
				int idx = timeIntegratorTerm.evalStartIdx_TimeIntegrator_QMatrixI + i;
				_timeTreeNodes[idx] = timeIntegratorTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_INTEGRATION, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck += M;

			timeIntegratorTerm.tmp_TimeTreeNode.reset();
		}


		/*********************************************************************
		 * Setup data containers
		 *********************************************************************/

		/*
		 * We need data containers for
		 *
		 *  - 2 iterations (k0 and k1)
		 *  - quadrature nodes M
		 *  - linear and explicit terms
		 */

		SWEET_ASSERT(M >= 1);

		int counterDataContainer = 0;
		if (implicitTerm.active)
			counterDataContainer += 2*M;

		if (explicitTerm.active)
			counterDataContainer += 2*M;

		if (timeIntegratorTerm.active)
		{
			counterDataContainer += 2*M;
			counterDataContainer += 2*M;
		}

		// For common data containers required for parallel executions
		if (runParallelDiagonalQDelta0Matrix)
			counterDataContainer += M;
		else
			counterDataContainer += 1;

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
		if (implicitTerm.active)
		{
			implicitTerm.dataContainer_Tendencies_k0 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			implicitTerm.dataContainer_Tendencies_k1 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;
		}

		if (explicitTerm.active)
		{
			explicitTerm.dataContainer_Tendencies_k0 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			explicitTerm.dataContainer_Tendencies_k1 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;
		}


		if (timeIntegratorTerm.active)
		{
			timeIntegratorTerm.dataContainer_Tendencies_k0 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			timeIntegratorTerm.dataContainer_Tendencies_k1 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			timeIntegratorTerm.dataContainer_Integration_U0 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			timeIntegratorTerm.dataContainer_Integration_U1 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;
		}

		if (runParallelDiagonalQDelta0Matrix)
		{
			/*
			 * For parallel evaluation we need independent data containers
			 * which can be accessed in parallel
			 */
			dataContainer_intermediateStates = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;
		}
		else
		{
			/*
			 * Without parallel evaluation we just need one container
			 */
			dataContainer_intermediateStates = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += 1;
		}

		SWEET_ASSERT(counterDataContainer == counterDataContainerCheck);

		return true;

	}


	void clear() override
	{
		TimeTree_Node_InteriorHelper::clear();
	}


	std::shared_ptr<TimeTree_Node_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<TimeTree_Node_Base>(new SDC_Martin(*this));
	}


	bool setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		/*
		 * Set time step size for linear terms
		 */
		if (implicitTerm.active)
		{
			for (int i = 0; i < M; i++)
			{
				double idx = implicitTerm.evalStartIdx_EulerBackward_QMatrix0+i;
				double dt = _dt*qMatDelta0(i, i);

				_timeTreeNodes[idx]->setTimeStepSize(dt);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}

			for (int i = 0; i < M; i++)
			{
				double idx = implicitTerm.evalStartIdx_EulerBackward_QMatrixI+i;
				double dt = _dt*qMatDeltaI(i, i);

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
		 * Set time step size for time integrator terms
		 */
		if (timeIntegratorTerm.active)
		{
			for (int i = 0; i < M; i++)
			{
				double idx = timeIntegratorTerm.evalStartIdx_TimeIntegrator_QMatrix0+i;
				double dt = _dt*qMatDelta0(i, i);

				_timeTreeNodes[idx]->setTimeStepSize(dt);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}

			for (int i = 0; i < M; i++)
			{
				double idx = timeIntegratorTerm.evalStartIdx_TimeIntegrator_QMatrixI+i;
				double dt = _dt*qMatDeltaI(i, i);

				_timeTreeNodes[idx]->setTimeStepSize(dt);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
		}

		return true;
	}


private:
	void _sweepPre(
			const sweet::Data::GenericContainer::Base& i_U,
			sweet::Data::GenericContainer::Base& o_U,
			double i_simulationTime
	)
	{
		if (initialSweepType == INITIAL_SWEEP_TYPE_QDELTA)
		{
			// Loop on nodes (can be parallelized if diagonal)
	#if SWEET_PARALLEL_SDC_OMP_MODEL
			#pragma omp parallel \
				num_threads(M) \
				if(runParallelDiagonalQDelta0Matrix)\
				default(none) \
				shared(i_U,o_U,i_simulationTime)
			#pragma omp for schedule(static,1)
	#endif
			for (size_t i = 0; i < M; i++) {

				sweet::Data::GenericContainer::Base *ts_state;
				if (runParallelDiagonalQDelta0Matrix)
					ts_state = dataContainer_intermediateStates[i];
				else
					ts_state = dataContainer_intermediateStates[0];

				// Initialize with U0
				ts_state->op_setVector(i_U);

				if (!diagonal)
				{
					// Add non-linear and linear terms from iteration k (already computed)
					for (size_t j = 0; j < i; j++) {
						if (explicitTerm.active)
							ts_state->op_addScalarMulVector(_dt*qMatDeltaE(i, j), *explicitTerm.dataContainer_Tendencies_k0[j]);
						if (implicitTerm.active)
							ts_state->op_addScalarMulVector(_dt*qMatDelta0(i, j), *implicitTerm.dataContainer_Tendencies_k0[j]);
					}
				}

				double t = i_simulationTime + _dt*tauNodes(i);

				if (implicitTerm.active)
				{
					// We abuse the buffer used by eval_explicit after this
					sweet::Data::GenericContainer::Base *tmp = implicitTerm.dataContainer_Tendencies_k0[i];

					// Step 1) **Integration**: Implicit solve with q0
					evalTimeStepper(implicitTerm.evalStartIdx_EulerBackward_QMatrix0+i, *ts_state, *tmp, t);
					tmp->swap(*ts_state);

					// Step 2) **Tendencies**: Evaluate and store linear term for k
					int idx = implicitTerm.evalStartIdx_Tendencies;
					if (runParallelDiagonalQDelta0Matrix)
						idx += i;
					evalTimeStepper(idx, *ts_state, *implicitTerm.dataContainer_Tendencies_k0[i], t);
				}

				if (explicitTerm.active)
				{
					// Step 1) **Tendencies**: Evaluate and store linear term for k
					int idx = explicitTerm.evalStartIdx_Tendencies;
					if (runParallelDiagonalQDelta0Matrix)
						idx += i;
					evalTimeStepper(idx, *ts_state, *explicitTerm.dataContainer_Tendencies_k0[i], t);
				}

				if (timeIntegratorTerm.active)
				{
					// We abuse the buffer used by eval_explicit after this
					sweet::Data::GenericContainer::Base *tmp = timeIntegratorTerm.dataContainer_Tendencies_k0[i];

					// Step 1) **Integration**: Integration with q0
					evalTimeStepper(timeIntegratorTerm.evalStartIdx_TimeIntegrator_QMatrix0+i, *ts_state, *tmp, t);
					tmp->swap(*ts_state);

					// Step 2) **Tendencies**: Evaluate and store linear term for k
					int idx = timeIntegratorTerm.evalStartIdx_Tendencies;
					if (runParallelDiagonalQDelta0Matrix)
						idx += i;
					evalTimeStepper(idx, *ts_state, *timeIntegratorTerm.dataContainer_Tendencies_k0[i], t);

					// **Backup U_0 Integration**
					timeIntegratorTerm.dataContainer_Integration_U0[i]->op_setVector(*ts_state);
				}

				if (!useEndUpdate && (i == M-1) && (nIters == 0))
					o_U.op_setVector(*ts_state);
			}
		}
		else if (initialSweepType == INITIAL_SWEEP_TYPE_COPY)
		{
			// Evaluate linear and non-linear with initial solution
			if (implicitTerm.active)
				evalTimeStepper(implicitTerm.evalStartIdx_Tendencies, i_U, *implicitTerm.dataContainer_Tendencies_k0[0], i_simulationTime);

			if (explicitTerm.active)
				evalTimeStepper(explicitTerm.evalStartIdx_Tendencies, i_U, *explicitTerm.dataContainer_Tendencies_k0[0], i_simulationTime);

			if (timeIntegratorTerm.active)
				evalTimeStepper(timeIntegratorTerm.evalStartIdx_Tendencies, i_U, *timeIntegratorTerm.dataContainer_Tendencies_k0[0], i_simulationTime);

			// Simply copy to each node as initial guess
	#if SWEET_PARALLEL_SDC_OMP_MODEL
			#pragma omp parallel \
				num_threads(M) \
				if(runParallelDiagonalQDelta0Matrix) \
				default(none)	\
				shared(i_U,o_U,i_simulationTime)
			#pragma omp for schedule(static,1)
	#endif
			for (size_t i = 0; i < M; i++)
			{
				// Include first node just to simulate same parallelization as for parallel updates
				if (i == 0)
					continue;

				if (implicitTerm.active)
					implicitTerm.dataContainer_Tendencies_k0[i]->op_setVector(*implicitTerm.dataContainer_Tendencies_k0[0]);

				if (explicitTerm.active)
					explicitTerm.dataContainer_Tendencies_k0[i]->op_setVector(*explicitTerm.dataContainer_Tendencies_k0[0]);

				if (timeIntegratorTerm.active)
				{
					timeIntegratorTerm.dataContainer_Tendencies_k0[i]->op_setVector(*timeIntegratorTerm.dataContainer_Tendencies_k0[0]);
					timeIntegratorTerm.dataContainer_Integration_U0[i]->op_setVector(i_U);
				}
			}

			if (!useEndUpdate && nIters == 0)
				SWEETErrorFatal("This doesn't make any sense");
		}
	}


	/*!
	 * Residual reconstruction (Q matrix update)
	 */
private:
	void _sweepMain_Residual(
			const sweet::Data::GenericContainer::Base& i_U,
			sweet::Data::GenericContainer::Base& o_U,
			double i_simulationTime,
			size_t k,	//! iteration number
			int i,		//! Quadrature node
			bool &o_residualBasedTermination	//!< Set to true if residual termination happened
	) {
		// Get buffer for state vector
		sweet::Data::GenericContainer::Base *ts_state;
		if (runParallelDiagonalQDelta0Matrix)
			ts_state = dataContainer_intermediateStates[i];
		else
			ts_state = dataContainer_intermediateStates[0];

		// Initialize state with U0
		ts_state->op_setVector(i_U);

		// Estimate residual by adding quadrature terms
		for (size_t j = 0; j < M; j++)
		{
			double a = _dt*qMat(i, j);

			if (explicitTerm.active)
				ts_state->op_addScalarMulVector(a, *explicitTerm.dataContainer_Tendencies_k0[j]);

			if (implicitTerm.active)
				ts_state->op_addScalarMulVector(a, *implicitTerm.dataContainer_Tendencies_k0[j]);

			if (timeIntegratorTerm.active)
				ts_state->op_addScalarMulVector(a, *timeIntegratorTerm.dataContainer_Tendencies_k0[j]);
		}

		/*
		 * Check for LMax norm on residual for stopping criteria
		 */
		if (residualStopActive && k >= 1)
		{
			/*
			 * We only do this for the last node.
			 * TODO: We could also do this on all nodes
			 * in parallel which wouldn't matter.
			 */
			if (i == M-1)
			{
				double residual = ts_state->reduceLMax();

				if (residual < residualStopThreshold)
					o_residualBasedTermination = true;
			}
		}
	}

	/*!
	 * Integrator corrections (QDelta matrix update)
	 */
private:
	void _sweepMain_Integrator(
			const sweet::Data::GenericContainer::Base& i_U,
			sweet::Data::GenericContainer::Base& o_U,
			double i_simulationTime,
			size_t k,	//! iteration number
			int i,		//! Quadrature node
			bool &o_residualBasedTermination	//!< Set to true if residual termination happened
	) {
		// Get buffer for state vector
		sweet::Data::GenericContainer::Base *ts_state;
		if (runParallelDiagonalQDelta0Matrix)
			ts_state = dataContainer_intermediateStates[i];
		else
			ts_state = dataContainer_intermediateStates[0];

		if (!diagonal) {

			// Add non-linear and linear terms from iteration k+1
			for (size_t j = 0; j < i; j++)
			{
				if (explicitTerm.active)
					ts_state->op_addScalarMulVector(_dt*qMatDeltaE(i, j), *explicitTerm.dataContainer_Tendencies_k1[j]);

				if (implicitTerm.active)
					ts_state->op_addScalarMulVector(_dt*qMatDeltaI(i, j), *implicitTerm.dataContainer_Tendencies_k1[j]);

				if (timeIntegratorTerm.active)
					ts_state->op_addScalarMulVector(1, *timeIntegratorTerm.dataContainer_Integration_U1[j]);
			}

			// Subtract non-linear and linear terms from iteration k
			for (size_t j = 0; j < i; j++)
			{
				if (explicitTerm.active)
					ts_state->op_addScalarMulVector(-_dt*qMatDeltaE(i, j), *explicitTerm.dataContainer_Tendencies_k0[j]);

				if (implicitTerm.active)
					ts_state->op_addScalarMulVector(-_dt*qMatDeltaI(i, j), *implicitTerm.dataContainer_Tendencies_k0[j]);

				if (timeIntegratorTerm.active)
					ts_state->op_addScalarMulVector(-1, *timeIntegratorTerm.dataContainer_Integration_U0[j]);
			}
		}

		double t = i_simulationTime + _dt*tauNodes(i);

		if (timeIntegratorTerm.active)
		{
			sweet::Data::GenericContainer::Base *tmp = timeIntegratorTerm.dataContainer_Tendencies_k1[i];

			// **Time Integration**: Compute R(U_1)
			evalTimeStepper(timeIntegratorTerm.evalStartIdx_TimeIntegrator_QMatrixI+i, *ts_state, *tmp, t);
			ts_state->swap(*tmp);

			timeIntegratorTerm.dataContainer_Integration_U1[i]->op_setVector(*ts_state);

			/*
			 * TODO: something's not working here.
			 */
			//ts_state->op_addScalarMulVector(-1, *timeIntegratorTerm.dataContainer_Integration_U0[i]);
		}

		if (implicitTerm.active)
		{
			// Last linear term from iteration k
			ts_state->op_addScalarMulVector(-_dt*qMatDeltaI(i, i), *implicitTerm.dataContainer_Tendencies_k0[i]);
		}

		if (implicitTerm.active)
		{
			// Implicit solve
			// we abuse the buffer used by eval_explicit after this
			sweet::Data::GenericContainer::Base *tmp = implicitTerm.dataContainer_Tendencies_k1[i];
			evalTimeStepper(implicitTerm.evalStartIdx_EulerBackward_QMatrixI+i, *ts_state, *tmp, t);
			ts_state->swap(*tmp);
		}


		if (!useEndUpdate && residualStopActive)
		{
			/*
			 * If residual stopping criteria is active, we always
			 * need to backup this state
			 */
			o_U.op_setVector(*ts_state);
			return;
		}

		if (!useEndUpdate && (i == M-1) && (k == nIters-1))
		{
			o_U.op_setVector(*ts_state);
		}
		else
		{
			// Evaluate and store linear term for k+1
			if (implicitTerm.active)
			{
				int idx = implicitTerm.evalStartIdx_Tendencies;
				if (runParallelDiagonalQDelta0Matrix)
					idx += i;
				evalTimeStepper(idx, *ts_state, *implicitTerm.dataContainer_Tendencies_k1[i], t);
			}

			if (explicitTerm.active)
			{
				int idx = explicitTerm.evalStartIdx_Tendencies;
				if (runParallelDiagonalQDelta0Matrix)
					idx += i;
				evalTimeStepper(idx, *ts_state, *explicitTerm.dataContainer_Tendencies_k1[i], t);
			}

			if (timeIntegratorTerm.active)
			{
				int idx = timeIntegratorTerm.evalStartIdx_Tendencies;
				if (runParallelDiagonalQDelta0Matrix)
					idx += i;
				evalTimeStepper(idx, *ts_state, *timeIntegratorTerm.dataContainer_Tendencies_k1[i], t);
			}

		}
	}


	/*!
	 * Main SDC sweeps
	 */
private:
	void _sweepMain(
			const sweet::Data::GenericContainer::Base& i_U,
			sweet::Data::GenericContainer::Base& o_U,
			double i_simulationTime,
			int k,	//! iteration number
			bool &o_residualBasedTermination	//!< Set to true if residual termination happened
	) {

		if (runParallelDiagonalQDelta0Matrix) {
			if (residualStopActive) {

				#if SWEET_PARALLEL_SDC_OMP_MODEL
					#pragma omp parallel \
						num_threads(M) \
						if(runParallelDiagonalQDelta0Matrix) \
						default(none) \
						shared(i_U,o_U,i_simulationTime,k,o_residualBasedTermination)
					#pragma omp for schedule(static,1)
				#endif
				for (size_t i = 0; i < M; i++)
					_sweepMain_Residual(i_U, o_U, i_simulationTime, k, i, o_residualBasedTermination);

				if (!o_residualBasedTermination)
				{
					#if SWEET_PARALLEL_SDC_OMP_MODEL
						#pragma omp parallel \
							num_threads(M) \
							if(runParallelDiagonalQDelta0Matrix) \
							default(none) \
							shared(i_U,o_U,i_simulationTime,k,o_residualBasedTermination)
						#pragma omp for schedule(static,1)
					#endif
					for (size_t i = 0; i < M; i++)
						_sweepMain_Integrator(i_U, o_U, i_simulationTime, k, i, o_residualBasedTermination);
				}
			}
			else
			{
				#if SWEET_PARALLEL_SDC_OMP_MODEL
					#pragma omp parallel \
						num_threads(M) \
						if(runParallelDiagonalQDelta0Matrix) \
						default(none) \
						shared(i_U,o_U,i_simulationTime,k,o_residualBasedTermination)
					#pragma omp for schedule(static,1)
				#endif
				for (size_t i = 0; i < M; i++) {
					_sweepMain_Residual(i_U, o_U, i_simulationTime, k, i, o_residualBasedTermination);
					_sweepMain_Integrator(i_U, o_U, i_simulationTime, k, i, o_residualBasedTermination);
				}
			}
		}
		else
		{
			for (size_t i = 0; i < M; i++) {
				_sweepMain_Residual(i_U, o_U, i_simulationTime, k, i, o_residualBasedTermination);
				_sweepMain_Integrator(i_U, o_U, i_simulationTime, k, i, o_residualBasedTermination);
			}
		}


		/*
		 * Swap k+1 and k values for next iteration (or end-point update)
		 *
		 * (Just swap the pointers into the vector structure)
		 */
		if (implicitTerm.active)
			std::swap(implicitTerm.dataContainer_Tendencies_k0, implicitTerm.dataContainer_Tendencies_k1);

		if (explicitTerm.active)
			std::swap(explicitTerm.dataContainer_Tendencies_k0, explicitTerm.dataContainer_Tendencies_k1);

		if (timeIntegratorTerm.active)
		{
			std::swap(timeIntegratorTerm.dataContainer_Tendencies_k0, timeIntegratorTerm.dataContainer_Tendencies_k1);
			std::swap(timeIntegratorTerm.dataContainer_Integration_U0, timeIntegratorTerm.dataContainer_Integration_U1);
		}
	}


private:
	void _sweepPost(
			const sweet::Data::GenericContainer::Base& i_U,
			sweet::Data::GenericContainer::Base& o_U,
			double i_timeStamp
	) {
		// Compute end point for last sweep
		if (useEndUpdate) {
			/*
			 * Use quadrature
			 */
			const Vec& w = endQuadWeights;

			// Get buffer for state vector
			sweet::Data::GenericContainer::Base *ts_state;
			ts_state = dataContainer_intermediateStates[0];


			// Compute collocation update

			// Initialize with U0
			ts_state->op_setVector(i_U);

			for (size_t j = 0; j < M; j++) {
				double a = _dt*w(j);

				if (explicitTerm.active)
					ts_state->op_addScalarMulVector(a, *explicitTerm.dataContainer_Tendencies_k0[j]);
					
				if (implicitTerm.active)
					ts_state->op_addScalarMulVector(a, *implicitTerm.dataContainer_Tendencies_k0[j]);

				if (timeIntegratorTerm.active)
					ts_state->op_addScalarMulVector(a, *timeIntegratorTerm.dataContainer_Tendencies_k0[j]);
			}

			o_U.op_setVector(*ts_state);
		}
	}


	bool _eval_integration(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)	override
	{
		// initialize nodes values and state
		_sweepPre(i_U, o_U, i_simulationTime);

		// perform sweeps

		bool residualBasedTermination = false;
		for (size_t k = 0; k < nIters && !residualBasedTermination; k++) {
			_sweepMain(i_U, o_U, i_simulationTime, k, residualBasedTermination);
		}

		if (residualBasedTermination)
			residualStopCounter++;

		_sweepPost(i_U, o_U, i_simulationTime);

		return true;
	}


	void printMule(const std::string &i_prefixStr = "")
	{
		std::cout << i_prefixStr << "SDC_IMEX:" << std::endl;
		printSDCInformationMule("sdc.");
	}

	void print(const std::string &i_prefixStr = "")
	{
		std::string newPrefix = i_prefixStr + "  ";
		std::cout << i_prefixStr << "SDC_IMEX(" << std::endl;
		std::cout << i_prefixStr << "  file: " << fileNameWithSDCCoefficients << std::endl;
		printSDCInformation(i_prefixStr);
	}
};

}}}

#endif
