/*
 * SDC_FP.hpp
 * 
 * 		SDC implementation using Fixed-Point formulation
 *
 *      Authors: Thibaut LUNET <thibaut.lunet@tuhh.de>, Martin SCHREIBER <schreiberx@gmail.com>
 *
 * Changelog:
 * 	2023-07-21: Martin converted Thibaut's code to TimeTree
 *  2024-02-01: Naming changed from SDC_Generic to SDC_FP
 */

#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_SDC_FP_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_SDC_FP_HPP
// #define SDC_DEBUG

#include <vector>
#include <string>
#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/InteriorNodes/SDC_Coefficients.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>


namespace sweet {
namespace TimeTree {
namespace InteriorNodes {

template <typename... Types>
void DEBUG(const Types&... args) {
	#ifdef SDC_DEBUG
	std::cout << "[SDC] ";
	((std::cout << args), ...) << std::endl;
	#endif
}

class SDC_FP	:
	public TimeTree_Node_InteriorHelper<SDC_FP>,
	public SDC_Coefficients
{
private:
	//! File name of SWEETDict with SDC coefficients
	std::string fileNameWithSDCCoefficients;

	/*!
	 * Explicit term
	 */
	class ExplicitTerm {
	public:
		bool active;

		///! temporary node which can be used during setup
		std::shared_ptr<TimeTree_Node_Base> tmp_TimeTreeNode;

		// Tree-Node index for the tendency evaluation function
		int tnIdx_eval;

		// Storage space for tendencies
		sweet::Data::GenericContainer::Base** tendencies_k0;
		sweet::Data::GenericContainer::Base** tendencies_k1;
		sweet::Data::GenericContainer::Base* tendency_U0;

		ExplicitTerm()	:
			active(false),
			tnIdx_eval(-1),
			tendencies_k0(nullptr),
			tendencies_k1(nullptr),
			tendency_U0(nullptr)
		{}
	};
	ExplicitTerm explicitTerm;

	/*!
	 * Implicit term
	 */
	class ImplicitTerm: public ExplicitTerm {
	public:
		// Tree-Node index for the Backward Euler solve used during pre sweep
		int tnIdx_preSweepSolve;
		// Tree-Node index for the Backward Euler solve used during sweep
		int tnIdx_sweepSolve;

		ImplicitTerm()	: 
			ExplicitTerm(),
			tnIdx_preSweepSolve(-1),
			tnIdx_sweepSolve(-1)
		{}
	};
	ImplicitTerm implicitTerm;

	//! Intermediate state(s) used during SDC sweep 
	// -- 1 data container required without parallelization
	// -- M data container required with parallelization 
	sweet::Data::GenericContainer::Base** stateContainers;

	// ! End-update state required for post sweep = INTERPOLATION
	sweet::Data::GenericContainer::Base* endUpdate_U1;


public:
	SDC_FP()	:
		stateContainers(nullptr)
	{
		setEvalAvailable(EVAL_INTEGRATION);
	}

	SDC_FP(
			const SDC_FP &i_src
	)	:
		TimeTree_Node_InteriorHelper<SDC_FP>(i_src),
		SDC_Coefficients(i_src),
		stateContainers(nullptr)
	{
		setEvalAvailable(EVAL_INTEGRATION);
	}

	~SDC_FP()
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
		retval.push_back("SDCFP");
		retval.push_back("SDC_FP");
		retval.push_back("SDC-FP");
		return retval;
	}

	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override
	{
		o_ostream << i_prefix << "InteriorNode: 'SDC-FP':" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Usage: SDC-FP(DeTermImplicit,DeTermExplicit,[parameter1],...)" << std::endl;
		o_ostream << i_prefix << "           Compute IMEX SDC time integration using the fixed-point formulation (ZeroToNodes)." << std::endl;
		o_ostream << i_prefix << "           DeTermStiff has to support 'tendencies' and 'backward_euler'" << std::endl;
		o_ostream << i_prefix << "           DeTermNonStiff has to support 'tendencies'" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Parameters:" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "    - file=[str]:" << std::endl;
		o_ostream << i_prefix << "        Filename with SDC coefficients" << std::endl;
		o_ostream << i_prefix << "    - explicit=[str]:" << std::endl;
		o_ostream << i_prefix << "    - e=[str]:" << std::endl;
		o_ostream << i_prefix << "        DETerm to treat explicitly" << std::endl;
		o_ostream << i_prefix << "    - implicit=[str]:" << std::endl;
		o_ostream << i_prefix << "    - i=[str]:" << std::endl;
		o_ostream << i_prefix << "        DETerm to treat implicitly" << std::endl;
		o_ostream << i_prefix << std::endl;

		return true;
	}

	bool _setupArgumentInternals()
	{
		if (fileNameWithSDCCoefficients == "")
			fileNameWithSDCCoefficients = shackSDC->fileName;

		parallel = shackSDC->runParallel_DiagonalQDelta0Matrix;

		if (fileNameWithSDCCoefficients == "")
			fileNameWithSDCCoefficients = "params_SDC.sweet";

		if (!implicitTerm.active && !explicitTerm.active)
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

				// if (a->key == "residualThreshold" || a->key == "r")
				// {
				// 	a->getValue(residualStopThreshold);
				// 	ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);

				// 	residualStopActive = true;
				// 	break;
				// }

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
		if (parallel && !diagonal)
			return error.set("Error! Only diagonal SDC can be run in parallel");

		if (!rightIsNode && postSweep == POST_SWEEP_LASTNODE)
			return error.set("Error! Cannot use post sweep LASTNODE if not rightIsNode");

		if (postSweep == POST_SWEEP_INTERPOLATION && parallel)
			return error.set("Error! Cannot use post sweep INTERPOLATION with parallel");

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
			if (parallel)
				counterTimeTreeNodes += M;
			else
				counterTimeTreeNodes += 1;

			/*
			 * TODO: avoid this if COPY is used or same implicit QDelta for each sweep
			 */
			// For first estimate and also (each) sweeps
			counterTimeTreeNodes += M*(nIters+1);
		}

		if (explicitTerm.active)
		{
			if (parallel)
				counterTimeTreeNodes += M;
			else
				counterTimeTreeNodes += 1;
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
			implicitTerm.tnIdx_eval = counterTimeTreeNodesCheck;

			/*
			 * Tendency evaluations
			 */
			if (parallel)
			{
				for (int i = 0; i < M; i++)
				{
					int idx = implicitTerm.tnIdx_eval + i;
					_timeTreeNodes[idx] = implicitTerm.tmp_TimeTreeNode->getInstanceCopy();
					_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
				}
				counterTimeTreeNodesCheck += M;
			}
			else
			{
				int idx = implicitTerm.tnIdx_eval;
				_timeTreeNodes[idx] = implicitTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);

				counterTimeTreeNodesCheck += 1;
			}


			/*
			 * Backward Euler time steps for Q0
			 */
			implicitTerm.tnIdx_preSweepSolve = counterTimeTreeNodesCheck;
			for (int m = 0; m < M; m++)
			{
				int idx = implicitTerm.tnIdx_preSweepSolve + m;
				_timeTreeNodes[idx] = implicitTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EULER_BACKWARD, &_evalFuns[idx]);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
			counterTimeTreeNodesCheck += M;

			/*
			 * Backward Euler time steps for QI
			 */
			implicitTerm.tnIdx_sweepSolve = counterTimeTreeNodesCheck;
			for (int k = 0; k < nIters; k++)
			{
				for (int m = 0; m < M; m++)
				{
					int idx = implicitTerm.tnIdx_sweepSolve + k*M+m;
					_timeTreeNodes[idx] = implicitTerm.tmp_TimeTreeNode->getInstanceCopy();
					_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EULER_BACKWARD, &_evalFuns[idx]);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
				}
				counterTimeTreeNodesCheck += M;
			}

			implicitTerm.tmp_TimeTreeNode.reset();
		}

		if (explicitTerm.active)
		{
			explicitTerm.tnIdx_eval = counterTimeTreeNodesCheck;

			if (parallel)
			{
				for (int i = 0; i < M; i++)
				{
					int idx = explicitTerm.tnIdx_eval + i;
					_timeTreeNodes[idx] = explicitTerm.tmp_TimeTreeNode->getInstanceCopy();
					_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
				}
				counterTimeTreeNodesCheck += M;
			}
			else
			{
				int idx = explicitTerm.tnIdx_eval;
				_timeTreeNodes[idx] = explicitTerm.tmp_TimeTreeNode->getInstanceCopy();
				_timeTreeNodes[idx]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[idx]);
				counterTimeTreeNodesCheck += 1;
			}

			explicitTerm.tmp_TimeTreeNode.reset();
		}

		SWEET_ASSERT(counterTimeTreeNodesCheck == counterTimeTreeNodes);

		/*********************************************************************
		 * Setup data containers
		 *********************************************************************/

		/*
		 * We need data containers for
		 *
		 * 	- Tendencies for each term (implicit and explicit)
		 *  	- 2 iterations on each quadrature nodes (k0 and k1) : 2*M
		 * 		- storage for the initial tendency : 1
		 *  - State vector used during sweep
		 *  	- M if parallel
		 * 		- 1 if not parallel
		 *  - End update (used only with POST_SWEEP_INTERPOLATION) : 1
		 */

		SWEET_ASSERT(M >= 1);

		int counterDataContainer = 0;
		if (implicitTerm.active)
			counterDataContainer += 2*M + 1;

		if (explicitTerm.active)
			counterDataContainer += 2*M + 1;

		// For common data containers required for parallel executions
		if (parallel)
			counterDataContainer += M;
		else
			counterDataContainer += 1;

		// End update
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
			implicitTerm.tendencies_k0 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			implicitTerm.tendencies_k1 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			implicitTerm.tendency_U0 = _tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += 1;
		}

		if (explicitTerm.active)
		{
			explicitTerm.tendencies_k0 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			explicitTerm.tendencies_k1 = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;

			explicitTerm.tendency_U0 = _tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += 1;
		}


		if (parallel)
		{
			/*
			 * For parallel evaluation we need independent data containers
			 * which can be accessed in parallel
			 */
			stateContainers = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += M;
		}
		else
		{
			/*
			 * Without parallel evaluation we just need one container
			 */
			stateContainers = &_tmpDataContainer[counterDataContainerCheck];
			counterDataContainerCheck += 1;
		}

		endUpdate_U1 = _tmpDataContainer[counterDataContainerCheck];
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
		return std::shared_ptr<TimeTree_Node_Base>(new SDC_FP(*this));
	}


	bool setTimeStepSize(double i_dt)	override
	{
		_dt = i_dt;

		/*
		 * Set time step size for linear terms
		 */
		if (implicitTerm.active)
		{
			for (int m = 0; m < M; m++)
			{
				double idx = implicitTerm.tnIdx_preSweepSolve + m;
				double dt = _dt*qMatDelta0(m, m);

				_timeTreeNodes[idx]->setTimeStepSize(dt);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
			}
		
			for (int k = 0; k < nIters; k++)
			{
				for (int m = 0; m < M; m++)
				{
					double idx = implicitTerm.tnIdx_sweepSolve + k*M+m;
					double dt = _dt*qMatDeltaI(k, m, m);

					_timeTreeNodes[idx]->setTimeStepSize(dt);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[idx]);
				}
			}
		}

		/*
		 * Set time step size for explicit terms
		 */
		if (explicitTerm.active)
		{
			// Nothing to do since it's just about tendencies
			//_timeTreeNodes[explicitTerm.tnIdx_eval]->setTimeStepSize(i_dt);
		}

		return true;
	}


private:

	sweet::Data::GenericContainer::Base*& _getStateContainer(int m) 
	{
		if (parallel)
			return stateContainers[m];
		else
			return stateContainers[0];
	}

	void _sweepPre(
			const sweet::Data::GenericContainer::Base& i_U,
			sweet::Data::GenericContainer::Base& o_U,
			double i_simulationTime
	)
	{
		DEBUG("starting pre sweep (evaluation of U0 tendency first ...)");

		// Evaluate linear and non-linear with initial solution
		// TODO : this should be optimized considering leftIsNode and rightIsNode
		if (implicitTerm.active)
			evalTimeStepper(implicitTerm.tnIdx_eval, i_U, *implicitTerm.tendency_U0, i_simulationTime);
		if (explicitTerm.active)
			evalTimeStepper(explicitTerm.tnIdx_eval, i_U, *explicitTerm.tendency_U0, i_simulationTime);

		if (preSweep == PRE_SWEEP_QDELTA)
		{
			DEBUG("using QDELTA initialization");
			// Loop on nodes (can be parallelized if diagonal)
	#if SWEET_PARALLEL_SDC_OMP_MODEL
			#pragma omp parallel \
				num_threads(M) \
				if(parallel)\
				default(none) \
				shared(i_U,o_U,i_simulationTime)
			#pragma omp for schedule(static,1)
	#endif
			for (size_t m = 0; m < M; m++) {

				DEBUG("update for m=", m);
				// Initialize state with U0
				sweet::Data::GenericContainer::Base *state = _getStateContainer(m);
				state->op_setVector(i_U);

				// Add dtau contributions (U0 tendencies)
				if (explicitTerm.active && dtauE(m) != 0.0) {
					DEBUG("adding dtau contribution for explicit term");
					state->op_addScalarMulVector(_dt*dtauE(m), *explicitTerm.tendency_U0);
				}
				if (implicitTerm.active && dtau0(m) != 0.0) {
					DEBUG("adding dtau contribution for implicit term");
					state->op_addScalarMulVector(_dt*dtau0(m), *implicitTerm.tendency_U0);
				}

				// Add terms tendency from previous nodes
				if (!diagonal)
				{
					DEBUG("adding tendency from previous nodes");
					for (size_t j = 0; j < m; j++) {
						if (explicitTerm.active)
							state->op_addScalarMulVector(_dt*qMatDeltaE(m, j), *explicitTerm.tendencies_k0[j]);
						if (implicitTerm.active)
							state->op_addScalarMulVector(_dt*qMatDelta0(m, j), *implicitTerm.tendencies_k0[j]);
					}
				}

				
				double t = i_simulationTime + _dt*tauNodes(m);

				// Implicit solve + evaluation of implicit tendency using updated state
				if (implicitTerm.active)
				{
					// Step 1) **Integration**: Implicit solve to get node solution
					sweet::Data::GenericContainer::Base *tmp = implicitTerm.tendencies_k0[m];
					// -- need to do this since evalTimeStepper does not overwrite input state 
					// 	  => use implicitTerm.tendencies_k0[m] as buffer
					// 	  -> clever idea of Martin, who may understood Z2N SDC more than I (viciously) thought ... 
					evalTimeStepper(implicitTerm.tnIdx_preSweepSolve + m, *state, *tmp, t);
					tmp->swap(*state);

					// Step 2) **Tendencies**: Evaluate and store implicit term for k
					int idx = implicitTerm.tnIdx_eval + m*((int) parallel);
					evalTimeStepper(idx, *state, *implicitTerm.tendencies_k0[m], t);
				}

				// Evaluation of explicit tendency using current state
				if (explicitTerm.active)
				{
					int idx = explicitTerm.tnIdx_eval + m*((int) parallel);
					evalTimeStepper(idx, *state, *explicitTerm.tendencies_k0[m], t);
				}

				// Prepare end-update solution for post sweep using INTERPOLATION (avoid data storage)
				if (postSweep == POST_SWEEP_INTERPOLATION && (nIters == 0))
				{
					// TODO : adapt for parallel !
					// Initialize with zero vector for first node
					if (m == 0)
						endUpdate_U1->op_setZero();
					// Add node contribution
					endUpdate_U1->op_addScalarMulVector(hCoeffs(m), *state);
				}

			}
		}

		else if (preSweep == PRE_SWEEP_COPY || preSweep == PRE_SWEEP_ZEROS)
		{
			// Initialize state with U0 (last state if parallel)
			sweet::Data::GenericContainer::Base *state = _getStateContainer(M-1);
			if (preSweep == PRE_SWEEP_COPY)
			{
				DEBUG("using COPY initialization");
				state->op_setVector(i_U);
			}
			else if (preSweep == PRE_SWEEP_ZEROS)
			{
				DEBUG("using ZEROS initialization");
				state->op_setZero();
			}
				

			// Copy tendency of U0 accross all nodes
	#if SWEET_PARALLEL_SDC_OMP_MODEL
			#pragma omp parallel \
				num_threads(M) \
				if(parallel) \
				default(none)	\
				shared(i_U,o_U,i_simulationTime)
			#pragma omp for schedule(static,1)
	#endif
			for (size_t m = 0; m < M; m++)
			{
				if (implicitTerm.active)
					implicitTerm.tendencies_k0[m]->op_setVector(*implicitTerm.tendency_U0);

				if (explicitTerm.active)
					explicitTerm.tendencies_k0[m]->op_setVector(*explicitTerm.tendency_U0);
			}
		}

	}

	/*!
	 * Sweep Update on one node
	 *  
	 * 1) Add quadrature terms (Q matrix update)
	 * 2) Add correction terms (QDelta matrix update)
	 * 3) Update node solution (Backward Euler solve)
	 * 4) Evaluate tendencies on the updated node solution
	 */
private:
	void _sweepNodeUpdate(
			const sweet::Data::GenericContainer::Base& i_U,
			sweet::Data::GenericContainer::Base& o_U,
			double i_simulationTime,
			size_t k,	//! iteration number
			int m		//! Quadrature node
	) {
		// Get buffer for state vector
		sweet::Data::GenericContainer::Base *state = _getStateContainer(m);

		// Simulation time on the node (for implicit solve and tendency evaluations)
		double t = i_simulationTime + _dt*tauNodes(m);

		// Initialize state with U0
		state->op_setVector(i_U);

		// 1) Add quadrature terms (Q matrix update)
		for (size_t j = 0; j < M; j++)
		{
			double a = _dt*qMat(m, j);

			if (explicitTerm.active)
				state->op_addScalarMulVector(a, *explicitTerm.tendencies_k0[j]);

			if (implicitTerm.active)
				state->op_addScalarMulVector(a, *implicitTerm.tendencies_k0[j]);
		}

		// 2) Add correction terms (QDelta matrix update)
		if (!diagonal) {

			// Add explicit and implicit terms from iteration k+1
			for (size_t j = 0; j < m; j++)
			{
				if (explicitTerm.active)
					state->op_addScalarMulVector(_dt*qMatDeltaE(m, j), *explicitTerm.tendencies_k1[j]);

				if (implicitTerm.active)
					state->op_addScalarMulVector(_dt*qMatDeltaI(k, m, j), *implicitTerm.tendencies_k1[j]);
			}

			// Subtract explicit and implicit terms from iteration k
			for (size_t j = 0; j < m; j++)
			{
				if (explicitTerm.active)
					state->op_addScalarMulVector(-_dt*qMatDeltaE(m, j), *explicitTerm.tendencies_k0[j]);

				if (implicitTerm.active)
					state->op_addScalarMulVector(-_dt*qMatDeltaI(k, m, j), *implicitTerm.tendencies_k0[j]);
			}
		}

		if (implicitTerm.active)
		{
			// Subtract last implicit term from iteration k
			if (implicitTerm.active)
				state->op_addScalarMulVector(-_dt*qMatDeltaI(k, m, m), *implicitTerm.tendencies_k0[m]);

			// 3) Update node solution (Backward Euler solve)
			sweet::Data::GenericContainer::Base *tmp = implicitTerm.tendencies_k1[m];
			evalTimeStepper(implicitTerm.tnIdx_sweepSolve + k*M+m, *state, *tmp, t);
			state->swap(*tmp);
			// -> need to do this since evalTimeStepper does not update input state
		}

		// Prepare end-update solution for post sweep using INTERPOLATION (avoid data storage)
		if (postSweep == POST_SWEEP_INTERPOLATION && (k == nIters-1))
		{
			// Initialize with zero vector for first node
			// TODO : adapt for parallel !
			if (m == 0) 
				endUpdate_U1->op_setZero();
			// Add node contribution
			endUpdate_U1->op_addScalarMulVector(hCoeffs(m), *state);
		}

		// Don't evaluate tendencies if :
		// - post sweep does not use QUADRATURE
		// - last iteration
		// - last node
		if (postSweep != POST_SWEEP_QUADRATURE && (m == M-1) && (k == nIters-1)) return;

		// 4) Evaluate tendencies on the updated node solution
		if (implicitTerm.active)
		{
			int idx = implicitTerm.tnIdx_eval + m*((int) parallel);
			evalTimeStepper(idx, *state, *implicitTerm.tendencies_k1[m], t);
		}
		if (explicitTerm.active)
		{
			int idx = explicitTerm.tnIdx_eval + m*((int) parallel);
			evalTimeStepper(idx, *state, *explicitTerm.tendencies_k1[m], t);
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
			int k	//! iteration number
	) {

		DEBUG("starting sweep for k=", k);
		if (parallel) {

			{
				#if SWEET_PARALLEL_SDC_OMP_MODEL
					#pragma omp parallel \
						num_threads(M) \
						if(parallel) \
						default(none) \
						shared(i_U,o_U,i_simulationTime,k)
					#pragma omp for schedule(static,1)
				#endif
				for (size_t i = 0; i < M; i++) {
					_sweepNodeUpdate(i_U, o_U, i_simulationTime, k, i);
				}
			}
		}
		else
		{
			for (size_t i = 0; i < M; i++) {
				_sweepNodeUpdate(i_U, o_U, i_simulationTime, k, i);
			}
		}

		/*
		 * Swap k+1 and k values for next iteration (or end-point update)
		 *
		 * (Just swap the pointers into the vector structure)
		 */
		if (implicitTerm.active)
			std::swap(implicitTerm.tendencies_k0, implicitTerm.tendencies_k1);

		if (explicitTerm.active)
			std::swap(explicitTerm.tendencies_k0, explicitTerm.tendencies_k1);
	}


private:
	void _sweepPost(
			const sweet::Data::GenericContainer::Base& i_U,
			sweet::Data::GenericContainer::Base& o_U,
			double i_timeStamp
	) {
		DEBUG("Starting post sweep");
		if (postSweep == POST_SWEEP_LASTNODE)
		{
			DEBUG("Using LASTNODE end-update");
			o_U.op_setVector(*_getStateContainer(M-1));
		}

		else if (postSweep == POST_SWEEP_QUADRATURE)
		{
			DEBUG("Using QUADRATURE end-update");
			/*
			 * Use quadrature formula to compute the end point
			 */
			const Vec& w = weights;

			// Get buffer for state vector
			sweet::Data::GenericContainer::Base *state = _getStateContainer(0);

			// Initialize with U0
			state->op_setVector(i_U);

			// Add quadrature terms
			for (size_t j = 0; j < M; j++) {
				double a = _dt*w(j);

				if (explicitTerm.active)
					state->op_addScalarMulVector(a, *explicitTerm.tendencies_k0[j]);

				if (implicitTerm.active)
					state->op_addScalarMulVector(a, *implicitTerm.tendencies_k0[j]);
			}

			o_U.op_setVector(*state);
		}

		else if (postSweep == POST_SWEEP_INTERPOLATION) 
		{
			DEBUG("Using INTERPOLATION end-update");
			// -> see _sweepNodeUpdate for pre-computation of endUpdate_U1
			o_U.op_setVector(*endUpdate_U1);
		}
	}


	bool _eval_integration(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)	override
	{
		DEBUG("starting time-step for t=", i_simulationTime);

		// initialize nodes values and state
		_sweepPre(i_U, o_U, i_simulationTime);

		// perform sweeps
		for (size_t k = 0; k < nIters; k++) {
			_sweepMain(i_U, o_U, i_simulationTime, k);
		}

		// compute end-update solution
		_sweepPost(i_U, o_U, i_simulationTime);

		return true;
	}


	void printMule(const std::string &i_prefixStr = "")
	{
		std::cout << i_prefixStr << "SDC_FP:" << std::endl;
		printSDCInformation("[MULE] sdc.");
	}

	void print(const std::string &i_prefixStr = "")
	{
		std::string newPrefix = i_prefixStr + "  ";
		std::cout << newPrefix << "SDC_FP(" << std::endl;
		std::cout << newPrefix << "  file: " << fileNameWithSDCCoefficients << std::endl;
		printSDCInformation(newPrefix);
	}
};

}}}

#endif
