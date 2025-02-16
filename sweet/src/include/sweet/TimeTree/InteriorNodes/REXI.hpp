#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_REXI_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_REXI_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <vector>
#include <string>
#include <sweet/ExpIntegration/ExpFunction.hpp>
#include <sweet/ExpIntegration/REXI/REXI.hpp>
#include <sweet/ExpIntegration/REXI/REXICoefficients.hpp>
#include <sweet/ExpIntegration/Shack.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>


namespace sweet {
namespace TimeTree {
namespace InteriorNodes {


class REXI	:
	public TimeTree_Node_InteriorHelper<REXI>
{
private:
	/*!
	 * Shack to exponential integration information including REXI information
	 */
	sweet::ExpIntegration::Shack *_shackExpIntegration;

	/*!
	 * String describing the exponential function ("phi0", "phi1" or ...)
	 */
	std::string _expFunctionString;

	/*!
	 * REXI coefficients
	 */
	std::vector<std::complex<double>> _rexiAlphas;
	std::vector<std::complex<double>> _rexiBetas;
	std::complex<double> _rexiGamma;
	bool _rexi_gamma_active;

	/*!
	 * Number of REXI terms processed by this MPI rank
	 */
	std::size_t _num_local_rexi_terms;

#if SWEET_MPI

	/*!
	 * MPI communicator
	 */
	MPI_Comm _mpi_comm;

	/*!
	 * Number of mpi ranks to be used
	 */
	int _mpi_comm_rank;

	/*!
	 * MPI ranks
	 */
	int _mpi_comm_size;
#endif

public:
	REXI()	:
		_shackExpIntegration(nullptr),
		_rexi_gamma_active(false),
		_num_local_rexi_terms(-1)
	{
		setEvalAvailable(EVAL_EXPONENTIAL);
		setEvalAvailable(EVAL_INTEGRATION);
	}


	~REXI()
	{
		clear();
	}


	/*!
	 * This constructor is used for creating copies of this time stepper.
	 *
	 * This is used in two different ways:
	 *
	 * 1) Creating new instances from the registry
	 *
	 * 2) Duplicate existing tree nodes with its parameters:
	 *    - to evaluate different phi functions
	 *    - to evaluate different REXI terms
	 *    - to support different time step sizes
	 */
	REXI(
			const REXI &i_src
	)	:
		TimeTree_Node_InteriorHelper(i_src)
	{
		_shackExpIntegration = i_src._shackExpIntegration;

		SWEET_ASSERT(_shackExpIntegration != nullptr);

		_rexiAlphas = i_src._rexiAlphas;
		_rexiBetas = i_src._rexiBetas;
		_rexiGamma = i_src._rexiGamma;
		_rexi_gamma_active = i_src._rexi_gamma_active;

		_num_local_rexi_terms = i_src._num_local_rexi_terms;
	}


	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("rexi");
		retval.push_back("REXI");
		return retval;
	}


	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)	override
	{
		_shackExpIntegration = io_shackDict->getAutoRegistration<sweet::ExpIntegration::Shack>();

		return TimeTree_Node_InteriorHelper::shackRegistration(io_shackDict);
	}


	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override
	{
		o_ostream << i_prefix << "InteriorNode: 'REXI':" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Usage: REXI(DeTerm,[parameter1],...)" << std::endl;
		o_ostream << i_prefix << "           Compute exponential integration based on REXI of DeTerm." << std::endl;
		o_ostream << i_prefix << "           Requires term to support 'rexiTerm'" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Parameters:" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "    - fun=[str]:" << std::endl;
		o_ostream << i_prefix << "        Specify the phi-function as a string ('phi0', 'phi1', ...)" << std::endl;
		o_ostream << i_prefix << std::endl;

		return true;
	}

	bool _setupArgumentInternals()
	{
		if (_timeTreeNodes.size() != 1)
			return error.set("One time node term needs to be given"+getNewLineDebugMessage());

		if (_expFunctionString == "")
			_expFunctionString = "phi0";

		return true;
	}


	bool setupTreeNodeByFunction(
			std::shared_ptr<sweet::TimeTree::TimeTreeIR::Function> &i_function,
			sweet::TimeTree::TimeTreeIR_2_TimeTreeNodes &i_tsAssemblation
	)	override
	{
		for (auto iter = i_function->arguments.begin(); iter != i_function->arguments.end(); iter++)
		{
			sweet::TimeTree::TimeTreeIR::Argument *a = iter->get();

			switch(a->argType)
			{
			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_FUNCTION:
			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION:
				if (_timeTreeNodes.size() == 1)
					return error.set("Only one term for exponential integration supported");

				_timeTreeNodes.push_back(std::shared_ptr<sweet::TimeTree::TimeTree_Node_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNodes.back()
					);

				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_KEY_VALUE:

				return error.set("Key not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_VALUE:
				if (_timeTreeNodes.size() != 0)
					return error.set("Only one term for exponential integration supported");

				_timeTreeNodes.push_back(std::shared_ptr<sweet::TimeTree::TimeTree_Node_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByName(
						a->value,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			default:
				SWEETErrorFatal("Internal error");
				return error.set("Internal error");
			}
		}

		// Provide debug message in case that something goes wrong with the arguments
		setDebugMessage(i_function->getDebugMessage());
		return _setupArgumentInternals();
	}


	/*!
	 * Special key/value setup
	 */
	bool setupByKeyValue(
			const std::string &i_key,
			const std::string &i_value
	) override
	{
		if (i_key == "ExpIntegrationFunction")
		{
			_expFunctionString = i_value;
			return true;
		}

		return false;
	}


	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper
	) override
	{
		SWEET_ASSERT(_timeTreeNodes.size() == 1);

#if SWEET_MPI
		_mpi_comm = MPI_COMM_WORLD;	// TODO: Make me more flexible in future versions
		MPI_Comm_rank(_mpi_comm, &_mpi_comm_rank);
		MPI_Comm_size(_mpi_comm, &_mpi_comm_size);
#endif

		/*
		 * Load REXI coefficients
		 */
		sweet::ExpIntegration::REXI::REXICoefficients<double> rexiCoefficients;
		sweet::ExpIntegration::REXI::REXI<> rexi;

		if (_shackExpIntegration->exp_method == "direct")
			return error.set("REXI: Direct exponential method is available with EXP() time tree function. REXI requires providing coefficients.");

		rexi.load(
				_shackExpIntegration,
				_expFunctionString,
				rexiCoefficients,
				_shackExpIntegration->verbosity
		);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(rexi);

		_rexiAlphas = rexiCoefficients.alphas;
		_rexiBetas = rexiCoefficients.betas;

		_rexiGamma = rexiCoefficients.gamma;
		if (_rexiGamma.real() != 0 || _rexiGamma.imag() != 0)
			_rexi_gamma_active = true;
		else
			_rexi_gamma_active = false;

		SWEET_ASSERT(_rexiAlphas.size() > 0);

		/*
		 * Compute for which REXI coefficients we're responsible for
		 */
		std::size_t N = _rexiAlphas.size();

#if SWEET_MPI
		_num_local_rexi_terms = N/_mpi_comm_size;
		if (_num_local_rexi_terms*_mpi_comm_size != N)
			_num_local_rexi_terms++;

		int offset = _mpi_comm_rank*_num_local_rexi_terms;

		_num_local_rexi_terms = std::min(offset+_num_local_rexi_terms, N)-offset;

#else
		_num_local_rexi_terms = N;
		int offset = 0;
#endif


		/*
		 * Allocate data structures, time steppers, etc.
		 */
		_timeTreeNodes.resize(_num_local_rexi_terms);
		_evalFuns.resize(_num_local_rexi_terms);
		_tmpDataContainer.resize(_num_local_rexi_terms);

		/*
		 * Initialize time tree node
		 */
#if SWEET_THREADING_TIME_REXI
		#pragma omp parallel for	\
				schedule(static)	\
				default(none)	\
				shared(_num_local_rexi_terms, _timeTreeNodes, i_deTermConfig, _tmpDataContainer)
#endif
		for (std::size_t i = 0; i < _num_local_rexi_terms; i++)
		{
			// Create instances from first timeTreeNode. We keep it in for equal split of for loop
			if (i != 0)
				_timeTreeNodes[i] = _timeTreeNodes[0]->getInstanceCopy();
		}


#if SWEET_THREADING_TIME_REXI
		#pragma omp parallel for schedule(static) \
				default(none)	\
				shared(_num_local_rexi_terms, _timeTreeNodes, i_deTermConfig, _tmpDataContainer, offset, N)
#endif
		for (std::size_t i = 0; i < _num_local_rexi_terms; i++)
		{
			/*
			 * Setup evaluation
			 */
			_timeTreeNodes[i]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_REXI_TERM, &_evalFuns[i]);
#if SWEET_THREADING_TIME_REXI
			error.forward(_timeTreeNodes[i]->error);
#else
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[i]);
#endif

			/*
			 * REXI alpha and beta coefficients
			 */
			// Don't use this assert since it will break with GCC8
			//SWEET_ASSERT(offset+i < N);
			_timeTreeNodes[i]->setupByKeyValue(
					"rexiTermAlpha",
					_rexiAlphas[offset+i]
				);

			_timeTreeNodes[i]->setupByKeyValue(
					"rexiTermBeta",
					_rexiBetas[offset+i]
				);

			if (i == 0 && offset == 0)
			{
				/*
				 * Add the gamma term with the workload handling the 1st REXI term
				 */
				if (_rexi_gamma_active)
				{
					_timeTreeNodes[i]->setupByKeyValue(
							"rexiTermGamma",
							_rexiGamma
						);
				}
			}

			/*
			 * setup data container for output, use argument "1" to request complex data
			 */
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();


#if SWEET_THREADING_TIME_REXI
			error.forward(_timeTreeNodes[i]->error);
#else
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[i]);
#endif
		}

		if (error.exists())
			return false;

		// default setup
		TimeTree_Node_Base::registerTimeStepperEval(
				i_evalType,
				o_timeStepper
			);

		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

		return true;
	}



	std::shared_ptr<TimeTree_Node_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<TimeTree_Node_Base>(new REXI(*this));
	}

	bool setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		SWEET_ASSERT(_timeTreeNodes.size() == _num_local_rexi_terms);

		std::size_t N = _timeTreeNodes.size();

#if SWEET_THREADING_TIME_REXI
		#pragma omp parallel for	\
			schedule(static)	\
			default(none)		\
			shared(N, _timeTreeNodes, i_dt)
#endif
		for (std::size_t i = 0; i < N; i++)
		{
			_timeTreeNodes[i]->setTimeStepSize(i_dt);
		}
		return true;
	}




private:
	bool _eval_integration(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)	override
	{
		return _eval_exponential(i_U, o_U, i_simulationTime);
	}


	/*!
	 * We also provide an exponential time integration for this one
	 * in order to transparently support 'exponential' time integration for
	 * either DE terms themselves, EXP and also REXI evaluations.
	 */
private:
	bool _eval_exponential(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)	override
	{
		const sweet::Data::GenericContainer::Base *i_U_real;

#if SWEET_MPI
		o_U.op_setVector(i_U);
		o_U.mpiBcast(_mpi_comm);
		i_U_real = &o_U;
#else
		i_U_real = &i_U;
#endif


		std::size_t N = _timeTreeNodes.size();

#if SWEET_THREADING_TIME_REXI
		#pragma omp parallel for	\
		schedule(static)	\
		default(none)								\
		shared(N, _timeTreeNodes,i_U_real,_tmpDataContainer,i_simulationTime)
#endif
		for (std::size_t i = 0; i < N; i++)
		{
			evalTimeStepper(
					i,
					*i_U_real,
					*_tmpDataContainer[i],
					i_simulationTime
				);
		}

		/*
		 * REDUCE operation
		 */
#if SWEET_MPI
		/*
		 * Step 1) Reduce to first tmpDataContainer.
		 *         Container #0 already contains data from the 1st REXI term
		 *
		 * Step 2) Call MPIReduce
		 */
		for (std::size_t i = 1; i < _timeTreeNodes.size(); i++)
		{
			_tmpDataContainer[0]->op_addVector(*_tmpDataContainer[i]);
		}

		o_U.mpiReduce(
				*_tmpDataContainer[0],
				_mpi_comm
			);

#else
		o_U.op_setZero();
		for (std::size_t i = 0; i < _timeTreeNodes.size(); i++)
		{
			o_U.op_addVector(*_tmpDataContainer[i]);
		}
#endif

		return true;
	}


	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "REXI(" << std::endl;
		std::cout << newPrefix << "  expFunctionString: '" << _expFunctionString << "'" << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}}}

#endif
