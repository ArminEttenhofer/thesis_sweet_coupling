#ifndef INCLUDE_SWEET_TIMETREE_TIMETREE_NODE_INTERIORHELPER_HPP
#define INCLUDE_SWEET_TIMETREE_TIMETREE_NODE_INTERIORHELPER_HPP

#include <sweet/TimeTree/TimeTree_Node_Base.hpp>
#include <sweet/TimeTree/TimeTreeIR_2_TimeTreeNodes.hpp>

namespace sweet {
namespace TimeTree {


/*!
 * Helper class for interior tree node
 *
 * It provides default member variables which are often used
 */
template <typename MainInteriorNodeClass>
class TimeTree_Node_InteriorHelper	:
	public TimeTree_Node_Base
{
protected:
	double _timestepSize;
	double &_dt = _timestepSize;

	// DE term to evaluate
	std::vector<std::shared_ptr<sweet::TimeTree::TimeTree_Node_Base>> _timeTreeNodes;
	std::vector<TimeTree_Node_Base::EvalFun> _evalFuns;

	// Number of stages to allocate buffers
	std::vector<sweet::Data::GenericContainer::Base*> _tmpDataContainer;

public:
	TimeTree_Node_InteriorHelper()	:
		_timestepSize(-1)
	{
	}

	virtual
	~TimeTree_Node_InteriorHelper()
	{
		clear();
	}


	/*
	 * Copy constructor
	 *
	 * Simply copy the raw data over here
	 */
	TimeTree_Node_InteriorHelper(
		const TimeTree_Node_InteriorHelper &i_src
	):
		TimeTree_Node_Base(i_src)
	{
		_timestepSize = i_src._timestepSize;

		// registered evaluation functions
		_registeredEvalTypes = i_src._registeredEvalTypes;

		// Data containers
		_tmpDataContainer.resize(i_src._tmpDataContainer.size());
		for (std::size_t i = 0; i < i_src._tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_src._tmpDataContainer[i]->getNewDataContainer();

		// Time tree nodes
		_timeTreeNodes.resize(i_src._timeTreeNodes.size());
		for (std::size_t i = 0; i < _timeTreeNodes.size(); i++)
			_timeTreeNodes[i] = i_src._timeTreeNodes[i]->getInstanceCopy();

		// _evalFuns is handled later (hopefully)
	}

	bool _helperSetupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,	//!< User-defined configuration
		TIME_STEPPER_TYPES i_thisTimeStepperEvalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper,
		TIME_STEPPER_TYPES i_childEvalType
	)
	{
		_evalFuns.resize(_timeTreeNodes.size());
		for (std::size_t i = 0; i < _evalFuns.size(); i++)
		{
			_timeTreeNodes[i]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, i_childEvalType, &_evalFuns[i]);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[i]);
		}

#if SWEET_XBRAID
		U_prev_solution = i_deTermConfig.getNewDataContainerInstance(Data::GenericContainer::Base::DATA_SIMULATION);
#endif

		// default setup
		TimeTree_Node_Base::registerTimeStepperEval(
				i_thisTimeStepperEvalType,
				o_timeStepper
			);
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

		return true;
	}

	bool setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		for (auto &i : _timeTreeNodes)
		{
			i->setTimeStepSize(_timestepSize);
		}

		return true;
	}


	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	)	override
	{
		for (auto &i : _timeTreeNodes)
		{
			i->shackRegistration(io_shackDict);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*i);
		}

		return true;
	}

	void clear() override
	{
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
		{
			SWEET_ASSERT(_tmpDataContainer[i] != nullptr);
			delete _tmpDataContainer[i];
		}

		_tmpDataContainer.clear();
	}


	/*
	 * Helper function to call the functions provided by pointers
	 *
	 * This calls the 1st one in the vector
	 */
	inline
	bool evalTimeStepper(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		SWEET_ASSERT(_timeTreeNodes.size() > 0);
		SWEET_ASSERT(_evalFuns.size() > 0);

		(_timeTreeNodes[0].get()->*_evalFuns[0])(i_U, o_U, i_simulationTime);
#if SWEET_DEBUG
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);
#endif
		return true;
	}

	/*
	 * Helper function to call the functions provided by pointers
	 */
	inline
	bool evalTimeStepper(
			std::size_t i_id,
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		SWEET_ASSERT(_timeTreeNodes.size() > i_id);
		SWEET_ASSERT(_evalFuns.size() > i_id);

		SWEET_ASSERT(_timeTreeNodes[i_id] != nullptr);
		SWEET_ASSERT(_evalFuns[i_id] != nullptr);

		(_timeTreeNodes[i_id].get()->*_evalFuns[i_id])(i_U, o_U, i_simulationTime);
#if SWEET_DEBUG
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[i_id]);
#endif
		return true;
	}

	std::shared_ptr<TimeTree_Node_Base> getInstanceCopy()	override
	{
		MainInteriorNodeClass *m = new MainInteriorNodeClass(static_cast<MainInteriorNodeClass&>(*this));
		return std::shared_ptr<TimeTree_Node_Base>(m);
	}

#if SWEET_XBRAID
public:
	bool storePrevSolution(
			sweet::Data::GenericContainer::Base* i_U		//!< Input simulation data
	) override
	{
		U_prev_solution->op_setVector(*i_U);
	}
#endif




};

}}

#endif
