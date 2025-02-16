#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_STRANGSPLITTING_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_STRANGSPLITTING_HPP

#include <sweet/Data/GenericContainer/ConfigBase.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/Error/Fatal.hpp>
#include <sweet/TimeTree/TimeTree_Node_Base.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>
#include <sweet/TimeTree/TimeTreeIR.hpp>
#include <sweet/TimeTree/TimeTreeIR_2_TimeTreeNodes.hpp>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <vector>


namespace sweet {
namespace TimeTree {
namespace InteriorNodes {


class StrangSplitting	:
		public TimeTree_Node_InteriorHelper<StrangSplitting>
{
private:
	// Order of Strang splitting
	int _order;

public:
	StrangSplitting()
	{
		setEvalAvailable(EVAL_INTEGRATION);
	}

	~StrangSplitting()
	{
		clear();
	}

	StrangSplitting(
			const StrangSplitting &i_src
	)	:
		TimeTree_Node_InteriorHelper<StrangSplitting>(i_src)
	{
		_order = i_src._order;
	}


	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("ss");
		retval.push_back("SS");
		return retval;
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


	bool _setupArgumentInternals()
	{
		if (_timeTreeNodes.size() != 2)
			return error.set("Only two time steppers supported in this strang splitting"+getNewLineDebugMessage());

		if (_order < 1 || _order > 2)
			return error.set("Only order 1 or 2 allowed for SS method"+getNewLineDebugMessage());

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
				_timeTreeNodes.push_back(std::shared_ptr<sweet::TimeTree::TimeTree_Node_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_VALUE:
				_timeTreeNodes.push_back(std::shared_ptr<sweet::TimeTree::TimeTree_Node_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByName(
						a->value,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION:
				return error.set("Key with functions not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_KEY_VALUE:
				if (a->key == "order" || a->key == "o")
				{
					a->getValue(_order);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
					break;
				}

				return error.set("Key not supported"+a->getNewLineDebugMessage());
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



	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override
	{
		o_ostream << i_prefix << "InteriorNode: 'Strang Splitting (SS)':" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Usage: SS(f1,f2,[parameter1],...)" << std::endl;
		o_ostream << i_prefix << "           f1: First Strang splitting term" << std::endl;
		o_ostream << i_prefix << "           f2: Second Strang splitting term" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Parameters:" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "    - order=[int]:" << std::endl;
		o_ostream << i_prefix << "        Specify the order of the Strang splitting" << std::endl;
		o_ostream << i_prefix << std::endl;

		return true;
	}

	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper
	) override
	{
		assert(_timeTreeNodes.size() == 2);
		_helperSetupConfigAndForwardTimeStepperEval(
				i_deTermConfig,
				i_evalType,
				o_timeStepper,
				EVAL_INTEGRATION
			);

		if (_order == 1)
			_tmpDataContainer.resize(1);
		else if (_order == 2)
			_tmpDataContainer.resize(2);

		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();

		return true;
	}

	bool setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		if (_order == 1)
		{
			for (auto &i : _timeTreeNodes)
			{
				i->setTimeStepSize(i_dt);
			}
			return true;
		}

		if (_order == 2)
		{
			_timeTreeNodes[0]->setTimeStepSize(i_dt*0.5);
			_timeTreeNodes[1]->setTimeStepSize(i_dt);
			return true;
		}

		SWEETErrorFatal("Internal error");
		return false;
	}




	bool _eval_integration(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)	override
	{
		if (_order == 1)
		{
			evalTimeStepper(0, i_U, *_tmpDataContainer[0], i_simulationTime);
#if SWEET_DEBUG
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);
#endif
			evalTimeStepper(1, *_tmpDataContainer[0], o_U, i_simulationTime);
#if SWEET_DEBUG
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[1]);
#endif
			return true;
		}
		
		if (_order == 2)
		{
			evalTimeStepper(0, i_U, *_tmpDataContainer[0], i_simulationTime);
#if SWEET_DEBUG
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);
#endif
			evalTimeStepper(1, *_tmpDataContainer[0], *_tmpDataContainer[1], i_simulationTime);
#if SWEET_DEBUG
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[1]);
#endif
			evalTimeStepper(0, *_tmpDataContainer[1], o_U, i_simulationTime+0.5*_dt);
#if SWEET_DEBUG
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);
#endif
			return true;
		}

		return error.set("Internal error: Wrong order");
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "SS(" << std::endl;
		std::cout << newPrefix << "  order: " << _order << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}}}

#endif
