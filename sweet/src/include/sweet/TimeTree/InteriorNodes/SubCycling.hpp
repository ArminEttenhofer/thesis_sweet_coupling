#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_SUBCYCLING_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_SUBCYCLING_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>
#include <vector>
#include <string>


namespace sweet {
namespace TimeTree {
namespace InteriorNodes {


class SubCycling	:
	public TimeTree_Node_InteriorHelper<SubCycling>
{
private:
	// Number of subcycling intervals
	int _subCyclingIntervals;


public:
	SubCycling()
	{
		setEvalAvailable(EVAL_INTEGRATION);
	}

	~SubCycling()
	{
		clear();
	}

	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("SUB");
		retval.push_back("subC");
		retval.push_back("subCycling");
		retval.push_back("SUBC");
		retval.push_back("SUBCYCLING");
		return retval;
	}

	bool _setupArgumentInternals()
	{
		if (_subCyclingIntervals < 1)
			return error.set("At least one interval for subcycling required"+getNewLineDebugMessage());

		if (_timeTreeNodes.size() == 0)
			return error.set("Subcycling requires one term/function"+getNewLineDebugMessage());

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

				if (_timeTreeNodes.size() >= 1)
					return error.set("Subcycling only supports a single function/DETerm"+a->getNewLineDebugMessage());
					
				_timeTreeNodes.push_back(std::shared_ptr<sweet::TimeTree::TimeTree_Node_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION:
				return error.set("Key with functions not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_KEY_VALUE:
				if (a->key == "intervals" || a->key == "i" || a->key == "n")
				{
					a->getValue(_subCyclingIntervals);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
					break;
				}

				return error.set("Key not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_VALUE:
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

		// provide debug message in case that something goes wrong with the arguments
		setDebugMessage(i_function->getDebugMessage());
		return _setupArgumentInternals();
	}


	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override
	{
		o_ostream << i_prefix << "InteriorNode: 'Sub cycling (SUBC)':" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Usage: SUBC(f,[parameter1],...)" << std::endl;
		o_ostream << i_prefix << "           f: Term to be time integrated multiple times for a smaller time step size" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Parameters:" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "    - n=[int]:" << std::endl;
		o_ostream << i_prefix << "        Specify number of subintervals to split time interval into" << std::endl;
		o_ostream << i_prefix << std::endl;

		return true;
	}


	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper
	) override
	{
		_helperSetupConfigAndForwardTimeStepperEval(
				i_deTermConfig,
				i_evalType,
				o_timeStepper,
				EVAL_INTEGRATION
			);
		
		_tmpDataContainer.resize(1);
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();
		
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

		return true;
	}

	bool setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		for (auto &i : _timeTreeNodes)
		{
			i->setTimeStepSize(_timestepSize/_subCyclingIntervals);
		}

		return true;
	}

	bool _eval_integration(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)	override
	{
		bool retval = true;
		if (_subCyclingIntervals == 1)
		{
			retval &= evalTimeStepper(i_U, o_U, i_simulationTime);
			return retval;
		}

		int stepsTodo = _subCyclingIntervals;

		// 1st step to use temporary buffer
		retval &= evalTimeStepper(i_U, *_tmpDataContainer[0], i_simulationTime);
		stepsTodo--;

		if (_subCyclingIntervals & 1)
		{
			/*
			 * Odd number of time steps
			 */

			for (int i = 1; i < _subCyclingIntervals-1; i+=2)
			{
				retval &= evalTimeStepper(*_tmpDataContainer[0], o_U, i_simulationTime);
				retval &= evalTimeStepper(o_U, *_tmpDataContainer[0], i_simulationTime);
				stepsTodo -= 2;
			}

			_tmpDataContainer[0]->swap(o_U);
		}
		else
		{
			/*
			 * Even number of time steps
			 */
			for (int i = 1; i < _subCyclingIntervals-1; i+=2)
			{
				retval &= evalTimeStepper(*_tmpDataContainer[0], o_U, i_simulationTime);
				retval &= evalTimeStepper(o_U, *_tmpDataContainer[0], i_simulationTime);
				stepsTodo -= 2;
			}

			retval &= evalTimeStepper(*_tmpDataContainer[0], o_U, i_simulationTime);
			stepsTodo--;
		}

		SWEET_ASSERT(stepsTodo == 0);
		return retval;
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "SUBCYCLING(" << std::endl;
		std::cout << newPrefix << "  subCyclingIntervals: " << _subCyclingIntervals << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}}}

#endif
