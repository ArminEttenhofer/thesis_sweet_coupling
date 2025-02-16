#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_EXPONENTIAL_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_EXPONENTIAL_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>
#include <vector>
#include <string>


namespace sweet {
namespace TimeTree {
namespace InteriorNodes {


class Exponential	:
	public TimeTree_Node_InteriorHelper<Exponential>
{
private:
	std::string _expFunctionString;

public:
	Exponential()
	{
		setEvalAvailable(EVAL_EXPONENTIAL);
		setEvalAvailable(EVAL_INTEGRATION);
	}


	~Exponential()
	{
		clear();
	}

	Exponential(
			const Exponential &i_src
	)	:
		TimeTree_Node_InteriorHelper<Exponential>(i_src)
	{
		_expFunctionString = i_src._expFunctionString;
	}


	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("exp");
		retval.push_back("EXP");
		return retval;
	}



	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override
	{
		o_ostream << i_prefix << "InteriorNode: 'Exponential':" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Usage: EXP(DeTerm,[parameter1],...)" << std::endl;
		o_ostream << i_prefix << "           Compute exponential integration of DeTerm." << std::endl;
		o_ostream << i_prefix << "           Requires term to support 'exponential'" << std::endl;
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
		if (_timeTreeNodes.size() == 0)
			return error.set("Some time node term needs to be given"+getNewLineDebugMessage());

		if (_expFunctionString != "")
			_timeTreeNodes[0]->setupByKeyValue("ExpIntegrationFunction", _expFunctionString);

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);

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

				if (a->key == "fun")
				{
					_expFunctionString = a->value;
					break;
				}

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

	bool setupByKeyValue(
			const std::string &i_key,
			const std::string &i_value
	) override
	{
		if (i_key == "ExpIntegrationFunction")
		{
			SWEET_ASSERT(_timeTreeNodes[0] != nullptr);
			_timeTreeNodes[0]->setupByKeyValue(i_key, i_value);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);

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
		_helperSetupConfigAndForwardTimeStepperEval(
				i_deTermConfig,
				i_evalType,
				o_timeStepper,
				EVAL_EXPONENTIAL
			);

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

		return true;
	}


	std::shared_ptr<TimeTree_Node_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<TimeTree_Node_Base>(new Exponential(*this));
	}


private:
	bool _eval_integration(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)	override
	{
		SWEET_ASSERT(_timeTreeNodes[0] != nullptr);
		return evalTimeStepper(i_U, o_U, i_simulationTime);
	}



	/*
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
		SWEET_ASSERT(_timeTreeNodes[0] != nullptr);
		return evalTimeStepper(i_U, o_U, i_simulationTime);
	}


	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "EXP(" << std::endl;
		std::cout << newPrefix << "  expFunctionString: '" << _expFunctionString << "'" << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}}}

#endif
