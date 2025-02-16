#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_ADDTENDENCIES_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_ADDTENDENCIES_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>
#include <vector>
#include <string>


namespace sweet {
namespace TimeTree {
namespace InteriorNodes {


/*!
 * Add tendencies of different tree nodes
 */
class AddTendencies	:
	public TimeTree_Node_InteriorHelper<AddTendencies>
{
public:
	AddTendencies()
	{
		setEvalAvailable(EVAL_TENDENCIES);
	}

	~AddTendencies()
	{
		clear();
	}

	AddTendencies(
			const AddTendencies &i_src
	)	:
		TimeTree_Node_InteriorHelper<AddTendencies>(i_src)
	{
	}


	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("addTendencies");
		retval.push_back("addT");
		retval.push_back("ADDTENDENCIES");
		retval.push_back("ADDT");
		return retval;
	}

	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override
	{
		o_ostream << i_prefix << "InteriorNode: 'AddTendencies':" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Usage: ADDT(d1,d2,d3,...)" << std::endl;
		o_ostream << i_prefix << "           Adds together evaluation of time tendencies" << std::endl;
		o_ostream << i_prefix << "           (In other words, it computes 'd1+d2+d3+...')" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Parameters: none" << std::endl;

		return true;
	}

	bool _setupArgumentInternals()
	{
		if (_timeTreeNodes.size() == 0)
			return error.set("No DE terms specified for time stepper"+getNewLineDebugMessage());

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
			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION:
			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_FUNCTION:
				_timeTreeNodes.push_back(std::shared_ptr<sweet::TimeTree::TimeTree_Node_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_KEY_VALUE:
				return error.set("Key-value not supported"+a->getNewLineDebugMessage());
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
				EVAL_TENDENCIES
			);
		
		_tmpDataContainer.resize(1);
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();
		
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

		return true;
	}

	std::shared_ptr<TimeTree_Node_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<TimeTree_Node_Base>(new AddTendencies(*this));
	}

private:
	bool _eval_tendencies(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)	override
	{
		o_U.op_setZero();
		SWEET_ASSERT(_timeTreeNodes.size() == _evalFuns.size());

		for (std::size_t i = 0; i < _evalFuns.size(); i++)
		{
			SWEET_ASSERT(_evalFuns[i] != nullptr);
		
			evalTimeStepper(
					i,
					i_U,
					*_tmpDataContainer[0],
					i_simulationTime
				);
			o_U.op_addVector(*_tmpDataContainer[0]);
		}

#if SWEET_DEBUG
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);
#endif

		return true;
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "Add(" << std::endl;
		std::cout << newPrefix << "  numDETerms: " << _timeTreeNodes.size() << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}}}

#endif
