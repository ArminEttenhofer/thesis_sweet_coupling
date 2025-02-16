#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_IMPLICITRUNGEKUTTA_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_IMPLICITRUNGEKUTTA_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>
#include <vector>
#include <string>


namespace sweet {
namespace TimeTree {
namespace InteriorNodes {


class ImplicitRungeKutta	:
		public TimeTree_Node_InteriorHelper<ImplicitRungeKutta>
{
private:
	/*
	 * We use an enum to identify different RK implementations
	 */
	enum IRKMethod
	{
		INVALID = -1,
		IRK1 = 1,	// 1st order forward Euler
		IRK2,
	};
	IRKMethod _rkMethodID;


	// Order of Runge-Kutta method
	int _order;

	// Particular RK method
	std::string _method;

	// Runge-Kutta stage storages
	//int _rkNumStages;

	// Damping factor for 2nd order IRK CN method. 0.5 means no damping
	double _crank_nicolson_damping_factor = 0.5;

	double _dt_explicit = -1;
	double _dt_implicit = -1;

public:
	ImplicitRungeKutta()	:
		_rkMethodID(INVALID),
		_order(-1),
		_method("std")
	{
		setEvalAvailable(EVAL_INTEGRATION);
	}

	~ImplicitRungeKutta()
	{
		clear();
	}


	ImplicitRungeKutta(
			const ImplicitRungeKutta &i_src
	)	:
		TimeTree_Node_InteriorHelper<ImplicitRungeKutta>(i_src)
	{
		_rkMethodID = i_src._rkMethodID;
		_order = i_src._order;
		_method = i_src._method;

		_crank_nicolson_damping_factor = i_src._crank_nicolson_damping_factor;

		_dt_explicit = i_src._dt_explicit;
		_dt_implicit = i_src._dt_implicit;
	}


	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("irk");
		retval.push_back("IRK");
		retval.push_back("implicitRungeKutta");
		return retval;
	}


	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override
	{
		o_ostream << i_prefix << "InteriorNode: 'ImplicitRungeKutta':" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Usage: IRK(DeTerm,[parameter1],...)" << std::endl;
		o_ostream << i_prefix << "           Compute implicit Runge-Kutta time integration." << std::endl;
		o_ostream << i_prefix << "           For 1st order, requires term to support 'backwardEuler'" << std::endl;
		o_ostream << i_prefix << "           For 2nd order, requires term to support 'backwardEuler' and 'timeTendencies'" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Parameters:" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "    - order=[int]:" << std::endl;
		o_ostream << i_prefix << "        1: Simply backward Euler" << std::endl;
		o_ostream << i_prefix << "        2: Crank-Nicolson scheme" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "    - damping=[float] [DEFAULT=0.5]:" << std::endl;
		o_ostream << i_prefix << "        Damping factor for Crank-Nicolson scheme" << std::endl;
		o_ostream << i_prefix << std::endl;

		return true;
	}


	bool _setupArgumentInternals()
	{
		_rkMethodID = INVALID;

		if (_order != -1)
		{
			if (_order == 1)
			{
				_method = "std";
			}
			else if (_order == 2)
			{
				_method = "crank_nicolson";
			}
			else
			{
				return error.set("Order of IRK method provided, only 1 or 2 supported"+getNewLineDebugMessage());
			}
		}

		if (_method == "cn" || _method == "crank_nicolson")
		{
			if (_order != -1 && _order != 2)
				return error.set("Order of Crank-Nicolson's method must be 2"+getNewLineDebugMessage());

			_rkMethodID = IRK2;
			_order = 2;
		}
		else if (_method == "std")
		{
			if (_order < 1 || _order > 2)
				return error.set("Order of IRK method needs to be 1 or 2"+getNewLineDebugMessage());

			_rkMethodID = IRK1;
			_order = 1;
		}
		else
		{
			return error.set("Unknown method '"+_method+"'");
		}

		if (_rkMethodID == INVALID)
			return error.set("Invalid time stepping method");

		if (_timeTreeNodes.size() != 1)
			return error.set("DE Term not specified for time stepper"+getNewLineDebugMessage());

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
				return error.set("key_function not supported!"+a->getNewLineDebugMessage());

			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_FUNCTION:
				if (_timeTreeNodes.size() > 0)
					return error.set("Only one DETerm is supported"+a->getNewLineDebugMessage());

				_timeTreeNodes.push_back(std::shared_ptr<sweet::TimeTree::TimeTree_Node_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;


			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_KEY_VALUE:
				if (a->key == "order" || a->key == "o")
				{
					a->getValue(_order);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
					break;
				}

				if (a->key == "method" || a->key == "m")
				{
					_method = a->value;
					break;
				}

				if (a->key == "damping" || a->key == "cn_damping")
				{
					a->getValue(_crank_nicolson_damping_factor);
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*a);
					break;
				}

				return error.set("Key not supported"+a->getNewLineDebugMessage());
				break;

			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_VALUE:
				if (_timeTreeNodes.size() > 0)
					return error.set("Only one DETerm is supported"+a->getNewLineDebugMessage());

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
		/*
		 * Manually setup the backward Euler
		 */
		_timeTreeNodes.resize(_order);
		_evalFuns.resize(_order);

		/*
		 * Setup temporary buffers
		 */
		if (_order == 1)
		{
			_timeTreeNodes[0]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EULER_BACKWARD, &_evalFuns[0]);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);
		}
		else if (_order == 2)
		{
			// Create new instance for explicit evaluation
			_timeTreeNodes[1] = _timeTreeNodes[0]->getInstanceCopy();

			// Crank-Nicolson
			_tmpDataContainer.resize(2);

			// Setup
			_timeTreeNodes[0]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EULER_BACKWARD, &_evalFuns[0]);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);

			_timeTreeNodes[1]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[1]);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[1]);
		}

		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();

		// Return time stepper for this routine
		TimeTree_Node_Base::registerTimeStepperEval(
				i_evalType,
				o_timeStepper
			);

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

		return true;

	}

	void clear() override
	{
		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			delete _tmpDataContainer[i];
		_tmpDataContainer.clear();
	}


	std::shared_ptr<TimeTree_Node_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<TimeTree_Node_Base>(new ImplicitRungeKutta(*this));
	}

	bool setTimeStepSize(double i_dt)	override
	{
		// DO NOT use this convenience function since we need to set different time step sizes
		//TimeTree_Node_InteriorHelper<ImplicitRungeKutta>::setTimeStepSize(i_dt)

		_timestepSize = i_dt;

		// Not required for explicit time stepper, but we do it
		if (_order == 1)
		{
			_timeTreeNodes[0]->setTimeStepSize(i_dt);
			_dt_explicit = i_dt;
		}
		else if (_order == 2)
		{
			// Crank-Nicolson: Use damping factor for explicit time stepper
			_dt_explicit = i_dt*(1.0-_crank_nicolson_damping_factor);
			_dt_implicit = i_dt*_crank_nicolson_damping_factor;

			// Implicit one
			_timeTreeNodes[0]->setTimeStepSize(_dt_implicit);

			// Explicit one
			_timeTreeNodes[1]->setTimeStepSize(_dt_explicit);


		}
		else
		{
			SWEETErrorFatal("Internal error");
		}

		return true;
	}

	bool _eval_integration(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)	override
	{
		switch(_rkMethodID)
		{
		case IRK1:	return _eval_timeIntegration_IRK1(i_U, o_U, i_simulationTime);
		case IRK2:	return _eval_timeIntegration_IRK2(i_U, o_U, i_simulationTime);
		default: return error.set("Internal error: Wrong IRK Method");
		}
	}

private:
	/*
	 * Backward Euler method:
	 *
	 * (U(t+1) - dt F(U(t+1))) = U(t)
	 */
	bool _eval_timeIntegration_IRK1(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		return evalTimeStepper(0, i_U, o_U, i_simulationTime);
	}

private:
	/*
	 * Crank-Nicolson method:
	 *
	 * (U(t+1) - q dt F(U(t+1))) = (U(t) + q dt F(U(t)))
	 *
	 * With q the CN damping factor with no damping for q=0.5
	 */
	bool _eval_timeIntegration_IRK2(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		// Forward Euler time step
		evalTimeStepper(1, i_U, *_tmpDataContainer[0], i_simulationTime);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

		_tmpDataContainer[1]->op_setVectorPlusScalarMulVector(
				i_U,
				_dt_explicit,
				*_tmpDataContainer[0]
			);

		// Backward Euler step
		evalTimeStepper(0, *_tmpDataContainer[1], o_U, i_simulationTime+_dt*_crank_nicolson_damping_factor);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*this);

		return true;
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "IRK(" << std::endl;
		std::cout << newPrefix << "  order: " << _order << std::endl;
		std::cout << newPrefix << "  method: " << _method << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}}}

#endif
