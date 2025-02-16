#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_LAWSONRUNGEKUTTA_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_LAWSONRUNGEKUTTA_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>
#include <sweet/TimeTree/InteriorNodes/Exponential.hpp>
#include <vector>
#include <string>

# define LRK_PHI_IDX 2


namespace sweet {
namespace TimeTree {
namespace InteriorNodes {

/*!
  \brief Lawson Runge-Kutta (Lawson 1967' paper)
 */
class LawsonRungeKutta	:
	public TimeTree_Node_InteriorHelper<LawsonRungeKutta>
{
private:
	// Different underlying RK implementations
	enum LRKMethod
	{
		INVALID = -1,
		LRK1 = 1,		// forward Euler
		LRK2_MIDPOINT,	// midpoint
		LRK3_HEUN,		// Heun's method
		LRK4_CLASSICAL,	// Runge-Kutta 4(5)

	};
	int _order;
	int _rkNumStages;

	LRKMethod _rkMethodID;
	std::string _method;

	// Coefficients of underlying RK method substeps c_i
	std::vector<double> _c;

public:
	// temporary container
	sweet::Data::GenericContainer::Base* tmp;

	// container for stages: ki = e(-cihL) f(tn+cih,pi)
	sweet::Data::GenericContainer::Base** k;
	// container for stages: pi = e(cihL) [yn + h{aijkj}]
	sweet::Data::GenericContainer::Base** p;

public:
	LawsonRungeKutta()	:
		_rkMethodID(INVALID),
		_order(-1),
		_method("std"),
		_rkNumStages(-1),
		tmp(nullptr)
	{
		setEvalAvailable(EVAL_INTEGRATION);
	}

	~LawsonRungeKutta()
	{
		clear();
	}

	LawsonRungeKutta(
			const LawsonRungeKutta &i_src
	)	:
		TimeTree_Node_InteriorHelper<LawsonRungeKutta>(i_src),
		tmp(nullptr)
	{
		_rkMethodID = i_src._rkMethodID;
		_order = i_src._order;
		_method = i_src._method;
		_rkNumStages = i_src._rkNumStages;
	}

	const std::vector<std::string> getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("lrk");
		retval.push_back("LRK");
		retval.push_back("LawsonRungeKutta");
		return retval;
	}

	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override
	{
		o_ostream << i_prefix << "InteriorNode: 'LawsonRungeKutta':" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Usage: LRK(DeTerm,[parameter1],...)" << std::endl;
		o_ostream << i_prefix << "           Compute explicit Runge-Kutta based time integrations." << std::endl;
		o_ostream << i_prefix << "           DeTerm has to support 'tendencies'" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Parameters:" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "    - order=[int]:" << std::endl;
		o_ostream << i_prefix << "        Specify the order of the time integration method" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "    - method=[str]:" << std::endl;
		o_ostream << i_prefix << "      - with order=2:" << std::endl;
		o_ostream << i_prefix << "        'midpoint': Midpoint's rule" << std::endl;
		o_ostream << i_prefix << "      - with order=3:" << std::endl;
		o_ostream << i_prefix << "        'heun3': Heun's 3rd order method" << std::endl;
		o_ostream << i_prefix << "      - with order=4:" << std::endl;
		o_ostream << i_prefix << "        'classical': Classical RK4 method" << std::endl;
		o_ostream << i_prefix << std::endl;

		return true;
	}

	bool _setupArgumentInternals()
	{
		if (_method == "midpoint")
		{
			_rkMethodID = LRK2_MIDPOINT;
			_order = 2;
			_rkNumStages = 2;
			_c.assign({0.5});
		}
		else if (_method == "heun3")
		{
			_rkMethodID = LRK3_HEUN;
			_order = 3;
			_rkNumStages = 3;
			_c.assign({1./3., 2./3.});
		}
		else if (_method == "classical")
		{
			_rkMethodID = LRK4_CLASSICAL;
			_order = 4;
			_rkNumStages = 4;
			_c.assign({0.5, 0.5, 1.});
		}
		else if (_method == "std" || _method == "default")
		{
			if (_order < 1 || _order > 4)
				return error.set("Order of LRK method needs to be 1, 2, 3 or 4"+getNewLineDebugMessage());

			switch(_order)
			{
				case 1:	_rkMethodID = LRK1; _rkNumStages = 1; break;
				case 2:	_rkMethodID = LRK2_MIDPOINT; _rkNumStages = 2; _c.assign({0.5}); break;
				case 3:	_rkMethodID = LRK3_HEUN; _rkNumStages = 3; _c.assign({1./3., 2./3.}); break;
				case 4:	_rkMethodID = LRK4_CLASSICAL; _rkNumStages = 4; _c.assign({0.5, 0.5, 1.}); break;
			}
		}
		else
		{
			return error.set("Unknown method '"+_method+"'");
		}
		SWEET_ASSERT(_rkNumStages >= 1);

		if (_rkMethodID == INVALID)
			return error.set("Invalid time stepping method");

		return true;
	}

	bool setupTreeNodeByFunction(
			std::shared_ptr<sweet::TimeTree::TimeTreeIR::Function> &i_function,
			sweet::TimeTree::TimeTreeIR_2_TimeTreeNodes &i_tsAssemblation
	)	override
	{
		_timeTreeNodes.push_back(std::shared_ptr<sweet::TimeTree::TimeTree_Node_Base>());
		_timeTreeNodes.push_back(std::shared_ptr<sweet::TimeTree::TimeTree_Node_Base>());

		for (auto iter = i_function->arguments.begin(); iter != i_function->arguments.end(); iter++)
		{
			TimeTreeIR::Argument *a = iter->get();

			switch(a->argType)
			{
			case TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION:
			case TimeTreeIR::Argument::ARG_TYPE_KEY_VALUE:
				if (a->key == "l" || a->key == "linear")
				{
					if (a->argType == TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION)
					{
						i_tsAssemblation.assembleTimeTreeNodeByFunction(
								a->function,
								_timeTreeNodes[0]
							);
					}
					else
					{
						i_tsAssemblation.assembleTimeTreeNodeByName(
								a->value,
								_timeTreeNodes[0]
							);
					}
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
					break;
				}
				if (a->key == "r" || a->key == "remainder")
				{
					if (a->argType == TimeTreeIR::Argument::ARG_TYPE_KEY_FUNCTION)
					{
						i_tsAssemblation.assembleTimeTreeNodeByFunction(
								a->function,
								_timeTreeNodes[1]
							);
					}
					else
					{
						i_tsAssemblation.assembleTimeTreeNodeByName(
								a->value,
								_timeTreeNodes[1]
							);
					}
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
					break;
				}
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
				return error.set("key not supported!"+a->getNewLineDebugMessage());

			case TimeTreeIR::Argument::ARG_TYPE_FUNCTION:
			case TimeTreeIR::Argument::ARG_TYPE_VALUE:
				// linear term comes first
				if (iter == i_function->arguments.begin())
				{
					if (a->argType == TimeTreeIR::Argument::ARG_TYPE_FUNCTION)
					{
						i_tsAssemblation.assembleTimeTreeNodeByFunction(
								a->function,
								_timeTreeNodes[0]
							);
					}
					else
					{
						i_tsAssemblation.assembleTimeTreeNodeByName(
								a->value,
								_timeTreeNodes[0]
							);
					}
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
					break;
				}
				// remainder comes second
				if (iter == i_function->arguments.begin() + 1)
				{
					if (a->argType == TimeTreeIR::Argument::ARG_TYPE_FUNCTION)
					{
						i_tsAssemblation.assembleTimeTreeNodeByFunction(
								a->function,
								_timeTreeNodes[1]
							);
					}
					else
					{
						i_tsAssemblation.assembleTimeTreeNodeByName(
								a->value,
								_timeTreeNodes[1]
							);
					}
					ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
					break;
				}
				return error.set("Too many DE terms, only 2 are supported!"+a->getNewLineDebugMessage());

			default:
				SWEETErrorFatal("Internal error");
				return error.set("Internal error");
			}
		}

		// check L() and N()
		if (!_timeTreeNodes[0] || !_timeTreeNodes[1])
			return error.set("Linear or remainder term was not provided"+getNewLineDebugMessage());

		if (!_timeTreeNodes[0]->isEvalAvailable(EVAL_EXPONENTIAL))
			return error.set("Linear term does not allow exponential evaluation"+getNewLineDebugMessage());

		if (!_timeTreeNodes[1]->isEvalAvailable(EVAL_TENDENCIES))
			return error.set("Remainder term does not allow tendencies evaluation"+getNewLineDebugMessage());
		
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
		TimeTree_Node_Base::registerTimeStepperEval(
				i_evalType,
				o_timeStepper
			);
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

		// Copy L() into phi0(), assuming c1 is always 0
		for (size_t i = 0; i < _rkNumStages - 1; i++) {
			_timeTreeNodes.push_back(_timeTreeNodes[0]->getInstanceCopy());
			_timeTreeNodes.push_back(_timeTreeNodes[0]->getInstanceCopy());
		}
		// Closure phi
		_timeTreeNodes.push_back(_timeTreeNodes[0]->getInstanceCopy());
		_evalFuns.resize(_timeTreeNodes.size());

		// Setup L and N
		_timeTreeNodes[0]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[0]);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);
		
		_timeTreeNodes[1]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[1]);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[1]);
		
		// Add phi functions
		for (size_t i = LRK_PHI_IDX; i < _timeTreeNodes.size(); i++) {
			_timeTreeNodes[i]->setupByKeyValue("ExpIntegrationFunction", "phi0");
			_timeTreeNodes[i]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EXPONENTIAL, &_evalFuns[i]);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[i]);
		}			

		// Setup temporary buffers: tmp, pi and ki
		_tmpDataContainer.resize(1 + _rkNumStages * 2);
		
		for (size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();

		tmp = _tmpDataContainer[0];
		p = &_tmpDataContainer[1];
		k = &_tmpDataContainer[1 + _rkNumStages];

		return true;
	}

	void clear() override
	{
		TimeTree_Node_InteriorHelper::clear();
	}

	std::shared_ptr<TimeTree_Node_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<TimeTree_Node_Base>(new LawsonRungeKutta(*this));
	}

	bool setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		// L() and N()
		_timeTreeNodes[0]->setTimeStepSize(i_dt);
		_timeTreeNodes[1]->setTimeStepSize(i_dt);

		// phi0(+/-cihL)
		for (size_t i = 0; i < _rkNumStages - 1; i++) {
			_timeTreeNodes[LRK_PHI_IDX + 2*i]->setTimeStepSize(_c[i]*i_dt);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[LRK_PHI_IDX + 2*i]);

			_timeTreeNodes[LRK_PHI_IDX + 2*i + 1]->setTimeStepSize(-1.0*_c[i]*i_dt);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[LRK_PHI_IDX + 2*i + 1]);
		}

		// Closure phi
		_timeTreeNodes.back()->setTimeStepSize(i_dt);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes.back());

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
			case LRK1:	return _eval_timeIntegration_LRK1(i_U, o_U, i_simulationTime);
			case LRK2_MIDPOINT:	return _eval_timeIntegration_LRK2_Midpoint(i_U, o_U, i_simulationTime);
			case LRK3_HEUN:	return _eval_timeIntegration_LRK3_Heun3(i_U, o_U, i_simulationTime);
			case LRK4_CLASSICAL:	return _eval_timeIntegration_LRK4_Classical(i_U, o_U, i_simulationTime);
			default: return error.set("Wrong RK evaluation");
		}
	}

private:
	/*
	 *	Lawson-Euler: u_{n+1} = phi0(hL) (u_n + h * N(u_n))
	 */
	bool _eval_timeIntegration_LRK1(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		/* 
		 * 0 | 0
		 * ------
		 *   | 1 
		 */

		// u_n + h * N(u_n)
		evalTimeStepper(1, i_U, *k[0], i_simulationTime);
		tmp->op_setVectorPlusScalarMulVector(i_U, _timestepSize, *k[0]);

		// phi0(hL) [u_n + h * N(u_n)]
		evalTimeStepper(LRK_PHI_IDX, *tmp, o_U, i_simulationTime);

		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);
		return true;
	}

private:
	bool _eval_timeIntegration_LRK2_Midpoint(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		/*
		 *  0  |
		 * 1/2 | 1/2
		 * --------------
		 *     |  0    1 
		 */

		double a2[1] = {0.5};
		double b[2] = {0.0, 1.0};
		/*
		 *	STAGE 1
		 *	p1 = phi(c1 == 0) yn
		 *	k1 = phi(-c1 == 0) N(p1)
		 */
		p[0]->op_setVector(i_U);
		evalTimeStepper(1, *p[0], *k[0], i_simulationTime);
		/*
		 *	STAGE 2
		 *	p2 = phi(c2hL) {yn + h*sum(aij)k1}
		 *	k2 = phi(-c2hL) N(p2)
		 */
		tmp->op_setVectorPlusScalarMulVector(i_U, _dt*a2[0], *k[0]);
		evalTimeStepper(LRK_PHI_IDX, *tmp, *p[1], i_simulationTime + _c[0]*_dt);
		evalTimeStepper(1, *p[1], *tmp, i_simulationTime + _c[0]*_dt);
		evalTimeStepper(LRK_PHI_IDX + 1, *tmp, *k[1], i_simulationTime + _c[0]*_dt);
		/*
		 *	Closure
		 */
		tmp->op_setVectorPlusScalarMulVector(i_U, _timestepSize*b[0], *k[0]);
		tmp->op_addScalarMulVector(_timestepSize*b[1], *k[1]);

		evalTimeStepper(LRK_PHI_IDX + 2, *tmp, o_U, i_simulationTime + _dt);

		return true;
	}

private:
	bool _eval_timeIntegration_LRK3_Heun3(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		/*
		 * c     a
		 * 0   |
		 * 1/3 | 1/3
		 * 2/3 | 0    2/3
		 * --------------
		 *     | 1/4  0   3/4
		 */
		double a2[1] = {1.0/3.0};
		double a3[2] = {0.0, 2.0/3.0};
		double b[3] = {1.0/4.0, 0.0, 3.0/4.0};
		/*
		 *	STAGE 1
		 *	p1 = phi(c1 == 0) yn
		 *	k1 = phi(-c1 == 0) N(p1)
		 */
		p[0]->op_setVector(i_U);
		evalTimeStepper(1, *p[0], *k[0], i_simulationTime);
		/*
		 *	STAGE 2
		 *	p2 = phi(c2hL) {yn + h*sum(aij)k1}
		 *	k2 = phi(-c2hL) N(p2)
		 */
		tmp->op_setVectorPlusScalarMulVector(i_U, _dt*a2[0], *k[0]);
		evalTimeStepper(LRK_PHI_IDX, *tmp, *p[1], i_simulationTime + _c[0]*_dt);
		evalTimeStepper(1, *p[1], *tmp, i_simulationTime + _c[0]*_dt);
		evalTimeStepper(LRK_PHI_IDX + 1, *tmp, *k[1], i_simulationTime + _c[0]*_dt);
		/*
		 *	STAGE 3
		 *	p3 = phi(c3hL) {yn + h*sum(aij)k2}
		 *	k3 = phi(-c3hL) N(p3)
		 */
		tmp->op_setVectorPlusScalarMulVector(i_U, _dt*a3[1], *k[1]);
		evalTimeStepper(LRK_PHI_IDX + 2, *tmp, *p[2], i_simulationTime + _c[1]*_dt);
		evalTimeStepper(1, *p[2], *tmp, i_simulationTime + _c[1]*_dt);
		evalTimeStepper(LRK_PHI_IDX + 3, *tmp, *k[2], i_simulationTime + _c[1]*_dt);
		/*
		 *	Closure
		 */
		tmp->op_setVectorPlusScalarMulVector(i_U, _timestepSize*b[0], *k[0]);
		tmp->op_addScalarMulVector(_timestepSize*b[1], *k[1]);
		tmp->op_addScalarMulVector(_timestepSize*b[2], *k[2]);

		evalTimeStepper(LRK_PHI_IDX + 4, *tmp, o_U, i_simulationTime + _dt);

		return true;
	}

private:
	bool _eval_timeIntegration_LRK4_Classical(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		/*
		 * c     a
		 * 0   |
		 * 1/2 | 1/2
		 * 1/2 | 0    1/2
		 * 1   | 0    0    1
		 * --------------
		 *     | 1/6  1/3  1/3  1/6
		 */
		double a2[1] = {0.5};
		double a3[2] = {0.0, 0.5};
		double a4[3] = {0.0, 0.0, 1.0};
		double b[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
		/*
		 *	STAGE 1
		 *	p1 = phi(c1 == 0) yn
		 *	k1 = phi(-c1 == 0) N(p1)
		 */
		p[0]->op_setVector(i_U);
		evalTimeStepper(1, *p[0], *k[0], i_simulationTime);
		/*
		 *	STAGE 2
		 *	p2 = phi(c2hL) {yn + h*sum(aij)k1}
		 *	k2 = phi(-c2hL) N(p2)
		 */
		tmp->op_setVectorPlusScalarMulVector(i_U, _dt*a2[0], *k[0]);
		evalTimeStepper(LRK_PHI_IDX, *tmp, *p[1], i_simulationTime + _c[0]*_dt);
		evalTimeStepper(1, *p[1], *tmp, i_simulationTime + _c[0]*_dt);
		evalTimeStepper(LRK_PHI_IDX + 1, *tmp, *k[1], i_simulationTime + _c[0]*_dt);
		/*
		 *	STAGE 3
		 *	p3 = phi(c3hL) {yn + h*sum(aij)k2}
		 *	k3 = phi(-c3hL) N(p3)
		 */
		tmp->op_setVectorPlusScalarMulVector(i_U, _dt*a3[1], *k[1]);
		evalTimeStepper(LRK_PHI_IDX + 2, *tmp, *p[2], i_simulationTime + _c[1]*_dt);
		evalTimeStepper(1, *p[2], *tmp, i_simulationTime + _c[1]*_dt);
		evalTimeStepper(LRK_PHI_IDX + 3, *tmp, *k[2], i_simulationTime + _c[1]*_dt);
		/*
		 *	STAGE 4
		 *	p4 = phi(c4hL) {yn + h*sum(aij)k3}
		 *	k4 = phi(c4hL) N(p4)
		 */
		tmp->op_setVectorPlusScalarMulVector(i_U, _dt*a4[2], *k[2]);
		evalTimeStepper(LRK_PHI_IDX + 4, *tmp, *p[3], i_simulationTime + _c[2]*_dt);
		evalTimeStepper(1, *p[3], *tmp, i_simulationTime + _c[2]*_dt);
		evalTimeStepper(LRK_PHI_IDX + 5, *tmp, *k[3], i_simulationTime + _c[2]*_dt);
		/*
		 *	Closure
		 */
		tmp->op_setVectorPlusScalarMulVector(i_U, _timestepSize*b[0], *k[0]);
		tmp->op_addScalarMulVector(_timestepSize*b[1], *k[1]);
		tmp->op_addScalarMulVector(_timestepSize*b[2], *k[2]);
		tmp->op_addScalarMulVector(_timestepSize*b[3], *k[3]);

		evalTimeStepper(LRK_PHI_IDX + 6, *tmp, o_U, i_simulationTime + _dt);

		return true;
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "LRK(" << std::endl;
		std::cout << newPrefix << "  order: " << _order << std::endl;
		std::cout << newPrefix << "  method: " << _method << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}}}

#endif
