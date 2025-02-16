#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_SETTLS_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_SETTLS_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>
#include <vector>
#include <string>



namespace sweet {
namespace TimeTree {
namespace InteriorNodes {


/*!
 * \brief SETTLS Semi-Lagrangian two level time integration scheme
 *
 * Following "Hortal2002", we write for the variables \f$X\f$
 *
 * \f[
 *	\frac{dX}{dt}=L+N
 * \f]
 *
 * with \f$d/dt\f$ the total derivative, \f$L\f$ the linearized part and \f$N\f$
 * all terms treated as nonlinearities.
 *
 * The 1997 version discretized this equation with an averaging as
 *
 * \f[
 * \frac{X_{A}^{t+\Delta t}-X_{D}^{t}}{\Delta t}=\frac{1}{2}\left(L_{A}^{t+\Delta t}+L_{D}^{t}\right)+\frac{1}{2}\left(N_{D}^{t+\frac{\Delta t}{2}}+N_{A}^{t+\frac{\Delta t}{2}}\right).
 * \f]
 *
 * Subscript \f$A\f$ marks the arrival points which are aligned at grid
 * cells and subscript \f$D\f$ marks the departure points which typically
 * require an interpolation.
 */
class SETTLS	:
	public TimeTree_Node_InteriorHelper<SETTLS>
{
private:
	//! Order of SETTLS method
	int _order;

	//! Number of states required for this method
	int numMultiStates;

	/*!
	 * Constructor
	 */
public:
	SETTLS()	:
		_order(-1),
		numMultiStates(-1)
	{
		setEvalAvailable(EVAL_INTEGRATION);
	}

	~SETTLS()
	{
		clear();
	}

	SETTLS(
			const SETTLS &i_src
	)	:
		TimeTree_Node_InteriorHelper<SETTLS>(i_src)
	{
		_order = i_src._order;
	}

	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("settls");
		retval.push_back("SETTLS");
		return retval;
	}

	bool _setupArgumentInternals()
	{

		if (_timeTreeNodes.size() != 3)
			return error.set("We require exactly 3 DE terms for SETTLS (a linear, nonlinear advection and a nonlinear remainder one)");

		//if (_order != 1 && _order != 2)
		if (_order != 2 && _order != -2)
			return error.set("Only order 1 or 2 supported for SETTLS");

		/*
		 * Time Tree node 1: For forward Euler
		 * Time Tree node 2: For backward Euler
		 * Time Tree node 3: Nonlinear advection / semi-Lagrangian
		 * Time Tree node 4: Nonlinear remainder term
		 */
		_timeTreeNodes.resize(4);

		// NR term -> 4th term
		_timeTreeNodes[2].swap(_timeTreeNodes[3]);

		// NA term -> 3rd term
		_timeTreeNodes[1].swap(_timeTreeNodes[2]);

		// Linear term for backward Euler -> 2nd term
		_timeTreeNodes[1] = _timeTreeNodes[0]->getInstanceCopy();


		for (auto &i : _timeTreeNodes)
		{
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*i);
		}

		return true;
	}

	bool setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		SWEET_ASSERT(_timeTreeNodes.size() == 4);
		if (_order == 1)
		{
			_timeTreeNodes[0]->setTimeStepSize(_timestepSize);
			//_timeTreeNodes[1]->setTimeStepSize(_timestepSize);
			_timeTreeNodes[2]->setTimeStepSize(_timestepSize);
			//_timeTreeNodes[3]->setTimeStepSize(_timestepSize);
		}
		else if (_order == 2 || _order == -2)
		{
			_timeTreeNodes[0]->setTimeStepSize(_timestepSize*0.5);
			_timeTreeNodes[1]->setTimeStepSize(_timestepSize*0.5);
			_timeTreeNodes[2]->setTimeStepSize(_timestepSize);
			//_timeTreeNodes[3]->setTimeStepSize(_timestepSize);
		}

		return true;
	}


	/*!
	 * SETTLS splits a PDE into three parts:
	 *
	 * 1. The linear part,
	 * 2. the nonlinear advection part
	 * 3. the remainder
	 *
	 * Since the nonlinear advection is treated with a Semi-Lagrangian scheme,
	 * we do not need to specify it.
	 *
	 * Hence, we can call it by
	 *
	 * 	SETTLS(linear part, nonlinear advection, nonlinear remainder part, options...)
	 */
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
				if (_timeTreeNodes.size() >= 3)
					return error.set("A 4th timestepper was provided, but only 3 allowed!"+a->getNewLineDebugMessage());

				_timeTreeNodes.push_back(std::shared_ptr<sweet::TimeTree::TimeTree_Node_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_VALUE:
				if (_timeTreeNodes.size() >= 3)
					return error.set("A 4th timestepper was provided, but only 3 allowed!"+a->getNewLineDebugMessage());

				_timeTreeNodes.push_back(std::shared_ptr<sweet::TimeTree::TimeTree_Node_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByName(
						a->value,
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
		o_ostream << i_prefix << "InteriorNode: 'Stable Extrapolation Two-Time-Level Scheme (SETTLS)':" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Usage: SETTLS(l,na,nr,[parameter1],...)" << std::endl;
		o_ostream << i_prefix << "           l: linear term treated implicitly (Crank-Nicolson style)." << std::endl;
		o_ostream << i_prefix << "           na: nonlinear advection term (requires support of Semi-Lagrangian evaluations)." << std::endl;
		o_ostream << i_prefix << "           nr: nonlinear remainder term (requires support of tendency evaluations)." << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Parameters:" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "    - order=[int]:" << std::endl;
		o_ostream << i_prefix << "        Specify the order of the time integration method (only order 2 currently supported)" << std::endl;
		o_ostream << i_prefix << std::endl;

		return true;
	}

	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper
	) override
	{
		/*
		 * We have 4 timeTreeNodes:
		 * L tendencies, L backward Euler, NA Semi-Lagrangian, NR term
		 */
		assert(_timeTreeNodes.size() == 4);
		_evalFuns.resize(_timeTreeNodes.size());

		_timeTreeNodes[0]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[0]);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);

		_timeTreeNodes[1]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EULER_BACKWARD, &_evalFuns[1]);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[1]);

		_timeTreeNodes[2]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_SEMI_LAGRANGIAN);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[2]);

		_timeTreeNodes[3]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[3]);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[3]);

		// default setup
		TimeTree_Node_Base::registerTimeStepperEval(
				i_evalType,
				o_timeStepper
			);
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

		_timeTreeNodes[2]->evalNA_getNumStates(&numMultiStates);
		ERROR_CHECK_COND_RETURN_BOOLEAN(*_timeTreeNodes[2]);

		SWEET_ASSERT(numMultiStates > 0);

		if (_order == 2 || _order == -2)
		{
			/*
			 * 0: Positions of SL method
			 * 1: Simulation data for previous time step
			 * 2,3,... Temporary simulation data
			 */
			_tmpDataContainer.resize(10);

			_tmpDataContainer[0] = i_deTermConfig.getNewDataContainerInstance(Data::GenericContainer::Base::DATA_SEMI_LAGRANGIAN_POSITIONS);

			for (std::size_t i = 1; i < _tmpDataContainer.size(); i++)
				_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance(Data::GenericContainer::Base::DATA_SIMULATION);
		}

#if SWEET_XBRAID
		U_prev_solution = i_deTermConfig.getNewDataContainerInstance(Data::GenericContainer::Base::DATA_SIMULATION);
#endif

		return true;

	}

	void clear() override
	{
		TimeTree_Node_InteriorHelper::clear();
	}


	std::shared_ptr<TimeTree_Node_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<TimeTree_Node_Base>(new SETTLS(*this));
	}


public:
	bool _eval_integration(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)	override
	{
		switch(_order)
		{
		case 1:	return _eval_timeIntegration_SETTLS_order1(i_U, o_U, i_simulationTime);
		case 2:	return _eval_timeIntegration_SETTLS_order2(i_U, o_U, i_simulationTime);
		case -2:	return _eval_timeIntegration_SETTLS_order2_modified(i_U, o_U, i_simulationTime);
		default: return error.set("Invalid order provided");
		}
		return false;
	}

private:
	bool _eval_timeIntegration_SETTLS_order1(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		SWEETErrorFatal("TODO");

		return true;
	}

private:
	inline
	void _eval_L_tendencies(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(0, i_U, o_U, i_simulationTime);
	}
private:
	inline
	void _eval_L_backwardEuler(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(1, i_U, o_U, i_simulationTime);
	}

private:
	inline
	void _eval_NR_tendencies(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(3, i_U, o_U, i_simulationTime);
	}

private:
	bool _eval_timeIntegration_SETTLS_order2(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		// Departure positions
		sweet::Data::GenericContainer::Base &pos_D = *_tmpDataContainer[0];

		// current state
		const sweet::Data::GenericContainer::Base &X1 = i_U;

		// state at previous time step
		sweet::Data::GenericContainer::Base &X0 = *_tmpDataContainer[1];

#if SWEET_XBRAID
		X0.op_setVector(*U_prev_solution);
#endif

		if (i_simulationTime == 0)
		{
			X0.op_setVector(X1);
		}

		const sweet::Data::GenericContainer::Base* states[2];
		states[0] = &X1;
		states[1] = &X0;

		_timeTreeNodes[2]->evalNA_departurePoints(
			states,	// IN: States
			_timestepSize,	// IN: timestep size
			pos_D	// OUT: departure position
		);

		/*
		 * Step 2) Midpoint rule
		 * Put everything together with midpoint rule and solve resulting Helmholtz problem
		 */

		/*
		 * Step 2a) Compute RHS
		 * R = X1_D + 1/2 dt L1_D + dt N*
		 */

		// Compute X1_D
		sweet::Data::GenericContainer::Base &X1_D = *_tmpDataContainer[2];
		_timeTreeNodes[2]->evalNA_interpolate(X1, pos_D, X1_D);

		// Compute L1_D
		sweet::Data::GenericContainer::Base &L1 = *_tmpDataContainer[3];
		_eval_L_tendencies(X1, L1, i_simulationTime);

		// Interpolate at departure points
		sweet::Data::GenericContainer::Base &L1_D = *_tmpDataContainer[4];
		_timeTreeNodes[2]->evalNA_interpolate(L1, pos_D, L1_D);


		/*
		 * Compute N*
		 * N*(t+0.5dt) = 1/2 ([ 2*N(t) - N(t-dt) ]_D + N(t))
		 */
		sweet::Data::GenericContainer::Base &N0 = *_tmpDataContainer[5];
		sweet::Data::GenericContainer::Base &N1 = *_tmpDataContainer[6];
		_eval_NR_tendencies(X0, N0, i_simulationTime);
		_eval_NR_tendencies(X1, N1, i_simulationTime);

		// 2*N(t) - N(t-dt)
		sweet::Data::GenericContainer::Base &tmp = *_tmpDataContainer[7];
		tmp.op_setVector(N0);
		tmp.op_mulScalar(-1);
		tmp.op_addScalarMulVector(2.0, N1);

		// [ 2*N(t) - N(t-dt) ]_D
		sweet::Data::GenericContainer::Base &N_star = *_tmpDataContainer[8];
		_timeTreeNodes[2]->evalNA_interpolate(tmp, pos_D, N_star);

		N_star.op_addVector(N1);
		N_star.op_mulScalar(0.5);

		/*
		 * Compute R = X_D + 1/2 dt L_D + dt N*
		 */
		sweet::Data::GenericContainer::Base &R = *_tmpDataContainer[9];
		R.op_setVectorPlusScalarMulVector(X1_D, 0.5*_dt, L1_D);
		R.op_addScalarMulVector(_dt, N_star);

		/*
		 * Step 2b) Solve Helmholtz problem
		 * X - 1/2 dt LX = R
		 */
		sweet::Data::GenericContainer::Base &X2 = o_U;
		_eval_L_backwardEuler(R, X2, i_simulationTime);

		/*
		 * Backup solution for next time step
		 */
		X0.op_setVector(X1);

		return true;
	}


private:
	bool _eval_timeIntegration_SETTLS_order2_modified(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		// Departure positions
		sweet::Data::GenericContainer::Base &pos_D = *_tmpDataContainer[0];

		// current state
		const sweet::Data::GenericContainer::Base &X1 = i_U;

		// state at previous time step
		sweet::Data::GenericContainer::Base &X0 = *_tmpDataContainer[1];

#if SWEET_XBRAID
		X0.op_setVector(*U_prev_solution);
#endif

		if (i_simulationTime == 0)
		{
			X0.op_setVector(X1);
		}

		const sweet::Data::GenericContainer::Base* states[2];
		states[0] = &X1;
		states[1] = &X0;

		_timeTreeNodes[2]->evalNA_departurePoints(
			states,	// IN: States
			_timestepSize,	// IN: timestep size
			pos_D	// OUT: departure position
		);

		/*
		 * Step 2) Midpoint rule
		 * Put everything together with midpoint rule and solve resulting Helmholtz problem
		 */

		/*
		 * Step 2a) Compute RHS
		 * R = X1_D + 1/2 dt L1_D + dt N*
		 */

		// Compute X1_D
		sweet::Data::GenericContainer::Base &X1_D = *_tmpDataContainer[2];
		_timeTreeNodes[2]->evalNA_interpolate(X1, pos_D, X1_D);

		// Compute L1_D
		sweet::Data::GenericContainer::Base &L1 = *_tmpDataContainer[3];
		_eval_L_tendencies(X1, L1, i_simulationTime);

		// Interpolate at departure points
		sweet::Data::GenericContainer::Base &L1_D = *_tmpDataContainer[4];
		_timeTreeNodes[2]->evalNA_interpolate(L1, pos_D, L1_D);


		/*
		 * Compute N*
		 * N*(t+0.5dt) = 1/2 ([ 2*N(t) - N(t-dt) ]_D + N(t))
		 */
		sweet::Data::GenericContainer::Base &N0 = *_tmpDataContainer[5];
		sweet::Data::GenericContainer::Base &N1 = *_tmpDataContainer[6];
		_eval_NR_tendencies(X0, N0, i_simulationTime);
		_eval_NR_tendencies(X1, N1, i_simulationTime);

		// 2*N(t) - N(t-dt)
		sweet::Data::GenericContainer::Base &tmp = *_tmpDataContainer[7];
		tmp.op_setVector(N0);
		tmp.op_mulScalar(-1);
		tmp.op_addScalarMulVector(1.0, N1);

		// [ 2*N(t) - N(t-dt) ]_D
		sweet::Data::GenericContainer::Base &N_star = *_tmpDataContainer[8];
		_timeTreeNodes[2]->evalNA_interpolate(tmp, pos_D, N_star);

		//////N_star.op_addVector(N1);
		N_star.op_mulScalar(0.5);

		/*
		 * Compute R = X_D + 1/2 dt L_D + dt N*
		 */
		sweet::Data::GenericContainer::Base &R = *_tmpDataContainer[9];
		R.op_setVectorPlusScalarMulVector(X1_D, 0.5*_dt, L1_D);
		R.op_addScalarMulVector(_dt, N_star);

		/*
		 * Step 2b) Solve Helmholtz problem
		 * X - 1/2 dt LX = R
		 */
		sweet::Data::GenericContainer::Base &X2 = o_U;
		_eval_L_backwardEuler(R, X2, i_simulationTime);

		/*
		 * Backup solution for next time step
		 */
		X0.op_setVector(X1);

		return true;
	}




	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "SETTLS(" << std::endl;
		std::cout << newPrefix << "  order: " << _order << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}}}

#endif
