#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_SLETDRK_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_SLETDRK_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>
#include <vector>
#include <string>


namespace sweet {
namespace TimeTree {
namespace InteriorNodes {


class SLETDRK	:
	public TimeTree_Node_InteriorHelper<SLETDRK>
{
private:
	// Order of Runge-Kutta method
	int _order;

	/*!
	 * Particular version of the method
	 * "default": phi0(DtL/2)[phi0(DtL/2)U]_* (real 2nd order method)
	 * "deprecated": phi0(DtL)[U]_*  ("2nd" order method)
	 */
	std::string _version;

	// Number of different phi variants we require (phi0, phi1, ...)
	int _numPhiVariants;
	int _numPsiVariants;

	sweet::Data::GenericContainer::Base* _pos_D;

public:
	SLETDRK()	:
		_order(-1),
		_version("default"),
		_numPhiVariants(-1),
		_numPsiVariants(-1)
	{
		setEvalAvailable(EVAL_INTEGRATION);
	}

	~SLETDRK()
	{
		clear();
	}

	SLETDRK(
			const SLETDRK &i_src
	)	:
		TimeTree_Node_InteriorHelper<SLETDRK>(i_src)
	{
		_order = i_src._order;
		_version = i_src._version;
		_numPhiVariants = i_src._numPhiVariants;
		_numPsiVariants = i_src._numPsiVariants;
	}

	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("sletdrk");
		retval.push_back("SLETDRK");
		return retval;
	}

	bool _setupArgumentInternals()
	{
		if (_timeTreeNodes.size() != 3)
			return error.set("We require three terms (linear + nonlinear advection + nonlinear remainder) for this time stepper!"+getNewLineDebugMessage());

		if (_order != 1 && _order != 2)
			return error.set("Only order 1 and 2 supported for SLETDRK");

		if (_version != "default" && _version != "deprecated")
			return error.set("Unknown method's version '"+_version+"'");

		if (_order == 1)
		{
			_numPhiVariants = 2; // phi0(DtL) and phi0(DtL/2)
			_numPsiVariants = 1; // psi1(DtL)
		}
		else if (_order == 2)
		{
			_numPhiVariants = 2; // phi0(DtL) and phi0(DtL/2)
			_numPsiVariants = 2; // psi1(DtL) and psi2(DtL)
		}

		// L,NA,NR --> NA,NR,L
		_timeTreeNodes[0].swap(_timeTreeNodes[1]);
		_timeTreeNodes[1].swap(_timeTreeNodes[2]);

		// We abuse the _timeTreeNodes to store also the different variants
		_timeTreeNodes.resize(2+_numPhiVariants+_numPsiVariants);

		/*
		 * _timeTreeNodes[0]:	 Nonlinear advection
		 * _timeTreeNodes[1]:	 Nonlinear remainder
		 * _timeTreeNodes[2]:	 phi0(...*L)
		 * _timeTreeNodes[3]:	 psi1(...*L)
		 * _timeTreeNodes[4]:	 psi2(...*L)
		 */

		for (int i = 3; i < 2+_numPhiVariants+_numPsiVariants; i++)
		{
			_timeTreeNodes[i] = _timeTreeNodes[2]->getInstanceCopy();
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[i]);
		}

		if (_order == 1)
		{
			SWEET_ASSERT(_timeTreeNodes.size() == 5);
			_timeTreeNodes[2]->setupByKeyValue("ExpIntegrationFunction", "phi0");
			_timeTreeNodes[3]->setupByKeyValue("ExpIntegrationFunction", "phi0");
			_timeTreeNodes[4]->setupByKeyValue("ExpIntegrationFunction", "psi1");
		}
		else if (_order == 2)
		{
			SWEET_ASSERT(_timeTreeNodes.size() == 6);
			_timeTreeNodes[2]->setupByKeyValue("ExpIntegrationFunction", "phi0");
			_timeTreeNodes[3]->setupByKeyValue("ExpIntegrationFunction", "phi0");
			_timeTreeNodes[4]->setupByKeyValue("ExpIntegrationFunction", "psi1");
			_timeTreeNodes[5]->setupByKeyValue("ExpIntegrationFunction", "psi2");
		}

		for (auto &i : _timeTreeNodes)
		{
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*i);
		}

		return true;
	}

	bool setTimeStepSize(double i_dt)	override
	{
		_timestepSize = i_dt;

		_timeTreeNodes[0]->setTimeStepSize(_timestepSize);

		if (_order == 1)
		{
			_timeTreeNodes[2]->setTimeStepSize(_dt);
			_timeTreeNodes[3]->setTimeStepSize(.5*_dt);
			_timeTreeNodes[4]->setTimeStepSize(_dt);
			return true;
		}
		else if (_order == 2)
		{
			_timeTreeNodes[2]->setTimeStepSize(_dt);
			_timeTreeNodes[3]->setTimeStepSize(.5*_dt);
			_timeTreeNodes[4]->setTimeStepSize(_dt);
			_timeTreeNodes[5]->setTimeStepSize(_dt);
			return true;
		}

		SWEETErrorFatal("Internal error");
		return false;
	}


	virtual
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
				if (a->key == "version" || a->key == "v")
				{
					_version = a->value;
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


	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper
	) override
	{
		_evalFuns.resize(_timeTreeNodes.size());

		_timeTreeNodes[0]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_SEMI_LAGRANGIAN);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);

		_timeTreeNodes[1]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[1]);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[1]);

		for (std::size_t i = 2; i < _timeTreeNodes.size(); i++)
		{
			_timeTreeNodes[i]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_EXPONENTIAL, &_evalFuns[i]);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[i]);
		}

		// default setup
		TimeTree_Node_Base::registerTimeStepperEval(
				i_evalType,
				o_timeStepper
			);
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

		/*
		 * Setup temporary buffers
		 */
		if (_order == 1)
		{
			if (_version == "deprecated")
				_tmpDataContainer.resize(6);
			else
				_tmpDataContainer.resize(9);
		}
		else if (_order == 2)
			_tmpDataContainer.resize(9);

		_tmpDataContainer[0] = i_deTermConfig.getNewDataContainerInstance(Data::GenericContainer::Base::DATA_SEMI_LAGRANGIAN_POSITIONS);
		for (std::size_t i = 1; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();

		// Departure points
		_pos_D = i_deTermConfig.getNewDataContainerInstance(Data::GenericContainer::Base::DATA_SEMI_LAGRANGIAN_POSITIONS);

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
		return std::shared_ptr<TimeTree_Node_Base>(new SLETDRK(*this));
	}


	bool _eval_integration(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)	override
	{

		bool split_phi0 = (_version == "default");

		switch(_order)
		{
		case 1:		return _eval_timeIntegration_SLETDRK1(i_U, o_U, i_simulationTime, split_phi0);
		case 2:		return _eval_timeIntegration_SLETDRK2(i_U, o_U, i_simulationTime, split_phi0);
		default: return error.set("Invalid order provided");
		}
		return false;
	}

private:
	inline
	void _evalNL(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(1, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalPhi0(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(2, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalPhi0_halfTimestep(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(3, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalPsi1(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(4, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalPsi2(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(5, i_U, o_U, i_simulationTime);
	}

private:
	bool _compute_departurePoints(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		////// Departure positions
		////sweet::Data::GenericContainer::Base &pos_D = *_tmpDataContainer[0];

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

		sweet::Data::GenericContainer::Base &aux = *_tmpDataContainer[0];
		_timeTreeNodes[0]->evalNA_departurePoints(
			states,	// IN: States
			_timestepSize,	// IN: timestep size
			aux	// OUT: departure position
		);

		_pos_D->op_setVector(aux);

		return true;
	}

private:
	bool _eval_timeIntegration_SLETDRK1(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime,
			bool split_phi0 = true
	)
	{

		/*************************************************************************************************
		 * Step 1) Compute departure points
		 *************************************************************************************************/
		_compute_departurePoints(i_U, o_U, i_simulationTime);
		ERROR_CHECK_COND_RETURN_BOOLEAN(*this);

		/*
		 * Backup solution for next time step
		 */
		_tmpDataContainer[1]->op_setVector(i_U);

		if (split_phi0)
		{
			/*
			 * U_{1} = \phi_{0}( \Delta t L / 2 ) [ \phi_{0}( \Delta t L / 2 ) U_{0} ]_* + \Delta t \phi_{0}( \Delta t L )[ \psi_{1}(\Delta tL) N(U_{0}) ]_*
			 *                                                                                                                                 NNNNNNNN
			 */
			sweet::Data::GenericContainer::Base &FU = *_tmpDataContainer[2];
			_evalNL(i_U, FU, i_simulationTime);

			/*
			 * U_{1} = \phi_{0}( \Delta t L / 2 ) [ \phi_{0}( \Delta t L / 2 ) U_{0} ]_* + \Delta t \phi_{0}( \Delta t L )[ \psi_{1}(\Delta tL) N(U_{0}) ]_*
			 *                                                                                                              NNNNNNNNNNNNNNNNNNN =======
			 */
			sweet::Data::GenericContainer::Base &psi1_FU = *_tmpDataContainer[3];
			_evalPsi1(FU, psi1_FU, i_simulationTime);

			/*
			 * U_{1} = \phi_{0}( \Delta t L / 2 ) [ \phi_{0}( \Delta t L / 2 ) U_{0} ]_* + \Delta t \phi_{0}( \Delta t L )[ \psi_{1}(\Delta tL) N(U_{0}) ]_*
			 *                                      ================================                 =====================   =================== =======  NN
			 */
			sweet::Data::GenericContainer::Base &psi1_FU_D = *_tmpDataContainer[4];
			_timeTreeNodes[0]->evalNA_interpolate(psi1_FU, *_pos_D, psi1_FU_D);

			/*
			 * U_{1} = \phi_{0}( \Delta t L / 2 ) [ \phi_{0}( \Delta t L / 2 ) U_{0} ]_* + \Delta t \phi_{0}( \Delta t L )[ \psi_{1}(\Delta tL) N(U_{0}) ]_*
			 *                                                                                      NNNNNNNNNNNNNNNNNNNNN   =================== =======
			 */
			sweet::Data::GenericContainer::Base &phi0_psi1_FU_D = *_tmpDataContainer[5];
			_evalPhi0(psi1_FU_D, phi0_psi1_FU_D, i_simulationTime);

			/*
			 * U_{1} = \phi_{0}( \Delta t L / 2 ) [ \phi_{0}( \Delta t L / 2 ) U_{0} ]_* + \Delta t \phi_{0}( \Delta t L )[ \psi_{1}(\Delta tL) N(U_{0}) ]_*
			 *                                      NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN                  =====================   =================== =======  ==
			 */
			sweet::Data::GenericContainer::Base &phi0_U = *_tmpDataContainer[6];
			_evalPhi0_halfTimestep(i_U, phi0_U, i_simulationTime);

			/*
			 * U_{1} = \phi_{0}( \Delta t L / 2 ) [ \phi_{0}( \Delta t L / 2 ) U_{0} ]_* + \Delta t \phi_{0}( \Delta t L )[ \psi_{1}(\Delta tL) N(U_{0}) ]_*
			 *                                      ================================  NN             =====================   =================== ======= ==
			 */
			sweet::Data::GenericContainer::Base &phi0_U_D = *_tmpDataContainer[7];
			_timeTreeNodes[0]->evalNA_interpolate(phi0_U, *_pos_D, phi0_U_D);

			/*
			 * U_{1} = \phi_{0}( \Delta t L / 2 ) [ \phi_{0}( \Delta t L / 2 ) U_{0} ]_* + \Delta t \phi_{0}( \Delta t L )[ \psi_{1}(\Delta tL) N(U_{0}) ]_*
			 *         NNNNNNNNNNNNNNNNNNNNNNNNN    ================================  ==            =====================   =================== =======   ==
			 */
			sweet::Data::GenericContainer::Base &phi0_phi0_U_D = *_tmpDataContainer[8];
			_evalPhi0_halfTimestep(phi0_U_D, phi0_phi0_U_D, i_simulationTime);

			/*
			 * U_{1} = \phi_{0}( \Delta t L / 2 ) [ \phi_{0}( \Delta t L / 2 ) U_{0} ]_* + \Delta t \phi_{0}( \Delta t L )[ \psi_{1}(\Delta tL) N(U_{0}) ]_*
			 *         =========================    ================================  == NNNNNNNNNN =====================   =================== =======   ==
			 */
			o_U.op_setVectorPlusScalarMulVector(phi0_phi0_U_D, _timestepSize, phi0_psi1_FU_D);
		}
		else
		{
			/*
			 * U_{1} = \phi_{0}( \Delta t L ) [ U_{0} + \Delta t \psi_{1}(\Delta tL) N(U_{0}) ]_*
			 *                                                                       NNNNNNNN
			 */
			sweet::Data::GenericContainer::Base &FU = *_tmpDataContainer[2];
			_evalNL(i_U, FU, i_simulationTime);

			/*
			 * U_{1} = \phi_{0}( \Delta t L ) [ U_{0} + \Delta t \psi_{1}(\Delta tL) N(U_{0}) ]_*
			 *                                                   NNNNNNNNNNNNNNNNNNN ========
			 */
			sweet::Data::GenericContainer::Base &psi1_FU = *_tmpDataContainer[3];
			_evalPsi1(FU, psi1_FU, i_simulationTime);

			/*
			 * U_{1} = \phi_{0}( \Delta t L ) [ U_{0} + \Delta t \psi_{1}(\Delta tL) N(U_{0}) ]_*
			 *                                  NNNNNNNNNNNNNNNN =================== ========
			 */
			sweet::Data::GenericContainer::Base &U_NL = *_tmpDataContainer[4];
			U_NL.op_setVectorPlusScalarMulVector(i_U, _timestepSize, psi1_FU);

			/*
			 * U_{1} = \phi_{0}( \Delta t L ) [ U_{0} + \Delta t \psi_{1}(\Delta tL) N(U_{0}) ]_*
			 *                                  ================ =================== ========  NN
			 */
			sweet::Data::GenericContainer::Base &U_NL_D = *_tmpDataContainer[5];
			_timeTreeNodes[0]->evalNA_interpolate(U_NL, *_pos_D, U_NL_D);

			/*
			 * U_{1} = \phi_{0}( \Delta t L ) [ U_{0} + \Delta t \psi_{1}(\Delta tL) N(U_{0}) ]_*
			 *         NNNNNNNNNNNNNNNNNNNNNN   ================ =================== ========  ==
			 */
			_evalPhi0(U_NL_D, o_U, i_simulationTime);
		}

		return true;
	}


private:
	bool _eval_timeIntegration_SLETDRK2(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime,
			bool split_phi0 = true
	)
	{

		/*
		 * U_{1} = U_{SL-ETD1RK} + \Delta t \phi_{0}( \Delta t L ) [ \psi_{2}(\Delta tL) N(U_{SL-ETD1RK}) - ( \psi_{2}(\Delta tL) N(U_{0}) )_* ]
		 *         NNNNNNNNNNNNN
		 */
		sweet::Data::GenericContainer::Base &U1 = *_tmpDataContainer[2];
		_eval_timeIntegration_SLETDRK1(i_U, U1, i_simulationTime, split_phi0);

		/*
		 * U_{1} = U_{SL-ETD1RK} + \Delta t \phi_{0}( \Delta t L ) [ \psi_{2}(\Delta tL) N(U_{SL-ETD1RK}) - ( \psi_{2}(\Delta tL) N(U_{0}) )_* ]
		 *         =============                                                                                                  NNNNNNNN
		 */
		sweet::Data::GenericContainer::Base &FU = *_tmpDataContainer[3];
		_evalNL(i_U, FU, i_simulationTime);

		/*
		 * U_{1} = U_{SL-ETD1RK} + \Delta t \phi_{0}( \Delta t L ) [ \psi_{2}(\Delta tL) N(U_{SL-ETD1RK}) - ( \psi_{2}(\Delta tL) N(U_{0}) )_* ]
		 *         =============                                                         NNNNNNNNNNNNNNNN                         ========
		 */
		sweet::Data::GenericContainer::Base &FU1 = *_tmpDataContainer[4];
		_evalNL(U1, FU1, i_simulationTime);

		/*
		 * U_{1} = U_{SL-ETD1RK} + \Delta t \phi_{0}( \Delta t L ) [ \psi_{2}(\Delta tL) N(U_{SL-ETD1RK}) - ( \psi_{2}(\Delta tL) N(U_{0}) )_* ]
		 *         =============                                                         NNNNNNNNNNNNNNNN     =================== NNNNNNNN
		 */
		sweet::Data::GenericContainer::Base &psi2_FU = *_tmpDataContainer[5];
		_evalPsi2(FU, psi2_FU, i_simulationTime);

		/*
		 * U_{1} = U_{SL-ETD1RK} + \Delta t \phi_{0}( \Delta t L ) [ \psi_{2}(\Delta tL) N(U_{SL-ETD1RK}) - ( \psi_{2}(\Delta tL) N(U_{0}) )_* ]
		 *         =============                                     NNNNNNNNNNNNNNNNNNN ================     =================== ========
		 */
		sweet::Data::GenericContainer::Base &psi2_FU1 = *_tmpDataContainer[6];
		_evalPsi2(FU1, psi2_FU1, i_simulationTime);

		/*
		 * U_{1} = U_{SL-ETD1RK} + \Delta t \phi_{0}( \Delta t L ) [ \psi_{2}(\Delta tL) N(U_{SL-ETD1RK}) - ( \psi_{2}(\Delta tL) N(U_{0}) )_* ]
		 *         =============                                     =================== ================     =================== ========  NN
		 */
		sweet::Data::GenericContainer::Base &psi2_FU_D = *_tmpDataContainer[7];
		_timeTreeNodes[0]->evalNA_interpolate(psi2_FU, *_pos_D, psi2_FU_D);

		/*
		 * U_{1} = U_{SL-ETD1RK} + \Delta t \phi_{0}( \Delta t L ) [ \psi_{2}(\Delta tL) N(U_{SL-ETD1RK}) - ( \psi_{2}(\Delta tL) N(U_{0}) )_* ]
		 *         =============            NNNNNNNNNNNNNNNNNNNNNN   =================== ================     =================== ========  ==
		 */
		sweet::Data::GenericContainer::Base &phi0_U = *_tmpDataContainer[8];
		psi2_FU1.op_subVector(psi2_FU_D);
		_evalPhi0(psi2_FU1, phi0_U, i_simulationTime);

		/*
		 * U_{1} = U_{SL-ETD1RK} + \Delta t \phi_{0}( \Delta t L ) [ \psi_{2}(\Delta tL) N(U_{SL-ETD1RK}) - ( \psi_{2}(\Delta tL) N(U_{0}) )_* ]
		 *         =============   NNNNNNNN ======================   =================== ================     =================== ========  ==
		 */
		o_U.op_setVectorPlusScalarMulVector(U1, _dt, phi0_U);

		return true;
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "SLETDRK(" << std::endl;
		std::cout << newPrefix << "  order: " << _order << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}}}

#endif
