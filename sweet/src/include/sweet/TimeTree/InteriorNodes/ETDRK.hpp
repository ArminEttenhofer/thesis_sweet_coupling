#ifndef INCLUDE_SWEET_TIMETREE_INTERIORNODES_ETDRK_HPP
#define INCLUDE_SWEET_TIMETREE_INTERIORNODES_ETDRK_HPP

#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/TimeTree/TimeTree_Node_InteriorHelper.hpp>
#include <vector>
#include <string>


namespace sweet {
namespace TimeTree {
namespace InteriorNodes {


class ETDRK	:
	public TimeTree_Node_InteriorHelper<ETDRK>
{
private:
	// Order of Runge-Kutta method
	int _order;

	// Number of different phi variants we require (phi0, phi1, ...)
	int _numPhiVariants;

public:
	ETDRK()	:
		_order(-1),
		_numPhiVariants(-1)
	{
		setEvalAvailable(EVAL_INTEGRATION);
	}

	~ETDRK()
	{
		clear();
	}

	ETDRK(
			const ETDRK &i_src
	)	:
		TimeTree_Node_InteriorHelper<ETDRK>(i_src)
	{
		_order = i_src._order;
		_numPhiVariants = i_src._numPhiVariants;
	}

	const std::vector<std::string>
	getNodeNames()	override
	{
		std::vector<std::string> retval;
		retval.push_back("etdrk");
		retval.push_back("ETDRK");
		return retval;
	}

	bool _setupArgumentInternals()
	{
		if (_timeTreeNodes.size() != 2)
			return error.set("We require two terms (linear + nonlinear) for this time stepper!"+getNewLineDebugMessage());

		if (_order != 1 && _order != 2 && _order != 4)
			return error.set("Only order 1, 2, and 4 supported for ETDRK");

		if (_order == 1)
		{
			_numPhiVariants = 2;
		}
		else if (_order == 2)
		{
			_numPhiVariants = 3;
		}
		else if (_order == 4)
		{
			_numPhiVariants = 7;
		}

		_timeTreeNodes[0].swap(_timeTreeNodes[1]);

		// We abuse the _timeTreeNodes to store also the different variants
		_timeTreeNodes.resize(1+_numPhiVariants);

		/*
		 * _timeTreeNodes[0]:	 Nonlinear part
		 * _timeTreeNodes[1]:	 phi0(...*L)
		 * _timeTreeNodes[2]:	 phi1(...*L)
		 * _timeTreeNodes[3]:	 phi2(...*L)
		 * _timeTreeNodes[4]:	 phi0(...*L)
		 * _timeTreeNodes[5]:	 ups1(...*L)
		 * _timeTreeNodes[6]:	 ups2(...*L)
		 * _timeTreeNodes[7]:	 ups3(...*L)
		 */

		for (int i = 2; i < 1+_numPhiVariants; i++)
		{
			_timeTreeNodes[i] = _timeTreeNodes[1]->getInstanceCopy();
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[i]);
		}

		if (_order == 1)
		{
			SWEET_ASSERT(_timeTreeNodes.size() == 3);
			_timeTreeNodes[1]->setupByKeyValue("ExpIntegrationFunction", "phi0");
			_timeTreeNodes[2]->setupByKeyValue("ExpIntegrationFunction", "phi1");
		}
		else if (_order == 2)
		{
			SWEET_ASSERT(_timeTreeNodes.size() == 4);
			_timeTreeNodes[1]->setupByKeyValue("ExpIntegrationFunction", "phi0");
			_timeTreeNodes[2]->setupByKeyValue("ExpIntegrationFunction", "phi1");
			_timeTreeNodes[3]->setupByKeyValue("ExpIntegrationFunction", "phi2");
		}
		else if (_order == 4)
		{
			SWEET_ASSERT(_timeTreeNodes.size() == 8);
			_timeTreeNodes[1]->setupByKeyValue("ExpIntegrationFunction", "phi0");
			_timeTreeNodes[2]->setupByKeyValue("ExpIntegrationFunction", "phi1");
			_timeTreeNodes[3]->setupByKeyValue("ExpIntegrationFunction", "phi2");

			_timeTreeNodes[4]->setupByKeyValue("ExpIntegrationFunction", "phi0");
			_timeTreeNodes[5]->setupByKeyValue("ExpIntegrationFunction", "ups1");
			_timeTreeNodes[6]->setupByKeyValue("ExpIntegrationFunction", "ups2");
			_timeTreeNodes[7]->setupByKeyValue("ExpIntegrationFunction", "ups3");
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
			_timeTreeNodes[1]->setTimeStepSize(_dt);
			_timeTreeNodes[2]->setTimeStepSize(_dt);
			return true;
		}
		else if (_order == 2)
		{
			_timeTreeNodes[1]->setTimeStepSize(_dt);
			_timeTreeNodes[2]->setTimeStepSize(_dt);
			_timeTreeNodes[3]->setTimeStepSize(_dt);
			return true;
		}
		else if (_order == 4)
		{
			_timeTreeNodes[1]->setTimeStepSize(0.5*_dt);
			_timeTreeNodes[2]->setTimeStepSize(0.5*_dt);
			_timeTreeNodes[3]->setTimeStepSize(0.5*_dt);

			_timeTreeNodes[4]->setTimeStepSize(_dt);
			_timeTreeNodes[5]->setTimeStepSize(_dt);
			_timeTreeNodes[6]->setTimeStepSize(_dt);
			_timeTreeNodes[7]->setTimeStepSize(_dt);
			return true;
		}

		SWEETErrorFatal("Internal error");
		return false;
	}



	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override
	{
		o_ostream << i_prefix << "InteriorNode: 'ETDRK':" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Node name & aliases: " << _getNodeNamesAsString() << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Usage: ETDRK(l,n,[parameter1],...)" << std::endl;
		o_ostream << i_prefix << "           Exponential time differencing for Runge-Kutta." << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "  - Parameters:" << std::endl;
		o_ostream << i_prefix << std::endl;
		o_ostream << i_prefix << "    - order=[int]:" << std::endl;
		o_ostream << i_prefix << "        Specify the order of the time integration method" << std::endl;
		o_ostream << i_prefix << std::endl;

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
				if (_timeTreeNodes.size() >= 2)
					return error.set("A 3rd timestepper was provided, but only 2 allowed!"+a->getNewLineDebugMessage());

				_timeTreeNodes.push_back(std::shared_ptr<sweet::TimeTree::TimeTree_Node_Base>());

				i_tsAssemblation.assembleTimeTreeNodeByFunction(
						a->function,
						_timeTreeNodes.back()
					);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(i_tsAssemblation);
				break;

			case sweet::TimeTree::TimeTreeIR::Argument::ARG_TYPE_VALUE:
				if (_timeTreeNodes.size() >= 2)
					return error.set("A 3rd timestepper was provided, but only 2 allowed!"+a->getNewLineDebugMessage());

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


	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper
	) override
	{
		_evalFuns.resize(_timeTreeNodes.size());

		_timeTreeNodes[0]->setupConfigAndForwardTimeStepperEval(i_deTermConfig, EVAL_TENDENCIES, &_evalFuns[0]);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*_timeTreeNodes[0]);

		for (std::size_t i = 1; i < _timeTreeNodes.size(); i++)
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
			_tmpDataContainer.resize(3);
		else if (_order == 2)
			_tmpDataContainer.resize(6);
		else if (_order == 4)
			_tmpDataContainer.resize(16);

		for (std::size_t i = 0; i < _tmpDataContainer.size(); i++)
			_tmpDataContainer[i] = i_deTermConfig.getNewDataContainerInstance();

		return true;
	}

	void clear() override
	{
		TimeTree_Node_InteriorHelper::clear();
	}


	std::shared_ptr<TimeTree_Node_Base> getInstanceCopy()	override
	{
		return std::shared_ptr<TimeTree_Node_Base>(new ETDRK(*this));
	}


	bool _eval_integration(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)	override
	{
		switch(_order)
		{
		case 1:	return _eval_timeIntegration_ETDRK1(i_U, o_U, i_simulationTime);
		case 2:	return _eval_timeIntegration_ETDRK2(i_U, o_U, i_simulationTime);
		case 4:	return _eval_timeIntegration_ETDRK4(i_U, o_U, i_simulationTime);
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
		evalTimeStepper(0, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalPhi0(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(1, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalPhi1(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(2, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalPhi2(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(3, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalUps0(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(4, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalUps1(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(5, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalUps2(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(6, i_U, o_U, i_simulationTime);
	}
	inline
	void _evalUps3(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		evalTimeStepper(7, i_U, o_U, i_simulationTime);
	}


private:
	bool _eval_timeIntegration_ETDRK1(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}	+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 *         NNNNNNNNNNNNNNNNNNNNNNNNNNNN
		 */
		sweet::Data::GenericContainer::Base &phi0_U = *_tmpDataContainer[0];
		_evalPhi0(i_U, phi0_U, i_simulationTime);

		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}	+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 *         ============================                               NNNNNNNN
		 */
		sweet::Data::GenericContainer::Base &FU = *_tmpDataContainer[1];
		_evalNL(i_U, FU, i_simulationTime);

		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}	+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 *         ============================           NNNNNNNNNNNNNNNNNNN =========
		 */
		sweet::Data::GenericContainer::Base &phi1_FU = *_tmpDataContainer[2];
		_evalPhi1(FU, phi1_FU, i_simulationTime);

		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}	+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 * NNNNN   ============================           =================== =========
		 */
		o_U.op_setVectorPlusScalarMulVector(phi0_U, _timestepSize, phi1_FU);

		return true;
	}


private:
	bool _eval_timeIntegration_ETDRK2(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		/*
		 * A_{n} = \psi_{0}(\Delta tL)U_{n} + \Delta t \psi_{1}(\Delta tL) F(U_{n})
		 *         NNNNNNNNNNNNNNNNNNNNNNNN
		 */
		sweet::Data::GenericContainer::Base &phi0_U = *_tmpDataContainer[0];
		_evalPhi0(i_U, phi0_U, i_simulationTime);

		/*
		 * A_{n} = \psi_{0}(\Delta tL)U_{n} + \Delta t \psi_{1}(\Delta tL) F(U_{n})
		 *         ========================                                NNNNNNNN
		 */
		sweet::Data::GenericContainer::Base &FU = *_tmpDataContainer[1];
		_evalNL(i_U, FU, i_simulationTime);

		/*
		 * A_{n} = \psi_{0}(\Delta tL)U_{n} + \Delta t \psi_{1}(\Delta tL) F(U_{n})
		 *         ========================            NNNNNNNNNNNNNNNNNNN ========
		 */
		sweet::Data::GenericContainer::Base &phi1_FU = *_tmpDataContainer[2];
		_evalPhi1(FU, phi1_FU, i_simulationTime);

		/*
		 * A_{n} = \psi_{0}(\Delta tL)U_{n} + \Delta t \psi_{1}(\Delta tL) F(U_{n})
		 * NNNNN   ========================            =================== ========
		 */
		sweet::Data::GenericContainer::Base &A = *_tmpDataContainer[3];
		A.op_setVectorPlusScalarMulVector(phi0_U, _dt, phi1_FU);

		/*
		 * U_{n+1} = A_{n} + \Delta t \psi_{2}(\Delta tL) [ (F(A_{n},t_{n}+\Delta t)-F(U_{n}) ]
		 *           =====                                   NNNNNNNNNNNNNNNNNNNNNNN
		 */
		sweet::Data::GenericContainer::Base &FA = *_tmpDataContainer[4];
		_evalNL(A, FA, i_simulationTime + _dt);

		/*
		 * U_{n+1} = A_{n} + \Delta t \psi_{2}(\Delta tL) [ (F(A_{n},t_{n}+\Delta t)-F(U_{n}) ]
		 *           =====                                ========================== NNNNNNNN
		 */
		FA.op_subVector(FU);

		/*
		 * U_{n+1} = A_{n} + \Delta t \psi_{2}(\Delta tL) [ (F(A_{n},t_{n}+\Delta t)-F(U_{n}) ]
		 *           =====            NNNNNNNNNNNNNNNNNNN =====================================
		 */
		sweet::Data::GenericContainer::Base &phi2_X = *_tmpDataContainer[5];
		_evalPhi2(FA, phi2_X, i_simulationTime);

		/*
		 * U_{n+1} = A_{n} + \Delta t \psi_{2}(\Delta tL) [ (F(A_{n},t_{n}+\Delta t)-F(U_{n}) ]
		 * =======   =====            =================== =====================================
		 */
		o_U.op_setVectorPlusScalarMulVector(A, _dt, phi2_X);

		return true;
	}


private:
	bool _eval_timeIntegration_ETDRK4(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	)
	{
		double dt_half = _dt*0.5;

		/*
		 * Precompute commonly used terms
		 */

		/*
		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + \Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
		 *         NNNNNNNNNNNNNNNNNNNNNNNNNNNN
		 */
		sweet::Data::GenericContainer::Base &phi0_U = *_tmpDataContainer[0];
		_evalPhi0(i_U, phi0_U, i_simulationTime);

		/*
		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + \Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
		 *         ============================                                   NNNNNNNN
		 */
		sweet::Data::GenericContainer::Base &FU = *_tmpDataContainer[1];
		_evalNL(i_U, FU, i_simulationTime);

		/*
		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + \Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
		 *         ============================           NNNNNNNNNNNNNNNNNNNNNNN ========
		 */
		sweet::Data::GenericContainer::Base &phi1 = *_tmpDataContainer[2];
		_evalPhi1(FU, phi1, i_simulationTime);


		/*
		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
		 * NNNNN   ============================               ======================= ========
		 */
		sweet::Data::GenericContainer::Base &A = *_tmpDataContainer[3];
		A.op_setVectorPlusScalarMulVector(phi0_U, dt_half, phi1);

		/*
		 * B_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(A_{n}, t_{n} + 0.5*\Delta t)
		 *         ============================                                       NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
		 */
		sweet::Data::GenericContainer::Base &FA = *_tmpDataContainer[4];
		_evalNL(A, FA, i_simulationTime+dt_half);

		/*
		 * B_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(A_{n}, t_{n} + 0.5*\Delta t)
		 *         ============================               ======================= ==============================
		 */
		_evalPhi1(FA, phi1, i_simulationTime);

		/*
		 * B_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(A_{n}, t_{n} + 0.5*\Delta t)
		 * =====   ============================               ======================= ==============================
		 */
		sweet::Data::GenericContainer::Base &B = *_tmpDataContainer[5];
		B.op_setVectorPlusScalarMulVector(phi0_U, dt_half, phi1);

		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)A_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 *         NNNNNNNNNNNNNNNNNNNNNNNNNNNN                                                                          ==============
		 */
		sweet::Data::GenericContainer::Base &phi0_A = *_tmpDataContainer[6];
		_evalPhi0(A, phi0_A, i_simulationTime);

		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)A_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 *         ============================                                            NNNNNNNNNNNNNNNNNNNNNNNNNNNNN ==============
		 */
		sweet::Data::GenericContainer::Base &FB = *_tmpDataContainer[7];
		_evalNL(B, FB, i_simulationTime+dt_half);

		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)A_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 *         ============================                                        NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
		 */
		sweet::Data::GenericContainer::Base &tmp1 = *_tmpDataContainer[8];
#if 1
		tmp1.op_setVector(FB);
		tmp1.op_mulScalar(2.0);
		tmp1.op_subVector(FU);
#else
		tmp1.op_setVectorPlusScalarMulVector(FU, -2.0, FB);
		tmp1.op_mulScalar(-1);
#endif


		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)A_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 *         ============================               NNNNNNNNNNNNNNNNNNNNNNNN =================================================
		 */
		_evalPhi1(tmp1, phi1, i_simulationTime);

		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)A_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 * NNNNN   ============================               ======================== =================================================
		 */
		sweet::Data::GenericContainer::Base &C = *_tmpDataContainer[9];
		C.op_setVectorPlusScalarMulVector(phi0_A, dt_half, phi1);


		/*
		 * R0 - R3
		 */
		sweet::Data::GenericContainer::Base &FC = *_tmpDataContainer[10];
		_evalNL(C, FC, i_simulationTime+_dt);

		const sweet::Data::GenericContainer::Base &R0 = i_U;
		sweet::Data::GenericContainer::Base &R1 = FU;

		sweet::Data::GenericContainer::Base &R2 = *_tmpDataContainer[11];
		R2.op_setVectorPlusScalarMulVector(FA, 1.0, FB);

		sweet::Data::GenericContainer::Base &R3 = FC;

		/*
		 * U_{n+1} =
		 * 		\psi_{0}(\Delta tL)R_{0}
		 * 			+ \Delta t
		 * 			(
		 * 				  \upsilon_{1}(\Delta tL) R_{1} +
		 * 				2*\upsilon_{2}(\Delta tL) R_{2} +
		 * 				  \upsilon_{3}(\Delta tL) R_{3}
		 * 			)
		 */

		sweet::Data::GenericContainer::Base &R0_ = *_tmpDataContainer[12];
		_evalUps0(R0, R0_, i_simulationTime);

		sweet::Data::GenericContainer::Base &R1_ = *_tmpDataContainer[13];
		_evalUps1(R1, R1_, i_simulationTime);

		sweet::Data::GenericContainer::Base &R2_ = *_tmpDataContainer[14];
		_evalUps2(R2, R2_, i_simulationTime);

		sweet::Data::GenericContainer::Base &R3_ = *_tmpDataContainer[15];
		_evalUps3(R3, R3_, i_simulationTime);


		o_U.op_setVectorPlusScalarMulVector(R0_, _dt, R1_);
		o_U.op_addScalarMulVector(2.0*_dt, R2_);
		o_U.op_addScalarMulVector(_dt, R3_);

		return true;
	}

	void print(const std::string &i_prefix = "")
	{
		std::string newPrefix = i_prefix + "  ";
		std::cout << i_prefix << "ETDRK(" << std::endl;
		std::cout << newPrefix << "  order: " << _order << std::endl;
		std::cout << i_prefix << ")" << std::endl;
	}
};

}}}

#endif
