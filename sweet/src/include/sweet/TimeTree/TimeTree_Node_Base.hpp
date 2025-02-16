/*
 * Make sure that we always try to include these classes to due
 * forward declarations
 */

#ifndef INCLUDE_SWEET_TIMETREE_TIMETREE_NODE_BASE_HPP
#define INCLUDE_SWEET_TIMETREE_TIMETREE_NODE_BASE_HPP

#include <vector>
#include <string>
#include <complex>
#include <sstream>
#include <sweet/Error/Base.hpp>
#include <sweet/Shacks/Dictionary.hpp>
#include <sweet/Data/GenericContainer/Base.hpp>
#include <sweet/Data/GenericContainer/ConfigBase.hpp>
#include <sweet/TimeTree/TimeTreeIR.hpp>


namespace sweet {
namespace TimeTree {

class TimeTreeIR_2_TimeTreeNodes;


/*!
 * \brief Base for node of time tree
 *
 * Node in time tree which can either represent a time integrator or a terminal DE node
 */
class TimeTree_Node_Base
{
public:
	sweet::Error::Base error;

#if SWEET_XBRAID
	sweet::Data::GenericContainer::Base* U_prev_solution = nullptr;
#endif

	/*!
	 * Different evaluations for time integration
	 */
protected:
	enum TIME_STEPPER_TYPES
	{
		EVAL_NONE = 0,
		EVAL_INTEGRATION = 1,
		EVAL_TENDENCIES,
		EVAL_EULER_BACKWARD,
		EVAL_REXI_TERM,
		EVAL_EXPONENTIAL,
		EVAL_SEMI_LAGRANGIAN,
	};


	/*!
	 * This will be set to the requested eval type.
	 * Only one kind of evaluation is allowed per tree node!
	 *
	 * This should help to setup an optimized evaluation.
	 */
protected:
	TIME_STEPPER_TYPES _evalTypeRequested;

	/*!
	 * Return human-readable string of eval type
	 */
public:
	const std::string evalTypeToString(TIME_STEPPER_TYPES i_evalType);


	/*!
	 * Type define "EvalFun" which is required to time step
	 */
public:
	typedef bool (TimeTree_Node_Base::*EvalFun)(
					const sweet::Data::GenericContainer::Base &i_U,
					sweet::Data::GenericContainer::Base &o_U,
					double i_simulationTime
			);

	/*!
	 * Keep track of a list of implemented evaluation functions
	 */
protected:
	std::vector<TIME_STEPPER_TYPES> _registeredEvalTypes;

public:
	TimeTree_Node_Base();

	virtual
	~TimeTree_Node_Base();

	/*!
	 * Copy constructor
	 */
public:
	TimeTree_Node_Base(
			const TimeTree_Node_Base &i_src
	);

	/*!
	 * Return vector of strings which are a unique IDs for this node.
	 *
	 * e.g.
	 * 	'ERK': fast gravity modes
	 * 	'SS': coriolis effect
	 * 	'l': lg+lc
	 */
	virtual
	const std::vector<std::string>
	getNodeNames() = 0;


	/*!
	 * Return a copy of this class
	 *
	 * This is a special case where we might like to reuse the already utilized shacks.
	 */
	virtual
	std::shared_ptr<TimeTree_Node_Base> getInstanceCopy() = 0;

	/*!
	 * Output a nicely formatted list of node names
	 */
protected:
public:
	std::string _getNodeNamesAsString()
	{
		std::ostringstream oss;
		const std::vector<std::string> nn = getNodeNames();
		for (std::size_t i = 0; i < nn.size(); i++)
		{
			oss << "'" << nn[i] << "'";
			if (i < nn.size()-1)
				oss << ", ";
		}
		return oss.str();
	}


	/*!
	 * Print out help for this tree node (including options, etc.)
	 */
public:
	virtual
	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	)
	{
		o_ostream << i_prefix << "TODO: Implement this for " << _getNodeNamesAsString() << std::endl;

		return true;
	}


	/*!
	 * Shack registration
	 */
	virtual bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) = 0;


	/*!
	 * Setup the treenode by a given function
	 */
	virtual
	bool setupTreeNodeByFunction(
			std::shared_ptr<sweet::TimeTree::TimeTreeIR::Function> &i_function,
			sweet::TimeTree::TimeTreeIR_2_TimeTreeNodes &i_tsAssemblation
	);



	/*!
	 * Final setup with
	 * - problem-specific configuration
	 * - evaluation type (which time interface to request)
	 * - callback handler to be initialized
	 *
	 * This should setup all internal storage spaces
	 */
	virtual
	bool setupConfigAndForwardTimeStepperEval(
			const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
			TIME_STEPPER_TYPES i_evalType,
			EvalFun *o_timeStepper = nullptr	//<! Setup eval function if not null.
		) = 0;



	/*!
	 * Default setup initializing EVAL_INTEGRATION
	 */
	bool setupConfigAndForwardTimeStepperEval(
			const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
			EvalFun *o_timeStepper	//<! Setup eval function if not null.
		);


	/*!
	 * This allows customizing the time stepper by its parent class
	 *
	 * This can be helpful, e.g., for exponential integration to specify the particular phi function:
	 *
	 * "expIntegratorFunction" => "phi0"
	 * or
	 * "expIntegratorFunction" => "phi1"
	 * ...
	 */
	virtual
	bool setupByKeyValue(
			const std::string &i_key,
			const std::string &i_value
	);


	/*!
	 * Setup with a floating point value
	 */
	virtual
	bool setupByKeyValue(
			const std::string &i_key,
			const double &i_value
	);


	/*!
	 * Setup with a complex value
	 */
	virtual
	bool setupByKeyValue(
			const std::string &i_key,
			const std::complex<double> &i_value
	);


	/*!
	 * Return the callback to the time stepper
	 */
protected:
	bool registerTimeStepperEval(
			TIME_STEPPER_TYPES i_evalType	//!< Type of time integrator
	);

	/*!
	 * Return the callback to the time stepper
	 */
protected:
	bool registerTimeStepperEval(
			TIME_STEPPER_TYPES i_evalType,	//!< Type of time integrator
			EvalFun *o_timeStepper	//!< Callback of time integrator
	);

	/*!
	 * Cleanup internal data structures
	 */
public:
	virtual
	void clear();


	/*!
	 * Set the time step size \f$\Delta t\f$.
	 *
	 * This is required, e.g., to setup certain data structures for an implicit time steppers
	 */
	virtual
	bool setTimeStepSize(double i_dt) = 0;

	/*!
	 * Return true if an evaluation function exists.
	 */
public:
	bool isEvalAvailable(TIME_STEPPER_TYPES i_evalType);


	/*!
	 * Return the supported evals as a string with evals comma separated
	 */
public:
	std::string getSupportedEvalsAsString();


	/*!
	 * Set an evaluation function to be available
	 */
public:
	void setEvalAvailable(TIME_STEPPER_TYPES i_evalType);


	/*!
	 * Return the time integration
	 */
private:
	virtual
	bool _eval_integration(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_simulationTime
	);


	/*!
	 * Optional: Return the time tendencies of the DE term
	 */
public:
	virtual
	bool _eval_tendencies(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_timeStamp
	);



	/*!
	 * Optional: Return the backward Euler evaluation of the term:
	 *
	 * \f[
	 *
	 * \frac{( U^{n+1} - U^{n} )}{\Delta t} = \frac{d}{dt} U^{n+1}\\
	 *
	 * \Leftrightarrow U^{n+1} - U^{n} = \Delta t \frac{d}{dt} U^{n+1}\\
	 *
	 * \Leftrightarrow (I - \Delta t \frac{d}{dt}) U^{n+1} = U^{n}
	 *
	 * \f]
	 *
	 * For a linear operator \f$ \frac{d}{dt}U = LU \f$ we would obtain
	 *
	 * \f[
	 * \Leftrightarrow (I - \Delta t L) U^{n+1} = U^{n}
	 * \f]
	 */
	virtual
	bool _eval_eulerBackward(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_timeStamp
	);

	/*!
	 * Optional: Return the REXI term evaluation
	 *
	 * These are terms of the form
	 *
	 * \f[
	 * 	U^{n+1} = \beta (\Delta t L - \alpha)^{-1} U^{n}
	 * \f]
	 *
	 * The \f$\alpha\f$ and \f$\beta\f$ values are complex-valued.
	 *
	 * They are set by calling setByKeyValue(...) with a complex-valued argument.
	 */
	virtual
	bool _eval_rexiTerm(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_timeStamp
	);

	/*!
	 * Optional: Return an evaluation of the exponential term
	 *
	 * \f[
	 * 	U^{n+1} = exp(\Delta t L) U^n
	 * \f]
	 */
	virtual
	bool _eval_exponential(
			const sweet::Data::GenericContainer::Base &i_U,
			sweet::Data::GenericContainer::Base &o_U,
			double i_timeStamp
	);

	/*!
	 * EVAL_SEMI_LAGRANGIAN: Return the number of states required for this method
	 */
public:
	virtual
	bool evalNA_getNumStates(
			int *o_numStates
	);


	/*!
	 * EVAL_SEMI_LAGRANGIAN: Return the departure points for the given states
	 *
	 * This function computes the departure points. There are various algorithms
	 * to accomplish this and higher-order can require multiple states.
	 *
	 * Hence, we provide the states as a vector.
	 *
	 * The 1st element refers to U(t)
	 * The 2nd element refers to U(t-dt)
	 * The 3rd element refers to U(t-2*dt)
	 * etc.
	 */
public:
	virtual
	bool evalNA_departurePoints(
			const sweet::Data::GenericContainer::Base* i_states[],	//!< Vector of states
			double i_timestepSize,		//!< Time step size (how far to go backwards in time along the trajectory)
			sweet::Data::GenericContainer::Base &o_departurePositions		//!< Computed departure positions
	);

	/*!
	 * EVAL_SEMI_LAGRANGIAN: Return the departure points for the given states
	 *
	 * Interpolate the input data at the given sampling positions
	 */
public:
	virtual
	bool evalNA_interpolate(
			const sweet::Data::GenericContainer::Base &i_U_input,		//!< Input simulation data
			const sweet::Data::GenericContainer::Base &i_samplingPositions,	//!< Sampling positions (computed by _evalNA_departurePoints)
			sweet::Data::GenericContainer::Base &o_U_samples			//!< Output samples
	);

	/*!
	 * Message to help debugging errors in the time tree (e.g. about the parsing of the time tree string)
	 */
private:
	std::string _debugMessage;

	/*!
	 * Set the debug message
	 */
public:
	std::string setDebugMessage(
			const std::string &i_debugMessage	//!< Message for debugging
	);

	/*!
	 * Return the debug message including a newline.
	 * This is handy for figuring out wrong parameters in the time tree.
	 *
	 * \return String with hopefully helpful debug information
	 */
	std::string getNewLineDebugMessage();

#if SWEET_XBRAID
public:
	virtual 
	bool storePrevSolution(
			sweet::Data::GenericContainer::Base* i_U		//!< Input simulation data
			////const sweet::Data::GenericContainer::ConfigBase* i_deTermConfig
	)
	{
		return true;
	}
#endif

};

}}

#endif
