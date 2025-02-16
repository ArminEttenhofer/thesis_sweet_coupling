#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_L_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_L_HPP


#include <sweet/Data/GenericContainer/CastHelper.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/ExpIntegration/ExpFunction.hpp>
#include <sweet/TimeTree/TimeTree_Node_LeafHelper.hpp>

#include "../../DataContainer/Simulation.hpp"
#include "../../DataContainer/Config.hpp"
#include "../../Shack.hpp"
#include "../../Benchmarks/Shack.hpp"

/*
 * Time tree node related includes
 */
#include "../../TimeHelpers/SphBandedMatrix_GridComplex.hpp"
#include "../../TimeHelpers/SphBandedMatrix_GridReal.hpp"


namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {


class l	:
	public sweet::TimeTree::TimeTree_Node_LeafHelper<l>,
	public sweet::Data::GenericContainer::CastHelper<
			DataContainer::Simulation,
			DataContainer::Config
		>
{
private:
	Shack *_shackPDESWESphere2D;
	sweet::Data::Sphere2D::Shack *_shackSphere2DDataOps;
	PDE_SWESphere2D::Benchmarks::Shack *_shackPDESWESphere2DBenchmarks;
	const sweet::Data::Sphere2D::Operators *_ops;
	const sweet::Data::Sphere2DComplex::Operators *_opsComplex;

	/*!
	 * Storage locations for EVAL_TENDENCIES
	 */
	class _EvalTendencies{
	public:
		/*!
		 * Coriolis effect
		 */
		sweet::Data::Sphere2D::DataGrid _fg;

		/*!
		 * Velocities
		 */
		sweet::Data::Sphere2D::DataGrid _ug;
		sweet::Data::Sphere2D::DataGrid _vg;

		/*!
		 * Other temporary variables
		 */
		sweet::Data::Sphere2D::DataGrid _tmpg1;
		sweet::Data::Sphere2D::DataGrid _tmpg2;
	} _evalDataTendencies;


	/*!
	 * Storage locations for EVAL_BACKWARD_EULER
	 */
	class _EvalBackwardEuler{
	public:
		//! Temporary array
		sweet::Data::Sphere2D::DataSpectral _rhs;

		//! Solver
		SphBandedMatrix_GridReal sphBandedMatrix_GridReal;

	} _evalDataBackwardEuler;



	/*!
	 * Storage locations for EVAL_REXI_TERM
	 */
	class _EvalREXITerm{
	public:
		//! Temporary array
		sweet::Data::Sphere2DComplex::DataSpectral _rhs;

		/*!
		 * Matrix solver
		 */
		SphBandedMatrix_GridComplex sphBandedMatrix_GridComplex;

		/*!
		 * Complex-valued time step size for complex-valued backward Euler
		 *
		 * This is used for REXI solvers of the form
		 * 	U_1 = (I-dt*L)^{-1} U_0
		 */
		std::complex<double> _rexiTermAlpha;
		std::complex<double> _rexiTermBeta;
		std::complex<double> _rexiTermGamma;
		bool _rexiTermGammaActive;

		_EvalREXITerm()	:
			_rexiTermAlpha(666,-666),
			_rexiTermBeta(666,-666),
			_rexiTermGamma(666,-666),
			_rexiTermGammaActive(false)
		{
		}

		_EvalREXITerm(const _EvalREXITerm &i_src)
		{
			_rexiTermAlpha = i_src._rexiTermAlpha;
			_rexiTermBeta = i_src._rexiTermBeta;
			_rexiTermGamma = i_src._rexiTermGamma;
			_rexiTermGammaActive = i_src._rexiTermGammaActive;
		}
	} _evalDataREXITerm;


public:
	l();

public:
	~l();

public:
	l(const l &i_val);

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override;

public:
	const std::vector<std::string> getNodeNames() override;

public:
	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override;

public:
	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper
	) override;

	/*!
	 * Set string-based value with a key
	 */
	bool setupByKeyValue(
			const std::string &i_key,
			const std::string &i_value
	) override;


	/*!
	 * Set complex value with a key
	 */
	bool setupByKeyValue(
			const std::string &i_key,
			const std::complex<double> &i_value
	) override;

	/*!
	 * Clear all data structures
	 */
	void clear() override;

	/*!
	 * Update the time step size.
	 *
	 * This should be called only if the time step size changed!
	 */
	bool setTimeStepSize(double i_dt)	override;

	/*!
	 * Return the time tendencies of the PDE term
	 */
private:
	bool _eval_tendencies(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	)	override;

	/*!
	 * Return the time tendencies of the PDE term
	 */
private:
	bool _eval_eulerBackward(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	)	override;

	/*!
	 * Return the evaluation of a single REXI term
	 */
private:
	bool _eval_rexiTerm(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	)	override;
};

}}}

#endif
