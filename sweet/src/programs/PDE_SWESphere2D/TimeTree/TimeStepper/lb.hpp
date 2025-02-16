#ifndef PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_LB_HPP
#define PROGRAMS_PDE_SWESPHERE2D_TIMETREE_TIMESTEPPER_LB_HPP


#include <sweet/Data/GenericContainer/CastHelper.hpp>
#include <sweet/Error/Base.hpp>
#include <sweet/TimeTree/TimeTree_Node_Base.hpp>
#include <sweet/TimeTree/TimeTree_Node_LeafHelper.hpp>
#include "../../DataContainer/Simulation.hpp"
#include "../../DataContainer/Config.hpp"
#include "../../Shack.hpp"
#include "../../Benchmarks/Shack.hpp"

/*
 * Time tree node related includes
 */

namespace PDE_SWESphere2D {
namespace TimeTree {
namespace TimeStepper {

class lb	:
	public sweet::TimeTree::TimeTree_Node_LeafHelper<lb>,
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

public:
	lb();
	~lb();

	lb(const lb &i_val);

public:
	bool shackRegistration(
			sweet::Shacks::Dictionary *io_shackDict
	) override;

public:
	bool outputHelp(
			std::ostream &o_ostream,
			const std::string &i_prefix = ""
	) override;

	const std::vector<std::string> getNodeNames() override;

	bool setupConfigAndForwardTimeStepperEval(
		const sweet::Data::GenericContainer::ConfigBase &i_deTermConfig,
		TIME_STEPPER_TYPES i_evalType,
		TimeTree_Node_Base::EvalFun *o_timeStepper
	) override;

	bool setupByKeyValue(
			const std::string &i_key,
			const std::string &i_value
	) override;

	bool setupByKeyValue(
			const std::string &i_key,
			const std::complex<double> &i_value
	) override;

	void clear() override;


public:
	bool _eval_tendencies(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	) override;

	/*
	 * Return the backward Euler time step
	 */
public:
	bool _eval_eulerBackward(
			const sweet::Data::GenericContainer::Base &i_u,
			sweet::Data::GenericContainer::Base &o_u,
			double i_timeStamp
	) override;

};

}}}

#endif
