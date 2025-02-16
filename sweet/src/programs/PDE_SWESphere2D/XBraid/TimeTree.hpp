#ifndef PROGRAMS_PDE_SWESPHERE2D_XBRAID_TIMETREE_HPP
#define PROGRAMS_PDE_SWESPHERE2D_XBRAID_TIMETREE_HPP

#include <sweet/Error/Base.hpp>
#include <sweet/XBraid/TimeTree.hpp>
#include <sweet/XBraid/Vector.hpp>
#include <programs/PDE_SWESphere2D/TimeTree/TimeTree.hpp>
#include <programs/PDE_SWESphere2D/DataContainer/Simulation.hpp>
#include <programs/PDE_SWESphere2D/XBraid/DataContainer.hpp>


namespace PDE_SWESphere2D {
namespace XBraid {

class TimeTree :
	public sweet::XBraid::TimeTree
{

public:

	sweet::Error::Base error;

	PDE_SWESphere2D::TimeTree::TimeTree* tsm = nullptr;

public:

	TimeTree()
	{
		tsm = new PDE_SWESphere2D::TimeTree::TimeTree;
	}

public:
	~TimeTree()
	{
		clear();
	}

public:
	void clear()
	{
		if (tsm)
		{
			tsm->clear();
			tsm = nullptr;
		}
	}

public:
	void setup(
			std::string &i_timestepping_method,
			sweet::Shacks::Dictionary* i_shackDict,
			sweet::Data::Sphere2D::Operators* i_op,
			sweet::Data::Sphere2DComplex::Operators* i_op_complex,
			PDE_SWESphere2D::XBraid::DataContainer* i_u,
			double i_dt
	)
	{
		tsm->setup_1_registerAllTimesteppers();
		tsm->setup_2_shackRegistration(i_shackDict);
		ERROR_FORWARD(*tsm);

		sweet::Shacks::ProgramArgumentsDictionary* shackDict = (sweet::Shacks::ProgramArgumentsDictionary*) i_shackDict;
		bool retval = tsm->setup_3_timestepper(
						i_timestepping_method,
						shackDict,
						i_op,
						i_op_complex,
						*i_u->data
					);

		if (!retval)
		{
			SWEETErrorFatal("Invalid time stepping method: " + i_timestepping_method);
			////int helpVerbosity = 0;
			////if (shackTimeDisc->timestepping_method == "helpall")
			////	helpVerbosity = 1;

			////timeSteppersNewTS.outputHelp(std::cout, "", helpVerbosity);

			////error.clear();
			////std::cout << "Finishing now..." << std::endl;
		}

		tsm->timeIntegrator->setTimeStepSize(i_dt);
		ERROR_FORWARD(*tsm);
	}

public:

	virtual
	bool runIntegration(
				const sweet::XBraid::Vector* i_U,
				const sweet::XBraid::Vector* o_U,
				double i_simulationTime
	)
	{

		PDE_SWESphere2D::XBraid::DataContainer* U1 = (PDE_SWESphere2D::XBraid::DataContainer*) i_U;
		PDE_SWESphere2D::XBraid::DataContainer* U2 = (PDE_SWESphere2D::XBraid::DataContainer*) o_U;

		tsm->runIntegration(
				*U1->data,
				*U2->data,
				i_simulationTime
			);

		return true;
	}

public:
	bool storePreviousSolution(
				const sweet::XBraid::Vector* i_U
	) override
	{
		PDE_SWESphere2D::XBraid::DataContainer* U = (PDE_SWESphere2D::XBraid::DataContainer*) i_U;

		tsm->timeIntegrator->storePrevSolution(U->data);

		return true;
	}

};

}}

#endif
