/*
 * Author: Thibaut LUNET <thibaut.lunet@tuhh.de>
 */
#ifndef INCLUDE_SWEET_PARALLELIZATION_SHACK_HPP
#define INCLUDE_SWEET_PARALLELIZATION_SHACK_HPP

#include <sweet/Shacks/Base.hpp>
#include <string>
#include <iostream>
#include <sweet/Tools/ProgramArguments.hpp>

#if SWEET_MPI
#   include <mpi.h>
#endif

namespace sweet {
namespace Parallelization {


/*!
 * Shack with everything around parallelization (OpenMP + MPI)
 */
class Shack : 
	public sweet::Shacks::Base
{
public:
	/**
	 * Wether MPI is used or not
	 */
	bool useMPI = false;

	/**
	 * MPI size of the main communicator
	 */
	int mpiSize = 1;

	/**
	 * MPI rank of the current process
	 */
	int mpiRank = 0;

	/**
	 * Whether or not this process is root in MPI_COMM_WORLD
	 *
	 * TODO: Change this so that it's not must for MPI_COMM_WORLD
	 */
	bool isMPIRoot = true;

public:
	
	Shack() {
#if SWEET_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
		MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
		useMPI = true;
		isMPIRoot = (mpiRank == 0);
#endif
	}

	void printProgramArguments(const std::string &i_prefix = "")  override
	{}

	bool processProgramArguments(sweet::Tools::ProgramArguments &i_pa)  override
	{
		return true;
	}

	void printShack(
		const std::string& i_prefix = ""
	) override
	{
		std::cout << i_prefix << std::endl;
		std::cout << i_prefix << "PARALLELIZATION:" << std::endl;
		std::cout << i_prefix << " + useMPI: " << useMPI << std::endl;
		std::cout << i_prefix << " + mpiSize: " << mpiSize << std::endl;
		std::cout << i_prefix << " + mpiRank: " << mpiRank << std::endl;

		std::cout << i_prefix << std::endl;
	}

};

}}

#endif
