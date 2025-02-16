#ifndef PROGRAMS_LIBPFASST_INTERFACE_LEVELSINGLETON_HPP
#define PROGRAMS_LIBPFASST_INTERFACE_LEVELSINGLETON_HPP

#include <sweet/Data/Sphere2D/Sphere2D.hpp>
#include <sweet/Data/Sphere2D/Operators.hpp>
#include "../../PDE_SWESphere2D/BenchmarksCombined.hpp"

/*
 * Class to store the configurations and operators at each level
 */
class LevelSingleton
{
public:
	int level;
	sweet::Data::Sphere2D::Config sphere2DDataConfig;
	sweet::Data::Sphere2D::Operators ops;

	PDE_SWESphere2D::Benchmarks::BenchmarksCombined benchmarks;
};

#endif
