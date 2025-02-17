# SWEET-OpenFOAM Coupling

## Preparation
- Install OpenFOAM
- Install SWEET in ```/sweet/``` following the README
- Install the OpenFOAM adapter in ```/openfoam-adapter/```

## Running
- Before running anything go to ```/sweet/``` and execute ```source activate.sh```
- Go to ```/coupling/```

### OpenFOAM
- ```./interfoam_case/BlockMesh``` Creates empty OpenFOAM domain
- ```./interfoam_case/SnappyMesh``` Creates OpenFOAM domain with obstacle, to switch out obstacle, place corresponding ```.stl``` file in ```interfoam_case/constant/triSurface``` and rename to ```geometry.stl```
- ```./interfoam_case/Allrun``` Runs the OpenFOAM solver

### SWEET Cart2D
- ```./build.sh``` Compiles the Cart2D program
- ```./combined.sh -p both``` Runs the program in coupled mode, and runs postprocessing
- ```./combined.sh``` Runs the program in benchmark mode, and runs postprocessing
- ```./plot.sh``` Postprocesses the output and renders images and videos


### SWEET Sphere2D
Go to ```/coupling/sphere/```
- ```run_erk.sh``` Build and Run Sphere2D program
- ```run_sdc.sh``` Run SDC configuration of Sphere2D program
- ```./plot.sh``` Run postprocessing and render images and videos

## Main Contribution Index
Location of the main implementations in the code. There are smaller changes scattered in many other files but these are central features.

### preCICE
- ```/coupling/precice-config.xml```

### OpenFOAM
- ```/coupling/interfoam_case/```

### SWEET
- Cart2D Program: ```sweet/src/programs/PDE_SWECart2D.cpp```
- Sphere2D Program: ```sweet/src/programs/PDE_SWESphere2D.cpp```
- Test Scenario Cart2D: ```sweet/src/programs/PDE_SWECart2D/Benchmarks/column.hpp```
- Test Scenario Sphere2D: ```sweet/src/programs/PDE_SWESphere2D/Benchmarks/column.hpp```

### OpenFOAM Adapter:
- ```openfoam-adapter/FF/Alpha.C```
- ```openfoam-adapter/FF/Velocity.C```
- ```openfoam-adapter/FF/DimensionReduction.C```
