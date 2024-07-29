# FlipSolver2d
![Linux build](https://github.com/ArtNlk/FlipSolver2d/actions/workflows/ubuntu-build.yml/badge.svg)
[![CMake windows-latest build](https://github.com/ArtNlk/FlipSolver2d/actions/workflows/windows-latest%20build.yml/badge.svg)](https://github.com/ArtNlk/FlipSolver2d/actions/workflows/windows-latest%20build.yml)

A 2d testing ground for writing FLIP/PIC fluid solver

## Building
1. After cloning, get submodules to setup imgui:
```
git submodule init
git submodule update --recursive
```
2. Create and go to build directory
3. Run:
```
cmake PATH_TO_CMakeLists_IN_REPO_ROOT -DCMAKE_BUILD_TYPE=Release
```
4. Wait for cmake to fetch dependencies and finish
5. Run:
```
cmake --build . --config Release -j NUM_OF_CORES
```
6. Outputs will be in following dirs:
  - Liquid2dRender:
    - scenes directory contains all scenes to load
    - config.json specifies what scene to load
  - bench:
    - Test benchmark, chages time to time, used as quick test for various things
## Directory structure
- FlipSolver2dLib:  
  Main solver library sources
- FlipSolver2dLib/threading:  
  Sources for thread pool
- Liquid2dRender:  
  Main app sources
- Liquid2dRender/include:  
  glad and khr headers
- Liquid2dRender/src:  
  glad source
- Liquid2dRender/scenes:  
  All scenes that are watched by cmake and that will be copied to output scenes directory
- bench:  
  Testing benchmark sources, change randomly, used to quickly test stuff outside of full gui app
- devDocs:  
  Various files to not forget some aspects of the simulation.  
  Are not meant to replace documentation and may be confusing

## Scene file specification
Each scene file is a json with two top-level keys:
- settings:  
  Holds global scene settings. Scene 0,0 is in top left, (X,I,U) points down, (Y,J,V) points right
  - simType:  
    flip, nbflip, smoke, fire  
    Specifies simulation type
  - domainSizeI, domainSizeJ:  
  Domain size in units (read meters)
  - resolution:  
  Number of grid cells along the shortest dimension, other dimension is calculated automatically
  - fps:  
  Target number of steps per second
  - density:  
  Fluid desity
  - maxSubsteps:  
  Maximum number of substeps per one big step
  - pcgIterLimit:  
  Max number of PCG solver iterations
  - seed:  
  Seed for RNG
  - particlesPerCell:  
  Target number of particles to keep in each cell
  - particleScale:  
  Radius of particles for SDF update, relative to grid cell size (0.5 means particle radius is half of grid cell size length). Keep 0.8 for 2d cases
  - picRatio:  
  PIC/FLIP velocity update methods ratio. Lower ratio - More FLIP, more splash, but less stability. 0.03 is usually a good spot
  - cflNumber:  
  Solver will try to keep particles from traveling more than this many cells per time step, substepping if necessary. 5 is a good start
  - scale:  
  Scale to be applied for whole scene. So domain 1x1 meter in size with scale 0.01 becoms 1x1 centimeter
  - viscosityEnabled:  
  Toggle viscosity solver for water sims
  - heavyViscosity:  
  If true, use more expensive but potentially more accurate viscosity model
  - globalAcceleration:  
  2 number array, holds gravitational acceleration. Mind the coordinate system.
  - parameterHandlingMethod:  
    grid, particle (default - particle)  
    Advect secondary parameters (viscosity, soot, fuel) via grid or particles. May cause interesting results, experimental
  - ambientTemperature:  
    Gas solver parameter. Ambient air temperature
  - temperatureDecayRate:  
    Gas solver parameter. How quickly does the gas cool down?
  - concentrationDecayRate:  
    Gas solver parameter. How quickly does smoke dissolve?
  - buoyancyFactor:  
    Gas solver parameter. Temperature weight in buoyancy calculations. Keep at 1 for realistic math
  - sootFactor:  
    Gas solver parameter. Soot weight in buoyancy calculations. Higher number - heavier the smoke is.
  - ignitionTemp:  
    Fire solver parameter. Fuel needs to be hotter than this to burn.
  - burnRate:  
    Fire solver parameter. How quickly the fuel burns?
  - smokeEmission:  
    Fire solver parameter. How much smoke does fuel emit when burnt?
  - heatEmission:  
    Fire solver parameter. How much heat does fuel produce when burnt?
  - billowing:  
    Fire solver parameter. How much does the fuel expand when burnt? Experimental, may not work in many versions.
- solver:
  - objects:  
      Holds all objects in the scene. For now, the only key under ```solver``` section.  
      Each object has:
    - type:  
      solid, source, sink, fluid  
      What type of object this object is
    - enabled:  
      Allows toggling objects for quick tinkering
    - verts:  
      Array of 2-number arrays. Specifies verticies of the object. Winding order does not matter, last vertex is connected to first
    
    Solids have:
    - friction:  
      How much does the obstacle slow down the fluid?
    
    Fluid sources and initial fluids have:
    - viscosity:  
      Water solver parameter. Initial fluid viscosity
    - temperature:  
      Gas solver parameter. Initial gas temperature
    - concentration:  
      Gas solver parameter. Initial smoke concentration
    - fuel:  
      Fire solver parameter. Initial amount of fuel in the emitted gas
    - transferVelocity:  
      If true, initial fluid velocity is set to ```velocity``` value
    - velocity:  
      Array of two numbers. If ```transferVelocity``` is true, gives fluid this initial velocity.
