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
## Direcory structure
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
