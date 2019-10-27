# This Template Merges two frameworks.  

On one hand, we have a bunch of projects built on Etienne's physics starter code: https://github.com/evouga/libigl-example-physics-project

On the other Polyscope is really shiny. http://polyscope.run/

My goal here was to produce a minimal front end for porting over our projects to polyscope.  Perhaps it might also be useful for someone else as well!     

Here are the readme's of the two projects, should clean them up a bit at some point.

# libIGL & Polyscope Example Project

Demonstrates basic functionality and project setup for using Polyscope with libIGL.

IGL code mimicks examples from the IGL tutorials.

Need to update this line in the CMakeLists.txt
```
"${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl/external/eigen"
```


## To download and build

```
git clone --recurse-submodules https://github.com/the13fools/geometry-processing-starter-kit.git
cd geometry-processing-starter-kit
mkdir build
cd build
cmake ..
make -j4

./bin/libIGL-polyscope-example ../bunnyhead.obj

./bin/libIGL-polyscope-example ../spot.obj
```

## Usage

After running the demo appliation as above, try clicking around the command UI in the upper right corner.

Notice that as quantities are added to Polyscope, they appear in the selection window on the left.

Try clicking a vertex to see the quantities associated with that vertex. Clicking may be easier if you first show edge by dragging the "edge width" slider in the left panel.

Check out `src/main.cpp` to see how easy it is! There are only 12 lines of Polyscope code in the whole demo to generate this visualization.


# libigl example physics project

A blank project example showing how to set up a simple physical simulation using libigl and cmake. Based on Alec Jacobson's libigl example project. This project contains some boilerplate that sets up a physical simulation to run in its own thread, with rendering provided by libigl.

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `example_bin` binary.

## Run

From within the `build` directory just issue:

    ./example_bin

A glfw app should launch displaying a 3D cube.
