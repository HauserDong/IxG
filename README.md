# INSATxGCS (IxG): Accelerating Planning On Graphs of Convex Sets by Interleaving Graph Search and Trajectory Optimization

INterleaved Search And Trajectory optimization (INSAT) [[1]](https://arxiv.org/abs/2101.12548), [[2]](https://arxiv.org/abs/2210.08627) is an algorithmic framework that combines graph search with trajectory optimization for kinodynamic planning in continuous non-convex spaces. Graphs of Convex Sets (GCS) [[3]](https://arxiv.org/abs/2101.11565), [[4]](https://arxiv.org/abs/2205.04422) is an optimization-based planning algorithm that optimizes over a graph made of convex sets as nodes to find continuous trajectories. 

IxG is a GCS accelerator that proposes INSAT as a solver for planning on the graph of convex sets.

> Since [the original repository](https://github.com/nrkumar93/ixg_deprecated) of IxG is deprecated and can not run properly, I modified it to be more user-friendly.

# Get Started

## Install Drake
First, install Drake through [APT](https://drake.mit.edu/apt.html). 

> If you want to use it through Python, you can install it through [Pip](https://drake.mit.edu/pip.html). 

> C++ API: https://drake.mit.edu/doxygen_cxx/index.html

> Python APT: https://drake.mit.edu/pydrake/index.html

## Set your mosek license
Get a [mosek](https://www.mosek.com/) license first.

Add the following line to your `~/.bashrc` file, so your code can find the license:
```
export MOSEKLM_LICENSE_FILE=/YOUR_MOSEK_FOLDER/mosek/mosek.lic
```
where YOUR_MOSEK_FOLDER is your folder to place the mosek license.

## Build this project
```
mkdir build
cd build
cmake ..
make -j4
```
After building it, you can see executables in your `build` folder.

## Running an example
In your `build` folder, run `run_insatxgcs` to test IxG algorithm.

1.  Bimanual environment:
```
./run_insatxgcs insatxgcs bimanual
```

2.  2D maze environment:
```
./run_insatxgcs insatxgcs maze2d
```