# OctreeBH

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://huchiayu.github.io/OctreeBH.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://huchiayu.github.io/OctreeBH.jl/dev)
[![Build Status](https://github.com/huchiayu/OctreeBH.jl/workflows/CI/badge.svg)](https://github.com/huchiayu/OctreeBH.jl/actions)
[![Coverage](https://codecov.io/gh/huchiayu/OctreeBH.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/huchiayu/OctreeBH.jl)

```OctreeBH``` is an implementation of octree for solving N-body problems using the [Barnes–Hut approximation](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation). Namely, tree nodes carry information that can be summed or compared in an efficient way. The package provides two main functionalities: (1) finding particles within a given radius (neighbor search) and (2) calculating gravitational acceleration in an N-body system. Neighbor search can be done in either the "gather" or "scatter" approach, which is particularly useful in [smoothed-particle hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics). Gravitational acceleration adopts the monopole approximation, i.e., distant particles are clustered together and treated as a point mass located at their center of mass. Boundary conditions can be open, periodic, or mixed (e.g. periodic in x & y while open in z). Spatial dimension (N) can be any arbitrary positive integer.


# Quick Start
Finding particles within a searching radius ```hsml0``` centered at ```x```: 
```
using OctreeBH
using StaticArrays

N = 3 #spatial dimension
T = Float64

Npart = 100000 #number of particles
hsml0 = 0.04 #searching radius
boxsizes = @SVector(ones(N)) #for periodic B.C.
topnode_length = @SVector(ones(N)) #size of the top tree node
center = topnode_length .* 0.5 #center of the top tree ndoe

hsml = ones(Npart) .* hsml0 #particle smoothing length 
mass = ones(Npart) #particle mass

#randomly distributing particles
X = [@SVector rand(N) for _ in 1:Npart]

#store particle data into the struct "Data" 
part = [Data{N,T}(X[i], i, hsml[i], mass[i]) for i in 1:Npart]

#build the tree
tree = buildtree(part, center, topnode_length);

x = @SVector rand(N) #search center

#get indices of the gather neighbors
idx_ngbs_g = get_gather_ngb_tree(x, hsml0, tree, boxsizes)

#get indices of the scatter neighbors
idx_ngbs_s = get_scatter_ngb_tree(x, tree, boxsizes)

#for constant hsml, the two ngbs are identical
@show length(idx_ngbs_g), length(idx_ngbs_s)
```

# Example
The example code visualization.jl visualizes the (quad)tree structure for a given particle configuration in 2D. It also illustrates the difference between the gather neighbors (lower left) and scatter neighbors (upper right), which are identical in this particular case as the smoothing length (search radius) is constant.

![vis_tree](https://user-images.githubusercontent.com/23061774/113936917-c0f69080-97f8-11eb-880e-4fd8019b9f49.png)


The example code [nbody.jl](https://github.com/huchiayu/OctreeBH.jl/blob/main/test/nbody.jl) solves a gravitational N-body system using the leapfrog integration scheme.

![movie](https://user-images.githubusercontent.com/23061774/112749075-417aed00-8fc0-11eb-8f18-9793b1e82f57.gif)


# Scaling 
OctreeBH.jl follows the N*log(N) scaling as opposed to the naive method that scales as N^2

![scaling_octree](https://user-images.githubusercontent.com/23061774/178747692-e877d0e5-59a6-46aa-a494-4f53ae19eab2.png)


# Author
Chia-Yu Hu @ Max Planck Institute for Extraterrestrial Physics (cyhu.astro@gmail.com)


