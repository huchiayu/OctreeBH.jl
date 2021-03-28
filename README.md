# OctreeBH

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://huchiayu.github.io/OctreeBH.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://huchiayu.github.io/OctreeBH.jl/dev)
[![Build Status](https://github.com/huchiayu/OctreeBH.jl/workflows/CI/badge.svg)](https://github.com/huchiayu/OctreeBH.jl/actions)
[![Coverage](https://codecov.io/gh/huchiayu/OctreeBH.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/huchiayu/OctreeBH.jl)

```OctreeBH``` is an implementation of octree that is suitable for N-body problems using the Barnes–Hut approximation. Namely, tree nodes carry information which can be summed or compared in an efficient way. Neighbor search can be done in either "gather" or "scatter" appraoch, which is particularly useful in [smoothed-particle hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics). Boundary conditions can be open, periodic, or mixed (e.g. periodic in x & y while open in z). Spatial dimension (N) can be any arbitrary postive integer.


# Quick Start
```
using OctreeBH
using StaticArrays

const BOXSIZE = 1.0
const N = 3 #spatial dimension
const Npart = 100000 #number of particles

#desired number of neighbor particles (ngbs) within a search radius (hsml)
const Nngb0 = 32

#searching radius that should contain roughly Nngb0 neighboring particles
const hsml0 = BOXSIZE * (Nngb0/(4*pi/3*Npart))^(1/N)

#randomly distributing particles
X = [@SVector rand(N) for _ in 1:Npart]

topnode_length = @SVector(ones(N)) * BOXSIZE  #actual length of tree
boxsizes       = @SVector(ones(N)) * BOXSIZE  #for periodic B.C.
center = topnode_length .* 0.5

#individual smoothing length (only relevant for scatter ngbs)
δ = 0. #add a small perturbation
hsml = ones(Npart) .* hsml0 .* ( 1 .+ δ .* (rand(Npart) .- 0.5) )

mass = ones(Npart) #particle mass

#build the tree
tree = buildtree(X, hsml, mass, center, topnode_length);

x = @SVector rand(N) #search center

#get indices of the gather neighbors
idx_ngbs_g = get_gather_ngb_tree(x, hsml0, tree, boxsizes)

#get indices of the scatter neighbors
idx_ngbs_s = get_scatter_ngb_tree(x, tree, boxsizes)

#for constant hsml, the two ngbs are identical
print( all(idx_ngbs_g.==idx_ngbs_s) )
```

# Example
N body gravitational system
![movie](https://user-images.githubusercontent.com/23061774/112749075-417aed00-8fc0-11eb-8f18-9793b1e82f57.gif)


# Author
Chia-Yu Hu @ Max Planck Institute for Extraterrestrial Physics (cyhu.astro@gmail.com)


