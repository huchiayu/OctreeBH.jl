using OctreeBH
using StaticArrays

const BOXSIZE = 1.0

const N = 3 #spatial dimension
const T = Float64

const Npart = 1000 #number of particles

const ANGLE = 0.7

#randomly distributing particles
X = [@SVector rand(N) for _ in 1:Npart]

topnode_length = @SVector(ones(N)) * BOXSIZE  #actual length of tree

center = topnode_length .* 0.5

#individual smoothing length (only relevant for scatter ngbs)
hsml = zeros(Npart)

mass = ones(Npart) #particle mass

part = [Data{N,T}(SVector(X[i]), i, hsml[i], mass[i]) for i in eachindex(X)]

#build the tree
tree = buildtree(part, center, topnode_length);

#search center
x = @SVector rand(N)

boxsizes = @SVector(ones(N)) * BOXSIZE  #for periodic B.C.

softening = BOXSIZE / Npart^(1/3.)
ga = GravTreeGather{N,T}()
gravity_treewalk!(ga,x,tree,ANGLE,softening,boxsizes)
@show ga.pot

ga_bf = GravTreeGather{N,T}()
gravity_bruteforce!(ga_bf,x,X,mass,boxsizes)
@show ga_bf.pot
