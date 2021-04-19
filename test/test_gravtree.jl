using OctreeBH
using StaticArrays
using Test

const BOXSIZE = 1.0

const N = 3 #spatial dimension
const T = Float64

Npart = 1000 #number of particles

@testset "gravtree" begin

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
#x = @SVector rand(N)

boxsizes = @SVector(ones(N)) * BOXSIZE  #for periodic B.C.

#softening = BOXSIZE / Npart^(1/3.)
softening = 0.0
ANGLE = 0.0

pot = zeros(Npart)
pot_bf = zeros(Npart)
for i in eachindex(X)
    ga = GravTreeGather{N,T}()
    gravity_treewalk!(ga,X[i],tree,ANGLE,softening,boxsizes)
    pot[i] = ga.pot

    ga_bf = GravTreeGather{N,T}()
    gravity_bruteforce!(ga_bf,X[i],X,mass,boxsizes)
    pot_bf[i] = ga_bf.pot
end
@test all(pot .≈ pot_bf)

end
