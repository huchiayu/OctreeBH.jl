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

#search center
x = @SVector rand(N)

#get indices of the gather neighbors
idx_ngbs_g = get_gather_ngb_tree(x, hsml0, tree, boxsizes)

#get indices of the scatter neighbors
idx_ngbs_s = get_scatter_ngb_tree(x, tree, boxsizes)

#for constant hsml, the two ngbs are identical
print( all(idx_ngbs_g.==idx_ngbs_s) )
