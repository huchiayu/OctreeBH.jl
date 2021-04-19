using OctreeBH
import OctreeBH: AbstractData, assign_additional_node_data!

using StaticArrays

mutable struct Data2{N,T<:Real} <: AbstractData{N,T}
	pos::SVector{N,T}
	idx::Int64
	hsml::T
	mass::T
	#define additional members below
	massH2::T
end

Data2{N,T}() where {N,T} = Data2{N,T}(zero(SVector{N,T}),0,0,0,0)

function assign_additional_node_data!(n::Data2, old::Data2, new::Data2)
	n.massH2 = old.massH2 + new.massH2
end

const BOXSIZE = 1.0

const N = 3 #spatial dimension
const T = Float64

const Npart = 100000 #number of particles

#desired number of neighbor particles (ngbs) within a search radius (hsml)
const Nngb0 = 32

#searching radius that should contain roughly Nngb0 neighboring particles
const hsml0 = BOXSIZE * (Nngb0/(4*pi/3*Npart))^(1/N)

#randomly distributing particles
X = [@SVector rand(N) for _ in 1:Npart]

topnode_length = @SVector(ones(N)) * BOXSIZE  #actual length of tree

center = topnode_length .* 0.5

#individual smoothing length (only relevant for scatter ngbs)
δ = 0.0 #add a small perturbation
hsml = ones(Npart) .* hsml0 .* ( 1 .+ δ .* (rand(Npart) .- 0.5) )

mass = ones(Npart) #particle mass

part = [Data2{N,T}(SVector(X[i]), i, hsml[i], mass[i], 0.6*mass[i]) for i in eachindex(X)]

#build the tree
tree = buildtree(part, center, topnode_length);

#search center
x = @SVector rand(N)

boxsizes = @SVector(ones(N)) * BOXSIZE  #for periodic B.C.

#get indices of the gather neighbors
idx_ngbs_g = get_gather_ngb_tree(x, hsml0, tree, boxsizes)

#get indices of the scatter neighbors
idx_ngbs_s = get_scatter_ngb_tree(x, tree, boxsizes)

#for constant hsml, the two ngbs are identical
@show length(idx_ngbs_g), length(idx_ngbs_s)
