using OctreeBH
using StaticArrays

const BOXSIZE = 1.0

const T = Float64
const N = 2 #spatial dimension
Npart = 300 #number of particles
hsml0 = 0.15 #searching radius

@testset "balance" begin

using Random
Random.seed!(1114)

#randomly distributing particles
X = [@SVector rand(N) for _ in 1:Npart]

topnode_length = @SVector(ones(N)) * BOXSIZE  #actual length of tree
center = topnode_length .* 0.5

hsml = ones(Npart) .* hsml0
mass = ones(Npart) #particle mass
part = [Data{N,T}(SVector(X[i]), i, hsml[i], mass[i]) for i in eachindex(X)]

#build the tree
@time tree = buildtree(part, center, topnode_length);

max_depth_0 = get_max_tree_depth(tree)

#set_max_depth_AMR!(tree, max_depth)

boxsizes = @SVector(ones(N)) * BOXSIZE  #for periodic B.C.


for i in max_depth_0:-1:7
	set_max_depth_AMR!(tree, i)
	@test i == get_max_tree_depth(tree)

	balance_all_level!(tree)
	@test i == get_max_tree_depth(tree)
end

end
