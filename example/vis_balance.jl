using OctreeBH
using StaticArrays
using PyPlot
#plt.style.use("dark_background")
#plt.style.use("default")


function plot_tree_grid(node::Node{N,T}, ix, iy, cc, ax) where {N,T}
	@assert N>=2
	if isLeaf(node)
    	xmin = node.center[ix] - 0.5*node.length[ix]
    	xmax = node.center[ix] + 0.5*node.length[ix]
    	ymin = node.center[iy] - 0.5*node.length[iy]
    	ymax = node.center[iy] + 0.5*node.length[iy]
    	color=cc
    	ax.plot([xmin,xmin],[ymin,ymax], c=color)
    	ax.plot([xmin,xmax],[ymin,ymin], c=color)
    	ax.plot([xmax,xmax],[ymin,ymax], c=color)
    	ax.plot([xmin,xmax],[ymax,ymax], c=color)
	else
		@inbounds for i in 1:2^N
	        #println("open this node")
			plot_tree_grid(node.child[i], ix, iy, cc, ax)
	    end
	end
end

const BOXSIZE = 1.0

const T = Float64
const N = 2 #spatial dimension
const Npart = 300 #number of particles
const hsml0 = 0.15 #searching radius

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

@show max_depth = get_max_tree_depth(tree)

#set_max_depth_AMR!(tree, max_depth)

boxsizes = @SVector(ones(N)) * BOXSIZE  #for periodic B.C.


ms = 6
ix, iy = 1, 2


clf()
fig, ax = subplots(1, 1, figsize=(8, 8))
ax.axis([0,boxsizes[ix],0,boxsizes[iy]])
ax.set_aspect(1)
#ax.plot(getindex.(X,ix), getindex.(X,iy), ".", c="tab:blue", ms=ms)
@time plot_tree_grid(tree, ix, iy, "grey", ax)
savefig("balance_none.png")

#max_depth = 9
@time for i in max_depth:-1:1
	println("balancing the tree at level ", i)
	OctreeBH.balance!(tree, i, tree.length[1])

	clf()
	fig, ax = subplots(1, 1, figsize=(8, 8))
	ax.axis([0,boxsizes[ix],0,boxsizes[iy]])
	ax.set_aspect(1)
	@time plot_tree_grid(tree, ix, iy, "grey", ax)
	savefig("balance_" * string(i) * ".png")
end
