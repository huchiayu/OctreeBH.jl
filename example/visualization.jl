using OctreeBH
using StaticArrays
using PyPlot
#plt.style.use("dark_background")
plt.style.use("default")


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

function draw_circle_nonperiodic(x, radius, ix, iy, cc, ax)
    circle = plt.Circle((x[ix]            , x[iy])            , radius, fill=false, color="tab:red")
    ax.add_artist(circle)
end

function draw_circle_periodic(x, radius, ix, iy, boxsizes, cc, ax)
	#draw a circle at point x and the four nearest periodic images
    circle = plt.Circle((x[ix]            , x[iy])            , radius, fill=false, color=cc)
    ax.add_artist(circle)
    circle = plt.Circle((x[ix]+boxsizes[1], x[iy])            , radius, fill=false, color=cc)
    ax.add_artist(circle)
    circle = plt.Circle((x[ix]-boxsizes[1], x[iy])            , radius, fill=false, color=cc)
    ax.add_artist(circle)
    circle = plt.Circle((x[ix]            , x[iy]+boxsizes[2]), radius, fill=false, color=cc)
    ax.add_artist(circle)
    circle = plt.Circle((x[ix]            , x[iy]-boxsizes[2]), radius, fill=false, color=cc)
    ax.add_artist(circle)
end


const BOXSIZE = 1.0

const T = Float64
const N = 2 #spatial dimension
const Npart = 100 #number of particles
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
tree = buildtree(part, center, topnode_length);

boxsizes = @SVector(ones(N)) * BOXSIZE  #for periodic B.C.


clf()
fig, ax = subplots(1, 1, figsize=(6, 6))
ms = 6
ix, iy = 1, 2

ax.plot(getindex.(X,ix), getindex.(X,iy), ".", c="tab:blue", ms=ms)

@time plot_tree_grid(tree, ix, iy, "grey", ax)

X0 = @SVector[0.1, 0.2]
ax.plot(X0[ix], X0[iy], "*", c="tab:green", ms=ms)
idx_ngbs = get_gather_ngb_tree(X0, hsml0, tree, boxsizes)
ax.plot(getindex.(X[idx_ngbs],1), getindex.(X[idx_ngbs],2), ".", c="tab:red", ms=ms)
draw_circle_periodic(X0, hsml0, ix, iy, boxsizes, "tab:green", ax)

X1 = @SVector[0.7, 0.8]
ax.plot(X1[ix], X1[iy], "*", c="tab:green", ms=ms)
idx_ngbs = get_scatter_ngb_tree(X1, tree, boxsizes)
ax.plot(getindex.(X[idx_ngbs],1), getindex.(X[idx_ngbs],2), ".", c="tab:red", ms=ms)
draw_circle_periodic(X1, hsml0, ix, iy, boxsizes, "tab:green", ax)
for i in eachindex(idx_ngbs)
	draw_circle_periodic(X[idx_ngbs[i]], hsml[idx_ngbs[i]], ix, iy, boxsizes, "tab:red", ax)
end

ax.axis([0,boxsizes[ix],0,boxsizes[iy]])
ax.set_aspect(1)
