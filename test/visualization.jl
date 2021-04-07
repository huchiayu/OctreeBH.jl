using OctreeBH
using PyPlot

function plot_quadtree(node::Node{N,T}, ix, iy, ax) where {N,T}
	@assert N>=2
	if isLeaf(node)
    	xmin = node.center[ix] - 0.5*node.length[ix]
    	xmax = node.center[ix] + 0.5*node.length[ix]
    	ymin = node.center[iy] - 0.5*node.length[iy]
    	ymax = node.center[iy] + 0.5*node.length[iy]
    	color="grey"
    	ax.plot([xmin,xmin],[ymin,ymax], c=color)
    	ax.plot([xmin,xmax],[ymin,ymin], c=color)
    	ax.plot([xmax,xmax],[ymin,ymax], c=color)
    	ax.plot([xmin,xmax],[ymax,ymax], c=color)
	else
		@inbounds for i in 1:2^N
	        #println("open this node")
			plot_quadtree(node.child[i], ix, iy, ax)
	    end
	end
end

function draw_circle_nonperiodic(x, radius, ix, iy, ax)
    circle = plt.Circle((x[ix]            , x[iy])            , radius, fill=false, color="tab:red")
    ax.add_artist(circle)
end

function draw_circle_periodic(x, radius, ix, iy, boxsizes, ax)
	#draw a circle at point x and the four nearest periodic images
    circle = plt.Circle((x[ix]            , x[iy])            , radius, fill=false, color="tab:red")
    ax.add_artist(circle)
    circle = plt.Circle((x[ix]+boxsizes[1], x[iy])            , radius, fill=false, color="tab:red")
    ax.add_artist(circle)
    circle = plt.Circle((x[ix]-boxsizes[1], x[iy])            , radius, fill=false, color="tab:red")
    ax.add_artist(circle)
    circle = plt.Circle((x[ix]            , x[iy]+boxsizes[2]), radius, fill=false, color="tab:red")
    ax.add_artist(circle)
    circle = plt.Circle((x[ix]            , x[iy]-boxsizes[2]), radius, fill=false, color="tab:red")
    ax.add_artist(circle)
end
