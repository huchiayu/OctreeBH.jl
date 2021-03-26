#ngb treewalk: recursively find gather ngbs for particle p within the searching radius
function get_ptl_in_box_gather(radius::T, p::SVector{N,T}, node::Node{N,T,D}, boxsizes::SVector{N,T}, idx_ngbs::Vector{Int64}) where {N,T,D<:AbstractData{N,T}}
    if isLeaf(node)
        #println("in a leaf node")
        if node.p != nothing
            #println("nonempty leaf")
	    	dist2 = get_distance2(node.p.pos, p, boxsizes, true)
            if dist2 < radius^2
                #println("push ", node.part)
				push!(idx_ngbs, node.p.idx)
            end
        end
    else
        #println("This is a node... check if node and radius overlap")
		@inbounds for j in 1:N
			if abs(nearest(node.center[j] - p[j], boxsizes[j])) > 0.5*node.length[j] + radius
				return  #no overlap, so node can be skipped!
			end
		end
		#=
		if get_distance2(node.center, p, boxsizes, true) + 0.75 * node.length[1]^2 < radius^2
			append!(idx_ngbs, node.n.idxs)
			return
		end
		=#
		#println("node and radius overlap, so we need to open the node and check its children")
        @inbounds for i in 1:2^N
			get_ptl_in_box_gather(radius, p, node.child[i], boxsizes, idx_ngbs)
        end
    end
end

#ngb treewalk: recursively find scatter ngbs for particle p
function get_ptl_in_box_scatter(p::SVector{N,T}, node::Node{N,T,D}, boxsizes::SVector{N,T}, idx_ngbs::Vector{Int64}) where {N,T,D<:AbstractData{N,T}}
    if isLeaf(node)
        #println("in a leaf node")
		if node.p != nothing
            #println("nonempty leaf")
	    	dist2 = get_distance2(node.p.pos, p, boxsizes, true)
            if dist2 < node.p.hsml^2
                #println("push ", node.part)
				push!(idx_ngbs, node.p.idx)
            end
        end
    else
		#println("This is a node... check if node and radius overlap")
		@inbounds for j in 1:N
			if abs(nearest(node.center[j] - p[j], boxsizes[j])) > 0.5*node.length[j] + node.n.hsml
				return  #no overlap, so node can be skipped!
			end
		end
		#println("node and radius overlap, so we need to open the node and check its children")
        @inbounds for i in 1:2^N
			get_ptl_in_box_scatter(p, node.child[i], boxsizes, idx_ngbs)
        end
    end
end

#driver function to get gather ngbs for particle x within the searching radius h
function get_gather_ngb_tree(x::SVector{N,T}, h::T, node::Node{N,T,D}, boxsizes::SVector{N,T}) where {N,T,D<:AbstractData{N,T}}
    idx_ngbs = Int64[]
    get_ptl_in_box_gather(h, x, node, boxsizes, idx_ngbs)
    return idx_ngbs
end

#driver function to get scatter ngbs for particle x; no need to specify a searching radius
function get_scatter_ngb_tree(x::SVector{N,T}, node::Node{N,T,D}, boxsizes::SVector{N,T}) where {N,T,D<:AbstractData{N,T}}
    idx_ngbs = Int64[]
    get_ptl_in_box_scatter(x, node, boxsizes, idx_ngbs)
    return idx_ngbs
end
