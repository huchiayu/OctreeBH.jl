mutable struct Node{N, T<:Real, D<:AbstractData{N,T}}
    center::SVector{N,T}
    length::SVector{N,T}
	p::Union{D, Nothing}           #initilize when an empty leaf becomes nonempty (0->1)
	n::Union{D, Nothing}           #initilize when a nonempty leaf becomes a node (1->2)
	child::Union{Array{Node{N,T,D},1}, Nothing}  #initilize when a nonempty leaf becomes a node (1->2)
end

Node{N,T,D}(center::SVector{N,T}, length::SVector{N,T}) where {N,T,D<:AbstractData{N,T}} =
Node{N,T,D}(center::SVector{N,T}, length::SVector{N,T}, nothing, nothing, nothing)

function isLeaf(node::Node{N,T,D}) where{N,T,D}
    return node.child == nothing ? true : false
end

function getChildIndex(pos::SVector{N,T}, node::Node{N,T,D}) where {N,T,D<:AbstractData{N,T}}
    idx = 0
    @inbounds for i in 1:N
        if pos[i] > node.center[i]
            idx |= 2^(i-1)
        end
        #println(i, 2^(i-1), idx)
    end
    return idx+1
end

function getoffset_precomp(N::T) where {T}
	offset_child = zeros(T,N,2^N)
	for i in 1:2^N
	    a = bitstring(i-1)
	    #offset = SVector(2*parse(Int, a[end])-1, 2*parse(Int, a[end-1])-1, 2*parse(Int, a[end-2])-1)
	    offset = T[]
	    @inbounds for j in 1:N
	        push!(offset,2*parse(T, a[end-(j-1)])-1)
	    end
	    #println(offset)
		offset_child[:,i] = offset
	end
	SMatrix{N,2^N,T}(offset_child)
end

const offset = getoffset_precomp(3) #precomputed offset array for 3D

getoffset(i::T, N::T) where {T<:Int} = ([2*parse(T, bitstring(i-1)[end-(j-1)])-1 for j in 1:N])

function insertpart!(p::D, node::Node{N,T,D}) where {N,T,D<:AbstractData{N,T}}
    if isLeaf(node)
        if node.p == nothing
            #println("at an empty leaf, insert particle and we're done")
			@assert node.n == nothing
			node.p = p
        else
            #println("at a leaf with a preexisting particle, so we split it")
			#initializing the child nodes and calculate their centers
            node.child = Array{Node{N,T,D},1}(undef,2^N)
            @inbounds for i in 1:2^N
                #println(getoffset(i,N))
				#childcenter = node.center + offset[:,i] * 0.25 .* node.length #precomputed offset. This leads to faster tree build but no flexible dimension
				childcenter = node.center + SVector{N}(getoffset(i,N)) * 0.25 .* node.length
                #println("childcenter=", childcenter, "  getoffset(i,N)=", getoffset(i,N))
				node.child[i] = Node{N,T,D}(childcenter, 0.5 .* node.length)
                #println("center = ", childcenter, "  edge_min = ", childcenter - 0.25*node.length,
                #        "  edge_max = ", childcenter + 0.25*node.length)
            end
			#println("insert back the preexsisting particle...")
			insertpart!(node.p, node.child[getChildIndex(node.p.pos, node)])
			#since it's not a leaf node anymore, initialize node data
			node.n = D()
			#push!(node.n.idxs, node.p.idx, p.idx)
			assign_node_data!(node.n, node.p, p)
            #println("insert the new particle...")
            insertpart!(p, node.child[getChildIndex(p.pos, node)])
        end
    else
        #println("open node")
		#push!(node.n.idxs, p.idx)
		assign_node_data!(node.n, node.n, p)
        insertpart!(p, node.child[getChildIndex(p.pos,node)])
    end
end

#build a tree for given particle locations X
#we distinguish between topnode_length (for the tree) and boxsizes (for PBC wrap)
function buildtree(data::Vector{D}, center::SVector{N,T}, topnode_length::SVector{N,T}) where {N,T,D<:AbstractData{N,T}}

    #the top tree
	tree = Node{N,T,D}(center, topnode_length)

	#insert particles one by one
    @inbounds for i in eachindex(data)
		insertpart!(data[i], tree)
    end
    return tree
end

#PBC wrap: typically used as nearest.(dx) where dx is a separation vector connecting two points
function nearest(x::T, boxsize::T) where {T}
	#can be negative as this is just a vector component
    return (x > 0.5*boxsize) ? (x - boxsize) : ( (x < -0.5*boxsize) ? (x + boxsize) : x )
end

function nearest_plain(x::T) where {T}
    if x < -boxHalf_X
        x = x + boxHalf_X*2
    else
        x = x
    end
    if x > boxHalf_X
        x = x - boxHalf_X*2
    else
        x = x
    end
    return x
end

function get_distance2(a::SVector{N,T}, b::SVector{N,T}, boxsizes::SVector{N,T}, periodic::Bool) where {N,T}
    c = a - b
    if periodic
        c = nearest.(c, boxsizes)
    end
    return sum(c .^ 2)
end
