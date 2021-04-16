#NE = [1,1]
#NW = [0,1]
#SE = [1,0]
#SW = [0,0]


#dim = 1: x-direction , or E & W
#dim = 2: y-direction , or N & S

#direction = 0 or 1 (control E or W if dim == 1; N or S if dim ==2)

#check balance level by level using a post-order treewalk (from deepest level up to root)

#return the bit of num at jth digit
getbit(num, j) = ( num >> (j-1) ) & 1

#flip jth digit of num and return
flipbit(num, j) = num ⊻ ( 1 << (j-1) )

const PERIODIC = true


#for j in 1:2
#    for i in 1:2^2
#        @show getbit(i-1,j), bitstring((i-1))[end-8:end], bitstring(flipbit(i-1,j))[end-8:end]
#    end
#end

function balance!(node::Node{N,T,D}, max_depth::Int, root_node_length::T) where{N,T,D}
	depth = Int(log2(round(root_node_length / node.length[1])))

    if isLeaf(node)
		#println("in a leaf node at level ", depth)
		if depth == max_depth #only do it at the deepest level
			if node.ID == 1 || node.ID == 2^N #we only need to check balance for two diagonal nodes
				#println("in a leaf node at level ", depth, "... check balance!")
				for i in 1:N
					#check both directions (0 & 1) along each dimension (i), which is N*2 directions in total
					check_ngb_node_and_refine_if_necessary(node, i, 0, node.length[1])
					check_ngb_node_and_refine_if_necessary(node, i, 1, node.length[1])
				end
			end
		end
    else
        #println("This is a node at level ", depth)
		#if depth >= max_depth and it's still a node,
		#then it's already checked during the previous iterations
		if depth < max_depth
	    	@inbounds for i in 1:2^N
	        	#println("open this node")
				balance!(node.child[i], max_depth, root_node_length)
	    	end
		end
    end
end


#search both directions for each dimension (dim)
function check_ngb_node_and_refine_if_necessary(node, dim, direction, length0)
    #if dim == 2 #find N & S ngbs
    #if direction == 1 #we're looking for the north ngb
    if node.parent == nothing
        #node is root! we're at the boundary! #TODO: peridoic B.C.?
        return nothing
    end

    #node is a south child, so it can give us a north ngb directly from its sibling (we need to distinguish between SW & SE)
    if getbit(node.ID-1, dim) != direction #meaning that node has a sibling along direction
        idx_ngb = flipbit(node.ID-1, dim) #flipbit = opposite direction in dim-th dimension
		#@show node.ID-1, dim, idx_ngb
        return node.parent.child[idx_ngb+1]
    end


    #Okay, node is a north child, so it can't give us a north ngb (it has no north siblings).
    #Therefore, we have to move up to the parent node. If that parent node is still not a south child,
    #then we keep moving up (to the parent's parent) until we find a parent who is a south child or until
    #we reach the root which means that we're at the boundary and there is no ngb along this direction
    uncle = check_ngb_node_and_refine_if_necessary(node.parent, dim, direction, length0) #an uncle is the sibling of a parent (duh)
    if uncle == nothing #root!
		if PERIODIC
			#@show "root!!!"
			uncle = node.parent #set uncle as node's parent and proceed (don't return!!!)
		else
			return uncle
		end
    end

    #if uncle is a leaf, then it's a least 2x bigger than node, so it doesn't matter if node is NW or NE
    #we check if the uncle needs to be refined
    if isLeaf(uncle)
        if uncle.length[1] / length0 > 2.5 #should be just 2 but use 2.5 to avoid round-off error
			#println("uncle.length = ", uncle.length[1], "  length0 = ", length0, "  refine the node!!!")
            refine_node!(uncle)
        else
            return uncle #only differ by a factor of 2 (one level)
        end
    end

    #uncle is not a leaf, so we return the ngb cousin
    idx_ngb = flipbit(node.ID-1, dim) #flipbit = opposite direction in dim-th dimension
    return uncle.child[idx_ngb+1] #SW (bit flip at the dim-th index)
end






function refine_node!(node::Node{N,T,D}) where{N,T,D}
    node.child = Array{Node{N,T,D},1}(undef,2^N)
    @inbounds for i in 1:2^N
        #println(getoffset(i,N))
        #childcenter = node.center + offset[:,i] * 0.25 .* node.length #precomputed offset. This leads to faster tree build but no flexible dimension
        childcenter = node.center + SVector{N}(getoffset(i,N)) * 0.25 .* node.length
        #println("childcenter=", childcenter, "  getoffset(i,N)=", getoffset(i,N))
        node.child[i] = Node{N,T,D}(childcenter, 0.5 .* node.length)
        #println("center = ", childcenter, "  edge_min = ", childcenter - 0.25*node.length,
        #        "  edge_max = ", childcenter + 0.25*node.length)
		node.child[i].parent = node
		node.child[i].ID = i
    end
    node.n = D() #optional as we just need the correct tree structure and the node data will be assign when we do the p2g mapping
end

#if node is a north child, then its southern ngb is its sibling, so the balance is fine, so we don't have to check its southern ngb
#we only need to check its northern ngb




#=

#search both directions for each dimension (dim)
check_ngb_node_and_refine_if_necessary(node, dim, direction, length0)
    if dim == 2 #find N & S ngbs
        if direction == 1 #we're looking for the north ngb
            if node.parent == nothing
                #node is root! we're at the boundary! #TODO: peridoic B.C.?
                return nothing
            end

            #node is a south child, so it can give us a north ngb directly from its sibling (we need to distinguish between SW & SE)
            if node.ID[dim] != direction
            if node.ID == [0,0] #SW
                return node.parent.child[0,1] #NW (bit flip at the dim-th index)
            elseif node.ID == [1,0] #SE
                return node.parent.child[1,1] #NE (bit flip at the dim-th index)
            end
            end

            #Okay, node is a north child, so it can't give us a north ngb (it has no north siblings).
            #Therefore, we have to move up to the parent node. If that parent node is still not a south child,
            #then we keep moving up (to the parent's parent) until we find a parent who is a south child or until
            #we reach the root which means that we're at the boundary and there is no ngb along this direction
            uncle = check_ngb_node_and_refine_if_necessary(node.parent, dim, direction, length0) #an uncle is the sibling of a parent (duh)
            if uncle == nothing #root!
                return uncle
            end

            #if uncle is a leaf, then it's a least 2x bigger than node, so it doesn't matter if node is NW or NE
            #we check if the uncle needs to be refined
            if isLeaf(uncle)
                if uncle.length / length0 > 2.5 #should be just 2 but use 2.5 to avoid round-off error
                    refine_node!(uncle)
                else
                    return uncle #only differ by a factor of 2 (one level)
                end
            end

            #uncle is not a leaf, so we return the ngb cousin
            if node.ID == [0,1] #NW
                return uncle.child[0,0] #SW (bit flip at the dim-th index)
            elseif node.ID == [1,1] #NE
                return uncle.child[1,0] #SE (bit flip at the dim-th index)
            end
        else  #direction == 0 #we're looking for the south ngb
            code...
        end #direction
    elseif dim == 1 #find E & W ngbs
        code...
    end #dim
end

=#
