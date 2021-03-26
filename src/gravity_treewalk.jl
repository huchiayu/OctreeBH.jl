mutable struct GravTreeGather{N,T}
	acc::SVector{N,T}
	pot::T
    mass::T
    #nodecenters::Vector{Vector{T}} #for debugging/visualization
    #nodelengths::Vector{Vector{T}} #for debugging/visualization
end

GravTreeGather{N,T}() where {N,T} = GravTreeGather{N,T}(zero(SVector{N,T}),0,0)

#when changing the argument of treewalk, remember to also change the recursive treewalk call inside the function body!
function gravity_treewalk!(ga::GravTreeGather{N,T}, p::SVector{N,T}, node::Node{N,T,D},
	openingangle::T, softening::T, boxsizes::SVector{N,T}) where {N,T,D}
    if isLeaf(node)
        #println("in a leaf node")
        if node.p != nothing
            #println("nonempty leaf")
            ga.mass += node.p.mass
            dx = nearest.(node.p.pos - p, boxsizes)
			r2 = sum(dx.^2)
			#@show r2, node.p.pos, p
			if r2 > 0.0 #exclude self contribution (r2=0)
				r_inv = sqrt(1.0 / (r2 + softening^2))
				r3_inv = r_inv^3
				ga.acc += node.p.mass * r3_inv * dx
				ga.pot += node.p.mass * r_inv
	            #push!(ga.nodecenters, node.p.pos)
	            #push!(ga.nodelengths, zeros(N))
			end
            #@show "it's a particle..." node.center, node.length
            #r2 = get_distance2(node.part, p, boxsizes, true)
        end
    else
        #println("This is a node... check the opening angle")
		#r2 = get_distance2(pos_c, p, boxsizes, true)
		#don't use get_distance2 as we need both r2 and dx (for vec2pix)
		dx = nearest.(node.n.pos - p, boxsizes)
		r2 = sum(dx.^2)
		if r2 > (node.length[1] / openingangle)^2
			#println("skip node ", i)
			ga.mass += node.n.mass
			r_inv = sqrt(1.0 / (r2 + softening^2))
			r3_inv = r_inv^3
			ga.acc += node.n.mass * r3_inv * dx
			ga.pot += node.n.mass * r_inv
			#push!(ga.nodecenters, node.center)
			#push!(ga.nodelengths, node.length)
			#@show "use this node!" node.center, node.length
		else
	        @inbounds for i in 1:2^N
	            #println("open this node")
	            gravity_treewalk!(ga, p, node.child[i], openingangle, softening, boxsizes)
	        end
		end
    end
end

function gravity_bruteforce!(ga::GravTreeGather{N,T}, p::SVector{N,T}, X::Vector{SVector{N,T}},
							 mass::Vector{T}, boxsizes::SVector{N,T}) where {N,T}
	for i in eachindex(X)
		ga.mass += mass[i]
		dx = nearest.(X[i] - p, boxsizes)
		r2 = sum(dx.^2)
		if r2 > 0.0 #exclude self contribution (r2=0)
			r_inv = sqrt(1.0 / r2)
			r3_inv = r_inv^3
			ga.acc += mass[i] * r3_inv * dx
			ga.pot += mass[i] * r_inv
		end
	end
end
