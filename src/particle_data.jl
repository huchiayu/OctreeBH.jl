abstract type AbstractData{N,T<:Real} end

mutable struct Data{N,T<:Real} <: AbstractData{N,T}
	pos::SVector{N,T}
	idx::Int64
	hsml::T
	mass::T
end

Data{N,T}() where {N,T} = Data{N,T}(zero(SVector{N,T}),0,0,0)

function assign_node_data!(n::D, old::D, new::D) where {D<:AbstractData}
	masssum = old.mass + new.mass
	inv_masssum = 1.0 / masssum
	n.pos  = (old.mass .* old.pos .+ new.mass .* new.pos) .* inv_masssum
	n.mass = masssum #this line has to be put after insertpart!, otherwise we're inserting back the wrong node.mass!!!
	n.hsml = old.hsml > new.hsml ? old.hsml : new.hsml
	assign_additional_node_data!(n, old, new)
end

function assign_additional_node_data!(n::Data, old::Data, new::Data)
end
