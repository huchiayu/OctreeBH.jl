module OctreeBH
using StaticArrays

export Node, Data
export buildtree, get_scatter_ngb_tree, get_gather_ngb_tree, nearest

include("particle_data.jl");
include("tree.jl");
include("ngb.jl");


end #module OctreeBH
