module OctreeBH
using StaticArrays

export Node, Data, GravTreeGather
export buildtree, get_scatter_ngb_tree, get_gather_ngb_tree
export gravity_treewalk!, gravity_bruteforce!
export nearest

include("particle_data.jl");
include("tree.jl");
include("ngb.jl");
include("gravity_treewalk.jl");


end #module OctreeBH
