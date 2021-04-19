module OctreeBH
using StaticArrays

export Node, Data, AbstractData, GravTreeGather
export buildtree, get_scatter_ngb_tree, get_gather_ngb_tree
export gravity_treewalk!, gravity_bruteforce!
export nearest
export isLeaf
export get_max_tree_depth, set_max_depth_AMR!, balance_all_level!

include("particle_data.jl");
include("tree.jl");
include("ngb.jl");
include("gravity_treewalk.jl");
include("balance.jl")

end #module OctreeBH
