module DiscreteOrderedMedian

    using DataStructures
    using JuMP
    using MathOptInterface
    using CPLEX
    using LinearAlgebra
    using Random
    using Graphs

    include("instancegen.jl")
    include("struct.jl")
    include("reader.jl")
    include("localsearch.jl")
    include("util.jl")
    include("setcover.jl")
    include("setpacking.jl")
    include("dompk.jl")
    include("domplb.jl")
    include("bnb.jl")
    include("domp2.jl")

    const MOI = MathOptInterface
    
    precompile(setcover, (Matrix{Bool}, Int64, Int64,))
    precompile(setpacking, (Matrix{Bool}, Int64, Int64,))
    precompile(dompk_pos, (DOMPData, BbNode, Int64, Int64, Int64, Vector{Int64},))
    precompile(dompk_neg, (DOMPData, BbNode, Int64, Int64, Int64, Vector{Int64},))
    precompile(build_coverage, (DOMPData, BbNode, Int64, Bool,))
    precompile(domp_lb!, (DOMPData, BbNode, Union{BbNode, Nothing}, Vector{Int64},))
    precompile(read_deleplanque, (String,))
    precompile(compute_sorted_distances, (DOMPData, Vector{Int64},))
    precompile(compute_weighted_cost, (DOMPData, Vector{Int64},))  
    precompile(modify_lambda, (DOMPData, Symbol,))
    precompile(strong_branching, (DOMPData, BbNode, Vector{Float64}, Vector{Float64},))
    precompile(can_recycle_solution, (BbNode, Vector{Int64},))
    precompile(iterated_local_search, (DOMPData, Vector{Vector{Int64}},))
    precompile(local_search!, (DOMPData, Int64, Vector{Int64}, Vector{Int64},))
    precompile(reduce_data_for_ls, (DOMPData, Vector{Vector{Int64}},))
    precompile(generate_euclidean, (Int64, Int64, Symbol,))
    precompile(generate_rand, (Int64, Int64, Symbol,))
end # module
