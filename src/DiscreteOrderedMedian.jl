module DiscreteOrderedMedian
    include("instancegen.jl")
    include("reader.jl")
    include("struct.jl")
    include("util.jl")
    include("setcover.jl")
    include("dompk.jl")
    include("domplb.jl")
    include("bnb.jl")

    precompile(setcover, (Matrix{Bool}, Int64, Int64,))
    precompile(dompk_pos, (DOMPData, BbNode, Int64, Int64, Int64,))
    precompile(build_coverage, (DOMPData, BbNode, Int64,))
    precompile(domp_lb!, (DOMPData, BbNode, Vector{Int64},))
    precompile(read_deleplanque, (String,))
    precompile(compute_sorted_distances, (DOMPData, Vector{Int64},))
    precompile(strong_branching, (DOMPData, BbNode, Vector{Float64}, Vector{Float64},))
end # module
