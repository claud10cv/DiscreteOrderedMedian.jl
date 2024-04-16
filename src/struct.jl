struct DOMPData
    D::Matrix{Int64}
    p::Int64
    lambda::Vector{Int64}
    uD::Vector{Int64} # unique values of D, sorted
end

struct BranchInfo
    j::Int64
    sense::Char
    bound::Int64
end

mutable struct BbNode
    branches::Vector{BranchInfo}
    lb::Int64
    ub::Int64
    xlb::Tuple{Vector{Float64}, Vector{Float64}}
    xub::Vector{Int64}
    xk::Vector{Vector{Int64}}
    ropt::Vector{Int64}
end

struct PseudoCost
    i::Int64
    down::Float64
    up::Float64
end

struct Result
    lb::Int64
    ub::Int64
    xlb::Tuple{Vector{Float64}, Vector{Float64}}
    xub::Vector{Int64}
end

struct OptimizerData
    optimizer::DataType
    set_attributes::Function
end

struct Parameters
    warm_starts::Bool
    var_fixing::Bool
    primal_heur::Bool
    symm_break::Bool
    time_limit::Float64
    optimizer_data::OptimizerData
end

BbNode() = BbNode(BranchInfo[], 0, typemax(Int64), (Float64[], Float64[]), Int64[], Vector{Vector{Int64}}(), Int64[])
BbNode(branches::Vector{BranchInfo}) = BbNode(branches, 0, typemax(Int64), (Float64[], Float64[]), Int64[], Vector{Vector{Int64}}(), Int64[])