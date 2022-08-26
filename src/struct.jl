struct DOMPData
    D::Matrix{Int64}
    p::Int64
    lambda::Vector{Int64}
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
    ropt::Vector{Int64}
end

BbNode() = BbNode(BranchInfo[], 0, typemax(Int64), (Float64[], Float64[]), Int64[], Int64[])
BbNode(branches::Vector{BranchInfo}) = BbNode(branches, 0, typemax(Int64), (Float64[], Float64[]), Int64[], Int64[])