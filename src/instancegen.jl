using LinearAlgebra
function generate_euclidean(n::Int64, p::Int64, itype::Symbol)::DOMPData
    lambda = rand(0 : 1000, n)
    if itype == :pcentre
        lambda = zeros(Int64, n)
        lambda[end] = 1
    elseif itype == :monotonic
        sort!(lambda)
    elseif itype == :pmedian
        lambda = ones(Int64, n)
    elseif itype == :kcentrum
        lambda = zeros(Int64, n)
        for j in ceil(Int64, n / 2) : n
            lambda[j] = 1
        end
    end
    pos = rand(1 : 10000, n, 2)
    D = zeros(Int64, n, n)
    for i in 1 : n, j in 1 : n
        diff = pos[i, :] - pos[j, :]
        D[i, j] = ceil(Int64, sqrt(dot(diff, diff)))
    end
    DOMPData(D, p, lambda)
end

function generate_rand(n::Int64, p::Int64, itype::Symbol)::DOMPData
    lambda = rand(0 : 100, n)
    if itype == :pcentre
        lambda = zeros(Int64, n)
        lambda[end] = 1
    elseif itype == :monotonic
        sort!(lambda)
    elseif itype == :pmedian
        lambda = ones(Int64, n)
    elseif itype == :kcentrum
        lambda = zeros(Int64, n)
        for j in ceil(Int64, n / 2) : n
            lambda[j] = 1
        end
    elseif itype == :random
    end
    D = rand(10000 : 100000, n, n)
    DOMPData(D, p, lambda)
end