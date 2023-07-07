function generate_euclidean(n::Int64, p::Int64, maxx::Int64 = 100, maxy::Int64 = 100, itype::Symbol = :pmedian; seed::Int64 = 0)::DOMPData
    rng = MersenneTwister(seed)
    lambda = rand(rng, 1 : 1000, n)
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
    pos = zeros(Int64, n, 2)
    for i in 1 : n
        pos[i, 1] = rand(rng, 0 : maxx)
        pos[i, 2] = rand(rng, 0 : maxy)
    end
    D = zeros(Int64, n, n)
    for i in 1 : n, j in 1 : n
        diff = pos[i, :] - pos[j, :]
        D[i, j] = ceil(Int64, sqrt(dot(diff, diff)))
        D[i, j] = max(1, D[i, j])
    end
    uD = unique(sort(vec(D)))
    DOMPData(D, p, lambda, uD)
end

function generate_rand(n::Int64, p::Int64, itype::Symbol; seed::Int64 = 0)::DOMPData
    rng = MersenneTwister(seed)
    lambda = rand(rng, 0 : 100, n)
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
    D = rand(rng, 10000 : 100000, n, n)
    uD = unique(sort(vec(D)))
    DOMPData(D, p, lambda, uD)
end