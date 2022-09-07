function compute_sorted_distances(data::DOMPData, x::Vector{Int64})::Vector{Int64}
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    open = [j for j in 1 : ncols if x[j] == 1]
    d = [(!isempty(open) ? minimum(data.D[i, open]) : 1000000000000) for i in 1 : nrows]
    return sort(d)
end

function compute_weighted_cost(data::DOMPData, d::Vector{Int64})
    return dot(data.lambda, d)
end

function find_most_fractional(x::Vector{Float64}, z::Vector{Float64})
    xmin, j = findmax([abs(y - round(y)) for y in z])
    if abs(xmin) < 1e-2
        xmin, j = findmax(abs.(x - round.(x)))
    end
    return xmin, j
end

function modify_lambda(data::DOMPData, ltype::Symbol)::DOMPData
    D = deepcopy(data.D)
    n = size(D, 1)
    p = data.p
    lambda = zeros(Int64, n)
    if ltype == :T0
        lambda = data.lambda
        # do nothing
    elseif ltype == :T1 #median
        for i in 1 : n
            lambda[i] = 1
        end
    elseif ltype == :T2 #centre
        lambda[end] = 1
    elseif ltype == :T3 #centrum
        for i in ceil(Int64, n / 2) : n
            lambda[i] = 1
        end
    elseif ltype == :T4 #k1-k2-trimmed mean, 25% each side
        for i in ceil(Int64, n / 4) : floor(Int64, 3 * n / 4)
            lambda[i] = 1
        end
    elseif ltype == :T5 # alternating ending in (0, 1)
        for i in n : -1 : 1
            if i % 2 == n % 2
                lambda[i] = 1
            else lambda[i] = 0
            end
        end
    elseif ltype == :T6 #alternating ending in (1, 0)
        for i in n : -1 : 1
            if i % 2 == n % 2
                lambda[i] = 0
            else lambda[i] = 1
            end
        end
    elseif ltype == :T7 # alternating ending in (0, 1, 1)
        for i in n : -1 : 1
            if i % 3 == n % 3
                lambda[i] = 1
            elseif i % 3 == (n - 1) % 3
                lambda[i] = 1
            else lambda[i] = 0
            end
        end
    elseif ltype == :T8 # alternating ending in (0, 0, 1)
        for i in n : -1 : 1
            if i % 3 == n % 3
                lambda[i] = 1
            else lambda[i] = 0
            end
        end
    elseif ltype == :T9 # (-1, 0,...,0, 1)
        lambda[1] = -1
        lambda[end] = 1
    elseif ltype == :T10 #(-1, a,...,a, 1)
        lambda[1] = -4
        lambda[end] = 4
        for i in 2 : n - 1
            lambda[i] = 1
        end
    end
    return DOMPData(D, p, lambda)
end

function strong_branching(data::DOMPData, parent::BbNode, pseudo::Matrix{Float64}, x::Vector{Float64}, z::Vector{Float64})::Vector{PseudoCost}
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    # println("branching with x = $([y for y in x if abs(y) > 1e-5])")
    # println("branching with z = $([y for y in z if abs(y) > 1e-5])")
    fract = [(i, abs(y - round(y)), y) for (i, y) in enumerate(z)]
    filter!(u -> u[2] > 1e-6, fract)
    # empty!(fract)
    if isempty(fract)
        fract = [(i, abs(y - round(y)), y) for (i, y) in enumerate(x)]
        filter!(u -> u[2] > 1e-6, fract)
    end        
    sort!(fract; lt = (u, v) -> u[2] > v[2])
    if size(fract, 1) > 5
        resize!(fract, 5)
    end
    pcosts = PseudoCost[]
    for (i, d, v) in fract
        if minimum(@view pseudo[i, :]) >= 1
            push!(pcosts, PseudoCost(i, pseudo[i, 1], pseudo[i, 2]))
        else
            lbound = rbound = 0
            for b in [BranchInfo(i, 'L', 0), BranchInfo(i, 'G', 1)]
                new_branch = deepcopy(parent.branches)
                new_ropt = deepcopy(parent.ropt)
                push!(new_branch, b)
                new_bbnode = BbNode(new_branch, 0, 0, ([], []), [], [Int64[] for k in 1 : nrows], new_ropt)
                new_xub = deepcopy(parent.xub)
                domp_lb!(data, new_bbnode, parent, new_xub)
                if b.sense == 'G'
                    rbound = new_bbnode.lb
                else
                    lbound = new_bbnode.lb
                end
            end
            pcdown = (lbound - parent.lb) / v
            pcup = (rbound - parent.lb) / (1 - v)
            pc = PseudoCost(i, pcdown, pcup)
            pseudo[i, 1] = (pseudo[i, 1] + pcdown) / 2
            pseudo[i, 2] = (pseudo[i, 2] + pcup) / 2
            push!(pcosts, pc)
        end
    end
    sort!(pcosts; lt = (u, v) -> min(u.down, u.up) > min(v.down, v.up) || (min(u.down, u.up) == min(v.down, v.up) && max(u.down, u.up) > max(v.down, v.up)))
    return pcosts
end

function can_recycle_solution(bbnode::BbNode, xk::Vector{Int64})::Bool
    for b in bbnode.branches
        j = b.j
        sense = b.sense
        bound = b.bound
        if sense == 'G' && bound == 1
            if xk[j] != 1 return false
            end
        elseif sense == 'L' && bound == 0
            if xk[j] != 0 return false
            end
        else @assert(false, ("unrecognized branching decision!"))
        end
    end
    return true
end