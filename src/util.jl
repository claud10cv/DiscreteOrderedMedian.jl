function compute_sorted_distances(data::DOMPData, x::Vector{Int64})::Vector{Int64}
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    open = [j for j in 1 : ncols if x[j] == 1]
    @assert(!isempty(open), "empty vector!")
    d = [(!isempty(open) ? minimum(data.D[i, open]) : round(Int64, 1e-4 * typemax(Int64))) for i in 1 : nrows]
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
    elseif ltype == :T2A
        lambda[end] = -1
    elseif ltype == :T2B
        lambda[p + 1] = 1
    elseif ltype == :T2C
        lambda[p + 1] = -1
    elseif ltype == :T2D
        lambda[round(Int64, n / 2)] = 1
    elseif ltype == :T2E
        lambda[round(Int64, n / 2)] = -1
    elseif ltype == :T2F
        lambda[ceil(Int64, 3 * n / 4)] = 1
    elseif ltype == :T2G
        lambda[ceil(Int64, 3 * n / 4)] = -1
    elseif ltype == :T3 #centrum
        for i in ceil(Int64, n / 2) + 1 : n
            lambda[i] = 1
        end
    elseif ltype == :T3A #complement centrum
        for i in 1 : ceil(Int64, n / 2)
            lambda[i] = 1
        end
    elseif ltype == :T4 #k1-k2-trimmed mean, 25% each side
        t0 = floor(Int64, n / 4) + 1
        tf = t0 + round(Int64, n / 2) - 1
        for i in t0 : tf
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
    elseif ltype == :T8A # alternating ending in (0, 1, 0)
        for i in n : -1 : 1
            if (i + 1) % 3 == n % 3
                lambda[i] = 1
            else lambda[i] = 0
            end
        end
    elseif ltype == :T8B # alternating ending in (1, 0, 0)
        for i in n : -1 : 1
            if (i + 2) % 3 == n % 3
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
    elseif ltype == :T11 # sort in ascending order
        lambda = sort(data.lambda)
    elseif ltype == :T12 # (-1, ..., -1)
        for i in 1 : n
            lambda[i] = -1
        end
    elseif ltype == :T13 # (1, -1, 1, -1, ...)
        for i in 1 : n
            lambda[i] = 2 * (i % 2) - 1
        end
    elseif ltype == :T14 # (-1, 1, -1, 1, ...)
        for i in 1 : n
            lambda[i] = 2 * ((i + 1)% 2) - 1
        end
    end
    uD = unique(sort(vec(D)))
    return DOMPData(D, p, lambda, uD)
end

function find_index_uD(data::DOMPData, d::Int64)::Int64
    l = 1
    u = size(data.uD, 1)
    while l < u
        m = ceil(Int64, (l + u) / 2)
        if data.uD[m] == d
            return m
        elseif data.uD[m] < d
            l = m
        else
            u = m - 1
        end
    end
    @assert(l == u, "binary search gone too far")
    return l
end

function strong_branching(data::DOMPData, params::Parameters, parent::BbNode, ub::Int64, pseudo::Matrix{Float64}, x::Vector{Float64}, z::Vector{Float64})::Vector{Tuple{PseudoCost, Union{Nothing, BbNode}, Union{Nothing, BbNode}}}
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    # println("branching with x = $([y for y in x if abs(y) > 1e-5])")
    # println("branching with z = $([y for y in z if abs(y) > 1e-5])")
    xy = [(i, x[i], z[i], abs(x[i] - round(x[i])), abs(z[i] - round(z[i]))) for i in 1 : nrows]
    filter!(u -> u[4] > 1e-6, xy)
    coefs = [1, 2]
    eps = 1e-5
    # sort!(xy; lt = (u, v) -> u[5] - v[5] > 1e-6 || (abs(u[5] - v[5]) <= 1e-6 && u[4] > v[4]))
    sort!(xy; lt = (u, v) -> (pseudo[u[1], 1] + pseudo[u[1], 2] + eps) * dot(coefs, u[4 : 5]) > (pseudo[v[1], 1] + pseudo[v[1], 2] + eps) * dot(coefs, v[4 : 5]))
    @assert(!isempty(xy), "empty candidate list!, x = $x, parent = $parent")
    # fract = [(i, abs(x[i] - round(x[i])), y) for (i, y) in enumerate(z) if y >= 1e-6 && abs(x[i] - round(x[i])) > 1e-6]
    # empty!(fract)
    # if isempty(fract)
    #     fract = [(i, abs(y - round(y)), y) for (i, y) in enumerate(x)]
    #     filter!(u -> u[2] > 1e-6, fract)
    # end        
    # sort!(fract; lt = (u, v) -> u[2] > v[2])
    candsize = 5
    if size(xy, 1) > candsize
        resize!(xy, candsize)
    end
    absgap = ub - parent.lb
    pcosts = Tuple{PseudoCost, Union{Nothing, BbNode}, Union{Nothing, BbNode}}[]
    fathomed = zeros(Int64, ncols)#falses(ncols)
    for (i, xi, zi, vi, wi) in xy
        minpc = minimum(@view pseudo[i, :])
        # println("min pc = $minpc")
        if minpc >= 1 * absgap
            push!(pcosts, (PseudoCost(i, pseudo[i, 1], pseudo[i, 2]), nothing, nothing))
        else
            lbound = rbound = 0
            newbbnode = BbNode[]
            for b in [BranchInfo(i, 'L', 0), BranchInfo(i, 'G', 1)]
                new_branch = deepcopy(parent.branches)
                new_ropt = deepcopy(parent.ropt)
                push!(new_branch, b)
                new_bbnode = BbNode(new_branch, 0, 0, ([], []), [], [Int64[] for k in 1 : nrows], new_ropt)
                new_xub = deepcopy(parent.xub)
                domp_lb!(data, params, new_bbnode, parent, new_xub)
                if b.sense == 'G'
                    rbound = new_bbnode.lb 
                else
                    lbound = new_bbnode.lb
                end
                push!(newbbnode, new_bbnode)
            end
            pcdown = (lbound - parent.lb) / xi
            pcup = (rbound - parent.lb) / (1 - xi)
            pc = PseudoCost(i, pcdown, pcup)
            pseudo[i, 1] = (pseudo[i, 1] + pcdown) / 2
            pseudo[i, 2] = (pseudo[i, 2] + pcup) / 2
            push!(pcosts, (pc, newbbnode[1], newbbnode[2]))
            if lbound >= ub fathomed[i] += 1 end
            if rbound >= ub fathomed[i] += 1 end
            # if max(lbound, rbound) >= ub 
            #     fathomed[i] = true 
            #     break
            # end
        end
    end
    sort!(pcosts; lt = (u, v) -> fathomed[u[1].i] > fathomed[v[1].i] || (fathomed[u[1].i] == fathomed[v[1].i] && u[1].down + u[1].up > v[1].down + v[1].up))
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

function reassign_vectors(data::DOMPData, x::Vector{Vector{Int64}}, d::Vector{Vector{Int64}}, ropt::Vector{Int64})::Vector{Vector{Int64}}
    return get_extreme_assignment_vectors(data, x, d, ropt)
end

function get_extreme_assignment_vectors(data::DOMPData, x::Vector{Vector{Int64}}, d::Vector{Vector{Int64}}, ropt::Vector{Int64})::Vector{Vector{Int64}}
    nrows = size(data.D, 1)
    rinds = [i for i in 1 : nrows if data.lambda[i] != 0]
    sort!(rinds; lt = (u, v) -> score(data, u, ropt) > score(data, v, ropt))
    y = [Int64[] for i in 1 : nrows]
    for (k, i) in enumerate(rinds)
        if !isempty(y[i]) continue
        else
            y[i] = x[i]
            for j in @view rinds[k + 1 : end]
                if isempty(y[j]) && ropt[j] == d[i][j]
                    y[j] = x[i]
                end
            end
        end
    end
    # println("about to return y = $y")
    return y
end

function get_assignment_minimum_support(data::DOMPData, x::Vector{Vector{Int64}}, d::Vector{Vector{Int64}}, ropt::Vector{Int64})::Vector{Vector{Int64}}
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    rinds = [i for i in 1 : nrows if data.lambda[i] != 0]
    nrinds = length(rinds)
    isfeas = falses(nrinds, nrinds)
    for i in 1 : nrinds, j in 1 : nrinds
        ii, jj = rinds[i], rinds[j]
        isfeas[i, j] = d[jj][ii] == ropt[ii] # the j-th sol provides the correct bound for the row i
    end
    nottaken = trues(nrinds)
    isassigned = falses(nrinds)
    support = falses(ncols)
    y = [Int64[] for i in 1 : nrows]
    let # take first vector as the one that covers the largest number of rows
        v, j = findmax([sum(@view isfeas[:, j]) for j in 1 : nrinds])
        jj = rinds[j]
        for i in 1 : nrinds
            if isfeas[i, j]
                ii = rinds[i]
                y[ii] = x[jj]
                isassigned[i] = true
            end
        end
        nottaken[j] = false
        support = support .| (x[jj] .> 0)
    end
    while !all(isassigned)
        vec = [(j, sum(support .| (x[rinds[j]] .> 0)), sum(isfeas[i, j] for i in 1 : nrinds if !isassigned[i])) for j in 1 : nrinds if nottaken[j]]
        sort!(vec; lt = (u, v) -> u[2] < v[2] || (u[2] == v[2] && u[3] > v[3]))
        for (j, s, f) in vec
            jj = rinds[j]
            if f > 0
                for i in 1 : nrinds
                    if isfeas[i, j] && !isassigned[i]
                        ii = rinds[i]
                        y[ii] = x[jj]
                        isassigned[i] = true
                    end
                end
                nottaken[j] = false
                support = support .| (x[jj] .> 0)
                break
            end
        end
    end
    return y
end

function score(data::DOMPData, k::Int64, d::Vector{Int64})::Int64
    nrows = size(data.D, 1)
    m = 1#abs(d[k] * data.lambda[k])
    if k <= data.p return 0
    elseif data.lambda[k] > 0 return m * (k - data.p + 1)
    elseif data.lambda[k] < 0 return m * (nrows - k + 1)
    else return 0
    end
end