function compute_sorted_distances(data::DOMPData, x::Vector{Int64})::Vector{Int64}
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    open = [j for j in 1 : ncols if x[j] == 1]
    d = 1000000000000 * ones(Int64, nrows)
    for i in 1 : nrows, j in open
        d[i] = min(d[i], data.D[i, j])
    end
    sort!(d)
    return d
end

function compute_weighted_cost(data::DOMPData, d::Vector{Int64})
    return dot(data.lambda, d)
end

function find_most_fractional(x::Vector{Float64}, z::Vector{Float64})
    xmin, j = findmax([abs(y - round(y)) for y in z])
    if abs(xmin) < 1e-2
        xmin, j = findmax([abs(y - round(y)) for y in x])
    end
    return xmin, j
end

function strong_branching(data::DOMPData, parent::BbNode, x::Vector{Float64}, z::Vector{Float64})
    fract = [(i, abs(y - round(y))) for (i, y) in enumerate(z)]
    filter!(u -> u[2] > 1e-6, fract)
    # empty!(fract)
    if isempty(fract)
        fract = [(i, abs(y - round(y))) for (i, y) in enumerate(x)]
        filter!(u -> u[2] > 1e-6, fract)
    end        
    sort!(fract; lt = (u, v) -> u[2] > v[2])
    if size(fract, 1) > 5
        resize!(fract, 5)
    end
    lb = 0
    besti = Int64[]
    bestd = Float64[]
    for (i, d) in fract
        new_branch = deepcopy(parent.branches)
        new_ropt = deepcopy(parent.ropt)
        push!(new_branch, BranchInfo(i, 'L', 0))
        new_bbnode = BbNode(new_branch, 0, 0, ([], []), [], new_ropt)
        new_xub = deepcopy(parent.xub)
        domp_lb!(data, new_bbnode, new_xub)
        if new_bbnode.lb > lb
            besti = [i]
            bestd = [d]
            lb = new_bbnode.lb
        elseif new_bbnode.lb >= lb
            push!(besti, i)
            push!(bestd, d)
        end
    end
    nbest = length(besti)
    fract = [(besti[k], bestd[k]) for k in 1 : nbest] 
    if nbest > 1
        lb = 0
        besti = []
        bestd = []
        for (i, d) in fract
            new_branch = deepcopy(parent.branches)
            new_ropt = deepcopy(parent.ropt)
            push!(new_branch, BranchInfo(i, 'G', 1))
            new_bbnode = BbNode(new_branch, 0, 0, ([], []), [], new_ropt)
            new_xub = deepcopy(parent.xub)
            domp_lb!(data, new_bbnode, new_xub)
            if new_bbnode.lb > lb
                besti = [i]
                bestd = [d]
                lb = new_bbnode.lb
            elseif new_bbnode.lb >= lb
                push!(besti, i)
                push!(bestd, d)
            end
        end
    end 
    return bestd[1], besti[1]
end