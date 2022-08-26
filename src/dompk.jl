function dompk_pos(data::DOMPData, bbnode::BbNode, k::Int64, lb::Int64, ub::Int64)::Tuple{Int64, Vector{Int64}}
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    open = [b.j for b in bbnode.branches if b.sense == 'G' && b.bound == 1]
    nopen = length(open)
    x = zeros(Int64, ncols)
    y = zeros(Int64, ncols)
    for j in open
        x[j] = 1
    end
    # if k == nrows
    #     println("starting with lb = $lb and ub = $ub")
    # end
    # println("there are $nopen open facilities already")
    while lb < ub
        r = floor(Int64, (lb + ub) / 2)
        cov, map, numcov = build_coverage(data, bbnode, r)
        # println("there are $numcov rows covered already") 
        # println("map = $map")
        sol = setcover(cov, k - numcov, data.p - nopen)
        if !isempty(sol)
            # if k == nrows
            #     println("feasible r = $r")
            # end
            ub = r
            for j in 1 : ncols
                y[j] = 0
            end
            for j in map[sol]
                y[j] = 1
            end
        else
            # if k == nrows
            #     println("infeasible r = $r")
            # end
            lb = r + 1
        end
    end
    return ub, (x + y)
end

function dompk_neg(data::DOMPData, bbnode::BbNode, k::Int64, lb::Int64, ub::Int64)::Int64
    open = [b.j for b in bbnode.branches if b.sense == 'G' && b.bound == 1]
    nopen = length(open)
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    x = zeros(Int64, ncols)
    y = zeros(Int64, ncols)
    for j in open
        x[j] = 1
    end
    while lb < ub
        r = ceil(Int64, (lb + ub) / 2)
        cov, map, numcov = build_coverage(data, bbnode, r)
        sol = setcover(cov, k - numcov, data.p - nopen)
        if !isempty(sol)
            lb = r
        else
            ub = r - 1
            for j in 1 : ncols
                y[j] = 0
            end
            for j in map[sol]
                y[j] = 1
            end
        end
    end
    return lb, (x + y)
end

function build_coverage(data::DOMPData, bbnode::BbNode, r::Int64)::Tuple{Matrix{Bool}, Vector{Int64}, Int64}
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    cov = falses(nrows, ncols)
    for i in 1 : nrows, j in 1 : ncols
        cov[i, j] = data.D[i, j] <= r
    end
    iscov = falses(nrows)
    for b in bbnode.branches
        j, s, bound = b.j, b.sense, b.bound
        if s == 'G' && bound == 1
            iscov = iscov .| cov[:, j]
        end
    end
    jinds = setdiff(collect(1 : ncols), [b.j for b in bbnode.branches])
    iinds = [i for i in 1 : nrows if !iscov[i]]
    newcov = cov[iinds, jinds]
    return newcov, jinds, sum(iscov)
end