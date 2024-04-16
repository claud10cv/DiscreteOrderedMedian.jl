function dompk_pos(data::DOMPData, 
                    params::Parameters,
                    bbnode::BbNode, 
                    k::Int64, 
                    lb::Int64, 
                    ub::Int64, 
                    xub::Vector{Int64},
                    supp::Vector{Int64})::Tuple{Int64, Vector{Int64}}
    if (lb == ub && !isempty(xub)) return ub, deepcopy(xub) end
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    open = [b.j for b in bbnode.branches if b.sense == 'G' && b.bound == 1]
    closed = [b.j for b in bbnode.branches if b.sense == 'L' && b.bound == 0]
    nopen = length(open)
    nclosed = length(closed)
    x = zeros(Int64, ncols)
    y = zeros(Int64, ncols)
    for j in open
        x[j] = 1
    end
    newxub = deepcopy(x)
    for i in 1 : ncols
        if sum(newxub) >= data.p break
        elseif !in(i, closed) newxub[i] = 1
        end
    end 
    dub = compute_sorted_distances(data, newxub)  
    if dub[k] < ub
        # println("improved ub from $ub to $(dub[k])")
        ub = dub[k]
        xub = newxub
    end
    lbidx = find_index_uD(data, lb)
    ubidx = find_index_uD(data, ub)
    while lbidx < ubidx
        ridx = floor(Int64, (lbidx + ubidx) / 2)
        r = data.uD[ridx]
        cov, map, numcov = build_coverage(data, bbnode, r, k, true)
        sol = setcover(params, cov, k - numcov, data.p - nopen, supp[map])
        if !isempty(sol)
            # println("feasible r = $r")
            ub = r
            ubidx = ridx
            for j in 1 : ncols
                y[j] = 0
            end
            for j in map[sol]
                y[j] = 1
            end
            xub = x + y
            if sum(xub) < data.p
                diff = data.p - sum(xub)
                others = [i for i in 1 : ncols if !in(i, closed) && xub[i] == 0]
                for i in others[1 : diff]
                    xub[i] = 1
                end
            end
        else
            # println("infeasible r = $r")
            lbidx = ridx + 1
            lb = data.uD[lbidx]
        end
    end
    if sum(y) == 0
        cov, map, numcov = build_coverage(data, bbnode, ub, k, true)
        sol = setcover(params, cov, k - numcov, data.p - nopen, supp[map])
        if isempty(sol) return ub, zeros(Int64, ncols) end
        for j in 1 : ncols
            y[j] = 0
        end
        for j in map[sol]
            y[j] = 1
        end
        xub = x + y
        if sum(xub) < data.p
            diff = data.p - sum(xub)
            others = [i for i in 1 : ncols if !in(i, closed) && xub[i] == 0]
            for i in others[1 : diff]
                xub[i] = 1
            end
        end
    end
    @assert(sum(xub) == data.p, "wrong number of columns in solution! x = $xub")
    return ub, deepcopy(xub)
end

function dompk_neg(data::DOMPData, 
                    params::Parameters,
                    bbnode::BbNode, 
                    k::Int64, 
                    lb::Int64, 
                    ub::Int64, 
                    xlb::Vector{Int64},
                    supp::Vector{Int64})::Tuple{Int64, Vector{Int64}}
    if (lb == ub && !isempty(xlb)) return lb, deepcopy(xlb) end
    # println("starting dompk_neg with lb = $lb, ub = $ub")
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    open = [b.j for b in bbnode.branches if b.sense == 'G' && b.bound == 1]
    closed = [b.j for b in bbnode.branches if b.sense == 'L' && b.bound == 0]
    nopen = length(open)
    nclosed = length(closed)
    x = zeros(Int64, ncols)
    y = zeros(Int64, ncols)
    for j in open
        x[j] = 1
    end
    newxlb = deepcopy(x)
    for i in 1 : ncols
        if sum(newxlb) >= data.p break
        elseif !in(i, closed) newxlb[i] = 1
        end
    end 
    dlb = compute_sorted_distances(data, newxlb)  
    if dlb[k] > lb
        lb = dlb[k]
        xlb = newxlb
    end
    lbidx = find_index_uD(data, lb)
    ubidx = find_index_uD(data, ub)
    while lbidx < ubidx
        ridx = ceil(Int64, (lbidx + ubidx) / 2)
        r = data.uD[ridx]
        rprev = ridx == 1 ? 0 : data.uD[ridx - 1]
        # println("lb = $lb, ub = $ub, r = $r")
        cov, map, numcov = build_coverage(data, bbnode, rprev, k - 1, false)
        sol = setpacking(params, cov, k - 1 - numcov, data.p - nopen, supp[map])
        if !isempty(sol)
            lbidx = ridx
            lb = r
            for j in 1 : ncols
                y[j] = 0
            end
            for j in map[sol]
                y[j] = 1
            end
            xlb = x + y# println(", I")
            @assert(sum(xlb) == data.p)
        else
            ubidx = ridx - 1
            ub = data.uD[ridx - 1]
        end
    end
    if sum(xlb) != data.p
        cov, map, numcov = build_coverage(data, bbnode, lb - 1, k - 1, false)
        sol = setpacking(params, cov, k - 1 - numcov, data.p - nopen, supp[map])
        if isempty(sol) return lb, zeros(Int64, ncols) end
        for j in 1 : ncols
            y[j] = 0
        end
        for j in map[sol]
            y[j] = 1
        end
        xlb = x + y
        @assert(sum(xlb) == data.p)
    end
    # println("result = $(lb)")
    # println("x = $xlb at line 218")
    @assert(sum(xlb) == data.p)
    return lb, deepcopy(xlb)
end

function build_coverage(data::DOMPData, bbnode::BbNode, r::Int64, k::Int64, iscovering::Bool = true)::Tuple{Matrix{Bool}, Vector{Int64}, Int64}
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
            iscov = iscov .| @view cov[:, j]
        end
    end
    noncov = [i for i in 1 : nrows if !iscov[i]]
    numcov = sum(iscov)
    jinds = setdiff(collect(1 : ncols), [b.j for b in bbnode.branches])
    iinds = [i for i in 1 : nrows if !iscov[i]]
    if iscovering
        isdom = falses(ncols)
        for j in jinds, j2 in jinds
            if j == j2 continue end
            if isdom[j] continue end
            if isdom[j2] continue end
            diff = cov[iinds, j2] .& .!cov[iinds, j]
            if sum(diff) == 0
                isdom[j2] = true
            end
        end
        dom = [j for j in 1 : ncols if isdom[j]]
        jinds = setdiff(jinds, dom)
    else #it is a packing instance
        dom = [j for j in 1 : ncols if sum(@view cov[noncov, j]) > k - numcov]
        jinds = setdiff(jinds, dom)
    end
    iinds = [i for i in iinds if sum(cov[i, jinds]) > 0]
    cov = cov[iinds, jinds]
    return cov, jinds, numcov
end