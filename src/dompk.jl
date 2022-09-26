function dompk_pos(data::DOMPData, 
                    bbnode::BbNode, 
                    k::Int64, 
                    lb::Int64, 
                    ub::Int64, 
                    xub::Vector{Int64},
                    supp::Vector{Int64})::Tuple{Int64, Vector{Int64}}
    if (lb == ub) return ub, deepcopy(xub) end
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
    newub = lb
    step = 1
    while newub < ub
        # println("lb, ub = $lb, $newub")
        r = newub
        cov, map, numcov = build_coverage(data, bbnode, r)
        sol = setcover(cov, k - numcov, data.p - nopen, supp[map])
        if !isempty(sol)
            ub = newub
            for j in 1 : ncols
                y[j] = 0
            end
            for j in map[sol]
                y[j] = 1
            end
            xub = x + y
            break
        else
            # println("infeasible r = $lb")
            lb = r + 1
            newub = min(ub, newub + step)
            step *= 2
        end
    end
    # println("starting BS from lb = $lb, ub = $ub")
    while lb < ub
        r = floor(Int64, (lb + ub) / 2)
        cov, map, numcov = build_coverage(data, bbnode, r)
        sol = setcover(cov, k - numcov, data.p - nopen, supp[map])
        if !isempty(sol)
            # println("feasible r = $r")
            ub = r
            for j in 1 : ncols
                y[j] = 0
            end
            for j in map[sol]
                y[j] = 1
            end
            xub = x + y
        else
            # println("infeasible r = $r")
            lb = r + 1
        end
    end
    return ub, deepcopy(xub)
end

function dompk_neg(data::DOMPData, 
                    bbnode::BbNode, 
                    k::Int64, 
                    lb::Int64, 
                    ub::Int64, 
                    xlb::Vector{Int64},
                    supp::Vector{Int64})::Tuple{Int64, Vector{Int64}}
    if (lb == ub) return lb, deepcopy(xlb) end
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
    newlb = ub
    step = 1
    while newlb > lb
        r = newlb
        cov, map, numcov = build_coverage(data, bbnode, r - 1, false)
        sol = setpacking(cov, k - 1 - numcov, data.p - nopen, supp[map])
        if !isempty(sol)
            # println("x = $xlb at line 94")
            lb = r
            for j in 1 : ncols
                y[j] = 0
            end
            for j in map[sol]
                y[j] = 1
            end
            xlb = x + y
            break
        else
            ub = r - 1
            newlb = max(lb, newlb - step)
            step *= 2
        end
    end
    # println("starting BS with lb = $lb, ub = $ub")
    while lb < ub
        r = ceil(Int64, (lb + ub) / 2)
        # println("lb = $lb, ub = $ub, r = $r")
        cov, map, numcov = build_coverage(data, bbnode, r - 1, false)
        sol = setpacking(cov, k - 1 - numcov, data.p - nopen, supp[map])
        if !isempty(sol)
            # println("feasible r = $r")
            # println("sol = $sol")
            lb = r
            for j in 1 : ncols
                y[j] = 0
            end
            for j in map[sol]
                y[j] = 1
            end
            xlb = x + y# println(", I")
        else
            # println("infeasible r = $r")
            ub = r - 1
            # println(", F")
        end
    end
    # println("result = $(lb)")
    # println("x = $xlb at line 133")
    return lb, deepcopy(xlb)
end

function build_coverage(data::DOMPData, bbnode::BbNode, r::Int64, applydom::Bool = true)::Tuple{Matrix{Bool}, Vector{Int64}, Int64}
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
    jinds = setdiff(collect(1 : ncols), [b.j for b in bbnode.branches])
    iinds = [i for i in 1 : nrows if !iscov[i]]
    if applydom
        dom = falses(ncols)
        for j in jinds, j2 in jinds
            if j == j2 continue end
            if dom[j] continue end
            if dom[j2] continue end
            diff = cov[iinds, j2] .& .!cov[iinds, j]
            if sum(diff) == 0
                # println("found dominance")
                dom[j2] = true
            end
        end
        dom = [j for j in 1 : ncols if dom[j]]
        jinds = setdiff(jinds, dom)
    end
    iinds = [i for i in iinds if sum(cov[i, jinds]) > 0]
    cov = cov[iinds, jinds]
    return cov, jinds, sum(iscov)
end