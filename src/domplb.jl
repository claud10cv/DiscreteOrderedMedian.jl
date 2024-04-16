function domp_lb!(data::DOMPData, params::Parameters, bbnode::BbNode, parent::Union{BbNode, Nothing}, global_xub::Vector{Int64})::Nothing#Tuple{Int64, Int64, Vector{Float64}, Vector{Int64}}
    lb = 0
    lk = ReentrantLock()
    global_dists = compute_sorted_distances(data, global_xub)
    global_ub = compute_weighted_cost(data, global_dists)
    ub = typemax(Int64)
    # println("ub domplb l6 = $ub")
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    mind = data.uD[1]
    maxd = data.uD[end]
    xlb_all = zeros(Float64, ncols)
    xlb_ub = zeros(Float64, ncols)
    xub = zeros(Int64, ncols)
    dk = [Int64[] for i in 1 : nrows]
    xk = [Int64[] for i in 1 : nrows]
    currgap = isnothing(parent) ? typemax(Int64) : global_ub - parent.lb
    currdelta = 0
    fixone = [b for b in bbnode.branches if b.sense == 'G' && b.bound == 1]
    nfixone = length(fixone)
    if nfixone >= data.p
        # println("fixed node")
        xfix = zeros(Int64, ncols)
        for b in fixone
            xfix[b.j] = 1
        end
        dfix = compute_sorted_distances(data, xfix)
        vfix = compute_weighted_cost(data, dfix)
        # println("value of node is $vfix")
        bbnode.xlb = deepcopy(xfix), deepcopy(xfix)
        bbnode.xub = deepcopy(xfix)
        bbnode.xk = [deepcopy(xfix) for i in 1 : nrows]
        bbnode.lb = bbnode.ub = vfix
        bbnode.ropt = deepcopy(dfix)
        return
    end

    support = zeros(Int64, ncols)
    if !isempty(global_xub) support = deepcopy(global_xub) end
    ordering = [i for i in 1 : nrows if data.lambda[i] != 0]
    if params.symm_break
        sort!(ordering; lt = (u, v) -> score(data, u, global_dists) > score(data, v, global_dists))
    end
    if !isnothing(parent) support = deepcopy(parent.xub) end
    fathomed = false
    Threads.@threads for k in ordering
        if currdelta >= currgap || fathomed continue end
        # println("starting lb for k = $k")
        if params.warm_starts && !isnothing(parent) && can_recycle_solution(bbnode, parent.xk[k])
            # println("can recycle solution $k = $(parent.xk[k])")
            # val, bbnode.xk[k] = dompk_pos(data, bbnode, k, bbnode.ropt[k], k == nrows ? maxd : bbnode.ropt[k + 1])
            # d = compute_sorted_distances(data, parent.xk[k])
            # @assert(val == d[k], ("distances do not coincide, with val = $val, d = $(d[k]) and parent.xk = $(parent.xk[k])"))
            # lk = ReentrantLock()
            lock(lk) do
                val = bbnode.ropt[k] = parent.ropt[k]
                xk[k] = deepcopy(parent.xk[k])
                support += xk[k]
            end
            dkk = compute_sorted_distances(data, xk[k])
            lock(lk) do
                dk[k] = dkk
            end
            # bbnode.ropt[k] = parent.ropt[k]
        elseif data.lambda[k] > 0
            valub = maxd + 1
            valubidx = length(data.uD) + 1
            xubk = Int64[]
            # lk = ReentrantLock()
            lock(lk) do
                for l in ordering
                    if !isempty(dk[l])
                        if dk[l][k] < valub
                            valub = dk[l][k]
                            valubidx = find_index_uD(data, valub)
                            xubk = xk[l]
                        end
                    end
                end
            end
            if !isnothing(parent)
                lock(lk) do
                    delta = abs(data.lambda[k] * (valub - parent.ropt[k]))
                    while currdelta + delta >= currgap
                        # println("reducing upper bound of dompk")
                        valubidx -= 1
                        valub = data.uD[valubidx]
                        xubk = Int64[]
                        delta = abs(data.lambda[k] * (valub - parent.ropt[k]))
                    end
                end
            end
            lbk = isnothing(parent) ? mind : parent.ropt[k]
            @assert(lbk <= valub, ("bounds inconsistent for lambda > 0, $k, $(lbk), $(valub), $xk, $dk"))
            val, solk = dompk_pos(data, params, bbnode, k, lbk, valub, xubk, support)
            if sum(solk) == 0
                lock(lk) do
                    # println("fathoming node")
                    fathomed = true
                end
            else
                @assert(sum(solk) == data.p, ("wrong number of columns: $(sum(solk)) instead of $(data.p)"))
                dkk = compute_sorted_distances(data, solk)
                @assert(sum(solk) > 0, ("empty solution vector!"))
                lock(lk) do
                    xk[k] = solk
                    dk[k] = dkk
                    support += solk
                end
            end
        elseif data.lambda[k] < 0
            vallb = mind - 1
            vallbidx = 0
            xlbk = Int64[]
            lock(lk) do
                for l in ordering
                    if !isempty(dk[l])
                        if dk[l][k] > vallb
                            vallb = dk[l][k]
                            vallbidx = find_index_uD(data, vallb)
                            xlbk = xk[l]
                        end
                    end
                end
            end
            if !isnothing(parent)
                lock(lk) do
                    delta = abs(data.lambda[k] * (vallb - parent.ropt[k]))
                    while currdelta + delta >= currgap
                        # println("reducing upper bound of dompk")
                        vallbidx += 1
                        vallb = data.uD[vallbidx]
                        xlbk = Int64[]
                        delta = abs(data.lambda[k] * (vallb - parent.ropt[k]))
                    end
                end
            end
            ubk = isnothing(parent) ? maxd : parent.ropt[k]
            @assert(ubk >= vallb, ("bounds inconsistent for lambda < 0, $k, $(vallb), $(ubk), $xk, $dk"))
            val, solk = dompk_neg(data, params, bbnode, k, vallb, ubk, xlbk, support)
            if sum(solk) == 0
                lock(lk) do
                    # println("fathoming node")
                    fathomed = true
                end
            else
                @assert(sum(solk) == data.p, ("wrong number of columns: $(sum(solk)) instead of $(data.p)"))
                dkk = compute_sorted_distances(data, solk)
                lock(lk) do
                    xk[k] = solk
                    dk[k] = dkk
                    support += solk
                end
                # println("solution of value $val")   
                @assert(val <= bbnode.ropt[k], ("discrepancy in computing val for lambda < 0, child = $val, parent = $(bbnode.ropt[k])"))
            end
        else
            @assert(false, "error doing stuff for lambda = 0")
        end
        lock(lk) do
            if !fathomed
                @assert(sum(xk[k]) == data.p, "xk is of the wrong size ($(sum(xk[k])) instead of $(data.p)), xk = $(xk[k])")
                dist = compute_sorted_distances(data, xk[k])
                # println("dk[$k] = $(dist[1 : 12])")
                # println("xk[$k] = $(xk[k])")
                ubk = compute_weighted_cost(data, dist)
                # println("ubk = $ubk")
                bbnode.ropt[k] = val
                currdelta += isnothing(parent) ? 0 : abs(data.lambda[k] * (val - parent.ropt[k]))
                lb += data.lambda[k] * val
                if ubk < ub
                    ub = ubk
                    xub = xk[k]
                end
            end
        end
    end
    # println("delta = $currdelta, gap = $currgap")
    if currdelta >= currgap || lb >= global_ub || fathomed
        bbnode.lb = global_ub + 1
        bbnode.ub = global_ub + 1
        # println("fathoming node")
        return
    end
    # println()
    bbnode.xk = reassign_vectors(data, xk, dk, bbnode.ropt)
    
    # println("ub domplb l147 = $ub")
    
    @assert(isnothing(parent) || lb >= parent.lb, ("lb is lower that that of parent, child.lb = $lb, parent.lb = $(parent.lb)"))
    dist = compute_sorted_distances(data, xub)
    # println("dist = $dist")
    # println("xub = $xub")
    ub = compute_weighted_cost(data, dist)

    # println("ub domplb = $ub")
    count_all = 0
    count_diff = 0
    # nit = 0
    delta = [(in(k, ordering) ? abs(data.lambda[k] * (dist[k] - bbnode.ropt[k])) + 1 : 0) for k in 1 : nrows]
    sdelta = sum(delta)
    for k in ordering
        # nit += 1
        pdelta = delta[k] / sdelta * 100
        sk = score(data, k, global_dists)
        @assert(!isempty(bbnode.xk[k]), "empty solution vector!")
        @assert(sum(bbnode.xk[k]) == data.p, "solution vector does not contain p points: $(bbnode.xk[k])")
        if !isempty(bbnode.xk[k]) && sum(bbnode.xk[k]) > 0
            xlb_all += sk * bbnode.xk[k]
            count_all += sk
            if delta[k] > 1 + 1e-5
                xlb_ub += sk * bbnode.xk[k]
                count_diff += sk
                # fractionality = maximum(xlb_ub / count_diff - round.(Int64, xlb_ub / count_diff))
                # println("fractionality = $fractionality")
                # if fractionality > 1e-6 break end
            end
        end
    end
    # println("aborted LB at iteration $nit of $(length(ordering))")
    xlb_all /= count_all
    xlb_ub /= count_diff
    for j in 1 : ncols
        if abs(xlb_all[j] - round(xlb_all[j])) < 1e-7
            xlb_all[j] = round(xlb_all[j])
        end
        if abs(xlb_ub[j] - round(xlb_ub[j])) < 1e-7
            xlb_ub[j] = round(xlb_ub[j])
        end
    end
    bbnode.xlb = xlb_all, xlb_ub
    bbnode.lb = lb
    bbnode.xub = xub
    bbnode.ub = ub
    # for r in 1 : nrows
    #     println("x[$r] = $(xk[r]), d[$r] = $(dk[r][r])")
    # end
    # println("lower bound of $(bbnode.lb)")
    return
    # return lb, ub, xlb, xub
end