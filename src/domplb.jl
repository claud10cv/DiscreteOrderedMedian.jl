function domp_lb!(data::DOMPData, bbnode::BbNode, parent::Union{BbNode, Nothing}, global_xub::Vector{Int64})::Nothing#Tuple{Int64, Int64, Vector{Float64}, Vector{Int64}}
    lb = 0
    lk = ReentrantLock()
    global_dists = compute_sorted_distances(data, global_xub)
    ub = typemax(Int64)#compute_weighted_cost(data, global_dists)
    # println("ub domplb l6 = $ub")
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    mind = minimum(data.D)
    maxd = maximum(data.D)
    xlb_all = zeros(Float64, ncols)
    xlb_ub = zeros(Float64, ncols)
    xub = zeros(Int64, ncols)
    dk = [Int64[] for i in 1 : nrows]
    xk = [Int64[] for i in 1 : nrows]
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
    sort!(ordering; lt = (u, v) -> score(data, u) > score(data, v))
    if !isnothing(parent) support = deepcopy(parent.xub) end
    Threads.@threads for k in ordering
        # println("starting lb for k = $k")
        if !isnothing(parent) && can_recycle_solution(bbnode, parent.xk[k])
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
            xubk = Int64[]
            # lk = ReentrantLock()
            lock(lk) do
                for l in ordering
                    if !isempty(dk[l])
                        if dk[l][k] < valub
                            valub = dk[l][k]
                            xubk = xk[l]
                        end
                    end
                end
            end
            lbk = isnothing(parent) ? mind : parent.ropt[k]
            @assert(lbk <= valub, ("bounds inconsistent for lambda > 0, $k, $(lbk), $(valub), $xk, $dk"))
            val, solk = dompk_pos(data, bbnode, k, lbk, valub, xubk, support)
            # if val == valub && !isempty(xubk)
            #     solk = xubk
            # end
            dkk = compute_sorted_distances(data, solk)
            @assert(sum(solk) > 0, ("empty solution vector!"))
            lock(lk) do
                xk[k] = solk
                dk[k] = dkk
                support += solk
            end
        elseif data.lambda[k] < 0
            vallb = mind - 1
            xlbk = Int64[]
            # lk = ReentrantLock()
            lock(lk) do
                for l in ordering
                    if !isempty(dk[l])
                        if dk[l][k] > vallb
                            vallb = dk[l][k]
                            xlbk = xk[l]
                        end
                    end
                end
            end
            ubk = isnothing(parent) ? maxd : parent.ropt[k]
            @assert(ubk >= vallb, ("bounds inconsistent for lambda < 0, $k, $(ubk), $(vallb), $xk, $dk"))
            # println("vallb = $vallb, bbnode.ropt[k] = $(bbnode.ropt[k])")
            # print("running dompk for k = $k, lambda < 0, lb = $vallb, ub = $ubk")
            val, solk = dompk_neg(data, bbnode, k, vallb, ubk, xlbk, support)
            # println(", val = $val")
            # if val == vallb && !isempty(xlbk)
            #     solk = xlbk
            # end
            @assert(sum(solk) > 0, ("empty solution vector!"))
            dkk = compute_sorted_distances(data, solk)
            lock(lk) do
                xk[k] = solk
                dk[k] = dkk
                support += solk
            end
            # println("solution of value $val")   
            @assert(val <= bbnode.ropt[k], ("discrepancy in computing val for lambda < 0, child = $val, parent = $(bbnode.ropt[k])"))
        else
            @assert(false, "error doing stuff for lambda = 0")
        end
        let
            # lk = ReentrantLock()
            @assert(sum(xk[k]) == data.p, "xk is of the wrong size ($(sum(xk[k])) instead of $(data.p)), xk = $(xk[k])")
            dist = compute_sorted_distances(data, xk[k])
            # println("dk[$k] = $(dist[1 : 12])")
            # println("xk[$k] = $(xk[k])")
            ubk = compute_weighted_cost(data, dist)
            # println("ubk = $ubk")
            lock(lk) do
                bbnode.ropt[k] = val
                lb += data.lambda[k] * val
                if ubk < ub
                    ub = ubk
                    xub = xk[k]
                end
            end
        end
    end
    
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
    for k in ordering
        # nit += 1
        sk = 1#log(2, score(data, k) + 1)
        if !isempty(bbnode.xk[k]) && sum(bbnode.xk[k]) > 0
            delta = abs(dist[k] - bbnode.ropt[k])
            xlb_all += sk * bbnode.xk[k]
            count_all += sk
            if delta > 0
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
    return
    # return lb, ub, xlb, xub
end