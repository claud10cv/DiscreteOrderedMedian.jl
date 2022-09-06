function domp_lb!(data::DOMPData, bbnode::BbNode, parent::Union{BbNode, Nothing}, global_xub::Vector{Int64})::Nothing#Tuple{Int64, Int64, Vector{Float64}, Vector{Int64}}
    lb = 0
    global_dists = compute_sorted_distances(data, global_xub)
    ub = compute_weighted_cost(data, global_dists)
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    mind = minimum(data.D)
    maxd = maximum(data.D)
    nneg_all = 0
    nneg_ub = 0
    xlb_all = zeros(Int64, ncols)
    xlb_ub = zeros(Int64, ncols)
    xub = zeros(Int64, ncols)
    dk = [zeros(Int64, nrows) for i in 1 : nrows]
    # resize!(bbnode.xk, nrows)
    for k in nrows : -1 : 1
        if data.lambda[k] != 0 && !isnothing(parent) && can_recycle_solution(bbnode, parent.xk[k])
            # println("can recycle solution $k = $(parent.xk[k])")
            # val, bbnode.xk[k] = dompk_pos(data, bbnode, k, bbnode.ropt[k], k == nrows ? maxd : bbnode.ropt[k + 1])
            # d = compute_sorted_distances(data, parent.xk[k])
            # @assert(val == d[k], ("distances do not coincide, with val = $val, d = $(d[k]) and parent.xk = $(parent.xk[k])"))
            val = bbnode.ropt[k] = parent.ropt[k]
            bbnode.xk[k] = deepcopy(parent.xk[k])
            dk[k] = compute_sorted_distances(data, bbnode.xk[k])
            # bbnode.ropt[k] = parent.ropt[k]
        elseif data.lambda[k] > 0
            valub = maxd + 1
            xubk = Int64[]
            for l in k + 1 : nrows
                if data.lambda[l] > 0
                    if dk[l][k] < valub
                        valub = dk[l][k]
                        xubk = bbnode.xk[l]
                    end
                end
            end
            @assert(bbnode.ropt[k] <= valub, ("bounds inconsistent, $(bbnode.ropt[k]), $(valub)"))
            val, solk = dompk_pos(data, bbnode, k, bbnode.ropt[k], valub, xubk)
            if val == valub && !isempty(xubk)
                solk = xubk
            end
            @assert(!isempty(solk), ("empty solution vector!"))
            bbnode.xk[k] = solk
            dk[k] = compute_sorted_distances(data, bbnode.xk[k])
        elseif data.lambda[k] < 0
            vallb = mind - 1
            xlbk = Int64[]
            for l in k + 1 : nrows
                if data.lambda[l] < 0
                    if dk[l][k] > vallb
                        vallb = dk[l][k]
                        xlbk = bbnode.xk[l]
                    end
                end
            end
            val, solk = dompk_neg(data, bbnode, k, vallb, bbnode.ropt[k], xlb[k])
            if val == vallb && !isempty(xlbk)
                solk = xlbk
            end
            @assert(!isempty(solk), ("empty solution vector!"))
            bbnode.xk[k] = solk
            dk[k] = compute_sorted_distances(data, bbnode.xk[k])    
            @assert(val <= bbnode.ropt[k], ("discrepancy in computing val for lambda < 0, child = $val, parent = $(bbnode.ropt[k])"))
        else
            if k == nrows
                bbnode.ropt[k] = maxd
            else
                bbnode.ropt[k] = bbnode.ropt[k + 1]
            end
        end
        if data.lambda[k] != 0
            bbnode.ropt[k] = val
            lb += data.lambda[k] * val
            dist = compute_sorted_distances(data, bbnode.xk[k])
            ubk = compute_weighted_cost(data, dist)
            if ubk < ub
                ub = ubk
                xub = bbnode.xk[k]
            end
            nneg_all += 1
            xlb_all += bbnode.xk[k]
        end
    end
    dist = compute_sorted_distances(data, xub)
    ub = compute_weighted_cost(data, dist)

    for k in 1 : nrows
        if !isempty(bbnode.xk[k]) && sum(bbnode.xk[k]) > 0 && bbnode.ropt[k] < dist[k]
            nneg_ub += 1
            xlb_ub += bbnode.xk[k]
        end
    end
    xlb_all /= nneg_all
    xlb_ub /= nneg_ub
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