function domp_lb!(data::DOMPData, bbnode::BbNode, global_xub::Vector{Int64})::Nothing#Tuple{Int64, Int64, Vector{Float64}, Vector{Int64}}
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
    xk = [Int64[] for k in 1 : nrows]
    for k in nrows : -1 : 1
        if data.lambda[k] > 0  
            val, xk[k] = dompk_pos(data, bbnode, k, bbnode.ropt[k], k == nrows ? maxd : bbnode.ropt[k + 1])
        elseif data.lambda[k] < 0
            val, xk[k] = dompk_neg(data, bbnode, k, mind, bbnode.ropt[k])
        end
        if data.lambda[k] != 0
            bbnode.ropt[k] = val
            lb += data.lambda[k] * val
            dist = compute_sorted_distances(data, xk[k])
            ubk = compute_weighted_cost(data, dist)
            if ubk < ub
                ub = ubk
                xub = xk[k]
            end
            nneg_all += 1
            xlb_all += xk[k]
        end
    end
    dist = compute_sorted_distances(data, xub)
    ub = compute_weighted_cost(data, dist)

    for k in 1 : nrows
        if !isempty(xk[k]) && bbnode.ropt[k] < dist[k]
            nneg_ub += 1
            xlb_ub += xk[k]
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