using DataStructures, Random
function local_search(data::DOMPData, xlb::Vector{Vector{Int64}})
    newdata, xlb, map = reduce_data_for_ls(data, xlb)
    nrows = size(newdata.D, 1)
    ncols = size(newdata.D, 2)
    best = zeros(Int64, 0)
    ub = typemax(Int64)
    for x in xlb
        if !isempty(x)
            d = compute_sorted_distances(newdata, x)
            val = compute_weighted_cost(newdata, d)
            if val < ub
                ub = val
                best = x
            end
        end
    end
    inset = [j for j in 1 : ncols if best[j] == 1]
    outset = [j for j in 1 : ncols if best[j] == 0]
    mod = true
    while mod
        mod = false
        shuffle!(inset)
        shuffle!(outset)
        for jin_ind in 1 : length(inset)
            for jout_ind in 1 : length(outset)
                jin = inset[jin_ind]
                jout = outset[jout_ind]
                best[jin] = 0
                best[jout] = 1
                d = compute_sorted_distances(newdata, best)
                val = compute_weighted_cost(newdata, d)
                if val < ub
                    ub = val
                    inset[jin_ind] = jout
                    outset[jout_ind] = jin
                    mod = true
                else
                    best[jin] = 1
                    best[jout] = 0
                end
            end
        end
    end
    ncols = size(data.D, 2)
    xsol = zeros(Int64, ncols)
    for j in map[inset]
        xsol[j] = 1
    end
    return xsol, ub
end
        
function reduce_data_for_ls(data::DOMPData, xlb::Vector{Vector{Int64}})::Tuple{DOMPData, Vector{Vector{Int64}}, Vector{Int64}}
    ncols = size(data.D, 2)
    map = [j for j in 1 : ncols if sum([x[j] for x in xlb if !isempty(x)]) > 0]
    D = data.D[:, map]
    data = DOMPData(D, data.p, data.lambda)
    xk = [isempty(x) ? Int64[] : x[map] for x in xlb]
    return data, xk, map
end 
