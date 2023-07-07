function local_search!(data::DOMPData, ub::Int64, inset::Vector{Int64}, outset::Vector{Int64})::Int64
    mod = true
    ncols = size(data.D, 2)
    x = zeros(Int64, ncols)
    for j in inset
        x[j] = 1
    end
    while mod
        mod = false
        shuffle!(inset)
        shuffle!(outset)
        for jin_ind in 1 : length(inset)
            for jout_ind in 1 : length(outset)
                jin = inset[jin_ind]
                jout = outset[jout_ind]
                x[jin] = 0
                x[jout] = 1
                d = compute_sorted_distances(data, x)
                val = compute_weighted_cost(data, d)
                if val < ub
                    ub = val
                    inset[jin_ind] = jout
                    outset[jout_ind] = jin
                    mod = true
                else
                    x[jin] = 1
                    x[jout] = 0
                end
            end
        end
    end
    return ub
end

function iterated_local_search(data::DOMPData, xlb::Vector{Vector{Int64}})
    # newdata, xlb, map = reduce_data_for_ls(data, xlb)
    newdata = data
    map = collect(1 : size(data.D, 2))
    nrows = size(newdata.D, 1)
    ncols = size(newdata.D, 2)
    best = zeros(Int64, ncols)
    ub = typemax(Int64)
    # println("x = $xlb")
    let
        for k in 1 : data.p
            best[k] = 1
        end
        d = compute_sorted_distances(newdata, best)
        val = compute_weighted_cost(newdata, d)
        ub = val
    end
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
    nonsuccess_ils = 0
    first = true
    max_shakes = 3
    num_shakes = 1
    max_unsuccess = 10
    while nonsuccess_ils < max_unsuccess
        inset_curr = deepcopy(inset)
        outset_curr = deepcopy(outset)
        if first
            first = false
        else
            for t in 1 : num_shakes
                rin = rand(1 : length(inset_curr))
                rout = rand(1 : length(outset_curr))
                aux = inset_curr[rin]
                inset_curr[rin] = outset_curr[rout]
                outset_curr[rout] = aux
            end
        end
        x = zeros(Int64, ncols)
        for j in inset_curr
            x[j] = 1
        end
        thisd = compute_sorted_distances(newdata, x)
        thisub = compute_weighted_cost(newdata, thisd)
        val = local_search!(newdata, thisub, inset_curr, outset_curr)
        if val < ub
            ub = val
            inset = inset_curr
            outset = outset_curr
            nonsuccess_ils = 0
            num_shakes = 1
            # println("successful LS")
        else
            # println("unsuccessful LS")
            nonsuccess_ils += 1
        end
        if nonsuccess_ils >= max_unsuccess
            if num_shakes < max_shakes
                num_shakes += 1
                nonsuccess_ils = 0
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
    uD = unique(sort(vec(D)))
    data = DOMPData(D, data.p, data.lambda, uD)
    xk = [isempty(x) ? Int64[] : x[map] for x in xlb]
    return data, xk, map
end 
