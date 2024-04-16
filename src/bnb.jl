function bnb(data::DOMPData, params::Parameters = default_parameters())::Result
    println("ITERATION,BEST_BOUND,BKS,GAP(%),T(s),FRACT,BR_VAR,BR_VAL,DEPTH,NODES_LEFT")
    t0 = time_ns()
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    maxd = maximum(data.D)
    mind = minimum(data.D)
    root = BbNode([], 0, 0, (zeros(ncols), zeros(ncols)), zeros(Int64, ncols), [Int64[] for k in 1 : nrows], zeros(Int64, nrows))
    for i in 1 : nrows
        root.ropt[i] = data.lambda[i] > 0 ? mind : (data.lambda[i] < 0 ? maxd : 0)
    end
    global_xub = zeros(Int64, ncols)
    for k in 1 : data.p
        global_xub[k] = 1
    end
    global_dists = compute_sorted_distances(data, global_xub)
    global_ub = compute_weighted_cost(data, global_dists)
    domp_lb!(data, params, root, nothing, global_xub)
    global_ub = root.ub
    global_lb = root.lb
    global_xlb = root.xlb
    global_xub = root.xub
    if params.primal_heur
        xls, lsub = iterated_local_search(data, root.xk)
        global_ub = lsub
        global_xub = xls
    end
    it = 0
    root_gap = ceil(100 * (global_ub - global_lb) / abs(global_ub + 1e-10) * 100) / 100
    t1 = time_ns()
    elapsed = ceil(100 * (t1 - t0) * 1e-9) / 100
    queue = PriorityQueue{BbNode, Int64}()
    push!(queue, root => root.lb)
    root_fract = length([i for i in 1 : ncols if abs(root.xlb[1][i] - round(root.xlb[1][i])) > 1e-6])
    println("$it,$global_lb,$global_ub,$root_gap,$elapsed,$root_fract,--,--,0,$(length(queue))")
    if root.lb >= root.ub return Result(root.lb, root.ub, root.xlb, root.xub) end
    pseudocosts = zeros(ncols, 2)
    next_restart = 20
    nrestarts = 0
    while !isempty(queue) && elapsed < params.time_limit
        bbnode = dequeue!(queue)
        global_lb = max(global_lb, bbnode.lb)
        global_xlb = bbnode.xlb
        gap = ceil(((global_ub - global_lb) / abs(global_ub + 1e-10) * 100) * 100) / 100
        if bbnode.lb >= global_ub continue end
        if it == next_restart && nrestarts < 2
            nrestarts += 1
            println("performing restart...")
            empty!(queue)
            enqueue!(queue, root => root.lb)
            next_restart *= 5
            continue
        end
        it += 1
        # find the most fractional variable
        # val, j = find_most_fractional(bbnode.xlb...)
        # find the best branch by strong branching
        branches = strong_branching(data, params, bbnode, global_ub, pseudocosts, bbnode.xlb...)
        j = branches[1][1].i
        left_child = branches[1][2]
        right_child = branches[1][3]
        # # println("branching on variable x[$j] with fraction $val")
        # # create branch xj <= 0
        if isnothing(left_child)
            left_branch = deepcopy(bbnode.branches)
            push!(left_branch, BranchInfo(j, 'L', 0))
            left_child = BbNode(left_branch, 0, 0, ([], []), [], [Int64[] for k in 1 : nrows], deepcopy(bbnode.ropt))
            domp_lb!(data, params, left_child, bbnode, global_xub)
        end
        if isnothing(right_child)
            right_branch = deepcopy(bbnode.branches)
            push!(right_branch, BranchInfo(j, 'G', 1))
            right_child = BbNode(right_branch, 0, 0, ([], []), [], [Int64[] for k in 1 : nrows], deepcopy(bbnode.ropt))
            domp_lb!(data, params, right_child, bbnode, global_xub)
        end
        
        children = [left_child, right_child]
        chdelta = Float64[]
        for child in children
            # domp_lb!(data, child, bbnode, global_xub)
            push!(chdelta, child.lb - bbnode.lb)
            # println("child bound = $(child.lb)")
            @assert(child.lb >= bbnode.lb, ("bound discrepancy of $(child.lb - bbnode.lb), child.lb = $(child.lb), parent.lb = $(bbnode.lb)"))
            if child.ub < global_ub
                global_ub = child.ub
                global_xub = child.xub
            end
            if child.lb < global_ub
                enqueue!(queue, child => child.lb)
            end
        end
        pseudocosts[j, 1] = (pseudocosts[j, 1] + chdelta[1] * bbnode.xlb[1][j]) / 2
        pseudocosts[j, 2] = (pseudocosts[j, 2] + chdelta[2] * (1 - bbnode.xlb[1][j])) / 2
        t1 = time_ns()
        elapsed = ceil(100 * (t1 - t0) * 1e-9) / 100
        bbnode_fract = length([i for i in 1 : ncols if abs(bbnode.xlb[1][i] - round(bbnode.xlb[1][i])) > 1e-6])
        println("$it,$global_lb,$global_ub,$gap,$elapsed,$bbnode_fract,x[$j],$(bbnode.xlb[1][j]),$(length(bbnode.branches)),$(length(queue))")
    end
    if isempty(queue)
        global_lb = global_ub
        depth = "--"
        fract = 0
    else
        bbnode, bestbound = peek(queue)
        global_lb = max(global_lb, bestbound)
        depth = length(bbnode.branches)
        fract = length([i for i in 1 : ncols if abs(bbnode.xlb[1][i] - round(bbnode.xlb[1][i])) > 1e-6])
    end 
    it += 1
    t1 = time_ns()
    elapsed = ceil(100 * (t1 - t0) * 1e-9) / 100
    gap = ceil(Int64, 100 * (global_ub - global_lb) / abs(global_ub + 1e-10) * 100) / 100
    println("$it,$global_lb,$global_ub,$(gap),$elapsed,$fract,--,--,$depth,$(length(queue))")
    return Result(global_lb, global_ub, global_xlb, global_xub)
end
        

