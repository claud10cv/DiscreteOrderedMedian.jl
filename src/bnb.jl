using DataStructures
function bnb(data::DOMPData; time_limit = 7200)::Tuple{Int64, Int64, Tuple{Vector{Float64}, Vector{Float64}}, Vector{Int64}}
    println("ITERATION\tBEST BOUND\tBKS\t\tGAP(%)\t\tT (s)\t\tNODES LEFT")
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
    domp_lb!(data, root, nothing, global_xub)
    global_ub = root.ub
    global_lb = root.lb
    global_xlb = root.xlb
    global_xub = root.xub
    xls, lsub = local_search(data, root.xk)
    global_ub = lsub
    global_xub = xls
    it = 0
    root_gap = ceil(100 * (global_ub - global_lb) / global_ub * 100) / 100
    t1 = time_ns()
    elapsed = ceil(100 * (t1 - t0) * 1e-9) / 100
    queue = PriorityQueue{BbNode, Int64}()
    push!(queue, root => root.lb)
    println("\t$it\t$global_lb\t\t$global_ub\t\t$root_gap\t\t$elapsed\t\t$(length(queue))")
    if root.lb >= root.ub return root.lb, root.ub, root.xlb, root.xub end
    while !isempty(queue) && elapsed < time_limit
        bbnode = dequeue!(queue)
        global_lb = bbnode.lb
        global_xlb = bbnode.xlb
        gap = ceil(((global_ub - global_lb) / global_ub * 100) * 100) / 100
        if bbnode.lb >= global_ub continue end
        it += 1
        # find the most fractional variable
        # val, j = find_most_fractional(bbnode.xlb...)
        # find the best branch by strong branching
        pcosts = strong_branching(data, bbnode, bbnode.xlb...)
        # println(pcosts[1])
        j = pcosts[1].i
        # println("branching on variable x[$j] with fraction $val")
        # create branch xj <= 0
        left_branch = deepcopy(bbnode.branches)
        push!(left_branch, BranchInfo(j, 'L', 0))
        left_child = BbNode(left_branch, 0, 0, ([], []), [], [Int64[] for k in 1 : nrows], deepcopy(bbnode.ropt))
        #create branch xj >= 1
        right_branch = deepcopy(bbnode.branches)
        push!(right_branch, BranchInfo(j, 'G', 1))
        right_child = BbNode(right_branch, 0, 0, ([], []), [], [Int64[] for k in 1 : nrows], deepcopy(bbnode.ropt))
        children = [left_child, right_child]
        # solve children
        # println("processing children")
        for child in children
            domp_lb!(data, child, bbnode, global_xub)
            # println("child bound = $(child.lb)")
            @assert(child.lb >= bbnode.lb, ("bound discrepancy of $(child.lb - bbnode.lb)"))
            if child.ub < global_ub
                global_ub = child.ub
                global_xub = child.xub
            end
            if child.lb < global_ub
                enqueue!(queue, child => child.lb)
            end
        end
        t1 = time_ns()
        elapsed = ceil(100 * (t1 - t0) * 1e-9) / 100
        println("\t$it\t$global_lb\t\t$global_ub\t\t$gap\t\t$elapsed\t\t$(length(queue))")
    end
    if isempty(queue)
        global_lb = global_ub
    else
        bbnode, bestbound = peek(queue)
        global_lb = bestbound
    end 
    it += 1
    t1 = time_ns()
    elapsed = ceil(100 * (t1 - t0) * 1e-9) / 100
    gap = ceil(Int64, 100 * (global_ub - global_lb) / global_ub * 100) / 100
    println("\t$it\t$global_lb\t\t$global_ub\t\t$(gap)\t\t$elapsed\t\t$(length(queue))")
    return global_lb, global_ub, global_xlb, global_xub
end
        

