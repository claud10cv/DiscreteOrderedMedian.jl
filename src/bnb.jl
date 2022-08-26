using DataStructures
function bnb(data::DOMPData)::Tuple{Int64, Int64, Tuple{Vector{Float64}, Vector{Float64}}, Vector{Int64}}
    t0 = time_ns()
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    root = BbNode([], 0, 0, (zeros(ncols), zeros(ncols)), zeros(Int64, ncols), zeros(Int64, nrows))
    global_xub = zeros(Int64, ncols)
    for k in 1 : data.p
        global_xub[k] = 1
    end
    global_dists = compute_sorted_distances(data, global_xub)
    global_ub = compute_weighted_cost(data, global_dists)
    domp_lb!(data, root, global_xub)
    if root.lb >= root.ub return root.lb, root.ub, root.xlb, root.xub end
    queue = PriorityQueue{BbNode, Int64}()
    push!(queue, root => root.lb)
    global_ub = root.ub
    global_lb = root.lb
    global_xlb = root.xlb
    global_xub = root.xub
    it = 0
    println("ITERATION\tBEST BOUND\tBKS\t\tGAP(%)\t\tT (s)\t\tNODES LEFT")
    while !isempty(queue)
        bbnode = dequeue!(queue)
        if bbnode.lb > global_ub continue end
        it += 1
        global_lb = bbnode.lb
        global_xlb = bbnode.xlb
        gap = ceil(((global_ub - global_lb) / global_ub * 100) * 100) / 100
        t1 = time_ns()
        elapsed = ceil(100 * (t1 - t0) * 1e-9) / 100
        println("\t$it\t$global_lb\t\t$global_ub\t\t$gap\t\t$elapsed\t\t$(length(queue))")
        if bbnode.lb >= global_ub continue end
        # println("processing node $bbnode")
        # create two branches using the most fractional criterion
        # val, j = find_most_fractional(bbnode.xlb...)
        val, j = strong_branching(data, bbnode, bbnode.xlb...)
        # println("branching on variable x[$j] with fraction $val")
        # create branch xj <= 0
        left_branch = deepcopy(bbnode.branches)
        push!(left_branch, BranchInfo(j, 'L', 0))
        left_child = BbNode(left_branch, 0, 0, ([], []), [], deepcopy(bbnode.ropt))
        #create branch xj >= 1
        right_branch = deepcopy(bbnode.branches)
        push!(right_branch, BranchInfo(j, 'G', 1))
        right_child = BbNode(right_branch, 0, 0, ([], []), [], deepcopy(bbnode.ropt))
        children = [left_child, right_child]
        # solve children
        for child in children
            domp_lb!(data, child, global_xub)
            @assert(child.lb >= bbnode.lb, ("bound discrepancy of $(child.lb - bbnode.lb)"))
            # println("created child with bound $(child.lb)")
            if child.ub < global_ub
                global_ub = child.ub
                global_xub = child.xub
            end
            if child.lb < global_ub
                enqueue!(queue, child => child.lb)
            end
            # println("created and solved child $child")
        end
    end
    if isempty(queue)
        global_lb = global_ub
    else
        bbnode = peek(queue)
        global_lb = bbnode.lb
    end 
    it += 1
    t1 = time_ns()
    elapsed = ceil(100 * (t1 - t0) * 1e-9) / 100
    println("\t$it\t$global_ub\t\t$global_ub\t\t0.0\t\t$elapsed\t\t0")
    return global_lb, global_ub, global_xlb, global_xub
end
        

