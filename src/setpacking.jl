function setpacking(params::Parameters, coverage::Matrix{Bool}, k::Int64, p::Int64, supp::Vector{Int64})::Vector{Int64}
    nrows = size(coverage, 1)
    ncols = size(coverage, 2)
    if k < 0 return [] end
    if p > ncols return [] end
    if p <= 0 return [1] end
    result = setpacking_localsearch(coverage, k, p, supp)
    if isempty(result)
        result = setpacking_exact(params, coverage, k, p, supp)
    end
    return result
end

function setpacking_localsearch(coverage::Matrix{Bool}, k::Int64, p::Int64, supp::Vector{Int64})::Vector{Int64}
    nrows = size(coverage, 1)
    ncols = size(coverage, 2)
    cols = collect(1 : ncols)
    result = Int64[]
    covered = falses(nrows)
    newcov = deepcopy(coverage)
    scov = zeros(Int64, ncols)
    for i in 1 : p
        # newcov = Matrix{Bool}(undef, nrows, ncols)
        for j in cols
            newcov[:, j] = coverage[:, j] .& .!covered
            scov[j] = sum(@view newcov[:, j])
        end
        sort!(cols; lt = (u, v) -> scov[u] < scov[v] || (scov[u] == scov[v] && supp[u] > supp[v]))
        jstar = popfirst!(cols)
        push!(result, jstar)
        covered = covered .| coverage[:, jstar]
    end
    if sum(covered) <= k 
        # println("feasible heuristic p = $p, k = $k, $result, coverage of $(sum(covered))")
        return sort(result)
    else 
        # println("infeasible heuristic")
        return Int64[]
    end
end

function setpacking_exact(params::Parameters, coverage::Matrix{Bool}, k::Int64, p::Int64, supp::Vector{Int64})::Vector{Int64}
    nrows = size(coverage, 1)
    ncols = size(coverage, 2)
    m = JuMP.direct_model(params.optimizer_data.optimizer())
    params.optimizer_data.set_attributes(m, params.time_limit + 1, 1)
    # println("solving setpacking with k = $k, p = $p, nrows = $nrows, ncols = $ncols")
    # println("cov = $coverage")
    # println("constructing model")
    # m = JuMP.direct_model(CPLEX.Optimizer())
    # set_optimizer_attributes(m, 
    #     # "CPXPARAM_MIP_Tolerances_LowerCutoff" => p - 1e-5, 
    #     # "CPXPARAM_MIP_Limits_UpperObjStop" => p - 1e-5,
    #     "CPXPARAM_ScreenOutput" => 0,
    #     "CPXPARAM_Threads" => 1,
    #     "CPXPARAM_MIP_Tolerances_MIPGap" => 0.1,
    #     "CPXPARAM_MIP_Limits_Solutions" => 1)
    conflict = Vector{Vector{Int64}}()
    let    
        g = SimpleGraph(ncols)
        for i in 1 : ncols, j in i + 1 : ncols
            s = sum(coverage[:, i] .| coverage[:, j])
            if s > k 
                add_edge!(g, (i, j))
            end
        end
        edg = collect(edges(g))
        for e in edg
            push!(conflict, [src(e), dst(e)])# = Graphs.edges(g)
        end
    end
    subset = Tuple{Int64, Int64}[]
    equiv = Tuple{Int64, Int64}[]
    for i in 1 : ncols, j in i + 1 : ncols
        s = .!coverage[:, i] .| coverage[:, j]
        t = .!coverage[:, j] .| coverage[:, i]
        alls = all(s)
        allt = all(t)
        if alls && allt
            push!(equiv, (i, j))
        elseif alls
            push!(subset, (i, j)) 
        elseif allt
            push!(subset, (j, i))
        end
    end
    
    # m = JuMP.Model(optimizer_with_attributes(Gurobi.Optimizer, 
    #     "Cutoff" => p + 1e-5, 
    #     "BestObjStop" => p + 1e-5,
    #     "OutputFlag" => 0))
    # println("hola")
    @variable(m, x[1 : ncols], Bin)
    @variable(m, y[1 : nrows], Bin)
    # println("supp = $supp")
    @objective(m, Max, dot(supp, x))
    nonempty_rows = [i for i in 1 : nrows if sum(@view coverage[i, :]) > 0]
    @constraint(m, [i in nonempty_rows], dot(coverage[i, :], x) - min(p, sum(@view coverage[i, :])) * y[i] <= 0)
    @constraint(m, sum(y) <= k)
    @constraint(m, sum(x) >= p)
    if !isempty(conflict)
        # println("adding $(length(cliques)) clique constraints")
        @constraint(m, [cl in conflict], sum(x[i] for i in cl) <= 1)
    end
    if !isempty(subset)
        # println("adding $(length(subset)) precedences")
        @constraint(m, [(i, j) in subset], x[i] - x[j] >= 0)
    end
    if false && !isempty(equiv)
        # println("adding $(length(subset)) precedences")
        @constraint(m, [(i, j) in equiv], x[i] - x[j] == 0)
    end
    optimize!(m)
    # println(termination_status(m))
    if has_values(m)#termination_status(m) != MOI.INFEASIBLE# && objective_bound(m) >= p - 1e-7
        # println("feasible setpacking")
        xval = round.(Int64, value.(x))
        yval = round.(Int64, value.(y))
        # println("xval = $([j for j in 1 : ncols if xval[j] == 1])")
        # println("yval = $([i for i in 1 : nrows if yval[i] == 1])")
        sol = [i for i in 1 : ncols if xval[i] == 1]
        return sol[1 : p]
    else 
        # println("infeasible setpacking")
        return Int64[]
    end
end
