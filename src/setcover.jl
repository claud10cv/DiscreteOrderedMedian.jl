function setcover(params::Parameters, coverage::Matrix{Bool}, k::Int64, p::Int64, supp::Vector{Int64})::Vector{Int64}
    nrows = size(coverage, 1)
    ncols = size(coverage, 2)
    if p < 0 return [] end
    if k <= 0 return collect(1 : min(p, ncols)) end
    result = setcover_localsearch(coverage, k, p, supp)
    if isempty(result)
        result = setcover_exact(params, coverage, k, p, supp)
    end
    return result
end

function setcover_localsearch(coverage::Matrix{Bool}, k::Int64, p::Int64, supp::Vector{Int64})::Vector{Int64}
    nrows = size(coverage, 1)
    ncols = size(coverage, 2)
    cols = collect(1 : ncols)
    result = Int64[]
    covered = falses(nrows)
    newcov = deepcopy(coverage)
    scov = zeros(Int64, ncols)
    for i in 1 : min(p, ncols)
        # newcov = Matrix{Bool}(undef, nrows, ncols)
        for j in cols
            newcov[:, j] = coverage[:, j] .& .!covered
            scov[j] = sum(@view newcov[:, j])
        end
        sort!(cols; lt = (u, v) -> scov[u] > scov[v] || (scov[u] == scov[v] && supp[u] > supp[v]))
        jstar = popfirst!(cols)
        push!(result, jstar)
        covered = covered .| coverage[:, jstar]
    end
    if sum(covered) >= k 
        # println("feasible heuristic p = $p, k = $k, $result, coverage of $(sum(covered))")
        return sort(result)
    else 
        # println("infeasible heuristic")
        return Int64[]
    end
end

function setcover_exact(params::Parameters, coverage::Matrix{Bool}, k::Int64, p::Int64, supp::Vector{Int64})::Vector{Int64}
    nrows = size(coverage, 1)
    ncols = size(coverage, 2)
    m = JuMP.direct_model(params.optimizer_data.optimizer())
    params.optimizer_data.set_attributes(m, params.time_limit + 1, 1)
    # m = JuMP.direct_model(CPLEX.Optimizer())
    # set_optimizer_attributes(m, 
    #     # "CPXPARAM_MIP_Tolerances_UpperCutoff" => p + 1e-5, 
    #     # "CPXPARAM_MIP_Limits_LowerObjStop" => p + 1e-5,
    #     "CPXPARAM_ScreenOutput" => 0,
    #     "CPXPARAM_Threads" => 1,
    #     "CPXPARAM_MIP_Tolerances_MIPGap" => 0.1,
    #     "CPXPARAM_MIP_Limits_Solutions" => 1)
    # m = JuMP.Model(optimizer_with_attributes(Gurobi.Optimizer, 
    #     "Cutoff" => p + 1e-5, 
    #     "BestObjStop" => p + 1e-5,
    #     "OutputFlag" => 0,
    #     "Threads" => 1,
    #     "SolutionLimit" => 1))
    @variable(m, x[1 : ncols], Bin)
    @variable(m, 0 <= y[1 : nrows] <= 1)
    nonempty_rows = [i for i in 1 : nrows if sum(@view coverage[i, :]) > 0]
    if length(nonempty_rows) != nrows
        println("some empty rows = $(nrows - length(nonempty_rows))")
    end
    @objective(m, Max, dot(supp, x))
    # @objective(m, Min, dot(rand(0 : 1, ncols), x))
    @constraint(m, [i in 1 : nrows], dot(coverage[i, :], x) - y[i] >= 0)
    @constraint(m, sum(y) >= k)
    @constraint(m, sum(x) <= p)
    optimize!(m)
    if has_values(m)#termination_status(m) != MOI.INFEASIBLE# && objective_bound(m) <= p + 1e-7
        xval = round.(Int64, value.(x))
        f = 1
        while f <= ncols && sum(xval) < p
            xval[f] = 1
            f += 1
        end
        return [i for i in 1 : ncols if xval[i] >= 1]
    else
        # JuMP.write_to_file(m, "error.lp") 
        return Int64[]
    end
end