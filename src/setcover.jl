using JuMP, CPLEX, MathOptInterface
const MOI = MathOptInterface

function setcover(coverage::Matrix{Bool}, k::Int64, p::Int64)::Vector{Int64}
    nrows = size(coverage, 1)
    ncols = size(coverage, 2)
    if p < 0 return [] end
    if k <= 0 return collect(1 : min(p, ncols)) end
    m = JuMP.Model(optimizer_with_attributes(CPLEX.Optimizer, 
        "CPXPARAM_MIP_Tolerances_UpperCutoff" => p + 1e-5, 
        "CPXPARAM_MIP_Limits_LowerObjStop" => p + 1e-5,
        "CPXPARAM_ScreenOutput" => 0,
        "CPXPARAM_Threads" => 1))
    # m = JuMP.Model(optimizer_with_attributes(Gurobi.Optimizer, 
    #     "Cutoff" => p + 1e-5, 
    #     "BestObjStop" => p + 1e-5,
    #     "OutputFlag" => 0))
    @variable(m, x[1 : ncols], Bin)
    @variable(m, 0 <= y[1 : nrows] <= 1)
    @objective(m, Min, sum(x))
    @constraint(m, [i in 1 : nrows], JuMP.dot(coverage[i, :], x) - y[i] >= 0)
    @constraint(m, sum(y) >= k)
    optimize!(m)
    if termination_status(m) != MOI.INFEASIBLE && objective_bound(m) <= p + 1e-7
        xval = value.(x)
        return [i for i in 1 : ncols if abs(xval[i]) > 1e-1]
    else return Int64[]
    end
end

function setpacking(coverage::Matrix{Bool}, k::Int64, p::Int64)::Vector{Int64}
    nrows = size(coverage, 1)
    ncols = size(coverage, 2)
    # println("solving setpacking with k = $k, p = $p, nrows = $nrows, ncols = $ncols")
    # println("cov = $coverage")
    if k < 0 return [] end
    if p > ncols return [] end
    if p <= 0 return [1] end
    # println("constructing model")
    m = JuMP.Model(optimizer_with_attributes(CPLEX.Optimizer, 
        "CPXPARAM_MIP_Tolerances_LowerCutoff" => p - 1e-5, 
        "CPXPARAM_MIP_Limits_UpperObjStop" => p - 1e-5,
        "CPXPARAM_ScreenOutput" => 0,
        "CPXPARAM_Threads" => 1))
    # m = JuMP.Model(optimizer_with_attributes(Gurobi.Optimizer, 
    #     "Cutoff" => p + 1e-5, 
    #     "BestObjStop" => p + 1e-5,
    #     "OutputFlag" => 0))
    @variable(m, x[1 : ncols], Bin)
    @variable(m, y[1 : nrows], Bin)
    @objective(m, Max, sum(x))
    nonempty_rows = [i for i in 1 : nrows if sum(coverage[i, :]) > 0]
    @constraint(m, [i in nonempty_rows], JuMP.dot(coverage[i, :], x) - sum(coverage[i, :]) * y[i] <= 0)
    @constraint(m, sum(y) <= k)
    optimize!(m)
    # println(termination_status(m))
    if termination_status(m) != MOI.INFEASIBLE && objective_bound(m) >= p - 1e-7
        xval = round.(Int64, value.(x))
        yval = round.(Int64, value.(y))
        # println("xval = $([j for j in 1 : ncols if xval[j] == 1])")
        # println("yval = $([i for i in 1 : nrows if yval[i] == 1])")
        sol = [i for i in 1 : ncols if xval[i] == 1]
        if length(sol) > p
            resize!(sol, p)
        end
    else return Int64[]
    end
end
