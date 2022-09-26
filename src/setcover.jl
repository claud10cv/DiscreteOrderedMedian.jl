using JuMP, CPLEX, MathOptInterface
const MOI = MathOptInterface

function setcover(coverage::Matrix{Bool}, k::Int64, p::Int64, supp::Vector{Int64})::Vector{Int64}
    nrows = size(coverage, 1)
    ncols = size(coverage, 2)
    if p < 0 return [] end
    if k <= 0 return collect(1 : min(p, ncols)) end
    m = JuMP.direct_model(CPLEX.Optimizer())
    set_optimizer_attributes(m, 
        # "CPXPARAM_MIP_Tolerances_UpperCutoff" => p + 1e-5, 
        # "CPXPARAM_MIP_Limits_LowerObjStop" => p + 1e-5,
        "CPXPARAM_ScreenOutput" => 0,
        "CPXPARAM_Threads" => 1,
        "CPXPARAM_MIP_Limits_Solutions" => 1)
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
    @constraint(m, [i in 1 : nrows], JuMP.dot(coverage[i, :], x) - y[i] >= 0)
    @constraint(m, sum(y) >= k)
    @constraint(m, sum(x) <= p)
    optimize!(m)
    if termination_status(m) != MOI.INFEASIBLE# && objective_bound(m) <= p + 1e-7
        xval = round.(Int64, value.(x))
        return [i for i in 1 : ncols if xval[i] >= 1]
    else
        # JuMP.write_to_file(m, "error.lp") 
        return Int64[]
    end
end

function setpacking(coverage::Matrix{Bool}, k::Int64, p::Int64, supp::Vector{Int64})::Vector{Int64}
    nrows = size(coverage, 1)
    ncols = size(coverage, 2)
    # println("solving setpacking with k = $k, p = $p, nrows = $nrows, ncols = $ncols")
    # println("cov = $coverage")
    if k < 0 return [] end
    if p > ncols return [] end
    if p <= 0 return [1] end
    # println("constructing model")
    m = JuMP.direct_model(CPLEX.Optimizer())
    set_optimizer_attributes(m, 
        # "CPXPARAM_MIP_Tolerances_LowerCutoff" => p - 1e-5, 
        # "CPXPARAM_MIP_Limits_UpperObjStop" => p - 1e-5,
        "CPXPARAM_ScreenOutput" => 0,
        "CPXPARAM_Threads" => 1,
        "CPXPARAM_MIP_Limits_Solutions" => 1)
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
    @constraint(m, [i in nonempty_rows], JuMP.dot(coverage[i, :], x) - min(p, sum(@view coverage[i, :])) * y[i] <= 0)
    @constraint(m, sum(y) <= k)
    @constraint(m, sum(x) >= p)
    optimize!(m)
    # println(termination_status(m))
    if termination_status(m) != MOI.INFEASIBLE# && objective_bound(m) >= p - 1e-7
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
