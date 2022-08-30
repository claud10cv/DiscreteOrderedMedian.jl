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
    # @constraint(m, sum(x) <= p)
    optimize!(m)
    if termination_status(m) != MOI.INFEASIBLE && objective_bound(m) <= p + 1e-7
        # println("termination status = $(termination_status(m))")
        # println("best bound = $(objective_bound(m)), p = $p")
        xval = value.(x)
        # println("xval = $(sum(xval)), p = $p")
        return [i for i in 1 : ncols if abs(xval[i]) > 1e-1]
    else return Int64[]
    end
end
