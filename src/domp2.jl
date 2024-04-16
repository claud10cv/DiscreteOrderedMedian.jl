function domp2(data::DOMPData, params::Parameters)::Tuple{Float64, Int64, Vector{Int64}}
    d = sort(unique(data.D))
    if d[1] != 0
        pushfirst!(d, 0)
    end
    nd = length(d)
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    m = JuMP.direct_model(params.optimizer_data.optimizer())
    params.optimizer_data.set_attributes(m, params.time_limit, 100000000)
    @variable(m, x[1 : nrows, 1 : ncols], Bin)
    @variable(m, y[1 : ncols], Bin)
    @variable(m, u[1 : nrows, 1 : nd], Bin)

    @objective(m, Min, sum(data.lambda[k] * (d[h] - d[h - 1]) * u[k, h] for k in 1 : nrows, h in 2 : nd))

    @constraint(m, [k in 1 : nrows], u[k, 1] == 0)
    @constraint(m, [k in 1 : nrows], u[k, nd] == 0)

    @constraint(m, sum(y) == data.p)
    @constraint(m, [i in 1 : nrows], sum(x[i, :]) == 1)
    @constraint(m, [i in 1 : nrows, j in 1 : ncols], x[i, j] - y[j] <= 0)
    @constraint(m, [k in 1 : nrows, h in 2 : nd - 1], u[k, h] - u[k, h + 1] >= 0)
    @constraint(m, [k in 1 : nrows - 1, h in 2 : nd], u[k + 1, h] - u[k, h] >= 0)
    @constraint(m, [h in 2 : nd], sum(x[i, j] for i in 1 : nrows for j in 1 : ncols if data.D[i, j] > d[h - 1]) - sum(u[k, h] for k in 1 : nrows) == 0)

    optimize!(m)
    bb = objective_bound(m)
    if has_values(m)
        yval = round.(Int64, value.(y))
        sol = [j for j in 1 : ncols if yval[j] == 1]
        dstar = compute_sorted_distances(data, yval)
        cost = compute_weighted_cost(data, dstar)
        return bb, cost, sol
    else
        return bb, typemax(Int64), Int64[]
    end
end