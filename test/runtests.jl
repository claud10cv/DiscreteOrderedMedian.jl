using DiscreteOrderedMedian
using Test

function feasible_result(res::DiscreteOrderedMedian.Result)::Bool
    return res.lb <= res.ub
end

function test_toy(seed::Int64, lambda::Symbol)::Bool
    data = DiscreteOrderedMedian.generate_euclidean(10, 5, 1000, 1000, :pmedian; seed = seed)
    data = DiscreteOrderedMedian.modify_lambda(data, lambda)
    params = DiscreteOrderedMedian.default_parameters()
    res = DiscreteOrderedMedian.bnb(data, params)
    return feasible_result(res)
end

@testset "Random Test" begin
    for seed in 1 : 10, lambda in [:T0, :T1, :T9]
        redirect_stdout(devnull) do
            @test test_toy(seed, lambda)
        end
    end
end