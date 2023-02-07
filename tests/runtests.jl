using DiscreteOrderedMedian
using Test

BASE_FOLDER = dirname(dirname(pathof(DiscreteOrderedMedian)))
# INSTANCE_FOLDER = "$BASE_FOLDER/instances"

function feasible_result(res::DiscreteOrderedMedian.Result)::Bool
    return res.lb <= res.ub
end

function test_toy(seed::Int64)::Bool
    data = DiscreteOrderedMedian.generate_euclidean(10, 5, 1000, 1000, :pmedian; seed = seed)
    res = DiscreteOrderedMedian.bnb(data)
    return feasible_result(res)
end

@testset "Random Test" begin
    for seed in 1 : 10
        print("running test for seed $seed...")
        redirect_stdout(devnull) do
            @test test_toy(seed)
        end
        println("done!")
    end
end