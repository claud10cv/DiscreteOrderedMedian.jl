module DiscreteOrderedMedian
    using DataStructures
    using JuMP
    using MathOptInterface
    using CPLEX
    using LinearAlgebra
    using Random
    using Graphs
    using PrecompileTools

    include("instancegen.jl")
    include("struct.jl")
    include("reader.jl")
    include("localsearch.jl")
    include("util.jl")
    include("setcover.jl")
    include("setpacking.jl")
    include("dompk.jl")
    include("domplb.jl")
    include("bnb.jl")
    include("domp2.jl")

    const MOI = MathOptInterface
    
    @setup_workload begin
        home = pkgdir(DiscreteOrderedMedian)
        @compile_workload begin
            redirect_stdout(devnull) do
                data0 = DiscreteOrderedMedian.read_deleplanque(joinpath(home, "data", "domp20p5v1.domp"))
                for lambda in [:T9, :T1, :T0]
                    data = DiscreteOrderedMedian.modify_lambda(data0, lambda)
                    params = DiscreteOrderedMedian.default_parameters()
                    DiscreteOrderedMedian.bnb(data, params)
                end
            end
        end
    end
    
    export bnb, modify_lambda, read_deleplanque, read_orlib, default_parameters
end # module
