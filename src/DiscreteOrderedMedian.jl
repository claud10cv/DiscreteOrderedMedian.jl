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
            data = DiscreteOrderedMedian.read_beasley(joinpath(home, "Beasley1.txt"))
            data = DiscreteOrderedMedian.modify_lambda(data, :T9)
            DiscreteOrderedMedian.bnb(data)
        end
    end
end # module
