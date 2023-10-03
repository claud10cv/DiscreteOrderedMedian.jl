# DiscreteOrderedMedian.jl
Branch-and-bound algorithm for the discrete ordered median problem

## Installation
From the REPL, execute the following line
```julia
julia> Pkg.add("https://github.com/claud10cv/DiscreteOrderedMedian.jl")
```

## Basic Usage
1. Read an instance from Deleplanque et al executing
```julia
julia> data = DiscreteOrderedMedian.read_deleplanque(filename)
```
2. Optionally, you can modify the lambda vector, executing for instance
```julia
julia> data = DiscreteOrderedMedian.modify_lambda(data, :T1)
```
3. Solve the problem using our B&B algorithm executing
```julia
julia> DiscreteOrderedMedian.bnb(data)
```

julia --project="/path/to/Rep" ./path/to/script_beas.jl "/path/to/data/file"
