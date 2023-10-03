using DiscreteOrderedMedian
data = DiscreteOrderedMedian.read_beasley(ARGS[1])
data = DiscreteOrderedMedian.modify_lambda(data, :T9)
DiscreteOrderedMedian.bnb(data)
