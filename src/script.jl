
data = DiscreteOrderedMedian.read_deleplanque("home/mattgru/orderedPMedian-MIP/Debug/JDD/domp20p5v1.domp")
data = DiscreteOrderedMedian.modify_lambda(data, :T1)
DiscreteOrderedMedian.bnb(data)