from persil.simplexchain import *
from persil.homology import *
from persil.graphical import *

# First step: create a list of simplices and filtration values (also called degrees)
list_simplex_degree = [([0],0), ([1],0), ([2],1), ([3],1), ([0, 1],1), ([1, 2],1), ([0, 3],2), ([2, 3],2), ([0, 2],3), ([0, 1, 2],4), ([0, 2, 3],5)]

# add them one by one in a fresh complex
fc = FilteredComplex()
for (simplex, value) in list_simplex_degree:
    fc.insert(simplex,value)

# compute persistent homology in Z/2Z
zc = ZomorodianCarlsson(fc)
zc.computeIntervals()
for i in range(2):
    print(zc.getIntervals(i))

