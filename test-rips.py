import random
from numpy import sqrt
from persil import *


def euclidianDistance(x,y):
    return sqrt( sum([(x[i] - y[i])**2 for i in range(len(x))]))


def randomPoints(n,D): # returns a list of n random points in the unit cube of dimension D
    l = [ tuple([random.random() for i in range(D)]) for j in range(n)]
    return l



l = randomPoints(1000,3)

r = RipsComplex(l,euclidianDistance,0.23,verbose = True)

r.compute_skeleton(2)

zc = ZomorodianCarlsson(r.complex, strict = True,verbose = True)

zc.computeIntervals()

persistence_diagram(zc.intervals[1])