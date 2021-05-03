from .simplexchain import *
from .homology import *

from numpy import sqrt
import matplotlib.pyplot as plt
import numpy as np


def euclidianDistance(x,y):
    res = 0
    for i in range(len(x)):
        res += (x[i]-y[i])**2
    return sqrt(res)




class RipsComplex:
    """
    Class used to process points in a metric space.

    Parameters:

    pointList : List of points.

    distance : function with signature point * point -> float

    threshold : float, maximum distance which will be considered
    when constructing the complex

    """
    def __init__(self, pointList, distance = euclidianDistance, threshold = None):
        self.points = pointList[:]
        self.distance = distance
        self.nPoints = len(pointList)
        self.compute_dist_matrix()
        if not threshold:
            self.threshold = max([max(a) for a in self.matrix])+1
        else:
            self.threshold = threshold

    def compute_dist_matrix(self):
        self.matrix = [[ self.distance(x,y) for x in self.points ] for y in self.points]



    def plot(self,threshold = None):
        """
        Plots the 1-skeleton of the points.
        It is possible to specify a threshold
        different from the threshold chosen
        during initialization.
        """
        if not threshold:
            threshold = self.threshold

        for x in self.points:
            plt.plot(x[0],x[1],'ro')
            for y in self.points:
                if x < y and self.distance(x,y) < threshold:
                    plt.plot([x[0],y[0]],[x[1],y[1]],'k-')
        plt.show()


    def compute_skeleton(self,maxDimension = None):
        """
        Computes the Rips-Vietoris complex of the points, with the
        chosen threshold, and up to a given maximum dimension. If
        no maximum dimension is specified, all simplices are computed.
        Algorithm comes from Afra Zomorodian :
        https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.210.426&rep=rep1&type=pdf
        """
        if not maxDimension:
            maxDimension = self.nPoints

        self.complex = FilteredComplex()
        prevSimplices = []
        currSimplices = []
        # initialize 1-skeleton
        for x in range(self.nPoints):
            for y in range(self.nPoints):
                self.complex.insert([x],0)
                if x < y and self.matrix[x][y] < self.threshold:
                    self.complex.insert([x,y],self.matrix[x][y])
                    currSimplices.append([x,y])


        # compute i+1-skeleton from i-skeleton
        for i in range(1,maxDimension):
            print("Computing {}-simplices".format(i))
            prevSimplices = currSimplices[:]
            currSimplices = []
            count = 0
            for l in prevSimplices:
                count +=1
                if count%100 ==0:
                    print(count)
                neighbours = self.lowerNeighbours(l)

                for v in neighbours:
                    s = Simplex(l+[v])
                    value = self.computeWeight(s)
                    self.complex.append(s,value)
                    currSimplices.append(l+[v])



    def lowerNeighbours(self,l):
        res = []
        for v in range(min(l)):
            inNbrs = True
            for u in l:
                if self.matrix[u][v] >= self.threshold:
                    inNbrs = False
                    break
            if inNbrs:
                res.append(v)

        return res

    def computeWeight(self,s):
        res = self.complex.degree(s)
        if res>=0:
            return res

        for f in s.faces():
            res = max(res,self.computeWeight(f))
        return res





