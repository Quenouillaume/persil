from .simplexchain import *
from numpy import inf

# FilteredComplex class (note that filtrations are non-decreasing), and ZomorodianCarlsson class, which is used to compute homology
# Maybe some other algos may be implemented later


class FilteredComplex:
    # the degree of a simplex is the lowest index for which it appears in the complex

    def __init__(self,warnings = False):
        self._simplices = [] # list of simplices
        self._degrees_dict = {} # contains the degrees. keys are simplices
        self._numSimplices = 0
        self._dimension = 0
        self._warnings = warnings
        self._maxDeg = 0


    def degree(self,s): # check if simplex s is already in the complex, returns the degree if it is, -1 otherwise
        if s in self._degrees_dict:
            return self._degrees_dict[s]
        else:
            return -1


    def append(self,s,d): #simplex as a list of vertices, degree. Insert so that the order is preserved
        # if the simplex is already in the complex, do nothing
        update = False
        if self.degree(s)>=0:
            if self._warnings:
                print("Face {} is already in the complex.".format(s))
            if self.degree(s) > d:
                if self._warnings:
                    print("However its degree is higher than {}: updating it to {}".format(self.degree(s),d))
                update = True
                self._degrees_dict.pop(s)
            else:
                if self._warnings:
                    print("Its degree is {} which is lower than {}: keeping it that way".format(self._degrees_dict[s],d))
                return

        # check that all faces are in the complex already. If not, warn the user and add faces (recursively)
        faces = s.faces()
        if s.dim>1:
            for f in faces:

                if self._warnings:
                    print("Inserting face {} as well".format(f))
                self.append(f,d)

        if not update:
            self._numSimplices += 1
            self._simplices.append(s)
        else:
            pass

        self._degrees_dict[s] = d
        self._dimension = max(self._dimension,s.dim)
        self._maxDeg = max(self._maxDeg,d)

    def insert(self,l,d):
        self.append(Simplex(l),d)

    def __str__(self):

        return  "\n".join(["{} : {}".format(s,self._degrees_dict[s]) for s in self._simplices])






class ZomorodianCarlsson:
    def __init__(self,filteredComplex,field = 2,strict = True,verbose = False):
        """
        Class for Zomorodian and Carlsson's algorithm for persistent homology.
        Initialization does not compute homology. Call self.computeIntervals
        for the actual computation.

        Arguments:
        - filteredComplex should be an instance of FilteredComplex, on which
        homology will be computed.
        - field: prime number, specifies the field over which homology is
        computed. Default value: 2
        - strict: Boolean. If set to True, homology elements of duration zero,
        i.e. intervals of the form (x,x), will be ignored. Default value: True
        - verbose: Boolean. If set to True, the computeIntervals method will
        output its progress throughout the algo. Default value: False

        """

        self.numSimplices = filteredComplex._numSimplices

        # first, order the simplices in lexico order on dimension, degree and then arbitrary order
        def key(s):
            d = filteredComplex.degree(s)
            return (s.dim,d,s)
        filteredComplex._simplices.sort(key = key)

        self.simplices = filteredComplex._simplices[:]
        self.dim = filteredComplex._dimension
        self.degrees = filteredComplex._degrees_dict.copy()
        self.field = field


        self.marked = [False for i in range(self.numSimplices)]
        self.T = [None for i in range(self.numSimplices)] # contains couples (index,chain)
        self.intervals = [[] for i in range(self.dim+1)] # contains homology intervals once the algo has finished
        self.pairs = []

        self._maxDeg = filteredComplex._maxDeg
        self._strict = strict
        self._verbose = verbose
        self._homologyComputed = False

    def addInterval(self,k,t,s):
        i = self.degrees[t]
        if not s:
            j = inf
        else:
            j = self.degrees[s]

        if i != j or (not self._strict):
            #if self._verbose:
                #print("Adding {}-interval ({},{})".format(k,i,j))
            self.intervals[k].append((i,j))
            self.pairs.append((t,s))
            #if i == 10:
            #    print(k,t,s)



    def computeIntervals(self):
        """
        Computes the homology group. After running this method,
        homology intervals become accessible with self.getIntervals(k).

        """
        if self._homologyComputed:
            print("Homology was already computed.")
            return
        if self._verbose:
            print("Beginning first pass")
        for j in range(self.numSimplices):
            if j%1000 == 0 and self._verbose:
                print('{}/{}'.format(j,self.numSimplices))
            s = self.simplices[j]
            #if self._verbose:
                #print("Examining {}. Removing pivot rows...".format(s))
            d = self.removePivotRows(s)
            #if self._verbose:
                #print("Done removing pivot rows")
            if d.isEmpty():
                #if self._verbose:
                    #print("Boundary is empty when pivots are removed: marking {}".format(s))
                self.marked[j] = True
            else:

                maxInd = self.maxIndex(d)
                t = self.simplices[maxInd]
                k = t.dim-1
                self.T[maxInd] = (s,d)
                self.addInterval(k,t,s)
                #if self._verbose:
                    #print("Boundary non-reducible: T{} is set to:".format(t))
                    #print(str(d))

        if self._verbose:
            print("First pass over, beginning second pass")
        for j in range(self.numSimplices):
            if j%1000 == 0 and self._verbose:
                print('{}/{}'.format(j,self.numSimplices))
            s = self.simplices[j]
            if self.marked[j] and not self.T[j]:
                k = s.dim -1
                #if self._verbose:
                    #print("Infinite interval found for {}.".format(s))
                self.addInterval(k,s,None)
        if self._verbose:
            print("Second pass over")
        self._homologyComputed = True



    def removePivotRows(self,s):
        d = simplexBoundary(s,self)
        for j in d.coeffs:
            if not self.marked[j]:
                d.coeffs[j] = 0
        d.purge()
        while not d.isEmpty():

            #if self._verbose:
                #print("Current chain d:")
                #print(str(d))

            maxInd = self.maxIndex(d)
            t = self.simplices[maxInd]
            #if self._verbose:
                #print("simplex with max index in d: {} with index {}".format(maxInd))

            if not self.T[maxInd]:
                #if self._verbose:
                    #print("{} is not in T: done removing pivot rows".format(t))
                break

            c = self.T[maxInd][1]
            q = c.getCoeff(t)
            #if self._verbose:
                #print("{} is in T with coeff {}: ".format(t,q),"##########",str(c),"##########",sep='\n'    )
            d = d - pow(q,self.field-2,self.field)*c
        return d


    def maxIndex(self,d):
        currmax = -1
        for j in d.coeffs:
            if j>currmax:
                currmax = j
        return currmax



    def getIntervals(self,d):
        """
        Returns the list of d-dimensional homology elements.
        Can only be run after computeIntervals has been run.
        """
        if not self._homologyComputed:
            print("Warning: homology has not yet been computed. This will return an empty list.")
        return self.intervals[d][:]

    def bettiNumber(self,k,l,p):
        res = 0
        if not self._homologyComputed:
            print("Warning: homology has not yet been computed. This will return 0.")
        for (i,j) in self.intervals[k]:
            if (i <= l and l + p < j ) and p>=0:
                    res+=1
        return res



