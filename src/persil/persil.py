import numpy as np

class Simplex:
    def __init__(self,l):
        self.vertices = l[:]
        self.vertices.sort()
        self.dim = len(l)



    def __eq__(self,other):
        return (self.vertices == other.vertices)

    def __lt__(self,other): # should only be used for simplices of same dimension
        for (x,y) in zip(self.vertices,other.vertices):
            if x>y:
                return False
        return True


    def __hash__(self):
        return hash(tuple(self.vertices))

    def __str__(self):
        return str(self.vertices)

    def __repr__(self):
        return "Simplex{}".format(str(self.vertices))


    def faces(self):
        res = []
        for i in range(self.dim):
            res.append(Simplex(self.vertices[:i]+self.vertices[i+1:]))
        return res



def simplexOrder(s1,s2): # returns True iff s1 <= s2
    for (x,y) in zip(s1.vertices,s2.vertices):
        if x>y:
            return False
    return True




def isFace(s1,s2):
    if s1.dim+1 != s2.dim:
        return False
    for i in s1.vertices:
        if i not in s2.vertices:
            return False
    return True



class SimplexChain:
    def __init__(self,simplexCoeffList,field = 2):
        self.coeffs = {}
        self.field = field
        for (s,c) in simplexCoeffList:
            self.coeffs[s] = c % self.field

    def getCoeff(self,s):
        if s in self.coeffs:
            return self.coeffs[s]
        else:
            return 0

    def purge(self):
        toBeRemoved= []
        for s in self.coeffs:
            if self.coeffs[s] == 0:
                toBeRemoved.append(s)
        for s in toBeRemoved:
            self.coeffs.pop(s)

    def isEmpty(self):
        for s in self.coeffs:
            if self.coeffs[s] != 0:
                return False
        return True

    def __add__(self,other):
        res = SimplexChain([],field = self.field)
        for s in self.coeffs:
            res.coeffs[s] = self.coeffs[s]
        for s in other.coeffs:
            if s in res.coeffs:
                res.coeffs[s] = (other.coeffs[s] + res.coeffs[s])%self.field
            else:
                res.coeffs[s] = other.coeffs[s]
            if res.coeffs[s] == 0:
                res.coeffs.pop(s)
        return res

    def __neg__(self):
        res = SimplexChain([],field = self.field)
        for s in self.coeffs:
            res.coeffs[s] = (-self.coeffs[s])%res.field
        return res


    def __sub__(self,other):
        return self + (-other)

    def __rmul__(self,other):
        res = SimplexChain([],field = self.field)
        for s in self.coeffs:
            res.coeffs[s] = (other*self.coeffs[s])%res.field
        return res

    def __str__(self):
        res = []
        for s in self.coeffs:
            res.append(str(self.coeffs[s]) + " * " + str(s))
        return "  " + "\n+ ".join(res)

    def __repr__(self):
        return str(self)




def boundary(schain):
    res = SimplexChain([],field = schain.field)

    for s in schain.coeffs:
        faces = s.faces()
        l = [(faces[i],(-1)**i) for i in range(s.dim) ]
        res += schain.coeffs[s] * SimplexChain(l,field = schain.field)
    return res

def simplexBoundary(s,sfield):
    res = SimplexChain([],field = sfield)
    if s.dim == 1:
        return res
    faces = s.faces()
    l = [(faces[i],(-1)**i) for i in range(s.dim) ]
    res += SimplexChain(l,field = sfield)
    return res



class FilteredComplex:
    # the degree of a simplex is the lowest index for which it appears in the complex

    def __init__(self,warnings = False):
        self._simplexes = [] # list of simplices, ordered by dimension, then by degree, then by an arbitrary predefined order to break ties
        self._degrees_dict = {} # also contains the degrees. keys are simplices
        self._numSimplexes = 0
        self._dimension = 0
        self._warnings = warnings
        self._maxDeg = 0


    def degree(self,s): # check if simplex s is already in the complex, returns the degree if it is, -1 otherwise
        if s in self._degrees_dict:
            return self._degrees_dict[s]
        else:
            return -1



    def order(self,s1,d1,s2,d2): # return True if s1 with degree d1 is greater than s2 with degree d2
        if s1.dim > s2.dim or (s1.dim == s2.dim and d1 > d2) or (s1.dim == s2.dim and d1 == d2 and s1 > s2 ):
            return True
        return False

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

        if update:
            j = self._simplexes.index(s)
            self._simplexes.pop(j)
            self._numSimplexes -= 1

        i = 0
        while i<self._numSimplexes:

            #if s.dim > self._simplexes[i].dim or (s.dim == self._simplexes[i].dim and d > self._degrees_dict[self._simplexes[i]]) or (s.dim == self._simplexes[i].dim and d == self._degrees_dict[self._simplexes[i]] and self._simplexes[i] < s ):
            if self.order(s,d,self._simplexes[i],self._degrees_dict[self._simplexes[i]]):
                i+=1
            else:
                break

        self._degrees_dict[s] = d
        self._simplexes = self._simplexes[:i] + [s] + self._simplexes[i:]
        self._numSimplexes += 1
        self._dimension = max(self._dimension,s.dim)
        self._maxDeg = max(self._maxDeg,d)

    def insert(self,l,d):
        self.append(Simplex(l),d)

    def __str__(self):

        return  "\n".join(["{} : {}".format(s,self._degrees_dict[s]) for s in self._simplexes])






class ZomorodianCarlsson:
    def __init__(self,filteredComplex,field = 2,strict = False,verbose = False):
        self.simplices = filteredComplex._simplexes[:]
        self.n = filteredComplex._numSimplexes
        self.dim = filteredComplex._dimension
        self.degrees = filteredComplex._degrees_dict.copy()
        self.field = field


        self.marked = set()
        self.T = {} # contains couples (simplex,chain)
        self.intervals = {}
        self.pairs = []

        self._maxDeg = filteredComplex._maxDeg
        self._strict = strict

        self.verbose = verbose

    def addInterval(self,k,t,s):
        i = self.degrees[t]
        if not s:
            j = np.inf
        else:
            j = self.degrees[s]

        if k not in self.intervals:
            self.intervals[k] = []
        if i != j or (not self._strict):
            if self.verbose:
                print("Adding {}-interval ({},{})".format(k,i,j))
            self.intervals[k].append((i,j))
            self.pairs.append((t,s))
            #if i == 10:
            #    print(k,t,s)



    def computeIntervals(self):
        for j in range(self.n):
            s = self.simplices[j]
            if self.verbose:
                print("Examining {}. Removing pivot rows...".format(s))
            d = self.removePivotRows(s)
            if self.verbose:
                print("Done removing pivot rows")
            if d.isEmpty():
                if self.verbose:
                    print("Boundary is empty when pivots are removed: marking {}".format(s))
                self.marked.add(s)
            else:

                t,maxInd = self.maxIndex(d)
                k = t.dim-1
                self.T[t] = (s,d)
                self.addInterval(k,t,s)
                if self.verbose:
                    print("Boundary non-reducible: T{} is set to:".format(t))
                    print(str(d))

        if self.verbose:
            print("First pass over, beginning second pass")
        for j in range(self.n):
            s = self.simplices[j]
            if s in self.marked and s not in self.T:
                k = s.dim -1
                if self.verbose:
                    print("Infinite interval found for {}.".format(s))
                self.addInterval(k,s,None)
        if self.verbose:
            print("Second pass over")



    def removePivotRows(self,s):
        d = simplexBoundary(s,self.field)
        for s in d.coeffs:
            if s not in self.marked:
                d.coeffs[s] = 0
        d.purge()
        while not d.isEmpty():

            if self.verbose:
                print("Current chain d:")
                print(str(d))

            t,maxInd = self.maxIndex(d)
            if self.verbose:print("simplex with max index in d: {} with index {}".format(t,maxInd))

            if t not in self.T:
                if self.verbose:
                    print("{} is not in T: done removing pivot rows".format(t))
                break

            c = self.T[t][1]
            q = c.getCoeff(t)
            if self.verbose:
                print("{} is in T with coeff {}: ".format(t,q),"##########",str(c),"##########",sep='\n'    )
            d = d - pow(q,d.field-2,d.field)*self.T[t][1]

        return d






    def maxIndex(self,d):
        currmax = -1
        res = None
        for s in d.coeffs:
            i = self.simplices.index(s)
            if i>currmax:
                currmax = i
                res = s
        return (res,currmax)



    def getIntervals(self,d):
        if d not in self.intervals:
            return []
        else:
            return self.intervals[d][:]

    def bettiNumber(self,k,l,p):
        res = 0
        if k not in self.intervals:
            return 0
        for (i,j) in self.intervals[k]:
            if (i <= l and l + p < j ) and p>=0:
                    res+=1
        return res












