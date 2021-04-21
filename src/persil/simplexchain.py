# Base classes for simplices and simplex chains


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

    def purge(self):   # removes simplices with coefficient 0
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