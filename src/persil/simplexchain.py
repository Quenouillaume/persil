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


    def faces(self,complex):
        res = []
        for i in range(self.dim):
            res.append(complex._indexBySimplex[Simplex(self.vertices[:i]+self.vertices[i+1:])])
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
    def __init__(self,simplexCoeffList,complex):
        self.coeffs = {}
        self.complex = complex
        for (j,c) in simplexCoeffList:
            self.coeffs[j] = c % self.complex.field

    def getCoeff(self,j):
        if j in self.coeffs:
            return self.coeffs[j]
        else:
            return 0

    def purge(self):   # removes simplices with coefficient 0
        toBeRemoved= []
        for j in self.coeffs:
            if self.coeffs[j] == 0:
                toBeRemoved.append(j)
        for j in toBeRemoved:
            self.coeffs.pop(j)

    def isEmpty(self):
        for j in self.coeffs:
            if self.coeffs[j] != 0:
                return False
        return True

    def __add__(self,other):
        res = SimplexChain([],self.complex)
        for j in self.coeffs:
            res.coeffs[j] = self.coeffs[j]
        for j in other.coeffs:
            if j in res.coeffs:
                res.coeffs[j] = (other.coeffs[j] + res.coeffs[j])%self.complex.field
            else:
                res.coeffs[j] = other.coeffs[j]
            if res.coeffs[j] == 0:
                res.coeffs.pop(j)
        return res

    def __neg__(self):
        res = SimplexChain([],self.complex)
        for j in self.coeffs:
            res.coeffs[j] = (-self.coeffs[j])%self.complex.field
        return res


    def __sub__(self,other):
        return self + (-other)

    def __rmul__(self,other):
        res = SimplexChain([],self.complex)
        for s in self.coeffs:
            res.coeffs[j] = (other*self.coeffs[j])%self.complex.field
        return res

    def __str__(self):
        res = []
        for j in self.coeffs:
            res.append(str(self.coeffs[j]) + " * " + str(self.complex.simplices[j]))
        return "  " + "\n+ ".join(res)

    def __repr__(self):
        return str(self)



def boundary(schain):
    res = SimplexChain([],schain.complex)

    for j in schain.coeffs:
        faces = j.faces()
        l = [(faces[i],(-1)**i) for i in range(shain.complex.simplices[j].dim) ]
        res += schain.coeffs[s] * SimplexChain(l,schain.complex)
    return res

def simplexBoundary(s,complex):
    res = SimplexChain([],complex)
    if s.dim == 1:
        return res
    faces = s.faces()
    l = [(faces[i],(-1)**i) for i in range(s.dim) ]
    res += SimplexChain(l,complex)
    return res