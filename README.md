# PersiL
A simple python project for computing **Persi**stent homo**L**ogy. My goal is to have some code nice enough to be able to add it to [SageMath](https://www.sagemath.org "Tiens tiens quelle surprise :)") !
The algorithm implemented is described in this [article by Afra Zomorodian and Gunnar Carlsson](https://geometry.stanford.edu/papers/zc-cph-05/zc-cph-05.pdf). Many thanks to them for writing such a clear article !

If you encounter any issues, please don't hesitate to contact me !

# Installation
## On Unix
Install directly from github with `python3 -m pip install git+https://github.com/quenouillaume/persil`
You can update the module with `python3 -m pip install --upgrade --user git+https://github.com/quenouillaume/persil`


## On Windows
Download repository as a ZIP file, then open a terminal in your downloads folder and enter:
`python /m pip install persil-main.zip`
To update, download the latest version as a ZIP file, then, as before, open a terminal in your downloads folder and enter:
`python /m pip install --upgrade --user persil-main.zip`

## In Sage
Load Sage with
```
module load gcc/8.3.0
module load sage/9.1
```
Then, to install:
```
sage -pip install --user git+https://github.com/quenouillaume/persil/
```
To upgrade, use:
```
sage -pip install --upgrade --user git+https://github.com/quenouillaume/persil/
```



# Guide

I have not yet written a guide, but it will come soon ! For now, you can use the following code as example. 
## Simple example
This code computes the peristent homology for the complex described in Zomorodian+Carlsson's article.

```python
from persil.simplexchain import *
from persil.homology import *
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
```
If everything works fine, it should output:

```python
[(0, 1), (1, 1), (1, 2), (0, inf)]
[(3, 4), (2, 5)]
```
## Graphical representation
On the same complex as above, shows the persistence diagram for dimension 0

```python
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

persistence_diagram(zc.intervals[0])
```
If everything works fine, this should draw a nice persistence diagram with one vertical line and 3 points !


## Rips-Vietoris Complex

The RipsVietoris class is used to construct Rips complexes from a cloud point. Initialize it as follows:

```python
r = RipsComplex(pointList,distance,threshold)
```
where:
* `pointList` is a list of points
* `distance` is a function taking two points as argument and returning their distance as a float
* `threshold` is the maximum distance to be considered when constructing the complex. It can be omitted, in which case the program will compute the entire Rips complex, but this can get quite long.
* Additionally, you can add the optional argument `verbose = True`, which will make the construction print info on the progress of the construction

Once the object is initialized, compute the Complex with:
```python
r.compute_skeleton(d)
```

where `d` is the maximum desired dimension. For instance if `d` is 2, the resulting complex will contain points, edges and triangles.

Finally, you can access the Rips complex with `r.complex`. You can then compute its homology like in the previous examples.
If you are working with points in the 2d plane, you can plot the point cloud and the neibourhood graph with `r.plot()`. This technically also works with points in higher dimension but it will only plot a projection on the two first axes.

Here is a full example:
```python
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


```




## TODO
* Organize files and code


## IN PROGRESS
* Persistence diagrams and barcodes
* Learn how to make python packages
* Comment code
* Add a guide
* Rips complex
