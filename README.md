# PersiL
A simple python project for computing **Persi**stent homo**L**ogy. My goal is to have some code nice enough to be able to add it to [SageMath](https://www.sagemath.org "Tiens tiens quelle surprise :)") !
The algorithm implemented is described in this [article by Afra Zomorodian and Gunnar Carlsson](https://geometry.stanford.edu/papers/zc-cph-05/zc-cph-05.pdf). Many thanks to them for writing such a clear article !

If you encounter any issues, please don't hesitate to contact me !

# Installation
## On Linux
Install directly from github with `python3 -m pip install https://github.com/quenouillaume/persil/archive/refs/tags/v0.3-beta.tar.gz`


## On Windows

Install directly from github with `python -m pip install https://github.com/quenouillaume/persil/archive/refs/tags/v0.3-beta.tar.gz`

## On your favorite python IDE
If you're using Pyzo, Jupyter Notebook or something similar, simply type `pip install https://github.com/quenouillaume/persil/archive/refs/tags/v0.3-beta.tar.gz` 


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




## TODO
* Organize files and code
* Persistence diagrams and barcodes

## IN PROGRESS
* Learn how to make python packages
* Comment code
* Add a guide
