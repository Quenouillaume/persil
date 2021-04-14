# PersiL
A simple python project for computing **Persi**stent homo**L**ogy. My goal is to have some code nice enough to be able to add it to [SageMath](https://www.sagemath.org "Tiens tiens quelle surprise :)") !
The algorithm implemented is described in this [article by Afra Zomorodian and Gunnar Carlsson](https://geometry.stanford.edu/papers/zc-cph-05/zc-cph-05.pdf). Many thanks to them for writing such a clear article !



# Installation
## On Linux
Install directly from github with `python3 -m pip install https://github.com/quenouillaume/persil/archive/refs/tags/v0.1-alpha.tar.gz`


## On Windows

Install directly from github with `python -m pip install https://github.com/quenouillaume/persil/archive/refs/tags/v0.1-alpha.tar.gz`

## On your favorite python IDE

Install directly from github with `pip install https://github.com/quenouillaume/persil/archive/refs/tags/v0.1-alpha.tar.gz`


# Guide

I have not yet written a guide, but it will come soon ! For now, you can use the following code as example. This code computes the peristent homology for the complex described in Zomorodian+Carlsson's article.

```python
from persil.persil import *
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




## TODO
* Organize files and code
* Persistence diagrams and barcodes

## IN PROGRESS
* Learn how to make python packages
* Comment code
* Add a guide
