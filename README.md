=======
# Graceful Graph Package

We provide here a sagemath implementation of the graceful graph package
**UPDATE 2016-05-03** 

The `Graceful Graph Package` is a combinatorial package designed to
investigate structural and combinatorial properties graphs with n-1
edges which admit graceful labelings.

# Installation 

A properly working install of [SageMath Version 7.1](http://sagemath.org/) is a prerequisite to using the 
graceful graph package. Download the sage file into your working directory. Load the package 
into a sage terminal session using the following command:

```python
sage: %runfile graceful_graph_package.sage
```

# Usage

The following goes through all the permutation of n>1 elements and outputs
the list of permutations which can be used for constructing
graceful graphs.

```python
sage: GracefulPermutations(4)
[[0, 1, 2, 3], [0, 2, 1, 3]]
```

The following goes through all the permutation of n>1 elements and outputs
the list of of functions derived from graceful permutations.

```python
sage: CountGracefulFunctions(4)
[[2, [0, 1, 2, 3]], [2, [0, 2, 1, 3]]]
```

The following goes through all the permutation of n>1 elements and outputs
the list of of functions on the vertices derived from graceful permutations.

```python
sage: GracefulFunctionsTuples(3)
[[[[(0, 0), (1, 2), (2, 0)], [(0, 0), (1, 0), (2, 0)]], [0, 1, 2]]]
```

The following converts a list of tuples associated with directed edges to
list of symbolic expression associated with the edges

```python
sage: Tuple2EdgeList([(0, 0), (1, 2), (2, 0)], 'x')
[0, x1 - x2, -x0 + x2]
```

The following saves in the working directory png files associated with
drawings of all gracefully labeled graphs having n-1 edges where n is the input

```python
sage: for l in GracefulDiGraphList(5):
....:     t=0
....:     for dg in l[0]:
....:         dg.plot().save(str(l[1]).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
....:         t=t+1
....:
```

The following saves only the ones which are trees.

```python
sage: for l in GracefulDiTreeList(5):
....:     t=0
....:     for dg in l[0]:
....:         dg.plot().save(str(l[1]).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
....:         t=t+1
....:
```

The following saves only the ones which are generator trees.

```python
sage: for l in GeneratorGracefulTreeDiGraphList(5):
....:     t=0
....:     for dg in l[0]:
....:         dg.plot().save(str(l[1]).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
....:         t=t+1
....:
```

(to be continued ...)

# Bug report

Please report any bugs, it is greatly appreciated.
