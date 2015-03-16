=======
# Graceful Graph Package

We provide here a sagemath implementation of the Integer formula encoding package
**UPDATE 2015-03-16** 

The `Graceful Graph Package` is a symbolic package designed to
investigate structural and combinatorial properties graphs with n-1
which admit graceful labelings.

# Installation 

A properly working install of [sage](http://sagemath.org/) is a prerequisite to using the 
graceful graph package. Download the sage file into your working directory. Load the package 
into a sage terminal session using the following command:

```python
sage: %runfile("Graceful_Trees_Package_code.sage")
```

# Usage

To create all connected graphs on n vertices which admit a graceful labeling we use
the following commands 

```python
sage: PlotGracefulGraphs(4)
```
The computation produces in the working directory png files depiciting the graphs obtain which
admit a graceful labeling. The permutation of the vertices which yield the graceful labeling is
given by the name of the file.
If instead we wanted a list of adjacency matrices for the graphs we ould use the following command

```python
sage: GracefulAdjacencyMatrices(4)
[
[0 1 1 1]  [0 0 1 1]  [0 0 0 1]  [0 0 0 1]  [0 1 1 1]  [0 0 1 1]
[1 0 0 0]  [0 0 1 0]  [0 0 1 1]  [0 0 0 1]  [1 0 0 0]  [0 0 0 1]
[1 0 0 0]  [1 1 0 0]  [0 1 0 0]  [0 0 0 1]  [1 0 0 0]  [1 0 0 0]
[1 0 0 0], [1 0 0 0], [1 1 0 0], [1 1 1 0], [1 0 0 0], [1 1 0 0],

[0 0 1 0]  [0 0 1 0]  [0 0 0 1]  [0 0 0 1]  [0 1 1 1]  [0 1 0 1]
[0 0 1 1]  [0 0 1 0]  [0 0 1 0]  [0 0 0 1]  [1 0 0 0]  [1 0 1 0]
[1 1 0 0]  [1 1 0 1]  [0 1 0 1]  [0 0 0 1]  [1 0 0 0]  [0 1 0 0]
[0 1 0 0], [0 0 1 0], [1 0 1 0], [1 1 1 0], [1 0 0 0], [1 0 0 0],

[0 0 1 0]  [0 0 1 0]  [0 1 1 1]  [0 1 1 0]  [0 1 1 1]  [0 1 0 1]
[0 0 0 1]  [0 0 1 0]  [1 0 0 0]  [1 0 0 1]  [1 0 0 0]  [1 0 0 0]
[1 0 0 1]  [1 1 0 1]  [1 0 0 0]  [1 0 0 0]  [1 0 0 0]  [0 0 0 1]
[0 1 1 0], [0 0 1 0], [1 0 0 0], [0 1 0 0], [1 0 0 0], [1 0 1 0],

[0 1 0 0]  [0 1 0 0]  [0 1 0 0]  [0 1 0 0]  [0 1 1 1]  [0 1 1 0]
[1 0 1 0]  [1 0 1 1]  [1 0 0 1]  [1 0 1 1]  [1 0 0 0]  [1 0 0 0]
[0 1 0 1]  [0 1 0 0]  [0 0 0 1]  [0 1 0 0]  [1 0 0 0]  [1 0 0 1]
[0 0 1 0], [0 1 0 0], [0 1 1 0], [0 1 0 0], [1 0 0 0], [0 0 1 0]
]
```


# Bug report

Please report any bugs, it is greatly appreciated.
