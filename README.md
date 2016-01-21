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

The following goes through all the permutation of n>1 elements and outputs
the list of permutations which can be used for constructing
graceful graphs with all the vertices have degree>0.

```python
sage: GracefulPermutations(4)
[[1, 2, 3], [2, 1, 3]]
```

The following goes through all the permutation of n>1 elements and outputs
the list of of functions derived from graceful permutations.

```python
sage: CountGracefulFunctions(4)
[[2, [1, 2, 3]], [2, [2, 1, 3]]]
```
The following goes through all the permutation of n>1 elements and outputs
the list of of functions on the vertices derived from graceful permutations.
```python
sage: GracefulFunctions(3)
[[[[0, 0], [2, 0]], [1, 2]]]
```

(to be continued ...)

# Bug report

Please report any bugs, it is greatly appreciated.
