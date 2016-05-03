#*************************************************************************#
#       Copyright (C) 2016 Edinah K. Gnang <kgnang@gmail.com>,            #
#                                                                         #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                                                                         #
#    This code is distributed in the hope that it will be useful,         #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU     #
#    General Public License for more details.                             #
#                                                                         #
#  The full text of the GPL is available at:                              #
#                                                                         #
#                  http://www.gnu.org/licenses/                           #
#*************************************************************************#


@cached_function
def GracefulPermutations(n):
    """
    Goes through all the permutation of n>1 elements and outputs
    the list of permutations which can be used for constructing
    graceful graphs with all the vertices have degree > 0.

    EXAMPLES:

    ::

        sage: GracefulPermutations(4)
        [[0, 1, 2, 3], [0, 2, 1, 3]]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Permutations of elements from 1 to (n-1)
    P=Permutations(n-1)
    # Intialization of the list collecting the graceful permutations
    L=[]
    # Loop collecting the Graceful Permutations
    for q in P:
        # Prepending 0 at the begining of the permutation
        p=[0]+[q[i] for i in range(n-1)]
        # Initialization of our boolean variable
        bl=True
        for i in range(1,n):
            # Verifying that the graceful criteria is met
            if not (p[i] < n-i or i >= p[i]):
                bl = False
                break
        if bl:
            L.append(p)
    return L

@cached_function
def CountGracefulFunctions(n):
    """
    Goes through all the list of permutations of n > 1 elements which allow 
    for the construction of at least one graceful graphs and outputs
    the list of permutations preceded by the count for the unmber of
    associated graceful functions.

    EXAMPLES:
    ::


        sage: CountGracefulFunctions(4)
        [[2, [0, 1, 2, 3]], [2, [0, 2, 1, 3]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Loop computing the counts.
    for p in GracefulPermutations(n):
        c=1
        # start from 1 because the loop edge is attached to the 0 vertex
        for i in range(1,n):
            if p[i] < n-i and i >= p[i]:
                c=2*c
        L.append([c,p])
    return L

@cached_function
def GracefulFunctionsTuples(n):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the list of functions on the vertices derived from graceful
    permutations.

    EXAMPLES:
    ::


        sage: GracefulFunctionsTuples(3)
        [[[[(0, 0), (1, 2), (2, 0)], [(0, 0), (1, 0), (2, 0)]], [0, 1, 2]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Loop going through the GracefulPermutations.
    for p in GracefulPermutations(n):
        # Initialization of the list of functions.
        c=[]
        for j in range(n):
            if j == 0:
                if len(c)==0:
                    c.append([(j, p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, p[j])]
            # testing that only the first criteria is met
            elif (j > 0) and (p[j]<(n-j)) and not (j>=p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
            # testing that only the second criteria is met
            elif j > 0 and not (p[j]<(n-j)) and (j>=p[j]):
                if len(c)==0:
                    c.append([(j, j-p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-p[j])]
            # testing that all two criterias are met
            elif j > 0 and (p[j]<(n-j)) and (j>=p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                    c.append([(j, j-p[j])])
                else:
                    d=copy(c)
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
                    for i in range(len(d)):
                        d[i]=d[i]+[(j, j-p[j])]
                    c=c+d
        L.append([c,p])
    return L

def GracefulPolynomialList(n, x):
    """
    Goes through graceful tuple list and produce a list of
    polynomials in the variables x associated with the tree.
    

    EXAMPLES:
    ::


        sage: GracefulPolynomialList(3, x)
        [[[-2*(x - 2)*x, 0], [0, 1, 2]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GracefulFunctionsTuples(n)
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        TmpL=[]
        for lt in l[0]:
            TmpL.append(BasicLagrangeInterpolation(lt,x))
        Rslt.append([TmpL,l[1]])
    return Rslt

def Tuple2EdgeList(Lt, c):
    """
    Goes through graceful tuple list and produce a list of
    polynomials in the variables x associated with the tree.
    

    EXAMPLES:
    ::


        sage: Tuple2EdgeList([(0, 0), (1, 2), (2, 0)], 'x')
        [0, x1 - x2, -x0 + x2]    


    AUTHORS:
    - Edinah K. Gnang
    """
    # returning of the list
    return [var(c+str(t[0]))-var(c+str(t[1])) for t in Lt]

def EdgeList2Tuples(Le, c):
    """
    Goes through graceful tuple list and produce a list of
    polynomials in the variables x associated with the tree.
    

    EXAMPLES:
    ::


        sage: EdgeList2Tuples([0, x1 - x2, -x0 + x2], 'x')
        [(0, 0), (1, 2), (2, 0)]    


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the temporary list
    TmpL=[var(c+str(i)) for i in range(len(Le))]
    # Initializing the list of tuples
    Lt=[(0, 0)]
    for p in Le[1:]:
        tmp=p.operands()
        Found=False
        for i in range(len(TmpL)):
            if TmpL[i]==tmp[0]:
                for j in range(len(TmpL)):
                    if -TmpL[j]==tmp[1]:
                        Lt.append((i, j))
                break
        if not Found:
            for i in range(len(TmpL)):
                if -TmpL[i]==tmp[0]:
                    for j in range(len(TmpL)):
                        if TmpL[j]==tmp[1]:
                            Lt.append((j, i))
                    break
    # returning of the list
    return Lt

def GracefulEdgeList(n, c):
    """
    Goes through graceful tuple list and produce a list of
    polynomials in the variables x associated with the tree.
    

    EXAMPLES:
    ::


        sage: GracefulEdgeList(3, 'x')
        [[[[0, x1 - x2, -x0 + x2], [0, -x0 + x1, -x0 + x2]], [0, 1, 2]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GracefulFunctionsTuples(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        TmpL=[]
        for lt in l[0]:
            TmpL.append(Tuple2EdgeList(lt, c))
        Rslt.append([TmpL, l[1]])
    return Rslt

@cached_function
def GeneratorGracefulFunctionsTuples(n):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the list of of generator graceful functions on the vertices 
    derived from graceful permutations.

    EXAMPLES:
    ::


        sage: GeneratorGracefulFunctionsTuples(5)
        [[[[(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],
           [(0, 0), (1, 0), (2, 4), (3, 0), (4, 0)]],
          [0, 1, 2, 3, 4]],
         [[], [0, 2, 1, 3, 4]],
         [[], [0, 3, 1, 2, 4]],
         [[[(0, 0), (1, 4), (2, 0), (3, 4), (4, 0)],
           [(0, 0), (1, 4), (2, 0), (3, 2), (4, 0)]],
          [0, 3, 2, 1, 4]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Loop going through the GracefulPermutations.
    for p in GracefulPermutations(n):
        # Initialization of the list of functions.
        c=[]
        for j in range(n):
            if j == 0:
                if len(c)==0:
                    c.append([(j, p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, p[j])]
            # testing that only the first criteria is met
            elif (j > 0) and (p[j]<(n-j)) and not (j>=p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
            # testing that only the second criteria is met
            elif j > 0 and not (p[j]<(n-j)) and (j>=p[j]):
                if len(c)==0:
                    c.append([(j, j-p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-p[j])]
            # testing that all two criterias are met
            elif j > 0 and (p[j]<(n-j)) and (j>=p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                    c.append([(j, j-p[j])])
                else:
                    d=copy(c)
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
                    for i in range(len(d)):
                        d[i]=d[i]+[(j, j-p[j])]
                    c=c+d
        # Initializing the list which will store the generators
        g=[]
        for il in c:
            # initializing the counter for the degree of 0 and n-1.
            cnt0=1; cnt1=1
            # computing the degrees of the two vertices.
            # the loop starts at 1 to avoid counting the loop edge
            for ix in range(1,n-1):
                if il[ix][1]==0:
                    cnt0=cnt0+1
                elif il[ix][1]==(n-1):
                    cnt1=cnt1+1
            if min(cnt0,cnt1)>1:
                g.append(il)
        L.append([g, p])
    return L

def is_Tree(A):
    """
    Returns an boolean value determining if the input unweighted adjacency matrix
    is associated with a tree. The implementation is based on a direct implementation
    of the matrix tree theorem.

    EXAMPLES:
    ::
        sage: is_Tree(Matrix([[0, 1, 0], [1, 0, 1], [0, 1, 0]]))
        True

    AUTHORS:
    - Edinah K. Gnang
    """
    if 1==((diagonal_matrix((A*ones_matrix(A.nrows(),1)).list())-A)[0:A.nrows()-1,0:A.ncols()-1]).det():
        return True
    else:
        return False

@cached_function
def GeneratorGracefulTreeFunctionsTuples(n):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the list of of generator graceful functions on the vertices 
    derived from graceful permutations.

    EXAMPLES:
    ::


        sage: GeneratorGracefulTreeFunctionsTuples(5)
        [[[[(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],
           [(0, 0), (1, 0), (2, 4), (3, 0), (4, 0)]],
          [0, 1, 2, 3, 4]],
         [[], [0, 2, 1, 3, 4]],
         [[], [0, 3, 1, 2, 4]],
         [[[(0, 0), (1, 4), (2, 0), (3, 4), (4, 0)],
           [(0, 0), (1, 4), (2, 0), (3, 2), (4, 0)]],
          [0, 3, 2, 1, 4]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Loop going through the GracefulPermutations.
    for p in GracefulPermutations(n):
        # Initialization of the list of functions.
        c=[]
        for j in range(n):
            if j == 0:
                if len(c)==0:
                    c.append([(j, p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, p[j])]
            # testing that only the first criteria is met
            elif (j > 0) and (p[j]<(n-j)) and not (j>=p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
            # testing that only the second criteria is met
            elif j > 0 and not (p[j] < (n-j)) and (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j-p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-p[j])]
            # testing that all two criterias are met
            elif j > 0 and (p[j] < (n-j)) and (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                    c.append([(j, j-p[j])])
                else:
                    d=copy(c)
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
                    for i in range(len(d)):
                        d[i]=d[i]+[(j, j-p[j])]
                    c=c+d
        # Initializing the list which will store the generators
        g=[]
        for il in c:
            # initializing the counter for the degree of 0 and n-1.
            cnt0=1; cnt1=1
            # computing the degrees of the two vertices.
            # the loop starts at 1 to avoid counting the loop edge
            for ix in range(1,n-1):
                if il[ix][1]==0:
                    cnt0=cnt0+1
                elif il[ix][1]==(n-1):
                    cnt1=cnt1+1
            # Initilization of the adjacency matrix
            if min(cnt0,cnt1)>1 and is_Tree( sum([Id[:,t[0]]*Id[t[1],:]+Id[:,t[1]]*Id[t[0],:] for t in il]) ):
                g.append(il)
        L.append([g, p])
    return L

@cached_function
def GracefulTreeFunctionsTuples(n):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the list of of generator graceful functions on the vertices 
    derived from graceful permutations.

    EXAMPLES:
    ::


        sage: GracefulTreeFunctionsTuples(4)
        [[[[(0, 0), (1, 2), (2, 0), (3, 0)], [(0, 0), (1, 0), (2, 0), (3, 0)]],
          [0, 1, 2, 3]],
        [[[(0, 0), (1, 3), (2, 3), (3, 0)], [(0, 0), (1, 3), (2, 1), (3, 0)]],
          [0, 2, 1, 3]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Loop going through the GracefulPermutations.
    for p in GracefulPermutations(n):
        # Initialization of the list of functions.
        c=[]
        for j in range(n):
            if j == 0:
                if len(c)==0:
                    c.append([(j, p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, p[j])]
            # testing that only the first criteria is met
            elif (j > 0) and (p[j]<(n-j)) and not (j>=p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
            # testing that only the second criteria is met
            elif j > 0 and not (p[j] < (n-j)) and (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j-p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-p[j])]
            # testing that all two criterias are met
            elif j > 0 and (p[j] < (n-j)) and (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                    c.append([(j, j-p[j])])
                else:
                    d=copy(c)
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
                    for i in range(len(d)):
                        d[i]=d[i]+[(j, j-p[j])]
                    c=c+d
        # Initializing the list which will store the generators
        g=[]
        for il in c:
            # Initilization of the adjacency matrix
            if is_Tree( sum([Id[:,t[0]]*Id[t[1],:]+Id[:,t[1]]*Id[t[0],:] for t in il]) ):
                g.append(il)
        L.append([g, p])
    return L

def GracefulDiGraphAdjacencyMatrixList(n):
    """
    Goes through all the permutation of graceful second index functions
    and outputs the list of undirected adjacency matrices of graceful 
    graphs.

    EXAMPLES:
    ::

        sage: GracefulDiGraphAdjacencyMatrixList(2)
        [[[
        [1 0]
        [1 0]
        ], [0, 1]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GracefulFunctionsTuples(n)
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        TmpL=[]
        for lt in l[0]:
            TmpL.append(sum([Id[:,t[0]]*Id[t[1],:] for t in lt]))
        Rslt.append([TmpL,l[1]])
    return Rslt

def GracefulDiGraphList(n):
    """
    Goes through all the permutation of graceful second index functions
    and outputs the list of undirected adjacency matrices of graceful 
    graphs.

    EXAMPLES:
    ::


        sage: for l in GracefulDiGraphList(5):
        ....:     t=0
        ....:     for dg in l[0]:
        ....:         dg.plot().save(str(l[1]).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
        ....:         t=t+1
        ....:


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GracefulFunctionsTuples(n)
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        TmpL=[]
        for lt in l[0]:
            TmpL.append(DiGraph(sum([Id[:,t[0]]*Id[t[1],:] for t in lt])))
        Rslt.append([TmpL,l[1]])
    return Rslt

def GracefulDiTreeList(n):
    """
    Goes through all the permutation of graceful second index functions
    and outputs the list of undirected adjacency matrices of graceful 
    graphs.

    EXAMPLES:
    ::


        sage: for l in GracefulDiTreeList(5):
        ....:     t=0
        ....:     for dg in l[0]:
        ....:         dg.plot().save(str(l[1]).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
        ....:         t=t+1
        ....:


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GracefulTreeFunctionsTuples(n)
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        TmpL=[]
        for lt in l[0]:
            TmpL.append(DiGraph(sum([Id[:,t[0]]*Id[t[1],:] for t in lt])))
        Rslt.append([TmpL,l[1]])
    return Rslt

def GeneratorGracefulDiGraphList(n):
    """
    Goes through all the permutation of graceful second index functions
    and outputs the list of undirected adjacency matrices of graceful 
    graphs.

    EXAMPLES:
    ::


        sage: for l in GeneratorGracefulDiGraphList(5):
        ....:     t=0
        ....:     for dg in l[0]:
        ....:         dg.plot().save(str(l[1]).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
        ....:         t=t+1
        ....:


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GeneratorGracefulFunctionsTuples(n)
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        TmpL=[]
        for lt in l[0]:
            TmpL.append(DiGraph(sum([Id[:,t[0]]*Id[t[1],:] for t in lt])))
        Rslt.append([TmpL,l[1]])
    return Rslt

def GeneratorGracefulTreeDiGraphList(n):
    """
    Goes through all the permutation of graceful second index functions
    and outputs the list of undirected adjacency matrices of graceful 
    graphs.

    EXAMPLES:
    ::


        sage: for l in GeneratorGracefulTreeDiGraphList(5):
        ....:     t=0
        ....:     for dg in l[0]:
        ....:         dg.plot().save(str(l[1]).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
        ....:         t=t+1
        ....:


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GeneratorGracefulTreeFunctionsTuples(n)
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        TmpL=[]
        for lt in l[0]:
            TmpL.append(DiGraph(sum([Id[:,t[0]]*Id[t[1],:] for t in lt])))
        Rslt.append([TmpL,l[1]])
    return Rslt

def BasicLagrangeInterpolation(L, x):
    """
    Implements the basic lagrange interpolation.
    The functions take as input a list of tuples
    and outputs a polynomial in the variable x


    EXAMPLES:
    ::


        sage: x=var('x')
        sage: BasicLagrangeInterpolation([(0,0), (1,1), (2,2), (3,3)], x)
        1/2*(x - 1)*(x - 2)*x - (x - 1)*(x - 3)*x + 1/2*(x - 2)*(x - 3)*x


    AUTHORS:
    - Edinah K. Gnang
    """
    # L is a list of tuples
    # Initialized the lenght of the list
    n = len(L)
    # Initialization of the function
    f = 0
    # Code for building the parts 
    for idx in range(len(L)):
        fk = 1
        for j in [i for i in range(len(L)) if i != idx]:
            fk = fk*((x-L[j][0])/(L[idx][0]-L[j][0]))
        f = f + L[idx][1]*fk
    return f

def GeneratorGracefulPolynomialList(n, x):
    """
    Goes through generator graceful tuple list and produce a list of
    polynomials.
    graphs.

    EXAMPLES:
    ::


        sage: GeneratorGracefulPolynomialList(5, x)
        [[[(x - 1)*(x - 3)*(x - 4)*x - 1/3*(x - 2)*(x - 3)*(x - 4)*x,
           (x - 1)*(x - 3)*(x - 4)*x],
          [0, 1, 2, 3, 4]],
         [[], [0, 2, 1, 3, 4]],
         [[], [0, 3, 1, 2, 4]],
         [[-2/3*(x - 1)*(x - 2)*(x - 4)*x - 2/3*(x - 2)*(x - 3)*(x - 4)*x,
           -1/3*(x - 1)*(x - 2)*(x - 4)*x - 2/3*(x - 2)*(x - 3)*(x - 4)*x],
          [0, 3, 2, 1, 4]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GeneratorGracefulFunctionsTuples(n)
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        TmpL=[]
        for lt in l[0]:
            TmpL.append(BasicLagrangeInterpolation(lt,x))
        Rslt.append([TmpL,l[1]])
    return Rslt

def EdgeListDegreeSequence(Le):
    """
    The method returns the degree sequence of a an input edge list.
    the inputs to this function is list of tuples specifying an
    edge list and a positive integer corresponding to the number
    of vertices in the graph. The loop edges is ingnored by this
    implementation.

    EXAMPLES:

    ::

        sage: EdgeListDegreeSequence([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)])
        [2, 1, 2, 1, 2]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list which will store the degree sequence
    L = [0 for v in range(len(Le))]
    # Main loop computing the degree sequence
    for t in Le:
        if t[0] != t[1]:
            L[t[0]]=L[t[0]]+1; L[t[1]]=L[t[1]]+1
    return L

def LeftChild(T):
    """
    The method returns list description of the left child by the first
    transformation rule.
    

    EXAMPLES:

    ::

        sage: LeftChild([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)])
        [(0, 0), (1, 2), (2, 4), (3, 0), (4, 0), (5, 0)]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    return T+[(len(T), 0)]

def MidChildrenI(T):
    """
    The method returns list description of the middle children by the second
    transformation rule. The implementation here splits the vertex with the
    smallest integer label

    EXAMPLES:

    ::

        sage: MidChildrenI([(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)])
        [[(0, 0), (1, 5), (2, 0), (3, 0), (4, 5), (5, 0)],
         [(0, 0), (1, 0), (2, 5), (3, 5), (4, 0), (5, 0)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the list and the vertex set
    L=[]; S=Set([t[0] for t in T[1:] if t[1] == 0])
    if S.cardinality() >= 2:
        for s in S.subsets():
            if (s.cardinality() >= 1) and (s.cardinality() < S.cardinality()):
                # Initializing the boolean variable asserting the splitting condition
                splittable = True
                for i in s:
                    if not (len(T)-i) in s:
                        splittable=False; break
                if splittable:
                    Tmp=[]
                    for t in T:
                        if t[0] in s and t[1] == 0:
                            Tmp.append((t[0], len(T)))
                        else:
                            Tmp.append(t)
                    Tmp.append((len(T), 0))
                    L.append(Tmp)
    return L

def MidChildrenII(T):
    """
    The method returns list description of the middle children by the second
    transformation rule. This implementation splits the vertex with the 
    largest integer label.

    EXAMPLES:

    ::

        sage: MidChildrenII([(0, 0), (1, 4), (2, 4), (3, 4), (4, 0)])
        [[(0, 0), (2, 5), (3, 5), (4, 0), (1, 0), (5, 0)],
         [(0, 0), (2, 0), (3, 0), (4, 5), (5, 1), (5, 0)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the list and the vertex set
    L=[]; S=Set([t[0]+1 for t in T if t[1] == len(T)-1]+[1])
    if S.cardinality() >= 2:
        for s in S.subsets():
            if (s.cardinality() >= 1) and (s.cardinality() < S.cardinality()):
                # Initializing the boolean variable asserting the splitting condition
                splittable = True
                for i in s:
                    if not (len(T)-i) in s:
                        splittable=False; break
                if splittable:
                    Tmp=[(0, 0)]
                    for t in T[1:]:
                        if t[0]+1 in s and t[1] == len(T)-1:
                            Tmp.append((t[0]+1, 0))
                        elif t[0] == len(T)-1 and t[1]==0 and 1 in s:
                            Tmp.append((1,0))
                        else:
                            if t[1]==0:
                                Tmp.append((t[1]+1, t[0]+1))
                            else:
                                Tmp.append((t[0]+1, t[1]+1))
                    Tmp.append((len(T), 0)); Tmp.sort()
                    L.append(copy(Tmp))
    return L

def RightChild(T):
    """
    The method returns list description of the left child by the second 
    transformation rule.

    EXAMPLES:

    ::

        sage: RightChild([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)])
        [(0, 0), (1, 5), (2, 3), (3, 5), (4, 1), (5, 0)]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    return [(0, 0), (1, len(T))] + [(T[i][0]+1, T[i][1]+1) for i in range(1,len(T)-1)] + [(len(T), 0)]

def BinaryTreeArrayList(T, DegSeq):
    """
    The method returns list description of an array implementation
    of a binary tree.

    EXAMPLES:

    ::

        sage: L=BinaryTreeArrayList([(0, 0), (1, 0)], [3, 1, 1, 1]); L
        [[(0, 0), (1, 0)],
         [(0, 0), (1, 0), (2, 0)],
         [(0, 0), (1, 2), (2, 0)],
         [(0, 0), (1, 0), (2, 0), (3, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0)],
         [(0, 0), (1, 2), (2, 0), (3, 0)],
         [(0, 0), (1, 3), (2, 3), (3, 0)],
         [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 1), (3, 1), (4, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 4), (3, 2), (4, 0)],
         [(0, 0), (1, 2), (2, 0), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 3), (3, 1), (4, 0)],
         [(0, 0), (1, 3), (2, 3), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 4), (3, 4), (4, 0)]]
        sage: for i in range(len(L)):
        ....:     T2DiGraph(L[i]).plot().save(str(i)+'.png')
        ....:
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list implementation
    # of binary tree storing the graceful graphs
    # whose degree sequence are lexicographically 
    # smaller then that of the second input tree T2.
    L=[T]
    # Initialization of the indexing variable
    indx=0
    while (sorted(EdgeListDegreeSequence(L[indx])))[::-1] <= (sorted(DegSeq))[::-1]:
        L.append(LeftChild(L[indx])); L.append(RightChild(L[indx]))
        indx=indx+1
    return L

def NodeSplitingList(T, DegSeq):
    """
    The method returns list description of an array implementation
    of a binary tree.

    EXAMPLES:

    ::

        sage: L=NodeSplitingList([(0, 0), (1, 0)], [3, 1, 1, 1]); L
        [[(0, 0), (1, 0)],
         [(0, 0), (1, 0), (2, 0)],
         [(0, 0), (1, 2), (2, 0)],
         [(0, 0), (1, 0), (2, 0), (3, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0)],
         [(0, 0), (1, 2), (2, 0), (3, 0)],
         [(0, 0), (1, 3), (2, 3), (3, 0)],
         [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)],
         [(0, 0), (1, 0), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 0), (3, 4), (4, 0)],
         [(0, 0), (1, 4), (2, 1), (3, 1), (4, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0), (4, 0)],
         [(0, 0), (2, 0), (3, 2), (4, 1), (4, 0)],
         [(0, 0), (1, 4), (2, 4), (3, 2), (4, 0)],
         [(0, 0), (1, 2), (2, 0), (3, 0), (4, 0)],
         [(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 3), (3, 1), (4, 0)],
         [(0, 0), (1, 3), (2, 3), (3, 0), (4, 0)],
         [(0, 0), (2, 0), (3, 4), (4, 1), (4, 0)],
         [(0, 0), (2, 4), (3, 0), (1, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 4), (3, 4), (4, 0)]]        
        sage: for i in range(len(L)):
        ....:     T2DiGraph(L[i]).plot().save(str(i)+'.png')
        ....:
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list implementation
    # of binary tree storing the graceful graphs
    # whose degree sequence are lexicographically 
    # smaller then that of the second input tree T2.
    L=[T]
    # Initialization of the indexing variable
    indx=0
    while (sorted(EdgeListDegreeSequence(L[indx])))[::-1] <= (sorted(DegSeq))[::-1]:
        L.append(LeftChild(L[indx]))
        TmpL=MidChildrenI(L[indx])
        for t in TmpL:
            if not t in L:
                L.append(t) 
        TmpL=MidChildrenII(L[indx])
        for t in TmpL:
            if not t in L:
                L.append(t) 
        L.append(RightChild(L[indx]))
        indx=indx+1
    return L

def NodeSplitingGeneratorList(T, DegSeq):
    """
    The method returns list description of an array implementation
    of a binary tree.

    EXAMPLES:

    ::

        sage: L=NodeSplitingGeneratorList([(0, 0), (1, 0)], [3, 1, 1, 1]); L
        [[(0, 0), (1, 0), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 0), (3, 4), (4, 0)],
         [(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (2, 0), (3, 4), (4, 1), (4, 0)],
         [(0, 0), (2, 4), (3, 0), (1, 0), (4, 0)]]
        sage: for i in range(len(L)):
        ....:     T2DiGraph(L[i]).plot().save(str(i)+'.png')
        ....:
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list implementation
    # of binary tree storing the graceful graphs
    # whose degree sequence are lexicographically 
    # smaller then that of the second input tree T2.
    L=NodeSplitingList(T, DegSeq)
    # Initializing the list which only stores generator trees.
    GnrtrLst=[]
    for t in L:
        # initializing the counter for the degree of 0 and n-1.
        cnt0=1; cnt1=1
        # computing the degrees of the two vertices.
        # the loop starts at 1 to avoid counting the loop edge
        for ix in range(1,len(t)-1):
            if t[ix][1]==0:
                cnt0=cnt0+1
            elif t[ix][1]==(len(t)-1):
                cnt1=cnt1+1
        # Initilization of the adjacency matrix
        if min(cnt0,cnt1) > 1: 
                GnrtrLst.append(t)
    return GnrtrLst

def T2DiGraph(T):
    """
    The method returns a directed graph object associated with 
    with the tuple list description of the directed graph

    EXAMPLES:

    ::

        sage: T2DiGraph([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)]).degree_sequence()
        [4, 2, 2, 1, 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(len(T))
    return DiGraph(sum([Id[:,t[0]]*Id[t[1],:] for t in T]))
