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
    Goes through all permutations of n > 1 elements and outputs
    the list of permutations which can be used to construct graceful
    graphs having no isolated vertices.

    EXAMPLES:

    ::

        sage: GracefulPermutations(4)
        [[0, 1, 2, 3], [0, 2, 1, 3]]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations of elements from 1 to (n-1).
    P=Permutations(n-1)
    # Initialization of the list of graceful permutations.
    L=[]
    # Loop collecting graceful permutations.
    for q in P:
        # Prepending 0 to the permutation.
        p=[0]+[q[i] for i in range(n-1)]
        # Initialization of the boolean variable.
        bl=True
        for i in range(1,n):
            # Checking the gracefulness condition.
            if not ( p[i] < n-i or p[i] <= i ):
                bl=False
                break
        if bl:
            L.append(p)
    return L

@cached_function
def CountGracefulFunctions(n):
    """
    Goes through the list of permutations of n > 1 elements
    which can be used to construct of at least one graceful
    graph having no isolated vertices. The method outputs
    the list of permutations preceded by the counts of 
    associated graceful graphs.

    EXAMPLES:

    ::

        sage: CountGracefulFunctions(4)
        [[2, [0, 1, 2, 3]], [2, [0, 2, 1, 3]]]
        sage: TpL=CountGracefulFunctions(5)
        sage: for l in TpL:
        ....:     L.append(l[0])
        ....:
        sage: L
        [4, 2, 2, 4]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list collecting graceful permutations.
    L=[]
    # Loop computing the counts.
    for p in GracefulPermutations(n):
        # Initialization of the count variable.
        c=1
        # Start from 1 because the loop edge is attached to 0.
        for i in range(1, n):
            if p[i] < n-i and p[i] <= i:
                c=2*c
        L.append([c, p])
    return L

def CountGracefulTrees(n):
    """
    Goes through all the permutation of n > 1 elements and enumerates 
    graceful trees on the vertices derived from graceful permutations.


    EXAMPLES:

    ::


        sage: CountGracefulTrees(4)
        [[2, [0, 1, 2, 3]], [2, [0, 2, 1, 3]]]
        sage: TpL=CountGracefulTrees(5)
        sage: for l in TpL:
        ....:     L.append(l[0])
        ....:
        sage: L
        [4, 2, 2, 4]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the list collecting graceful permutations.
    L=[]
    # Loop going through the graceful permutations.
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
            elif (j > 0) and (p[j] < n-j) and not (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
            # testing that only the second criteria is met
            elif (j > 0) and not (p[j] < n-j) and (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j-p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-p[j])]
            # testing that all two criterias are met
            elif (j > 0) and (p[j] < n-j) and (j >= p[j]):
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
        g = 0
        for il in c:
            # Initilization of the adjacency matrix
            if is_Tree( sum([Id[:,t[0]]*Id[t[1],:]+Id[:,t[1]]*Id[t[0],:] for t in il]) ):
                g = g + 1
        L.append([g, p])
    return L

@cached_function
def SignedPermutations(n):
    """
    Goes through all permutations of n > 1 elements and outputs
    the list of signed permutations associated with graceful graphs.

    EXAMPLES:

    ::

        sage: L=SignedPermutations(4)
        sage: L
        [[0, 1, -2, -3], [0, -1, -2, -3], [0, 2, 1, -3], [0, 2, -1, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of graceful permutations.
    L=[]
    # Loop going through graceful permutations.
    for p in GracefulPermutations(n):
        # Initialization of the list of functions.
        c=[]
        for j in range(n):
            if j == 0:
                if len(c) == 0:
                    c.append([p[j]])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[p[j]]
            # Testing that only the first gracefulness criteria is met.
            elif (j > 0) and (p[j] < (n-j)) and not (j >= p[j]):
                if len(c) == 0:
                    c.append([p[j]])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[p[j]]
            # Testing that only the second gracefulness criteria is met.
            elif (j > 0) and not (p[j] < (n-j)) and (j >= p[j]):
                if len(c) == 0:
                    c.append([-p[j]])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[-p[j]]
            # Testing that both gracefulness criterias are met.
            elif (j > 0) and (p[j] < (n-j)) and (j >= p[j]):
                if len(c) == 0:
                    c.append([p[j]])
                    c.append([-p[j]])
                else:
                    d=copy(c)
                    for i in range(len(c)):
                        c[i]=c[i]+[p[j]]
                    for i in range(len(d)):
                        d[i]=d[i]+[-p[j]]
                    c=c+d
        L.append([c,p])
    RsltL=[]
    for l in L:
        RsltL=RsltL+l[0]
    return RsltL

def LeftSignedPermutation(sigma):
    """
    The method returns a list description of the signed permutation
    associated with the left child. The method splits the node with
    the smallest label to create a leaf edge.

    EXAMPLES:

    ::

        sage: sigma=[0, 1, -2]; LeftSignedPermutation(sigma)
        [0, 1, -2, -3]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    return sigma+[-len(sigma)]

def MiddleSignedPermutationsI(sigma):
    """
    The method returns a list description of the first family of middle 
    children. The method splits the vertex with the smallest label to
    create a non leaf edge.

    EXAMPLES:

    ::

        sage: MiddleSignedPermutationsI([0, -1, -2, -3, -4])
        [[0, 4, -2, -3, 1, -5], [0, -1, 3, 2, -4, -5]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the list and the vertex set
    L=[]; S=Set([i for i in range(1,len(sigma)) if i+sigma[i]==0])
    # Going through all the subsets of vertices adjacent to 0
    if S.cardinality() >= 2:
        for s in S.subsets():
            if (s.cardinality() >= 1) and (s.cardinality() < S.cardinality()):
                # Initializing the boolean variable asserting the splitting condition
                splittable=True
                # for (i-0) in s
                for i in s:
                    if not (len(sigma)-i) in s:
                        splittable=False; break
                if splittable:
                    Tmp=[]
                    for t in range(len(sigma)):
                        #if t in s and t+sigma[t] == 0:
                        if t in s:
                            Tmp.append(sigma[t]+len(sigma))
                        else:
                            Tmp.append(sigma[t])
                    Tmp.append(-len(sigma))
                    L.append(Tmp)
    return L

def MiddleSignedPermutationsII(sigma):
    """
    The method returns list description of the middle children by the second
    transformation rule. This implementation splits the vertex with the 
    largest integer label.

    EXAMPLES:

    ::

        sage: MiddleSignedPermutationsII([0, 3, 2, -1, -4])
        [[0, 4, -2, -3, -1, -5]]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the list and the vertex set
    L=[]; S=Set([0]+[i for i in range(1,len(sigma)) if i+sigma[i]==len(sigma)-1])
    # Going through all the subsets of vertices adjacent to len(sigma)-1
    if S.cardinality() >= 2:
        for s in S.subsets():
            if (s.cardinality() >= 1) and (s.cardinality() < S.cardinality()):
                # Initializing the boolean variable asserting the splitting condition
                splittable=True
                # for (i-0) in s
                for i in s:
                    if not (len(sigma)-(i+1)) in [v+1 for v in s]:
                        splittable=False; break
                if splittable:
                    Tmp=[0]
                    for t in range(1,len(sigma)):
                        if t in [v+1 for v in s]:
                            Tmp.append(-t)
                        else:
                            if t-1==0:
                                Tmp.append(len(sigma)-1)
                            else:
                                Tmp.append(sigma[t-1])
                    Tmp.append(-len(sigma))
                    L.append(Tmp)
    return L

def RightSignedPermutation(sigma):
    """
    The method returns list description of the signed permutation
    associated with the right child by the first transformation rule.

    EXAMPLES:

    ::

        sage: sigma=[0, 1, -2]; RightSignedPermutation(sigma)
        [0, 2, 1, -3]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    return [0,len(sigma)-1]+sigma[1:-1]+[-len(sigma)]

def NodeSplittingSignedPermutationList(T, nbv):
    """
    The method returns list of signed permutations associated with
    gracefully labeled trees obtained via node splitting for trees
    having up to nbv vertices.
    Note that there is a symmetry relating all the descendants
    of [0,-1,-2] and the descendants of [0,1,-2]. 


    EXAMPLES:

    ::

        sage: L=NodeSplittingSignedPermutationList([0, -1, -2], 3); L
        [[0, -1, -2], [0, -1, -2, -3], [0, 2, -1, -3]]
       

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
    while len(L[indx]) <= nbv:
        L.append(LeftSignedPermutation(L[indx]))
        TmpL=MiddleSignedPermutationsI(L[indx])
        for t in TmpL:
            if not t in L:
                L.append(t) 
        TmpL=MiddleSignedPermutationsII(L[indx])
        for t in TmpL:
            if not t in L:
                L.append(t) 
        L.append(RightSignedPermutation(L[indx]))
        indx=indx+1
    return L

def GraphSignedPermutationOrbitsT(T):
    """
    Obtain the orbits of signed function defined by
    the isomorphism classes of associated trees. The tree
    isomorphism function is doing the heavy lifting here.
    The functions takes as input a tuple edge list.


    EXAMPLES:

    ::


        sage: GraphSignedPermutationOrbitsT([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)])
        [[(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 3), (3, 1), (4, 0)],
         [(0, 0), (1, 4), (2, 0), (3, 2), (4, 0)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=SignedPermutations(len(T))
    # Initialization of the list storing the equivalence class of trees.
    cL=[]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        if SP2Graph(s).is_isomorphic(T2Graph(T)) and not s in cL:
            cL.append(SP2T(s))
    return cL

def GraphSignedPermutationOrbitsTII(T):
    """
    Obtain the orbits of signed function defined by
    the isomorphism classes of associated directed graph.
    The tree isomorphism function is doing the heavy 
    lifting here. The functions takes as input a tuple
    edge list. The difference with the implementation
    above is the fact that the graph is treated as a
    directed graph when checking for isomorphism


    EXAMPLES:

    ::


        sage: GraphSignedPermutationOrbitsTII([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)])
        [[(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0), (4, 0)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=SignedPermutations(len(T))
    # Initialization of the list storing the equivalence class of trees.
    cL=[]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        if T2DiGraph(SP2T(s)).is_isomorphic(T2DiGraph(T)) and not s in cL:
            cL.append(SP2T(s))
    return cL

def GraphSignedPermutationOrbitsSP(sigma):
    """
    Obtain the orbits of signed function defined by
    the isomorphism classes of associated trees. The tree
    isomorphism function is doing the heavy lifting here.
    The functions takes as input a signed permutation.


    EXAMPLES:

    ::


        sage: GraphSignedPermutationOrbitsSP([0, -1, -2, -3])
        [[0, -1, -2, -3], [0, 2, 1, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=SignedPermutations(len(sigma))
    # Initialization of the list storing the equivalence class of trees.
    cL=[]
    # Loop perfomring the isomorphism binning.
    for s in L:
        if SP2Graph(s).is_isomorphic(SP2Graph(sigma)) and not s in cL:
            cL.append(s)
    return cL

def GraphSignedPermutationOrbitsSPII(sigma):
    """
    Obtain the orbits of signed function defined by
    the isomorphism classes of associated drected graph.
    The tree isomorphism function is doing the heavy 
    lifting here.The functions takes as input a signed
    permutation. The difference with the previous 
    implementation is that the graphs here are directed
    and rooted


    EXAMPLES:

    ::


        sage: GraphSignedPermutationOrbitsSPII([0, -1, -2, -3])
        [[0, -1, -2, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=SignedPermutations(len(sigma))
    # Initialization of the list storing the equivalence class of trees.
    cL=[]
    # Loop perfomring the isomorphism binning.
    for s in L:
        if T2DiGraph(SP2T(s)).is_isomorphic(T2DiGraph(SP2T(sigma))) and not s in cL:
            cL.append(s)
    return cL

@cached_function
def GracefulFunctionsTuples(n):
    """
    Goes through all permutations of n > 1 elements and outputs the
    list of directed edges as pairs of vertices derived from graceful
    permutations.

    EXAMPLES:

    ::


        sage: GracefulFunctionsTuples(3)
        [[[[(0, 0), (1, 2), (2, 0)], [(0, 0), (1, 0), (2, 0)]], [0, 1, 2]]]
        sage: L=[]
        sage: TpL=GracefulFunctionsTuples(4)
        sage: for l in TpL:
        ....:     L=L+l[0]
        ....:
        sage: L
        [[(0, 0), (1, 2), (2, 0), (3, 0)],
         [(0, 0), (1, 0), (2, 0), (3, 0)],
         [(0, 0), (1, 3), (2, 3), (3, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Loop going through the list of graceful permutations.
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
            # Testing that only the first criteria is met.
            elif (j > 0) and (p[j] < n-j) and not (j >= p[j]):
                if len(c) == 0:
                    c.append([(j, j+p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
            # Testing that only the second criteria is met.
            elif (j > 0) and not (p[j] < n-j) and (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j-p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-p[j])]
            # Testing that all two criterias are met
            elif (j > 0) and (p[j] < n-j) and (j >= p[j]):
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

@cached_function
def HiddenGracefulFunctionsTuples(n):
    """
    Goes through all permutations of n > 1 elements and outputs the
    list of directed edges as pairs of vertices derived from graceful
    permutations. The graceful labeling is hidden by the action of
    elements of the symmetric group on the vertices.

    EXAMPLES:

    ::


        sage: HiddenGracefulFunctionsTuples(3)
        [[[[(0, 0), (1, 2), (2, 0)], [(0, 0), (1, 0), (2, 0)]], [0, 1, 2], [0, 1, 2]],
         [[[(0, 0), (2, 1), (1, 0)], [(0, 0), (2, 0), (1, 0)]], [0, 1, 2], [0, 2, 1]],
         [[[(1, 1), (0, 2), (2, 1)], [(1, 1), (0, 1), (2, 1)]], [0, 1, 2], [1, 0, 2]],
         [[[(1, 1), (2, 0), (0, 1)], [(1, 1), (2, 1), (0, 1)]], [0, 1, 2], [1, 2, 0]],
         [[[(2, 2), (0, 1), (1, 2)], [(2, 2), (0, 2), (1, 2)]], [0, 1, 2], [2, 0, 1]],
         [[[(2, 2), (1, 0), (0, 2)], [(2, 2), (1, 2), (0, 2)]], [0, 1, 2], [2, 1, 0]]]
        sage: L=[]
        sage: TpL=HiddenGracefulFunctionsTuples(3)
        sage: for l in TpL:
        ....:     L=L+l[0]
        ....:
        sage: L
        [[(0, 0), (1, 0), (2, 0)],
         [(0, 0), (1, 2), (2, 0)],
         [(0, 0), (2, 0), (1, 0)],
         [(0, 0), (2, 1), (1, 0)],
         [(1, 1), (0, 1), (2, 1)],
         [(1, 1), (0, 2), (2, 1)],
         [(1, 1), (2, 0), (0, 1)],
         [(1, 1), (2, 1), (0, 1)],
         [(2, 2), (0, 1), (1, 2)],
         [(2, 2), (0, 2), (1, 2)],
         [(2, 2), (1, 0), (0, 2)],
         [(2, 2), (1, 2), (0, 2)]]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the list of graceful permutations
    GpL=GracefulPermutations(n)
    # Intialization of the list collecting the graceful permutations.
    L=[]
    for qm in Permutations(n):
        # fixing up the permutation
        q=[qm[i]-1 for i in range(n)]
        # Loop going through the list of graceful permutations.
        for p in GpL:
            # Initialization of the list of functions.
            c=[]
            for j in range(n):
                if j == 0:
                    if len(c)==0:
                        c.append([(q[j], q[p[j]])])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[(q[j], q[p[j]])]
                # Testing that only the first criteria is met.
                elif (j > 0) and (p[j] < n-j) and not (j >= p[j]):
                    if len(c) == 0:
                        c.append([(q[j], q[j+p[j]])])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[(q[j], q[j+p[j]])]
                # Testing that only the second criteria is met.
                elif (j > 0) and not (p[j] < n-j) and (j >= p[j]):
                    if len(c)==0:
                        c.append([(q[j], q[j-p[j]])])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[(q[j], q[j-p[j]])]
                # Testing that all two criterias are met
                elif (j > 0) and (p[j] < n-j) and (j >= p[j]):
                    if len(c)==0:
                        c.append([(q[j], q[j+p[j]])])
                        c.append([(q[j], q[j-p[j]])])
                    else:
                        d=copy(c)
                        for i in range(len(c)):
                            c[i]=c[i]+[(q[j], q[j+p[j]])]
                        for i in range(len(d)):
                            d[i]=d[i]+[(q[j], q[j-p[j]])]
                        c=c+d
            L.append([c,p,q])
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
    Goes through graceful tuple list and produces a list of
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

def Tuple2EdgeListII(Lt, c):
    """
    Goes through graceful tuple list and produces a list of
    polynomials in the variables x associated with the tree.
    

    EXAMPLES:
    ::


        sage: Tuple2EdgeList([(0, 0), (1, 2), (2, 0)], 'x')
        [1, x1/x2, x2/x0]    


    AUTHORS:
    - Edinah K. Gnang
    """
    # returning of the list
    return [1]+[var(c+str(t[0]))*var(c+str(t[1])) for t in Lt[1:]]


def SP2T(sp):
    """
    Returns the tuple encoding of the input signed permutation


    EXAMPLES:

    ::

        sage: sp=[0, 1, 2, -3, -4]
        sage: SP2T(sp)
        [(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)]


    AUTHORS:
    - Edinah K. Gnang
    """
    return [(i,i+sp[i]) for i in range(len(sp))]

def T2SP(T):
    """
    Returns the tuple encoding of the input signed permutation


    EXAMPLES:

    ::

        sage: T=[(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)]
        sage: T2SP(T)
        [0, 1, 2, -3, -4]


    AUTHORS:
    - Edinah K. Gnang
    """
    return [t[1]-t[0] for t in T]

def EdgeList2Tuples(Le, c):
    """
    Goes through graceful tuple list and produce a list of
    polynomials in the variables x associated with the tree.


    EXAMPLES:

    ::

        sage: x0, x1, x2 = var('x0, x1, x2')
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
    edge list expressed as subtraction of vertex variables.
    

    EXAMPLES:
    ::


        sage: GracefulEdgeList(3, 'x')
        [[[[0, x1 - x2, -x0 + x2], [0, -x0 + x1, -x0 + x2]], [0, 1, 2]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of graceful directed edge list.
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
        sage: L=[]
        sage: TpL=GeneratorGracefulFunctionsTuples(5)
        sage: for l in TpL:
        ....:     L=L+l[0]
        ....:
        sage: L
        [[(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (1, 0), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 0), (3, 4), (4, 0)],
         [(0, 0), (1, 4), (2, 0), (3, 2), (4, 0)]]


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
            # Testing that only the first criteria is met
            elif (j > 0) and (p[j]<(n-j)) and not (j>=p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
            # Testing that only the second criteria is met
            elif (j > 0) and not (p[j]<(n-j)) and (j>=p[j]):
                if len(c)==0:
                    c.append([(j, j-p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-p[j])]
            # Testing that both criterias are met
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
            # Initializing the counter for the degree of 0 and n-1.
            cnt0=1; cnt1=1
            # Computing the degrees of the two vertices.
            # The loop starts at 1 to avoid counting the loop edge
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
    Returns an boolean value determining if the input unweighted
    adjacency matrix is associated with a tree. The implementation
    is based on a implementation of the matrix tree theorem.


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
def GeneratorGracefulNonTreeFunctionsTuples(n):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the list of fuction associated with generator graceful non trees
    on the vertices derived from graceful permutations.


    EXAMPLES:

    ::


        sage: GeneratorGracefulNonTreeFunctionsTuples(7)
        [[[], [0, 1, 2, 3, 4, 5, 6]],
         [[], [0, 1, 3, 2, 4, 5, 6]],
         [[], [0, 1, 4, 2, 3, 5, 6]],
         [[], [0, 1, 4, 3, 2, 5, 6]],
         [[], [0, 2, 1, 3, 4, 5, 6]],
         [[], [0, 2, 3, 1, 4, 5, 6]],
         [[[(0, 0), (1, 3), (2, 6), (3, 4), (4, 1), (5, 0), (6, 0)]],
          [0, 2, 4, 1, 3, 5, 6]],
         [[], [0, 2, 4, 3, 1, 5, 6]],
         [[], [0, 3, 1, 2, 4, 5, 6]],
         [[], [0, 3, 2, 1, 4, 5, 6]],
         [[], [0, 3, 4, 1, 2, 5, 6]],
         [[[(0, 0), (1, 4), (2, 6), (3, 1), (4, 3), (5, 0), (6, 0)]],
          [0, 3, 4, 2, 1, 5, 6]],
         [[], [0, 4, 1, 2, 3, 5, 6]],
         [[], [0, 4, 1, 3, 2, 5, 6]],
         [[], [0, 4, 2, 1, 3, 5, 6]],
         [[], [0, 4, 2, 3, 1, 5, 6]],
         [[], [0, 4, 3, 1, 2, 5, 6]],
         [[], [0, 4, 3, 2, 1, 5, 6]],
         [[], [0, 5, 1, 2, 3, 4, 6]],
         [[[(0, 0), (1, 6), (2, 3), (3, 5), (4, 0), (5, 2), (6, 0)]],
          [0, 5, 1, 2, 4, 3, 6]],
         [[], [0, 5, 1, 3, 2, 4, 6]],
         [[], [0, 5, 1, 3, 4, 2, 6]],
         [[], [0, 5, 2, 1, 3, 4, 6]],
         [[], [0, 5, 2, 1, 4, 3, 6]],
         [[], [0, 5, 2, 3, 1, 4, 6]],
         [[], [0, 5, 2, 3, 4, 1, 6]],
         [[], [0, 5, 3, 1, 2, 4, 6]],
         [[[(0, 0), (1, 6), (2, 5), (3, 2), (4, 0), (5, 3), (6, 0)]],
          [0, 5, 3, 1, 4, 2, 6]],
         [[], [0, 5, 3, 2, 1, 4, 6]],
         [[], [0, 5, 3, 2, 4, 1, 6]],
         [[], [0, 5, 4, 1, 2, 3, 6]],
         [[], [0, 5, 4, 1, 3, 2, 6]],
         [[], [0, 5, 4, 2, 1, 3, 6]],
         [[], [0, 5, 4, 2, 3, 1, 6]],
         [[], [0, 5, 4, 3, 1, 2, 6]],
         [[], [0, 5, 4, 3, 2, 1, 6]]]
        sage: L=[]
        sage: TpL=GeneratorGracefulNonTreeFunctionsTuples(7)
        sage: for l in TpL:
        ....:     L=L+l[0]
        ....:
        sage: L
        [[(0, 0), (1, 3), (2, 6), (3, 4), (4, 1), (5, 0), (6, 0)],
         [(0, 0), (1, 4), (2, 6), (3, 1), (4, 3), (5, 0), (6, 0)],
         [(0, 0), (1, 6), (2, 3), (3, 5), (4, 0), (5, 2), (6, 0)],
         [(0, 0), (1, 6), (2, 5), (3, 2), (4, 0), (5, 3), (6, 0)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the list collecting graceful permutations.
    L=[]
    # Loop going through the graceful permutations.
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
            elif (j > 0) and (p[j] < n-j) and not (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
            # testing that only the second criteria is met
            elif (j > 0) and not (p[j] < n-j) and (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j-p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-p[j])]
            # testing that all two criterias are met
            elif (j > 0) and (p[j] < n-j) and (j >= p[j]):
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
            if min(cnt0,cnt1)>1 and not is_Tree( sum([Id[:,t[0]]*Id[t[1],:]+Id[:,t[1]]*Id[t[0],:] for t in il]) ):
                g.append(il)
        L.append([g, p])
    return L

@cached_function
def GeneratorGracefulTreeFunctionsTuples(n):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the list of fuction associated with generator graceful trees on
    the vertices derived from graceful permutations.


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
        sage: L=[]
        sage: TpL=GeneratorGracefulTreeFunctionsTuples(5)
        sage: for l in TpL:
        ....:     L=L+l[0]
        ....:
        sage: L
        [[(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (1, 0), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 0), (3, 4), (4, 0)],
         [(0, 0), (1, 4), (2, 0), (3, 2), (4, 0)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the list collecting graceful permutations.
    L=[]
    # Loop going through the graceful permutations.
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
            elif (j > 0) and (p[j] < n-j) and not (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
            # testing that only the second criteria is met
            elif (j > 0) and not (p[j] < n-j) and (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j-p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-p[j])]
            # testing that all two criterias are met
            elif (j > 0) and (p[j] < n-j) and (j >= p[j]):
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
def TreeSignedPermutations(n):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the list of of signed permutation associated with trees .


    EXAMPLES:

    ::


        sage: TreeSignedPermutations(4)
        [[0, 1, -2, -3], [0, -1, -2, -3], [0, 2, 1, -3], [0, 2, -1, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of function Tulples
    TpL=GracefulTreeFunctionsTuples(n)
    L=[]
    for l in TpL:
        L=L+l[0]
    # Initialization of the Result list
    RsltL=[]
    for l in L:
        RsltL.append([t[1]-t[0] for t in l])
    return RsltL

@cached_function
def NonTreeSignedPermutations(n):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the list of of signed permutation associated with non trees .


    EXAMPLES:

    ::


        sage: NonTreeSignedPermutations(7)
        [[0, 2, 4, 1, -3, -5, -6],
         [0, 3, 4, -2, -1, -5, -6],
         [0, 5, 1, 2, -4, -3, -6],
         [0, 5, 3, -1, -4, -2, -6]]        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of function Tulples
    TpL=GracefulNonTreeFunctionsTuples(n)
    L=[]
    for l in TpL:
        L=L+l[0]
    # Initialization of the Result list
    RsltL=[]
    for l in L:
        RsltL.append([t[1]-t[0] for t in l])
    return RsltL

def GraphSignedPermutationClasses(n):
    """
    Obtain the equivalence classes of signed function defined by
    the isomorphism classes of associated trees. The tree
    isomorphism function is doing the heavy lifting here.


    EXAMPLES:

    ::


        sage: GraphSignedPermutationClasses(4)
        [[[0, 1, -2, -3], [0, 2, -1, -3]], [[0, -1, -2, -3], [0, 2, 1, -3]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=SignedPermutations(n)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ [L[0]] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2Graph([(j,j+s[j]) for j in range(n)]).is_isomorphic(T2Graph([(j,j+cL[i][0][j]) for j in range(n)])):
                if not s in cL[i]:
                    cL[i].append(s)
                nwT=False
                break
        if nwT==True:
            cL.append([s])
    return cL

def GracefulGraphPermutationClasses(sz):
    """
    Obtain the equivalence classes of permutation graphs 
    defined by the isomorphism classes of associated trees.
    The tree isomorphism function is doing the heavy lifting here.


    EXAMPLES:

    ::


        sage: GracefulGraphPermutationClasses(4)
        [[[0, 1, 2, 3], [0, 2, 1, 3]], [[0, 1, 2, 3], [0, 2, 1, 3]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of signed permutations
    L=GraphSignedPermutationClasses(sz)
    # Obtaining the graceful permutation orbit
    RsLt=[]
    for l in L:
        Tmp=Set([sum(abs(sigma[i])*x^i for i in range(sz)) for sigma in l]).list()
        RsLt.append([[f.subs(x=0)]+[f.coefficient(x^i) for i in range(1,sz)] for f in Tmp])
    return RsLt

def GraphSignedPermutationClassesII(n):
    """
    Obtain the equivalence classes of signed function defined by
    the isomorphism classes of associated trees. The tree
    isomorphism function is doing the heavy lifting here.
    The difference with the function above is that it stores
    only one graph per equivalence classes


    EXAMPLES:

    ::


        sage: GraphSignedPermutationClassesII(4)
        [[0, 1, -2, -3], [0, -1, -2, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=SignedPermutations(n)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism test
    for s in L:
        nwG=True
        for i in range(len(cL)):
            if T2Graph([(j,j+s[j]) for j in range(n)]).is_isomorphic(T2Graph([(j,j+cL[i][j]) for j in range(n)])):
                nwG=False
                break
        if nwG==True:
            cL.append(s)
    return cL

def TreeSignedPermutationClasses(n):
    """
    Obtain the equivalence classes of signed function defined by
    the isomorphism classes of associated trees. The tree
    isomorphism function is doing the heavy lifting here.


    EXAMPLES:

    ::


        sage: TreeSignedPermutationClasses(4)
        [[[0, 1, -2, -3], [0, 2, -1, -3]], [[0, -1, -2, -3], [0, 2, 1, -3]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=TreeSignedPermutations(n)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ [L[0]] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2Graph([(j,j+s[j]) for j in range(n)]).is_isomorphic(T2Graph([(j,j+cL[i][0][j]) for j in range(n)])):
                if not s in cL[i]:
                    cL[i].append(s)
                nwT=False
                break
        if nwT==True:
            cL.append([s])
    return cL

def GracefulTreePermutationClasses(sz):
    """
    Obtain the equivalence classes of permutation graphs 
    defined by the isomorphism classes of associated trees.
    The tree isomorphism function is doing the heavy lifting here.


    EXAMPLES:

    ::


        sage: GracefulTreePermutationClasses(4)
        [[[0, 1, 2, 3], [0, 2, 1, 3]], [[0, 1, 2, 3], [0, 2, 1, 3]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of signed permutations
    L=TreeSignedPermutationClasses(sz)
    # Obtaining the graceful permutation orbit
    RsLt=[]
    for l in L:
        Tmp=Set([sum(abs(sigma[i])*x^i for i in range(sz)) for sigma in l]).list()
        RsLt.append([[f.subs(x=0)]+[f.coefficient(x^i) for i in range(1,sz)] for f in Tmp])
    return RsLt

def TreeSignedPermutationClassesII(n):
    """
    Obtain the equivalence classes of signed function defined by
    the isomorphism classes of associated trees. The tree
    isomorphism function is doing the heavy lifting here.
    the difference with the implementation above is that
    we only store one tree per equivalence class


    EXAMPLES:

    ::


        sage: TreeSignedPermutationClassesII(4)
        [[0, 1, -2, -3], [0, -1, -2, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=TreeSignedPermutations(n)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2Graph([(j,j+s[j]) for j in range(n)]).is_isomorphic(T2Graph([(j,j+cL[i][j]) for j in range(n)])):
                nwT=False
                break
        if nwT==True:
            cL.append(s)
    return cL

def NonTreeSignedPermutationClasses(n):
    """
    Obtain the equivalence classes of signed function defined by
    the isomorphism classes of associated non trees. The graph 
    isomorphism function is doing the heavy lifting here.


    EXAMPLES:

    ::


        sage: NonTreeSignedPermutationClasses(7)
        [[[0, 2, 4, 1, -3, -5, -6],
          [0, 3, 4, -2, -1, -5, -6],
          [0, 5, 1, 2, -4, -3, -6],
          [0, 5, 3, -1, -4, -2, -6]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=NonTreeSignedPermutations(n)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ [L[0]] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2Graph([(j,j+s[j]) for j in range(n)]).is_isomorphic(T2Graph([(j,j+cL[i][0][j]) for j in range(n)])):
                if not s in cL[i]:
                    cL[i].append(s)
                nwT=False
                break
        if nwT==True:
            cL.append([s])
    return cL

def NonTreeSignedPermutationClassesII(n):
    """
    Obtain the equivalence classes of signed function defined by
    the isomorphism classes of associated non trees. The graph 
    isomorphism function is doing the heavy lifting here.


    EXAMPLES:

    ::


        sage: NonTreeSignedPermutationClassesII(7)
        [[0, 2, 4, 1, -3, -5, -6]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=NonTreeSignedPermutations(n)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2Graph([(j,j+s[j]) for j in range(n)]).is_isomorphic(T2Graph([(j,j+cL[i][j]) for j in range(n)])):
                nwT=False
                break
        if nwT==True:
            cL.append(s)
    return cL

def NodeSplittingSignedPermutationClasses(nbv):
    """
    Obtain the equivalence classes of signed function defined by
    the isomorphism classes of associated trees obtained via node
    splittings. The tree isomorphism function is doing the heavy lifting here.


    EXAMPLES:

    ::

        sage: L=NodeSplittingSignedPermutationClasses(4); L
        [[[0, -1, -2, -3], [0, 2, 1, -3]], [[0, 2, -1, -3], [0, 1, -2, -3]]]        
       

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list implementation
    # of binary tree storing the graceful graphs
    # whose degree sequence are lexicographically 
    # smaller then that of the second input tree T2.
    L=[[0,-1]]
    # Initialization of the indexing variable
    indx=0
    while len(L[indx]) <= nbv:
        L.append(LeftSignedPermutation(L[indx]))
        TmpL=MiddleSignedPermutationsI(L[indx])
        for t in TmpL:
            if not t in L:
                L.append(t) 
        TmpL=MiddleSignedPermutationsII(L[indx])
        for t in TmpL:
            if not t in L:
                L.append(t) 
        L.append(RightSignedPermutation(L[indx]))
        indx=indx+1
    TmpL=[l for l in L if len(l)==nbv]; L=TmpL
    # Initialization of the list storing the equivalence class of trees.
    cL=[ [L[0]] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if SP2Graph(s).is_isomorphic(SP2Graph(cL[i][0])):
                if not s in cL[i]:
                    cL[i].append(s)
                nwT=False
                break
        if nwT==True:
            cL.append([s])
    return cL

def NodeSplittingSignedPermutationClassesII(nbv):
    """
    Obtain the equivalence classes of signed function defined by
    the isomorphism classes of associated trees obtained via node
    splittings. The tree isomorphism function is doing the heavy lifting here.
    The difference with the implementation above is that only one tree is stored
    per equivalence classes.


    EXAMPLES:

    ::

        sage: L=NodeSplittingSignedPermutationClassesII(4); L
        [[0, -1, -2, -3], [0, 2, -1, -3]]        
       

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list implementation
    # of binary tree storing the graceful graphs
    # whose degree sequence are lexicographically 
    # smaller then that of the second input tree T2.
    L=[[0,-1]]
    # Initialization of the indexing variable
    indx=0
    while len(L[indx]) <= nbv:
        L.append(LeftSignedPermutation(L[indx]))
        TmpL=MiddleSignedPermutationsI(L[indx])
        for t in TmpL:
            if not t in L:
                L.append(t) 
        TmpL=MiddleSignedPermutationsII(L[indx])
        for t in TmpL:
            if not t in L:
                L.append(t) 
        L.append(RightSignedPermutation(L[indx]))
        indx=indx+1
    TmpL=[l for l in L if len(l)==nbv]; L=TmpL
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if SP2Graph(s).is_isomorphic(SP2Graph(cL[i])):
                nwT=False
                break
        if nwT==True:
            cL.append(s)
    return cL

@cached_function
def GracefulNonTreeFunctionsTuples(n):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the list of graceful functions on the vertices derived from 
    graceful permutations.


    EXAMPLES:

    ::


        sage: GracefulTreeFunctionsTuples(7)
        [[[], [0, 1, 2, 3, 4, 5, 6]],
         [[], [0, 1, 3, 2, 4, 5, 6]],
         [[], [0, 1, 4, 2, 3, 5, 6]],
         [[], [0, 1, 4, 3, 2, 5, 6]],
         [[], [0, 2, 1, 3, 4, 5, 6]],
         [[], [0, 2, 3, 1, 4, 5, 6]],
         [[[(0, 0), (1, 3), (2, 6), (3, 4), (4, 1), (5, 0), (6, 0)]],
          [0, 2, 4, 1, 3, 5, 6]],
         [[], [0, 2, 4, 3, 1, 5, 6]],
         [[], [0, 3, 1, 2, 4, 5, 6]],
         [[], [0, 3, 2, 1, 4, 5, 6]],
         [[], [0, 3, 4, 1, 2, 5, 6]],
         [[[(0, 0), (1, 4), (2, 6), (3, 1), (4, 3), (5, 0), (6, 0)]],
          [0, 3, 4, 2, 1, 5, 6]],
         [[], [0, 4, 1, 2, 3, 5, 6]],
         [[], [0, 4, 1, 3, 2, 5, 6]],
         [[], [0, 4, 2, 1, 3, 5, 6]],
         [[], [0, 4, 2, 3, 1, 5, 6]],
         [[], [0, 4, 3, 1, 2, 5, 6]],
         [[], [0, 4, 3, 2, 1, 5, 6]],
         [[], [0, 5, 1, 2, 3, 4, 6]],
         [[[(0, 0), (1, 6), (2, 3), (3, 5), (4, 0), (5, 2), (6, 0)]],
          [0, 5, 1, 2, 4, 3, 6]],
         [[], [0, 5, 1, 3, 2, 4, 6]],
         [[], [0, 5, 1, 3, 4, 2, 6]],
         [[], [0, 5, 2, 1, 3, 4, 6]],
         [[], [0, 5, 2, 1, 4, 3, 6]],
         [[], [0, 5, 2, 3, 1, 4, 6]],
         [[], [0, 5, 2, 3, 4, 1, 6]],
         [[], [0, 5, 3, 1, 2, 4, 6]],
         [[[(0, 0), (1, 6), (2, 5), (3, 2), (4, 0), (5, 3), (6, 0)]],
          [0, 5, 3, 1, 4, 2, 6]],
         [[], [0, 5, 3, 2, 1, 4, 6]],
         [[], [0, 5, 3, 2, 4, 1, 6]],
         [[], [0, 5, 4, 1, 2, 3, 6]],
         [[], [0, 5, 4, 1, 3, 2, 6]],
         [[], [0, 5, 4, 2, 1, 3, 6]],
         [[], [0, 5, 4, 2, 3, 1, 6]],
         [[], [0, 5, 4, 3, 1, 2, 6]],
         [[], [0, 5, 4, 3, 2, 1, 6]]]
        sage: L=[]
        sage: TpL=GracefulTreeFunctionsTuples(7)
        sage: for l in TpL:
        ....:     L=L+l[0]
        ....:
        sage: L
        [[(0, 0), (1, 3), (2, 6), (3, 4), (4, 1), (5, 0), (6, 0)],
         [(0, 0), (1, 4), (2, 6), (3, 1), (4, 3), (5, 0), (6, 0)],
         [(0, 0), (1, 6), (2, 3), (3, 5), (4, 0), (5, 2), (6, 0)],
         [(0, 0), (1, 6), (2, 5), (3, 2), (4, 0), (5, 3), (6, 0)]]


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
            # Testing that only the first criteria is met
            elif (j > 0) and (p[j]<(n-j)) and not (j>=p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
            # Testing that only the second criteria is met
            elif (j > 0) and not (p[j] < (n-j)) and (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j-p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-p[j])]
            # Testing that all two criterias are met
            elif (j > 0) and (p[j] < (n-j)) and (j >= p[j]):
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
            if not is_Tree( sum([Id[:,t[0]]*Id[t[1],:]+Id[:,t[1]]*Id[t[0],:] for t in il]) ):
                g.append(il)
        L.append([g, p])
    return L

@cached_function
def GracefulTreeFunctionsTuples(n):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the list of graceful functions on the vertices 
    derived from graceful permutations.


    EXAMPLES:

    ::


        sage: GracefulTreeFunctionsTuples(4)
        [[[[(0, 0), (1, 2), (2, 0), (3, 0)], [(0, 0), (1, 0), (2, 0), (3, 0)]],
          [0, 1, 2, 3]],
        [[[(0, 0), (1, 3), (2, 3), (3, 0)], [(0, 0), (1, 3), (2, 1), (3, 0)]],
          [0, 2, 1, 3]]]
        sage: L=[]
        sage: TpL=GracefulTreeFunctionsTuples(4)
        sage: for l in TpL:
        ....:     L=L+l[0]
        ....:
        sage: L
        [[(0, 0), (1, 2), (2, 0), (3, 0)],
         [(0, 0), (1, 0), (2, 0), (3, 0)],
         [(0, 0), (1, 3), (2, 3), (3, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0)]]


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
            # Testing that only the first criteria is met
            elif (j > 0) and (p[j]<(n-j)) and not (j>=p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
            # Testing that only the second criteria is met
            elif (j > 0) and not (p[j] < (n-j)) and (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j-p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-p[j])]
            # Testing that all two criterias are met
            elif (j > 0) and (p[j] < (n-j)) and (j >= p[j]):
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
        [[(0, 0), (1, 0), (2, 5), (3, 5), (4, 0), (5, 0)],
         [(0, 0), (1, 5), (2, 0), (3, 0), (4, 5), (5, 0)]]


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

def NodeSplittingList(T, DegSeq):
    """
    The method returns list description of an array implementation
    of a binary tree.

    EXAMPLES:

    ::

        sage: L=NodeSplittingList([(0, 0), (1, 0)], [3, 1, 1, 1]); L
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
         [(0, 0), (1, 4), (2, 0), (3, 2), (4, 0)],
         [(0, 0), (1, 4), (2, 4), (3, 2), (4, 0)],
         [(0, 0), (1, 2), (2, 0), (3, 0), (4, 0)],
         [(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],
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

def NodeSplittingListII(T, nbv):
    """
    The method returns list description of an array implementation
    of a binary tree.

    EXAMPLES:

    ::

        sage: L=NodeSplittingListII([(0, 0), (1, 0)], 3); L
        [[(0, 0), (1, 0)],
         [(0, 0), (1, 0), (2, 0)],
         [(0, 0), (1, 2), (2, 0)],
         [(0, 0), (1, 0), (2, 0), (3, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0)],
         [(0, 0), (1, 2), (2, 0), (3, 0)],
         [(0, 0), (1, 3), (2, 3), (3, 0)]]
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
    while len(L[indx]) <= nbv:
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

def NodeSplittingGeneratorList(T, DegSeq):
    """
    The method returns list description of an array implementation
    of a binary tree.

    EXAMPLES:

    ::

        sage: L=NodeSplittingGeneratorList([(0, 0), (1, 0)], [3, 1, 1, 1]); L
        [[(0, 0), (1, 0), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 0), (3, 4), (4, 0)],
         [(0, 0), (1, 4), (2, 0), (3, 2), (4, 0)],
         [(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)]]
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
    L=NodeSplittingList(T, DegSeq)
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

def T2Graph(T):
    """
    The method returns an undirected graph object associated with 
    with the tuple list description of the directed graph

    EXAMPLES:

    ::

        sage: T2Graph([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)]).degree_sequence()
        [4, 2, 2, 1, 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(len(T))
    return Graph(sum(Id[:,t[0]]*Id[t[1],:] for t in T[1:])+sum(Id[:,t[0]]*Id[t[1],:] for t in T[1:]).transpose())

def SP2Graph(sigma):
    """
    The method returns an undirected graph object associated with 
    with the tuple list description of the directed graph

    EXAMPLES:

    ::

        sage: SP2Graph([(0-0), -(1-2), -(2-4), -(3-0), -(4-0)]).degree_sequence()
        [4, 2, 2, 1, 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the Tree tuple function
    T=[(i,i+sigma[i]) for i in range(len(sigma))]
    # Initialization of the identity matrix
    Id=identity_matrix(len(T))
    return Graph(sum(Id[:,t[0]]*Id[t[1],:] for t in T[1:])+sum(Id[:,t[0]]*Id[t[1],:] for t in T[1:]).transpose())

def TreeFunctionList(n):
    """
    Goes through all the functions and determines which ones
    are associated with trees.

    EXAMPLES:
    ::
        sage: TreeFunctionList(4)
        [[0, 0, 0],
        [2, 0, 0],
        [3, 0, 0],
        [0, 1, 0],
        [3, 1, 0],
        [0, 3, 0],
        [2, 3, 0],
        [3, 3, 0],
        [0, 0, 1],
        [2, 0, 1],
        [0, 1, 1],
        [0, 3, 1],
        [0, 0, 2],
        [2, 0, 2],
        [3, 0, 2],
        [0, 1, 2]]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[n for i in range(n-1)]; Lf=[]
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Appending the function to the list
        if is_Tree(sum([Id[:,j]*Id[f[j-1],:]+Id[:,f[j-1]]*Id[j,:] for j in range(1,n)])):
            Lf.append(f)
    return Lf

@cached_function
def GracefulK_arryTreeFunctionsTuples(n,K):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the list of graceful functions associated with K-arry trees 
    (degree of all vertices is bounded by K) derived from graceful
    permutations.


    EXAMPLES:

    ::


        sage: GracefulK_arryTreeFunctionsTuples(4,2)
        [[[[(0, 0), (1, 2), (2, 0), (3, 0)]], [0, 1, 2, 3]],
         [[[(0, 0), (1, 3), (2, 1), (3, 0)]], [0, 2, 1, 3]]]
        sage: L=[]
        sage: TpL=GracefulK_arryTreeFunctionsTuples(4,2)
        sage: for l in TpL:
        ....:     L=L+l[0]
        ....:
        sage: L
        [[(0, 0), (1, 2), (2, 0), (3, 0)], [(0, 0), (1, 3), (2, 1), (3, 0)]]


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
            # Testing that only the first criteria is met
            elif (j > 0) and (p[j]<(n-j)) and not (j>=p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
            # Testing that only the second criteria is met
            elif (j > 0) and not (p[j] < (n-j)) and (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j-p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-p[j])]
            # Testing that all two criterias are met
            elif (j > 0) and (p[j] < (n-j)) and (j >= p[j]):
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
            if Set(EdgeListDegreeSequence(il)).issubset(Set(range(1,K+1))) and is_Tree( sum([Id[:,t[0]]*Id[t[1],:]+Id[:,t[1]]*Id[t[0],:] for t in il]) ):
                g.append(il)
        L.append([g,p])
    return L

def GracefulK_arryDiTreeList(n,K):
    """
    Goes through all the permutation of graceful second index functions
    and outputs the list of undirected adjacency matrices of graceful 
    Karry tree (the degree of all vertices is bounded by K).

    EXAMPLES:

    ::


        sage: for l in GracefulK_arryDiTreeList(5,2):
        ....:     t=0
        ....:     for dg in l[0]:
        ....:         dg.plot().save(str(l[1]).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
        ....:         t=t+1
        ....:


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GracefulK_arryTreeFunctionsTuples(n,K)
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

def NodeSplittingK_arryTreesList(T, DegSeq, K):
    """
    The method returns list of functions
    tuples corresponding to gracefully labeled
    trees whose degree sequence is lexicographically
    smaller and have vertices of degree bounded by K.
    exceeds the length of DegSeq.


    EXAMPLES:

    ::

        sage: L=NodeSplittingK_arryTreesList([(0, 0), (1, 0)], [3, 1, 1, 1], 3); L
        [[(0, 0), (1, 0)],
         [(0, 0), (1, 0), (2, 0)],
         [(0, 0), (1, 2), (2, 0)],
         [(0, 0), (1, 0), (2, 0), (3, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0)],
         [(0, 0), (1, 2), (2, 0), (3, 0)],
         [(0, 0), (1, 3), (2, 3), (3, 0)],
         [(0, 0), (1, 0), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 0), (3, 4), (4, 0)],
         [(0, 0), (1, 4), (2, 1), (3, 1), (4, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 0), (3, 2), (4, 0)],
         [(0, 0), (1, 4), (2, 4), (3, 2), (4, 0)],
         [(0, 0), (1, 2), (2, 0), (3, 0), (4, 0)],
         [(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 3), (3, 1), (4, 0)],
         [(0, 0), (1, 3), (2, 3), (3, 0), (4, 0)]]
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
    while (sorted(EdgeListDegreeSequence(L[indx])))[::-1] <= (sorted(DegSeq))[::-1] and len(L[indx])<=len(DegSeq)+1:
        #L.append(LeftChild(L[indx]))
        t=LeftChild(L[indx])
        if Set(EdgeListDegreeSequence(t)).issubset(Set(range(1,K+1))):
            L.append(t) 
        TmpL=MidChildrenI(L[indx])
        for t in TmpL:
            if not t in L and Set(EdgeListDegreeSequence(t)).issubset(Set(range(1,K+1))):
                L.append(t)
        TmpL=MidChildrenII(L[indx])
        for t in TmpL:
            if not t in L and Set(EdgeListDegreeSequence(t)).issubset(Set(range(1,K+1))):
                L.append(t) 
        #L.append(RightChild(L[indx]))
        t=RightChild(L[indx])
        if Set(EdgeListDegreeSequence(t)).issubset(Set(range(1,K+1))):
            L.append(t) 
        indx=indx+1
    return L

def GracefulK_arryDiTreeNodeSplittingList(T, DegSeq, K):
    """
    Starts from the input tree T and performs the node splitting routines
    untill we encounter a K_arry tree (all vertices have degree bounded 
    by K) whose degree sequence is lexicographically greater then the 
    input degree sequence DegSeq.


    EXAMPLES:

    ::

        sage: t=0
        sage: for dg in GracefulK_arryDiTreeNodeSplittingList([(0, 0), (1, 0)], [1, 1, 2, 2], 2):
        ....:     dg.plot().save(str(t)+'.png')
        ....:     t=t+1
        ....:


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of directed graph as edge list.
    L=NodeSplittingK_arryTreesList(T, DegSeq, K)
    # Initialization of the list storing the result
    Rslt=[]
    for lt in L:
        # Initialization of the identity matrix
        Id=identity_matrix(len(lt))
        Rslt.append(DiGraph(sum([Id[:,t[0]]*Id[t[1],:] for t in lt])))
    return Rslt

def TupleSubgraphList(T,k):
    """
    Returns a list of tuple description of subgraph on
    k vertices.
    

    EXAMPLES:
    ::


        sage: TupleSubgraphList([(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)], 3)
        [[(1, 0), (2, 0), (3, 0)],
         [(2, 0), (3, 0), (4, 0)],
         [(1, 0), (3, 0), (4, 0)],
         [(1, 0), (2, 0), (4, 0)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the vertex set
    V=Set(T[1:])
    # Initialization of the list
    L=[]
    for s in V.subsets(k):
        L.append(s.list())
    # Sorting the list
    for i in range(len(L)):
        L[i].sort()
    return L

def T2GraphII(T,sz):
    """
    The method returns an undirected graph object associated with 
    with the tuple list description of the directed graph on sz
    vertices.

    EXAMPLES:

    ::

        sage: T2GraphII([(1, 2), (3, 0)],5).degree_sequence()
        [4, 2, 2, 1, 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    return Graph(sum(Id[:,t[0]]*Id[t[1],:] for t in T)+sum(Id[:,t[0]]*Id[t[1],:] for t in T).transpose())

def generate_script_fast_graceful_list(sz):
    """
    The produces an implementation an optimal sage script for
    displaying and listing gracefully labled graph on sz vertices.
    sz must be greater than 4. Make sure to include the
    Hypermatrix Package in the working directory

    EXAMPLES:

    ::

        sage: !touch fast_listing_of_graceful_trees_on_10_vertices.sage
        sage: !rm fast_listing_of_graceful_trees_on_10_vertices.sage
        sage: generate_script_fast_graceful_list(10)
        sage: !ls fast_listing_of_graceful_trees_on_10_vertices.sage
        fast_listing_of_graceful_trees_on_10_vertices.sage
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Creating the string corresponding to the file name
    filename = 'fast_listing_of_graceful_trees_on_'+str(sz)+'_vertices.sage'
    # Opening the file
    f = open(filename,'w')
    f.write('# Loading the Hypermatrix Package\n')
    f.write("load('./Hypermatrix_Algebra_tst.sage')\n\n")
    f.write('# Loading the Graceful Graph Package\n')
    f.write("load('./graceful_graph_package.sage')\n\n")
    f.write('# Initialization of the size parameter\n')
    f.write('sz='+str(sz)+'\n\n')
    f.write('# Initialization of the identity matrix\n')
    f.write('Id=identity_matrix(sz)\n\n')
    f.write('# Initialization of the Permutations\n')
    f.write('P=Permutations(floor((sz-1)/2))\n\n')
    f.write('# List storing the result\n')
    f.write('RsLt=[] \n\n')
    f.write('# Main loop\n')
    f.write('for n_2 in Set(range(1,2)+[sz-(i+1) for i in range(1,2)]):\n')
    # variable storing the spaces
    sp = '    '
    for i in range(3,floor(sz/2)+1):
        tmpString='for n_'+str(i)+' in Set(range(1,'+str(i)+')+[sz-(i+1) for i in range(1,'+str(i)+')]).difference(Set(['
        for j in range(2,i-1):
            tmpString=tmpString+'n_'+str(j)+','
        tmpString=tmpString+'n_'+str(i-1)+'])):\n'
        f.write(sp+tmpString)
        sp=sp+'    '
    f.write(sp+"Tmp=HM(sz,1,'zero'); Tmp[sz-1,0]=sz-1\n")
    for i in range(2,floor(sz/2)+1):
        f.write(sp+'Tmp[n_'+str(i)+',0]=sz-'+str(i)+'\n')
    f.write(sp+"# Indentifying the indices which are zero\n")
    f.write(sp+"IndX=[i for i in range(1,sz) if Tmp[i,0].is_zero()]\n")
    f.write(sp+"# Loop filling up the rest of the list\n")
    f.write(sp+"for q in P:\n")
    f.write(sp+"    for i in range(len(IndX)):\n")
    f.write(sp+"        Tmp[IndX[i],0]=q[i]\n")
    f.write(sp+"    p=Tmp.copy().list()\n")
    f.write(sp+"    # Initialization of the list of functions.\n")
    f.write(sp+"    c=[]\n")
    f.write(sp+"    for j in range(sz):\n")
    f.write(sp+"        if j == 0:\n")
    f.write(sp+"            if len(c)==0:\n")
    f.write(sp+"                c.append([(j, p[j])])\n")
    f.write(sp+"            else:\n")
    f.write(sp+"                for i in range(len(c)):\n")
    f.write(sp+"                    c[i]=c[i]+[(j, p[j])]\n")
    f.write(sp+"        # testing that only the first criteria is met\n")
    f.write(sp+"        elif (j > 0) and (p[j]<(sz-j)) and not (j>=p[j]):\n")
    f.write(sp+"            if len(c)==0:\n")
    f.write(sp+"                c.append([(j, j+p[j])])\n")
    f.write(sp+"            else:\n")
    f.write(sp+"                for i in range(len(c)):\n")
    f.write(sp+"                    c[i]=c[i]+[(j, j+p[j])]\n")
    f.write(sp+"        # testing that only the second criteria is met\n")
    f.write(sp+"        elif j > 0 and not (p[j] < (sz-j)) and (j >= p[j]):\n")
    f.write(sp+"            if len(c)==0:\n")
    f.write(sp+"                c.append([(j, j-p[j])])\n")
    f.write(sp+"            else:\n")
    f.write(sp+"                for i in range(len(c)):\n")
    f.write(sp+"                    c[i]=c[i]+[(j, j-p[j])]\n")
    f.write(sp+"        # testing that all two criterias are met\n")
    f.write(sp+"        elif j > 0 and (p[j] < (sz-j)) and (j >= p[j]):\n")
    f.write(sp+"            if len(c)==0:\n")
    f.write(sp+"                c.append([(j, j+p[j])])\n")
    f.write(sp+"                c.append([(j, j-p[j])])\n")
    f.write(sp+"            else:\n")
    f.write(sp+"                d=copy(c)\n")
    f.write(sp+"                for i in range(len(c)):\n")
    f.write(sp+"                    c[i]=c[i]+[(j, j+p[j])]\n")
    f.write(sp+"                for i in range(len(d)):\n")
    f.write(sp+"                    d[i]=d[i]+[(j, j-p[j])]\n")
    f.write(sp+"                c=c+d\n")
    f.write(sp+"    # Initializing the list which will store the generators\n")
    f.write(sp+"    RsLt=RsLt+c\n")
    f.write(sp+"    tc=0\n")
    f.write(sp+"    for il in c:\n")
    f.write(sp+"        # Initilization of the adjacency matrix\n")
    f.write(sp+"        if is_Tree( sum([Id[:,t[0]]*Id[t[1],:]+Id[:,t[1]]*Id[t[0],:] for t in il]) ):\n")
    f.write(sp+"            DiGraph(sum([Id[:,t[0]]*Id[t[1],:] for t in il])).plot().save(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(tc)+'.png')\n")
    f.write(sp+"            tc=tc+1\n")
    # Closing the file
    f.close()

def generate_labeling_script(T):
    """
    The produces an implementation an optimal sage script for
    displaying and listing of graceful lablings of the input
    tree specified as a function tuple including the self loop
    at the vertex labeled 0.
    Make sure to include the Hypermatrix Package in
    the working directory


    EXAMPLES:

    ::

        sage: !touch graceful_labeling_of_a_trees_on_6_vertices.sage
        sage: !rm graceful_labeling_of_a_trees_on_6_vertices.sage
        sage: generate_labeling_scriptII([(0,0),(1,2),(2,3),(3,4),(4,5),(5,0)])
        sage: !ls graceful_labeling_of_a_trees_on_6_vertices.sage
        graceful_labeling_of_a_trees_on_6_vertices.sage
        

    AUTHORS:
    - Edinah K. Gnang
    """
    sz=len(T)
    # Creating the string corresponding to the file name
    filename='graceful_labeling_of_a_trees_on_'+str(sz)+'_vertices.sage'
    # Opening the file
    f=open(filename,'w')
    f.write('# Loading the Hypermatrix Package\n')
    f.write("load('./Hypermatrix_Algebra_tst.sage')\n\n")
    f.write('# Loading the Graceful Graph Package\n')
    f.write("load('./graceful_graph_package.sage')\n\n")
    f.write('# Initialization of the Tuple function\n')
    f.write('T='+str(T)+'\n\n')
    f.write('# Initialization of the size parameter\n')
    f.write('sz=len(T)\n\n')
    f.write('# Initialization of the identity matrix\n')
    f.write('Id=identity_matrix(sz)\n\n')
    f.write('# Initialization of the Permutations\n')
    f.write('P=Permutations(floor((sz-1)/2))\n\n')
    f.write('# List storing the result\n')
    f.write('RsLt=[]\n\n')
    f.write('# Main loop\n')
    f.write('for n_2 in Set(range(1,2)+[sz-(i+1) for i in range(1,2)]):\n')
    # variable storing the spaces
    sp = '    '
    for i in range(3,floor(sz/2)+1):
        tmpString='for n_'+str(i)+' in Set(range(1,'+str(i)+')+[sz-(i+1) for i in range(1,'+str(i)+')]).difference(Set(['
        for j in range(2,i-1):
            tmpString=tmpString+'n_'+str(j)+','
        tmpString=tmpString+'n_'+str(i-1)+'])):\n'
        f.write(sp+tmpString)
        sp=sp+'    '
    f.write(sp+"Tmp=HM(sz,1,'zero'); Tmp[sz-1,0]=sz-1\n")
    for i in range(2,floor(sz/2)+1):
        f.write(sp+'Tmp[n_'+str(i)+',0]=sz-'+str(i)+'\n')
    f.write(sp+"# Indentifying the indices which are zero\n")
    f.write(sp+"IndX=[i for i in range(1,sz) if Tmp[i,0].is_zero()]\n")
    f.write(sp+"TmpTple=[(i,Tmp[i,0]) for i in range(1,sz) if i in IndX]\n")
    f.write(sp+"CompTmpTple=[]\n")
    f.write(sp+"for i in Set(range(1,sz)).difference(Set(IndX)):\n")
    f.write(sp+"    if Tmp[i,0]+i > sz-1:\n")
    f.write(sp+"        CompTmpTple.append((i,i-Tmp[i,0]))\n")
    f.write(sp+"    else:\n")
    f.write(sp+"        CompTmpTple.append((i,i+Tmp[i,0]))\n")
    f.write(sp+"subGrph=False\n")
    f.write(sp+"for t in TupleSubgraphList(T,sz-1-floor((sz-1)/2)):\n")
    f.write(sp+"    if T2GraphII(CompTmpTple,sz).is_isomorphic(T2GraphII(t,sz)):\n")
    f.write(sp+"        subGrph=True; break\n")
    f.write(sp+"if subGrph:\n")
    f.write(sp+"    # Loop filling up the rest of the list\n")
    f.write(sp+"    for q in P:\n")
    f.write(sp+"        for i in range(len(IndX)):\n")
    f.write(sp+"            Tmp[IndX[i],0]=q[i]\n")
    f.write(sp+"        p=Tmp.copy().list()\n")
    f.write(sp+"        # Initialization of the list of functions.\n")
    f.write(sp+"        c=[]\n")
    f.write(sp+"        for j in range(sz):\n")
    f.write(sp+"            if j == 0:\n")
    f.write(sp+"                if len(c)==0:\n")
    f.write(sp+"                    c.append([(j, p[j])])\n")
    f.write(sp+"                else:\n")
    f.write(sp+"                    for i in range(len(c)):\n")
    f.write(sp+"                        c[i]=c[i]+[(j, p[j])]\n")
    f.write(sp+"            # testing that only the first criteria is met\n")
    f.write(sp+"            elif (j > 0) and (p[j]<(sz-j)) and not (j>=p[j]):\n")
    f.write(sp+"                if len(c)==0:\n")
    f.write(sp+"                    c.append([(j, j+p[j])])\n")
    f.write(sp+"                else:\n")
    f.write(sp+"                    for i in range(len(c)):\n")
    f.write(sp+"                        c[i]=c[i]+[(j, j+p[j])]\n")
    f.write(sp+"            # testing that only the second criteria is met\n")
    f.write(sp+"            elif j > 0 and not (p[j] < (sz-j)) and (j >= p[j]):\n")
    f.write(sp+"                if len(c)==0:\n")
    f.write(sp+"                    c.append([(j, j-p[j])])\n")
    f.write(sp+"                else:\n")
    f.write(sp+"                    for i in range(len(c)):\n")
    f.write(sp+"                        c[i]=c[i]+[(j, j-p[j])]\n")
    f.write(sp+"            # testing that all two criterias are met\n")
    f.write(sp+"            elif j > 0 and (p[j] < (sz-j)) and (j >= p[j]):\n")
    f.write(sp+"                if len(c)==0:\n")
    f.write(sp+"                    c.append([(j, j+p[j])])\n")
    f.write(sp+"                    c.append([(j, j-p[j])])\n")
    f.write(sp+"                else:\n")
    f.write(sp+"                    d=copy(c)\n")
    f.write(sp+"                    for i in range(len(c)):\n")
    f.write(sp+"                        c[i]=c[i]+[(j, j+p[j])]\n")
    f.write(sp+"                    for i in range(len(d)):\n")
    f.write(sp+"                        d[i]=d[i]+[(j, j-p[j])]\n")
    f.write(sp+"                    c=c+d\n")
    f.write(sp+"        # Initializing the list which will store the trees\n")
    f.write(sp+"        TmpRsLt=[]\n")
    f.write(sp+"        for t in c:\n")
    f.write(sp+"            if T2Graph(t).is_isomorphic(T2Graph(T)):\n")
    f.write(sp+"                TmpRsLt.append(t)\n")
    f.write(sp+"        RsLt=RsLt+TmpRsLt\n")
    f.write(sp+"        tc=0\n")
    f.write(sp+"        for il in TmpRsLt:\n")
    f.write(sp+"            # Initilization of the adjacency matrix\n")
    f.write(sp+"            if is_Tree( sum([Id[:,t[0]]*Id[t[1],:]+Id[:,t[1]]*Id[t[0],:] for t in il]) ):\n")
    f.write(sp+"                DiGraph(sum([Id[:,t[0]]*Id[t[1],:] for t in il])).plot().save(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(tc)+'.png')\n")
    f.write(sp+"            tc=tc+1\n")
    # Closing the file
    f.close()

def generate_labeling_scriptII(T):
    """
    The produces an implementation an optimal sage script for
    displaying and listing of graceful lablings of the input
    tree specified as a function tuple including the self loop
    at the vertex labeled 0.
    Make sure to include the Hypermatrix Package in
    the working directory


    EXAMPLES:

    ::

        sage: !touch graceful_labeling_of_a_trees_on_6_vertices.sage
        sage: !rm graceful_labeling_of_a_trees_on_6_vertices.sage
        sage: generate_labeling_scriptII([(0,0),(1,2),(2,3),(3,4),(4,5),(5,0)])
        sage: !ls graceful_labeling_of_a_trees_on_6_vertices.sage
        graceful_labeling_of_a_trees_on_6_vertices.sage
        

    AUTHORS:
    - Edinah K. Gnang
    """
    sz=len(T)
    # Creating the string corresponding to the file name
    filename='graceful_labeling_of_a_trees_on_'+str(sz)+'_vertices.sage'
    # Opening the file
    f=open(filename,'w')
    f.write('# Loading the Hypermatrix Package\n')
    f.write("load('./Hypermatrix_Algebra_tst.sage')\n\n")
    f.write('# Loading the Graceful Graph Package\n')
    f.write("load('./graceful_graph_package.sage')\n\n")
    f.write('# Initialization of the Tuple function\n')
    f.write('T='+str(T)+'\n\n')
    f.write('def Label(T):\n')
    f.write('    # Initialization of the size parameter\n')
    f.write('    sz=len(T)\n\n')
    f.write('    # Initialization of the identity matrix\n')
    f.write('    Id=identity_matrix(sz)\n\n')
    f.write('    # Initialization of the Permutations\n')
    f.write('    P=Permutations(floor((sz-1)/2))\n\n')
    f.write('    # List storing the result\n')
    f.write('    RsLt=[]\n\n')
    f.write('    # Main loop\n')
    f.write('    for n_2 in Set(range(1,2)+[sz-(i+1) for i in range(1,2)]):\n')
    # variable storing the spaces
    sp = '        '
    for i in range(3,floor(sz/2)+1):
        tmpString='for n_'+str(i)+' in Set(range(1,'+str(i)+')+[sz-(i+1) for i in range(1,'+str(i)+')]).difference(Set(['
        for j in range(2,i-1):
            tmpString=tmpString+'n_'+str(j)+','
        tmpString=tmpString+'n_'+str(i-1)+'])):\n'
        f.write(sp+tmpString)
        sp=sp+'    '
    f.write(sp+"Tmp=HM(sz,1,'zero'); Tmp[sz-1,0]=sz-1\n")
    for i in range(2,floor(sz/2)+1):
        f.write(sp+'Tmp[n_'+str(i)+',0]=sz-'+str(i)+'\n')
    f.write(sp+"# Indentifying the indices which are zero\n")
    f.write(sp+"IndX=[i for i in range(1,sz) if Tmp[i,0].is_zero()]\n")
    f.write(sp+"TmpTple=[(i,Tmp[i,0]) for i in range(1,sz) if i in IndX]\n")
    f.write(sp+"CompTmpTple=[]\n")
    f.write(sp+"for i in Set(range(1,sz)).difference(Set(IndX)):\n")
    f.write(sp+"    if Tmp[i,0]+i > sz-1:\n")
    f.write(sp+"        CompTmpTple.append((i,i-Tmp[i,0]))\n")
    f.write(sp+"    else:\n")
    f.write(sp+"        CompTmpTple.append((i,i+Tmp[i,0]))\n")
    f.write(sp+"subGrph=False\n")
    f.write(sp+"for t in TupleSubgraphList(T,sz-1-floor((sz-1)/2)):\n")
    f.write(sp+"    if T2GraphII(CompTmpTple,sz).is_isomorphic(T2GraphII(t,sz)):\n")
    f.write(sp+"        subGrph=True; break\n")
    f.write(sp+"if subGrph:\n")
    f.write(sp+"    # Loop filling up the rest of the list\n")
    f.write(sp+"    for q in P:\n")
    f.write(sp+"        for i in range(len(IndX)):\n")
    f.write(sp+"            Tmp[IndX[i],0]=q[i]\n")
    f.write(sp+"        p=Tmp.copy().list()\n")
    f.write(sp+"        # Initialization of the list of functions.\n")
    f.write(sp+"        c=[]\n")
    f.write(sp+"        for j in range(sz):\n")
    f.write(sp+"            if j == 0:\n")
    f.write(sp+"                if len(c)==0:\n")
    f.write(sp+"                    c.append([(j, p[j])])\n")
    f.write(sp+"                else:\n")
    f.write(sp+"                    for i in range(len(c)):\n")
    f.write(sp+"                        c[i]=c[i]+[(j, p[j])]\n")
    f.write(sp+"            # testing that only the first criteria is met\n")
    f.write(sp+"            elif (j > 0) and (p[j]<(sz-j)) and not (j>=p[j]):\n")
    f.write(sp+"                if len(c)==0:\n")
    f.write(sp+"                    c.append([(j, j+p[j])])\n")
    f.write(sp+"                else:\n")
    f.write(sp+"                    for i in range(len(c)):\n")
    f.write(sp+"                        c[i]=c[i]+[(j, j+p[j])]\n")
    f.write(sp+"            # testing that only the second criteria is met\n")
    f.write(sp+"            elif j > 0 and not (p[j] < (sz-j)) and (j >= p[j]):\n")
    f.write(sp+"                if len(c)==0:\n")
    f.write(sp+"                    c.append([(j, j-p[j])])\n")
    f.write(sp+"                else:\n")
    f.write(sp+"                    for i in range(len(c)):\n")
    f.write(sp+"                        c[i]=c[i]+[(j, j-p[j])]\n")
    f.write(sp+"            # testing that all two criterias are met\n")
    f.write(sp+"            elif j > 0 and (p[j] < (sz-j)) and (j >= p[j]):\n")
    f.write(sp+"                if len(c)==0:\n")
    f.write(sp+"                    c.append([(j, j+p[j])])\n")
    f.write(sp+"                    c.append([(j, j-p[j])])\n")
    f.write(sp+"                else:\n")
    f.write(sp+"                    d=copy(c)\n")
    f.write(sp+"                    for i in range(len(c)):\n")
    f.write(sp+"                        c[i]=c[i]+[(j, j+p[j])]\n")
    f.write(sp+"                    for i in range(len(d)):\n")
    f.write(sp+"                        d[i]=d[i]+[(j, j-p[j])]\n")
    f.write(sp+"                    c=c+d\n")
    f.write(sp+"        # Initializing the list which will store the generators\n")
    f.write(sp+"        TmpRsLt=[]\n")
    f.write(sp+"        for t in c:\n")
    f.write(sp+"            if T2Graph(t).is_isomorphic(T2Graph(T)):\n")
    f.write(sp+"                return t\n")
    # Closing the file
    f.close()

def generate_script_fast_partial_graceful_list(sz):
    """
    The produces an implementation an optimal sage script for
    displaying and listing partial gracefully labled graph on
    sz vertices, sz must be greater than 4. Make sure to include the
    Hypermatrix Package in the working directory


    EXAMPLES:

    ::

        sage: !touch fast_listing_of_partial_graceful_trees_on_10_vertices.sage
        sage: !rm fast_listing_of_partial_graceful_trees_on_10_vertices.sage
        sage: generate_script_fast_partial_graceful_list(10)
        sage: !ls fast_listing_of_partial_graceful_trees_on_10_vertices.sage
        fast_listing_of_partial_graceful_trees_on_10_vertices.sage
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Creating the string corresponding to the file name
    filename = 'fast_listing_of_partial_graceful_trees_on_'+str(sz)+'_vertices.sage'
    # Opening the file
    f = open(filename,'w')
    f.write('# Loading the Hypermatrix Package\n')
    f.write("load('./Hypermatrix_Algebra_tst.sage')\n\n")
    f.write('# Loading the Graceful Graph Package\n')
    f.write("load('./graceful_graph_package.sage')\n\n")
    f.write('# Initialization of the size parameter\n')
    f.write('sz='+str(sz)+'\n\n')
    f.write('# Initialization of the identity matrix\n')
    f.write('Id=identity_matrix(sz)\n\n')
    f.write('# List storing the result\n')
    f.write('RsLt=[]; tc=0\n\n')
    f.write('# Main loop\n')
    f.write('for n_2 in Set(range(1,2)+[sz-(i+1) for i in range(1,2)]):\n')
    # variable storing the spaces
    sp='    '
    for i in range(3,floor(sz/2)+1):
        tmpString='for n_'+str(i)+' in Set(range(1,'+str(i)+')+[sz-(i+1) for i in range(1,'+str(i)+')]).difference(Set(['
        for j in range(2,i-1):
            tmpString=tmpString+'n_'+str(j)+','
        tmpString=tmpString+'n_'+str(i-1)+'])):\n'
        f.write(sp+tmpString)
        sp=sp+'    '
    f.write(sp+"Tmp=HM(sz,1,'zero'); Tmp[sz-1,0]=sz-1\n")
    for i in range(2,floor(sz/2)+1):
        f.write(sp+'Tmp[n_'+str(i)+',0]=sz-'+str(i)+'\n')
    f.write(sp+"# Initializing the list of tuples\n")
    f.write(sp+"c=[(0,0)]\n")
    f.write(sp+"for i in range(1,sz):\n")
    f.write(sp+"    if i-Tmp[i,0] >= 0 and not Tmp[i,0].is_zero():\n")
    f.write(sp+"        c.append((i,i-Tmp[i,0]))\n")
    f.write(sp+"    elif i-Tmp[i,0] < 0 and not Tmp[i,0].is_zero():\n")
    f.write(sp+"        c.append((i,i+Tmp[i,0]))\n")
    f.write(sp+"# Initializing the list which will store the partial trees\n")
    f.write(sp+"print c\n")
    f.write(sp+"RsLt=RsLt+c\n")
    f.write(sp+"T2GraphII(c,sz).plot().save(str(tc)+'.png')\n")
    f.write(sp+"tc=tc+1\n")
    # Closing the file
    f.close()

def generate_script_fast_partial_graceful_listII(sz):
    """
    The produces an implementation an optimal sage script for
    displaying and listing partial gracefully labled graph on
    sz vertices, sz must be greater than 4. Make sure to include the
    Hypermatrix Package in the working directory. The difference with
    inplementation above is that this function displays one partial
    tree per equivalence classe.


    EXAMPLES:

    ::

        sage: !touch fast_listing_of_partial_graceful_trees_on_10_verticesII.sage
        sage: !rm fast_listing_of_partial_graceful_trees_on_10_verticesII.sage
        sage: generate_script_fast_partial_graceful_listII(10)
        sage: !ls fast_listing_of_partial_graceful_trees_on_10_verticesII.sage
        fast_listing_of_partial_graceful_trees_on_10_verticesII.sage
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Creating the string corresponding to the file name
    filename = 'fast_listing_of_partial_graceful_trees_on_'+str(sz)+'_verticesII.sage'
    # Opening the file
    f = open(filename,'w')
    f.write('# Loading the Hypermatrix Package\n')
    f.write("load('./Hypermatrix_Algebra_tst.sage')\n\n")
    f.write('# Loading the Graceful Graph Package\n')
    f.write("load('./graceful_graph_package.sage')\n\n")
    f.write('# Initialization of the size parameter\n')
    f.write('sz='+str(sz)+'\n\n')
    f.write('# List storing the result\n')
    f.write('RsLt=[]; tc=0\n\n')
    f.write('# Main loop\n')
    f.write('for n_2 in Set(range(1,2)+[sz-(i+1) for i in range(1,2)]):\n')
    # variable storing the spaces
    sp='    '
    for i in range(3,floor(sz/2)+1):
        tmpString='for n_'+str(i)+' in Set(range(1,'+str(i)+')+[sz-(i+1) for i in range(1,'+str(i)+')]).difference(Set(['
        for j in range(2,i-1):
            tmpString=tmpString+'n_'+str(j)+','
        tmpString=tmpString+'n_'+str(i-1)+'])):\n'
        f.write(sp+tmpString)
        sp=sp+'    '
    f.write(sp+"Tmp=HM(sz,1,'zero'); Tmp[sz-1,0]=sz-1\n")
    for i in range(2,floor(sz/2)+1):
        f.write(sp+'Tmp[n_'+str(i)+',0]=sz-'+str(i)+'\n')
    f.write(sp+"# Initializing the list of tuples\n")
    f.write(sp+"c=[]\n")
    f.write(sp+"for i in range(1,sz):\n")
    f.write(sp+"    if i-Tmp[i,0] >= 0 and not Tmp[i,0].is_zero():\n")
    f.write(sp+"        c.append((i,i-Tmp[i,0]))\n")
    f.write(sp+"    elif i-Tmp[i,0] < 0 and not Tmp[i,0].is_zero():\n")
    f.write(sp+"        c.append((i,i+Tmp[i,0]))\n")
    f.write(sp+"# Initializing the list which will store the partial trees\n")
    f.write(sp+"if len(RsLt) == 0:\n")
    f.write(sp+"    RsLt=RsLt+[c]\n")
    f.write(sp+"    tc=tc+1\n")
    f.write(sp+"    print c\n")
    f.write(sp+"    T2GraphII(c,sz).plot().save(str(tc)+'.png')\n")
    f.write(sp+"elif len(RsLt) > 0:\n")
    f.write(sp+"    Is_Absent=True\n")
    f.write(sp+"    for Tple in RsLt:\n")
    f.write(sp+"        if T2GraphII(Tple,sz).is_isomorphic(T2GraphII(c,sz)):\n")
    f.write(sp+"            Is_Absent=False; break\n")
    f.write(sp+"    if Is_Absent:\n")
    f.write(sp+"        RsLt=RsLt+[c]\n")
    f.write(sp+"        tc=tc+1\n")
    f.write(sp+"        print c\n")
    f.write(sp+"        T2GraphII(c,sz).plot().save(str(tc)+'.png')\n")
    # Closing the file
    f.close()

def RandomGraphSignedPermutationClasses(n):
    """
    Returns a random graph selected uniformly among the
    unlabeled graphs which admits a graceful labelin.


    EXAMPLES:

    ::


        sage: RandomGraphSignedPermutationClasses(3)
        [0, -1, -2]


    AUTHORS:
    - Edinah K. Gnang
    """
    L=GraphSignedPermutationClassesII(n)
    return L[randint(0,len(L)-1)]

def RandomTreeSignedPermutationClasses(n):
    """
    Returns a random tress selected uniformly among the
    unlabeled graphs which admits a graceful labelin.


    EXAMPLES:

    ::


        sage: RandomTreeSignedPermutationClasses(3)
        [0, -1, -2]


    AUTHORS:
    - Edinah K. Gnang
    """
    L=TreeSignedPermutationClassesII(n)
    return L[randint(0,len(L)-1)]

def Formator(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs matrix
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    working over SR for both A and b.

    EXAMPLES:

    ::

        sage: x,y = var('x,y')
        sage: CnstrLst = [x+y==1, x-y==2]
        sage: VrbLst = [x, y]
        sage: [A,b] = Formator(CnstrLst, VrbLst)
        sage: A
        [ 1  1]
        [ 1 -1]
        sage: b
        [1]
        [2]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initializing the Matrix
    A=Matrix(SR,len(CnstrLst),len(VrbLst),zero_matrix(len(CnstrLst),len(VrbLst)))
    b=vector(SR, [eq.rhs() for eq in CnstrLst]).column()
    for r in range(len(CnstrLst)):
        for c in range(len(VrbLst)):
            A[r,c]=(CnstrLst[r]).lhs().coefficient(VrbLst[c])
    return [A,b]


def CountGracefulTreeLabelings(sigma):
    """
    Implements the counting formla derived via Cramer's rule
    The implementation is very slow because it computes a 
    sum over all permutations and all signed permutations.
    The function does not work for non trees. Because the
    incidence matrix inversion formula breaks for non
    trees due to the fact that the Nullsapce has dimension
    greater then one.


    EXAMPLES:

    ::

        sage: CountGracefulTreeLabelings([0, 1, -2])[0]
        2


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the tuple representation
    T=SP2T(sigma)
    # Initialization of the size parameter
    sz = len(T)
    # Initialization of the list of variables
    X = [var('x'+str(i)) for i in range(sz)]
    Y = [0]+[var('y'+str(i)) for i in range(1,sz)]
    # Initialization of the edge function
    F=[T[i][1] for i in range(sz)]
    # Initialization of the edge list
    Le = Tuple2EdgeList(T,'x')
    # Initialization of the list of equations
    Eq = [0*X[0] == Y[0]]+[Le[i] == Y[i] for i in range(1,sz)]
    # Formating the constraints to extract the incidence
    # matrix and the edge vector on the right hand side.
    [B,b] = Formator(Eq,X)
    # Alternative derivation of the solution via matrix inversion
    M = zero_matrix(SR,sz,sz)
    M[:sz-1,:sz-1] = B[1:,:sz-1].inverse()
    Lf = (M*Matrix(SR,sz,1,Y[1:]+[0])-min_symbolic((M*Matrix(SR,sz,1,Y[1:]+[0])).list())*ones_matrix(SR,sz,1)).list()
    # Initialization of list of signed permutations
    SP = SignedPermutations(sz)
    # Initialization of the list of permutations.
    P = Permutations(sz-1)
    # Initialization of the result list
    RsLt=[]
    # Initialization of the sum
    f = 0
    for tmq in P:
        q = [0]+[tmq[i] for i in range(sz-1)]
        for p in SP:
            # Initializing the vertex variable values
            Xn = [Lf[i].subs([Y[q[j]] == p[j] for j in range(1,sz)]) for i in range(sz)]
            tmpf1 = prod([prod([(j-Xn[i])/(j-i) for j in range(sz,2*sz-1)]) for i in range(sz)])^2
            tmpf2 = prod([(Xn[j]-Xn[i])/(j-i) for i in range(sz) for j in range(sz) if i<j])^2
            f=f+tmpf1*tmpf2
            if tmpf1 != 0 and tmpf2 != 0:
                RsLt.append([[X[i]==Xn[i] for i in range(sz)], [Le[i]==(B*Matrix(SR,sz,1,Xn)).list()[i] for i in range(sz)]])
                #print 'Xn = ', [X[i]==Xn[i] for i in range(sz)]
                #print 'B*Xn=', [Le[i]==(B*Matrix(SR,sz,1,Xn)).list()[i] for i in range(sz)]
                #print '\n'
    return [f,RsLt]

