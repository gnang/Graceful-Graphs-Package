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


# Loading the Hypermatrix Algebra Package.
load('./Hypermatrix_Algebra_tst.sage')

@cached_function
def GracefulPermutations(sz):
    """
    Goes through all permutations of sz > 1 elements and outputs
    the list of permutations which can be used to construct graceful
    graphs having no isolated vertices.

    EXAMPLES:

    ::

        sage: GracefulPermutations(4)
        [[0, 1, 2, 3], [0, 2, 1, 3]]
        sage: [factorial(floor((i+3-1)/2))*factorial(ceil((i+3-1)/2)) for i in range(7)]
        [1, 2, 4, 12, 36, 144, 576]
        sage: [len(GracefulPermutations(i+3)) for i in range(7)]
        [1, 2, 4, 12, 36, 144, 576]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations of elements from 1 to (n-1).
    P=Permutations(sz-1)
    # Initialization of the list of graceful permutations.
    L=[]
    # Loop collecting graceful permutations.
    for q in P:
        # Prepending 0 to the permutation.
        p=[0]+[q[i] for i in range(sz-1)]
        # Initialization of the boolean variable.
        bl=True
        for i in range(1,sz):
            # Checking the gracefulness condition.
            if not ( p[i] < sz-i or p[i] <= i ):
                bl=False
                break
        if bl:
            L.append(p)
    return L

@cached_function
def GracefulPermutationsII(sz,j):
    """
    Goes through all permutations of sz > 1 elements and outputs
    the list of permutations which can be used to construct graceful
    graphs having no isolated vertices. It is enought to cover values
    j from 0 to floor((sz-1)/2) inclusively.

    EXAMPLES:

    ::

        sage: GracefulPermutationsII(4,0)
        [[0, 1, 2, 3], [0, 2, 1, 3]]
        sage: [factorial(floor((i+3-1)/2))*factorial(ceil((i+3-1)/2)) for i in range(7)]
        [1, 2, 4, 12, 36, 144, 576]
        sage: [len(GracefulPermutationsII(i+3,0)) for i in range(7)]
        [1, 2, 4, 12, 36, 144, 576]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations of elements from 1 to (n-1).
    P=Permutations(sz)
    # Initialization of the list of graceful permutations.
    L=[]
    # Loop collecting graceful permutations.
    for q in P:
        # Prepending 0 to the permutation.
        p=[q[i]-1 for i in range(sz)]
        # Initialization of the boolean variable.
        bl=True
        for i in range(1,sz):
            # Checking the gracefulness condition.
            if not ( p[i]-j < sz-i or p[i]-j <= i ):
                bl=False
                break
        if bl:
            L.append(p)
    return L

@cached_function
def CountGracefulFunctions(sz):
    """
    Goes through the list of permutations of sz > 1 elements
    which can be used to construct of at least one graceful
    graph having no isolated vertices. The method outputs
    the list of permutations preceded by the counts of 
    associated graceful graphs.


    EXAMPLES:

    ::

        sage: CountGracefulFunctions(4)
        [[2, [0, 1, 2, 3]], [2, [0, 2, 1, 3]]]
        sage: TpL=CountGracefulFunctions(5); L=[]
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
    for p in GracefulPermutations(sz):
        # Initialization of the count variable.
        c=1
        # Start from 1 because the loop edge is attached to 0.
        for i in range(1,sz):
            if p[i] < sz-i and p[i] <= i:
                c=2*c
        L.append([c,p])
    return L

@cached_function
def CountGracefulFunctionsII(sz,j):
    """
    Goes through the list of permutations of sz > 1 elements
    which can be used to construct of at least one graceful
    graph having no isolated vertices. The method outputs
    the list of permutations preceded by the counts of 
    associated graceful graphs.


    EXAMPLES:

    ::

        sage: CountGracefulFunctionsII(4,0)
        [[2, [0, 1, 2, 3]], [2, [0, 2, 1, 3]]]
        sage: sz=5; TpL=CountGracefulFunctionsII(sz,0); L=[]
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
    for p in GracefulPermutationsII(sz,j):
        # Initialization of the count variable.
        c=1
        # Start from 1 because the loop edge is attached to 0.
        for i in range(1,sz):
            if p[i]-j < sz-i and p[i]-j <= i:
                c=2*c
        L.append([c, p])
    return L

def CountGracefulTrees(sz):
    """
    Goes through all the permutation of sz > 1 elements and enumerates 
    graceful trees on the vertices derived from graceful permutations.


    EXAMPLES:

    ::


        sage: CountGracefulTrees(4)
        [[2, [0, 1, 2, 3]], [2, [0, 2, 1, 3]]]
        sage: TpL=CountGracefulTrees(5); L=[]
        sage: for l in TpL:
        ....:     L.append(l[0])
        ....:
        sage: L
        [4, 2, 2, 4]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    # Initialization of the list collecting graceful permutations.
    L=[]
    # Loop going through the graceful permutations.
    for p in GracefulPermutations(sz):
        # Initialization of the list of functions.
        c=[]
        for j in range(sz):
            if j == 0:
                if len(c)==0:
                    c.append([(j, p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, p[j])]
            # testing that only the first criteria is met
            elif (j > 0) and (p[j] < sz-j) and not (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j+p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
            # testing that only the second criteria is met
            elif (j > 0) and not (p[j] < sz-j) and (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j-p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-p[j])]
            # testing that all two criterias are met
            elif (j > 0) and (p[j] < sz-j) and (j >= p[j]):
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

def CountGracefulTreesII(sz,k):
    """
    Goes through all the permutation of sz > 1 elements and enumerates 
    graceful trees on the vertices derived from graceful permutations.


    EXAMPLES:

    ::


        sage: CountGracefulTreesII(4,0)
        [[2, [0, 1, 2, 3]], [2, [0, 2, 1, 3]]]
        sage: TpL=CountGracefulTreesII(5,0); L=[]
        sage: for l in TpL:
        ....:     L.append(l[0])
        ....:
        sage: L
        [4, 2, 2, 4]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    # Initialization of the list collecting graceful permutations.
    L=[]
    # Loop going through the graceful permutations.
    for p in GracefulPermutationsII(sz,k):
        # Initialization of the list of functions.
        c=[]
        for j in range(sz):
            if j == 0:
                if len(c)==0:
                    c.append([(j, abs(p[j]-k))])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, abs(p[j]-k))]
            # testing that only the first criteria is met
            #elif (j > 0) and (p[j]-k < sz-j) and not (j >= p[j]-k):
            elif (j > 0) and (abs(p[j]-k) < sz-j) and not (j >= abs(p[j]-k)):
                if len(c)==0:
                    #c.append([(j, j+p[j]-k)])
                    c.append([(j, j+abs(p[j]-k))])
                    #print c
                else:
                    for i in range(len(c)):
                        #c[i]=c[i]+[(j, j+p[j]-k)]
                        c[i]=c[i]+[(j, j+abs(p[j]-k))]
                        #print c
            # testing that only the second criteria is met
            #elif (j > 0) and not (p[j]-k < sz-j) and (j >= p[j]-k):
            elif (j > 0) and not (abs(p[j]-k) < sz-j) and (j >= abs(p[j]-k)):
                if len(c)==0:
                    #c.append([(j, j-p[j]+k])
                    c.append([(j, j-abs(p[j]-k))])
                    #print c
                else:
                    for i in range(len(c)):
                        #c[i]=c[i]+[(j, j-p[j]+k)]
                        c[i]=c[i]+[(j, j-abs(p[j]-k))]
                        #print c
            # testing that all two criterias are met
            #elif (j > 0) and (p[j]-k < sz-j) and (j >= p[j]-k):
            elif (j > 0) and (abs(p[j]-k) < sz-j) and (j >= abs(p[j]-k)):
                if len(c)==0:
                    #c.append([(j, j+p[j]-k)])
                    c.append([(j, j+abs(p[j]-k))])
                    #c.append([(j, j-p[j]+k)])
                    c.append([(j, j-abs(p[j]-k))])
                else:
                    d=copy(c)
                    for i in range(len(c)):
                        #c[i]=c[i]+[(j, j+p[j]-k)]
                        c[i]=c[i]+[(j, j+abs(p[j]-k))]
                    for i in range(len(d)):
                        #d[i]=d[i]+[(j, j-p[j]+k)]
                        d[i]=d[i]+[(j, j-abs(p[j]-k))]
                    c=c+d
        # Initializing the list which will store the generators
        g = 0
        for il in c:
            #print 'c=',c
            #print 'il=',il
            #print sum([Id[:,t[0]]*Id[t[1],:]+Id[:,t[1]]*Id[t[0],:] for t in il])
            #print '\n'
            # Initilization of the adjacency matrix
            if is_Tree( sum([Id[:,t[0]]*Id[t[1],:]+Id[:,t[1]]*Id[t[0],:] for t in il]) ):
                g = g + 1
        L.append([g, p])
    return L

@cached_function
def SignedPermutations(sz):
    """
    Goes through all permutations of sz > 1 elements and outputs
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
    for p in GracefulPermutations(sz):
        # Initialization of the list of functions.
        c=[]
        for j in range(sz):
            if j == 0:
                if len(c) == 0:
                    c.append([p[j]])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[p[j]]
            # Testing that only the first gracefulness criteria is met.
            elif (j > 0) and (p[j] < (sz-j)) and not (j >= p[j]):
                if len(c) == 0:
                    c.append([p[j]])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[p[j]]
            # Testing that only the second gracefulness criteria is met.
            elif (j > 0) and not (p[j] < (sz-j)) and (j >= p[j]):
                if len(c) == 0:
                    c.append([-p[j]])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[-p[j]]
            # Testing that both gracefulness criterias are met.
            elif (j > 0) and (p[j] < (sz-j)) and (j >= p[j]):
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

def FunctionalMonomialGracefulTreeTerms(sz,A):
    """
    Outputs the sum of monomials associated with functional trees
    rooted at 0. The input is an HM


    EXAMPLES:

    ::

        sage: sz=4; f=FunctionalMonomialGracefulTreeTerms(sz,HM(sz,sz,'a')); f
        a00*a10*a20*a30 + a00*a12*a20*a30 + a00*a13*a21*a30 + a00*a13*a23*a30


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the Identity matrix
    Id=identity_matrix(sz); SP=SignedPermutations(sz)
    # Initialization of the matrix
    dt=0
    for sigma in SP:
        M=Matrix2HM(sum(A[i,i+sigma[i]]*Id[:,i]*Id[i+sigma[i],:] for i in rg(sz)))
        # Initialization of the directed Laplcian
        Ml=HM(2,(M*HM(sz,1,'one')).list(),'diag')-M
        dt=dt+A[0,0]*Deter(Ml.matrix()[1:,1:])
    return dt

def FunctionalMonomialGracefulGraphTerms(sz,A):
    """
    Outputs the sum of monomials associated with functional trees
    rooted at 0. The input is an HM


    EXAMPLES:

    ::

        sage: sz=4; FunctionalMonomialGracefulGraphTerms(sz,HM(sz,sz,'a'))
        a00*a10*a20*a30 + a00*a12*a20*a30 + a00*a13*a21*a30 + a00*a13*a23*a30


    AUTHORS:
    - Edinah K. Gnang
    """
    return sum(prod(A[i,i+sigma[i]] for i in rg(sz)) for sigma in SignedPermutations(sz))

@cached_function
def SignedPermutationsII(sz,k):
    """
    Goes through all permutations of sz > 1 elements and outputs
    the list of signed permutations associated with graceful graphs.
    The difference with the implementation above is that the gracefulness
    criteria is associated with the labeling of star graph where the label k
    is placed at the center of the star tree.


    EXAMPLES:

    ::

        sage: L=SignedPermutationsII(4,0)
        sage: L
        [[0, 1, -2, -3], [0, -1, -2, -3], [0, 2, 1, -3], [0, 2, -1, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of graceful permutations.
    L=[]
    # Loop going through graceful permutations.
    for p in GracefulPermutationsII(sz,k):
        # Initialization of the list of functions.
        c=[]
        for j in range(sz):
            if j == 0:
                if len(c) == 0:
                    c.append([abs(p[j]-k)])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[abs(p[j]-k)]
            # Testing that only the first gracefulness criteria is met.
            elif (j > 0) and (abs(p[j]-k) < (sz-j)) and not (j >= abs(p[j]-k)):
                if len(c) == 0:
                    c.append([abs(p[j]-k)])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[abs(p[j]-k)]
            # Testing that only the second gracefulness criteria is met.
            elif (j > 0) and not (abs(p[j]-k) < (sz-j)) and (j >= abs(p[j]-k)):
                if len(c) == 0:
                    c.append([-abs(p[j]-k)])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[-abs(p[j]-k)]
            # Testing that both gracefulness criterias are met.
            elif (j > 0) and (abs(p[j]-k) < (sz-j)) and (j >= abs(p[j]-k)):
                if len(c) == 0:
                    c.append([abs(p[j]-k)])
                    c.append([-abs(p[j]-k)])
                else:
                    d=copy(c)
                    for i in range(len(c)):
                        c[i]=c[i]+[abs(p[j]-k)]
                    for i in range(len(d)):
                        d[i]=d[i]+[-abs(p[j]-k)]
                    c=c+d
        L.append([c,p])
    RsltL=[]
    for l in L:
        RsltL=RsltL+l[0]
    return RsltL

@cached_function
def SignedPermutationsAlphaI(sz):
    """
    Goes through all permutations of sz > 1 elements and outputs
    the list of signed permutations associated with graceful graphs.

    EXAMPLES:

    ::

        sage: L=SignedPermutationsAlphaI(4)
        sage: L
        [[[0, 1, -2, -3], 1],
         [[0, -1, -2, -3], 1],
         [[0, 2, 1, -3], y^8],
         [[0, 2, -1, -3], y^8]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of variables
    x,y=var('x,y')
    return [[sigma, prod(prod(y^(sz^j) for j in rg(min(i, abs(sigma[i])), max(i, abs(sigma[i])))) for i in rg(sz))] for sigma in SignedPermutations(sz)]

def SignedPermutationsAlphaII(sz):
    """
    Goes through all permutations of sz > 1 elements and outputs
    the list of signed permutations associated with graceful graphs.

    EXAMPLES:

    ::

        sage: L=SignedPermutationsAlphaII(4)
        sage: L
        [[[0, 1, -2, -3], y0^2*y1^3*y2],
         [[0, -1, -2, -3], y0^3*y1^2*y2],
         [[0, 2, 1, -3], y0*y1^2*y2^3],
         [[0, 2, -1, -3], y0*y1^3*y2^2]]        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of variables
    X=var_list('x',sz); Y=var_list('y',sz)
    return [[sigma, prod(prod(Y[j] for j in rg(min(i, i+sigma[i]), max(i, i+sigma[i]))) for i in rg(sz))] for sigma in SignedPermutations(sz)]

def GraphSignedPermutationOrbitsSPIII(T):
    """
    Obtain the orbits of signed function defined by
    the isomorphism classes of associated trees. The tree
    isomorphism function is doing the heavy lifting here.
    The functions takes as input a tuple edge list.


    EXAMPLES:

    ::


        sage: GraphSignedPermutationOrbitsSPIII([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)])
        [[0, 1, 2, -3, -4], [0, 2, -1, -3, -4]]


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
        if SP2DiGraph(s).is_isomorphic(T2DiGraph(T)) and not s in cL:
            cL.append(s)
    return cL

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

def GeneralGraphSignedPermutationOrbitsT(T,j):
    """
    Obtain the orbits of signed function defined by
    the isomorphism classes of associated trees. The tree
    isomorphism function is doing the heavy lifting here.
    The functions takes as input a tuple edge list.


    EXAMPLES:

    ::


        sage: GeneralGraphSignedPermutationOrbitsT([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],0)
        [[(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0), (4, 0)],
         [(0, 0), (1, 4), (2, 3), (3, 1), (4, 0)],
         [(0, 0), (1, 4), (2, 0), (3, 2), (4, 0)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=SignedPermutationsII(len(T),j)
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

def GeneralGraphSignedPermutationOrbitsTII(T,j):
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


        sage: GeneralGraphSignedPermutationOrbitsTII([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],0)
        [[(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0), (4, 0)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=SignedPermutationsII(len(T),j)
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

def GeneralGraphSignedPermutationOrbitsSP(sigma,j):
    """
    Obtain the orbits of signed function defined by
    the isomorphism classes of associated trees. The tree
    isomorphism function is doing the heavy lifting here.
    The functions takes as input a signed permutation.


    EXAMPLES:

    ::


        sage: GeneralGraphSignedPermutationOrbitsSP([0, -1, -2, -3],0)
        [[0, -1, -2, -3], [0, 2, 1, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=SignedPermutationsII(len(sigma),j)
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

def GeneralGraphSignedPermutationOrbitsSPII(sigma):
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


        sage: GeneralGraphSignedPermutationOrbitsSPII([0, -1, -2, -3],0)
        [[0, -1, -2, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=SignedPermutationsII(len(sigma),j)
    # Initialization of the list storing the equivalence class of trees.
    cL=[]
    # Loop perfomring the isomorphism binning.
    for s in L:
        if T2DiGraph(SP2T(s)).is_isomorphic(T2DiGraph(SP2T(sigma))) and not s in cL:
            cL.append(s)
    return cL

@cached_function
def GracefulFunctionsTuples(sz):
    """
    Goes through all permutations of sz > 1 elements and outputs the
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
    for p in GracefulPermutations(sz):
        # Initialization of the list of functions.
        c=[]
        for j in range(sz):
            if j == 0:
                if len(c)==0:
                    c.append([(j, p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, p[j])]
            # Testing that only the first criteria is met.
            elif (j > 0) and (p[j] < sz-j) and not (j >= p[j]):
                if len(c) == 0:
                    c.append([(j, j+p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+p[j])]
            # Testing that only the second criteria is met.
            elif (j > 0) and not (p[j] < sz-j) and (j >= p[j]):
                if len(c)==0:
                    c.append([(j, j-p[j])])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-p[j])]
            # Testing that all two criterias are met
            elif (j > 0) and (p[j] < sz-j) and (j >= p[j]):
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
def GracefulFunctionsTuplesII(sz,k):
    """
    Goes through all permutations of sz > 1 elements and outputs the
    list of directed edges as pairs of vertices derived from graceful
    permutations. The difference with the previous implementations
    comes from the fact that we are allowing vertices other the 0
    to be placed at the center of the star tree model.


    EXAMPLES:

    ::


        sage: GracefulFunctionsTuplesII(3,0)
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
    for p in GracefulPermutationsII(sz,k):
        # Initialization of the list of functions.
        c=[]
        for j in range(sz):
            if j == 0:
                if len(c)==0:
                    c.append([(j, abs(p[j]-k))])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, abs(p[j]-k))]
            # Testing that only the first criteria is met.
            elif (j > 0) and (abs(p[j]-k) < sz-j) and not (j >= abs(p[j]-k)):
                if len(c) == 0:
                    c.append([(j, j+abs(p[j]-k))])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+abs(p[j]-k))]
            # Testing that only the second criteria is met.
            elif (j > 0) and not (abs(p[j]-k) < sz-j) and (j >= abs(p[j]-k)):
                if len(c)==0:
                    c.append([(j, j-abs(p[j]-k))])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-abs(p[j]-k))]
            # Testing that all two criterias are met
            elif (j > 0) and (abs(p[j]-k) < sz-j) and (j >= abs(p[j]-k)):
                if len(c)==0:
                    c.append([(j, j+abs(p[j]-k))])
                    c.append([(j, j-abs(p[j]-k))])
                else:
                    d=copy(c)
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+abs(p[j]-k))]
                    for i in range(len(d)):
                        d[i]=d[i]+[(j, j-abs(p[j]-k))]
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
        [[(0, 0), (1, 2), (2, 0)],
         [(0, 0), (1, 0), (2, 0)],
         [(0, 0), (2, 1), (1, 0)],
         [(0, 0), (2, 0), (1, 0)],
         [(1, 1), (0, 2), (2, 1)],
         [(1, 1), (0, 1), (2, 1)],
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
        [0, -x1 + x2, x0 - x2]


    AUTHORS:
    - Edinah K. Gnang
    """
    # returning of the list
    return [var(c+str(t[1]))-var(c+str(t[0])) for t in Lt]

def Tuple2EdgeListII(Lt, c):
    """
    Goes through graceful tuple list and produces a list of
    polynomials in the variables x associated with the tree.
    

    EXAMPLES:
    ::


        sage: Tuple2EdgeListII([(0, 0), (1, 2), (2, 0)], 'x')
        [1, x1/x2, x2/x0]    


    AUTHORS:
    - Edinah K. Gnang
    """
    # returning of the list
    return [1]+[var(c+str(t[0]))/var(c+str(t[1])) for t in Lt[1:]]

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

def is_TreeII(A):
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
    # Initialization of the dierected Laplacian
    lA=diagonal_matrix((A*ones_matrix(A.nrows(),1)).list())-A
    # Initialization of the list of sumbratrices
    Lmtr=[Matrix(ZZ,sz-1,sz-1,[lA[u,v] for u in range(sz) for v in range(sz) if u!=t if v!=t]) for t in range(sz)]
    if sum(A[t,t]*Lmtr[t].det() for t in range(sz)) == 1:
        return True
    else:
        return False

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
def TreeSignedPermutationsII(sz,j):
    """
    Goes through all the permutation of sz > 1 elements and outputs
    the list of of signed permutation associated with trees.
    the input j is associated with the shift to the permutation
    in the permutation encoding not to be confused with the root.


    EXAMPLES:

    ::


        sage: TreeSignedPermutationsII(4,0)
        [[0, 1, -2, -3], [0, -1, -2, -3], [0, 2, 1, -3], [0, 2, -1, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of function Tulples
    TpL=GracefulTreeFunctionsTuplesII(sz,j)
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
    the isomorphism classes of associated graphs. The tree
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

def GraphTupleClasses(n):
    """
    Obtain the equivalence classes of tuple function defined by
    the isomorphism classes of associated trees. The graph
    isomorphism function is doing the heavy lifting here.


    EXAMPLES:

    ::


        sage: GraphTupleClasses(4)
        [[[(0, 0), (1, 2), (2, 0), (3, 0)], [(0, 0), (1, 3), (2, 1), (3, 0)]],
         [[(0, 0), (1, 0), (2, 0), (3, 0)], [(0, 0), (1, 3), (2, 3), (3, 0)]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of signed permutations equivalence classes
    L=GraphSignedPermutationClasses(n)
    # Initialization the final list
    fL=[]
    for l in L:
        fL.append([SP2T(sigma) for sigma in l])
    return fL

def GracefulGraphPermutationClasses(sz):
    """
    Obtain the equivalence classes of permutation graphs 
    defined by the isomorphism classes of associated trees.
    The tree isomorphism function is doing the heavy lifting here.


    EXAMPLES:

    ::


        sage: len(GracefulGraphPermutationClasses(4))
        2


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
    Obtains representative equivalence classes of signed 
    function defined by the isomorphism classes of associated
    trees. The tree isomorphism function is doing the heavy
    lifting here. The difference with the function above is
    that it stores only one graph per equivalence classes


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

def RepresentativeGraphSignedPermutationClasses(sz,j):
    """
    Obtains representative equivalence classes of signed 
    function defined by the isomorphism classes of associated
    trees. The tree isomorphism function is doing the heavy
    lifting here. The difference with the function above is
    that it stores only one graph per equivalence classes


    EXAMPLES:

    ::


        sage: RepresentativeGraphSignedPermutationClasses(4,0)
        [[0, 1, -2, -3], [0, -1, -2, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=SignedPermutationsII(sz,j)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism test
    for s in L:
        nwG=True
        for i in range(len(cL)):
            if T2Graph([(j,j+s[j]) for j in range(sz)]).is_isomorphic(T2Graph([(j,j+cL[i][j]) for j in range(sz)])):
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


        sage: len(GracefulTreePermutationClasses(4))
        2


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
    Obtains representative equivalence classes of signed 
    function defined by the isomorphism classes of 
    associated trees. The tree isomorphism function is
    doing the heavy lifting here. The difference with
    the implementation above is that we only store one
    tree per equivalence class. The isomorphism class here
    is that of a functional graph. The loop edge is removed
    when testing isomorphism.


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

def FunctionalTreeSignedPermutationClassesII(sz):
    """
    Obtains representative equivalence classes of signed 
    function defined by the isomorphism classes of 
    associated trees. The tree isomorphism function is
    doing the heavy lifting here. The difference with
    the implementation above is that we only store one
    tree per equivalence class. The isomorphism class here
    is that of a functional trees. The loop edges is accounted
    for when considring the automorphism group.


    EXAMPLES:

    ::


        sage: FunctionalTreeSignedPermutationClassesII(4)
        [[0, 1, -2, -3], [0, -1, -2, -3], [0, 2, 1, -3], [0, 2, -1, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=TreeSignedPermutations(sz)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2DiGraphII([(j,j+s[j]) for j in range(sz)],sz).is_isomorphic(T2DiGraphII([(j,j+cL[i][j]) for j in range(sz)],sz)):
                nwT=False
                break
        if nwT==True:
            cL.append(s)
    return cL

def GeneralTreeSignedPermutationClassesII(sz,j):
    """
    Obtains representative equivalence classes of signed 
    function defined by the isomorphism classes of 
    associated trees. The tree isomorphism function is
    doing the heavy lifting here. The difference with
    the implementation above is that we only store one
    tree per equivalence class


    EXAMPLES:

    ::


        sage: GeneralTreeSignedPermutationClassesII(4,0)
        [[0, 1, -2, -3], [0, -1, -2, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=TreeSignedPermutationsII(sz,j)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2Graph([(j,j+s[j]) for j in range(sz)]).is_isomorphic(T2Graph([(j,j+cL[i][j]) for j in range(sz)])):
                nwT=False
                break
        if nwT==True:
            cL.append(s)
    return cL

def RepresentativeTreeSignedPermutationClasses(sz,j):
    """
    Obtains representative equivalence classes of signed 
    function defined by the isomorphism classes of 
    associated trees. The tree isomorphism function is
    doing the heavy lifting here. The difference with
    the implementation above is that we only store one
    tree per equivalence class


    EXAMPLES:

    ::


        sage: RepresentativeTreeSignedPermutationClasses(4,0)
        [[0, 1, -2, -3], [0, -1, -2, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=TreeSignedPermutationsII(sz,j)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2Graph([(j,j+s[j]) for j in range(sz)]).is_isomorphic(T2Graph([(j,j+cL[i][j]) for j in range(sz)])):
                nwT=False
                break
        if nwT==True:
            cL.append(s)
    return cL

def TreeSignedPermutationClassesIII(n):
    """
    Obtains representative equivalence classes of signed 
    function defined by the isomorphism classes of 
    associated trees. The tree isomorphism function is
    doing the heavy lifting here. The difference with
    the implementation above is that we only store one
    tree per equivalence class and ignores the loop edges.


    EXAMPLES:

    ::


        sage: TreeSignedPermutationClassesIII(4)
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
            if T2GraphII([(j,j+s[j]) for j in range(n) if j!=j+s[j]],n).is_isomorphic(T2GraphII([(j,j+cL[i][j]) for j in range(n) if j!=j+s[j]],n)):
                nwT=False
                break
        if nwT==True:
            cL.append(s)
    return cL

def RepresentativeTreeSignedPermutationClassesII(sz,k):
    """
    Obtains representative equivalence classes of signed 
    function defined by the isomorphism classes of 
    associated trees. The tree isomorphism function is
    doing the heavy lifting here. The difference with
    the implementation above is that we only store one
    tree per equivalence class and does not ignores the
    loop edges.


    EXAMPLES:

    ::


        sage: RepresentativeTreeSignedPermutationClassesIII(4,0)
        [[0, 1, -2, -3], [0, -1, -2, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=TreeSignedPermutations(sz)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2DiGraphII([(j,j+s[j]) for j in range(sz)],sz).is_isomorphic(T2DiGraphII([(j,j+cL[i][j]) for j in range(sz)],sz)):
                nwT=False
                break
        if nwT==True:
            cL.append(s)
    return [T2SP(switch_sink(SP2T(sg),k)) for sg in cL]

def RepresentativeTreeTupleClasses(sz,k):
    """
    Obtains representative equivalence classes of tuple
    function rooted at j defined by the isomorphism classes of 
    associated trees. The tree isomorphism function is
    doing the heavy lifting here. The difference with
    the implementation above is that we only store one
    tree per equivalence class and does not ignores the
    loop edges.


    EXAMPLES:

    ::


        sage: RepresentativeTreeTupleClasses(4,0)
        [[(0, 0), (1, 2), (2, 0), (3, 0)],
         [(0, 0), (1, 0), (2, 0), (3, 0)],
         [(0, 0), (1, 3), (2, 3), (3, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0)]]        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=TreeSignedPermutations(sz)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2DiGraphII([(j,j+s[j]) for j in range(sz)],sz).is_isomorphic(T2DiGraphII([(j,j+cL[i][j]) for j in range(sz)],sz)):
                nwT=False
                break
        if nwT==True:
            cL.append(s)
    return [switch_sink(SP2T(sg),k) for sg in cL]

def HiddenNonTreeFunctionsTupleClasses(n):
    """
    Obtains representative equivalence classes of 
    function tuples defined by the isomorphism classes of 
    associated trees. The tree isomorphism function is
    doing the heavy lifting here. The difference with
    the implementation above is that we only store one
    tree per equivalence class. The number of vertices
    must be greater or equal to 7


    EXAMPLES:

    ::


        sage: HiddenNonTreeFunctionsTupleClasses(7)
        [[(2, 2), (4, 1), (5, 3), (1, 0), (0, 4), (6, 2), (3, 2)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations of elements from 1 to (n-1).
    P=Permutations(n)
    # Initialization of the list of permutations
    L=NonTreeSignedPermutations(n)
    tmp=P.random_element()
    # fixing up the permutation
    q=[tmp[i]-1 for i in rg(n)]
    tmptp=SP2T(L[0])
    # Initialization of the list storing the equivalence class of trees.
    cL=[ [(q[t[0]], q[t[1]]) for t in tmptp]  ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2GraphII([(j,j+s[j]) for j in range(n)],n).is_isomorphic(T2GraphII(cL[i],n)):
                nwT=False
                break
        if nwT==True:
            # Initialization of the random permutation
            tmp=P.random_element()
            # fixing up the permutation
            q=[tmp[i]-1 for i in rg(n)]
            tmptp=SP2T(s)
            cL.append([(q[t[0]], q[t[1]]) for t in tmptp])
    return cL

def HiddenTreeFunctionsTupleClasses(n):
    """
    Obtains representative equivalence classes of 
    function tuples defined by the isomorphism classes of 
    associated trees. The tree isomorphism function is
    doing the heavy lifting here. The difference with
    the implementation above is that we only store one
    tree per equivalence class


    EXAMPLES:

    ::


        sage: HiddenTreeFunctionsTupleClasses(4)
        [[(2, 2), (1, 0), (0, 2), (3, 2)],
         [(3, 3), (2, 3), (0, 3), (1, 3)],
         [(1, 1), (2, 3), (0, 3), (3, 1)],
         [(3, 3), (0, 1), (2, 0), (1, 3)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations of elements from 1 to (n-1).
    P=Permutations(n)
    # Initialization of the list of permutations
    L=TreeSignedPermutations(n)
    tmp=P.random_element()
    # fixing up the permutation
    q=[tmp[i]-1 for i in rg(n)]
    tmptp=SP2T(L[0])
    # Initialization of the list storing the equivalence class of trees.
    cL=[ [(q[t[0]], q[t[1]]) for t in tmptp]  ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2GraphII([(j,j+s[j]) for j in range(n)],n).is_isomorphic(T2GraphII(cL[i],n)):
                nwT=False
                break
        if nwT==True:
            # Initialization of the random permutation
            tmp=P.random_element()
            # fixing up the permutation
            q=[tmp[i]-1 for i in rg(n)]
            tmptp=SP2T(s)
            cL.append([(q[t[0]], q[t[1]]) for t in tmptp])
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

@cached_function
def GracefulNonTreeFunctionsTuples(n):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the list of graceful functions on the vertices derived from 
    graceful permutations.


    EXAMPLES:

    ::


        sage: GracefulNonTreeFunctionsTuples(7)
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
        sage: TpL=GracefulNonTreeFunctionsTuples(7)
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

@cached_function
def GracefulTreeFunctionsTuplesII(sz,k):
    """
    Goes through all the permutation of sz > 1 elements and outputs
    the list of graceful functions on the vertices 
    derived from graceful permutations.


    EXAMPLES:

    ::


        sage: GracefulTreeFunctionsTuplesII(4,0)
        [[[[(0, 0), (1, 2), (2, 0), (3, 0)], [(0, 0), (1, 0), (2, 0), (3, 0)]],
          [0, 1, 2, 3]],
        [[[(0, 0), (1, 3), (2, 3), (3, 0)], [(0, 0), (1, 3), (2, 1), (3, 0)]],
          [0, 2, 1, 3]]]
        sage: L=[]
        sage: TpL=GracefulTreeFunctionsTuplesII(4,0)
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
    Id=identity_matrix(sz)
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Loop going through the GracefulPermutations.
    for p in GracefulPermutationsII(sz,k):
        # Initialization of the list of functions.
        c=[]
        for j in range(sz):
            if j == 0:
                if len(c)==0:
                    c.append([(j, abs(p[j]-k))])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, abs(p[j]-k))]
            # Testing that only the first criteria is met
            elif (j > 0) and (abs(p[j]-k)<(sz-j)) and not (j>=abs(p[j]-k)):
                if len(c)==0:
                    c.append([(j, j+abs(p[j]-k))])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+abs(p[j]-k))]
            # Testing that only the second criteria is met
            elif (j > 0) and not (abs(p[j]-k) < (sz-j)) and (j >= abs(p[j]-k)):
                if len(c)==0:
                    c.append([(j, j-abs(p[j]-k))])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j-abs(p[j]-k))]
            # Testing that all two criterias are met
            elif (j > 0) and (abs(p[j]-k) < (sz-j)) and (j >= abs(p[j]-k)):
                if len(c)==0:
                    c.append([(j, j+abs(p[j]-k))])
                    c.append([(j, j-abs(p[j]-k))])
                else:
                    d=copy(c)
                    for i in range(len(c)):
                        c[i]=c[i]+[(j, j+abs(p[j]-k))]
                    for i in range(len(d)):
                        d[i]=d[i]+[(j, j-abs(p[j]-k))]
                    c=c+d
        # Initializing the list which will store the generators
        g=[]
        for il in c:
            # Initilization of the adjacency matrix
            if is_Tree( sum([Id[:,t[0]]*Id[t[1],:]+Id[:,t[1]]*Id[t[0],:] for t in il]) ):
                g.append(il)
        L.append([g,p])
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
    and outputs the list of directed graceful graphs.


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
    and outputs the list of directed graceful graphs.


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

def GracefulDiNonTreeList(n):
    """
    Goes through all the permutation of graceful second index functions
    and outputs the list of directed adjacency matrices of graceful 
    graphs.

    EXAMPLES:
    ::


        sage: for l in GracefulDiNonTreeList(5):
        ....:     t=0
        ....:     for dg in l[0]:
        ....:         dg.plot().save(str(l[1]).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
        ....:         t=t+1
        ....:


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GracefulNonTreeFunctionsTuples(n)
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

def SignedPermutationListDegreeSequence(sp):
    """
    The method returns the degree sequence of a an input edge list.
    the inputs to this function is list of signed permutations 
    specifying an edge list and a positive integer corresponding to
    the number of vertices in the graph. The loop edges is ingnored
    by this implementation.


    EXAMPLES:

    ::

        sage: SignedPermutationListDegreeSequence([0, 1, 2, -3, -4])
        [2, 1, 2, 1, 2]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    return EdgeListDegreeSequence(SP2T(sp))

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

def T2DiGraphII(T,sz):
    """
    The method returns a directed graph object associated with 
    with the tuple list description of the directed graph

    EXAMPLES:

    ::

        sage: T2DiGraphII([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],5).degree_sequence()
        [4, 2, 2, 1, 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    return DiGraph(sum([Id[:,t[0]]*Id[t[1],:] for t in T]))

def T2Adjacency(T,sz):
    """
    The method returns a directed graph object associated with 
    with the tuple list description of the directed graph

    EXAMPLES:

    ::

        sage: T2Adjacency([(0, 1), (1, 2), (2, 0), (3, 3)],4) 
        [0 1 0 0]
        [0 0 1 0]
        [1 0 0 0]
        [0 0 0 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    return (sum([Id[:,t[0]]*Id[t[1],:] for t in T]))

def T2AdjacencyII(T,sz,c='a'):
    """
    The method returns a directed graph object associated with 
    with the tuple list description of the directed graph
    the 1s are replaced with symbolic variables.

    EXAMPLES:

    ::

        sage: T2AdjacencyII([(0, 1), (1, 2), (2, 0), (3, 3)],4)
        [  0 a01   0   0]
        [  0   0 a12   0]
        [a20   0   0   0]
        [  0   0   0 a33]        
        

    AUTHORS:
    - Edinah K. Gnang
    """
    A=HM(sz,sz,c)
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    return (sum(A[t[0],t[1]]*Id[:,t[0]]*Id[t[1],:] for t in T))

def T2Incidence(T,sz):
    """
    The method returns a directed graph object associated with 
    with the tuple list description of the directed graph

    EXAMPLES:

    ::

        sage: T2Incidence([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],5)
        [ 0  0  0  1  1]
        [ 0 -1  0  0  0]
        [ 0  1 -1  0  0]
        [ 0  0  0 -1  0]
        [ 0  0  1  0 -1] 
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    B=HM(sz,len(T), 'zero')
    T.sort()
    clindx=0
    for t in T:
        B[t[0],clindx]=-1; B[t[1],clindx]=B[t[1],clindx]+1
        clindx=clindx+1
    return B.matrix()

def T2IncidenceHM(T,sz):
    """
    The method returns a directed graph object associated with 
    with the tuple list description of the directed graph

    EXAMPLES:

    ::

        sage: T2Incidence([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],5).matrix()
        [ 0  0  0  1  1]
        [ 0 -1  0  0  0]
        [ 0  1 -1  0  0]
        [ 0  0  0 -1  0]
        [ 0  0  1  0 -1] 
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    B=HM(sz,len(T), 'zero')
    T.sort()
    clindx=0
    for t in T:
        B[t[0],clindx]=-1; B[t[1],clindx]=B[t[1],clindx]+1
        clindx=clindx+1
    return B

def T2IncidenceII(T,sz):
    """
    The method returns a directed graph object associated with 
    with the tuple list description of the directed graph

    EXAMPLES:

    ::

        sage: T2IncidenceII([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],5)
        [   0    0    0  a30  a40]
        [   0 -a12    0    0    0]
        [   0  a12 -a24    0    0]
        [   0    0    0 -a30    0]
        [   0    0  a24    0 -a40]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    A=HM(sz,sz,'a')
    # Initialization of the identity matrix
    B=HM(sz,len(T), 'zero')
    T.sort()
    clindx=0
    for t in T:
        B[t[0],clindx]=-A[t[0],t[1]]; B[t[1],clindx]=B[t[1],clindx]+A[t[0],t[1]]
        clindx=clindx+1
    return B.matrix()

def T2Graph(T):
    """
    The method returns an undirected graph object associated with 
    with the tuple list description of the directed graph

    EXAMPLES:

    ::

        sage: T2Graph([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)]).degree_sequence()
        [2, 2, 2, 1, 1]
        

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
        [2, 2, 2, 1, 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the Tree tuple function
    T=[(i,i+sigma[i]) for i in range(len(sigma))]
    # Initialization of the identity matrix
    Id=identity_matrix(len(T))
    return Graph(sum(Id[:,t[0]]*Id[t[1],:] for t in T[1:])+sum(Id[:,t[0]]*Id[t[1],:] for t in T[1:]).transpose())

def SP2DiGraph(sigma):
    """
    The method returns an undirected graph object associated with 
    with the tuple list description of the directed graph

    EXAMPLES:

    ::

        sage: SP2DiGraph([(0-0), -(1-2), -(2-4), -(3-0), -(4-0)]).degree_sequence()
        [4, 2, 2, 1, 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the Tree tuple function
    T=[(i, i+sigma[i]) for i in range(len(sigma))]
    # Initialization of the identity matrix
    Id=identity_matrix(len(T))
    return DiGraph(sum(Id[:,t[0]]*Id[t[1],:] for t in T))

def TupleTreeFunctionList(sz):
    """
    Goes through all the functions and determines which ones
    are associated with trees.

    EXAMPLES:
    ::
        sage: TupleTreeFunctionList(4)
        [[(0, 0), (1, 0), (2, 0), (3, 0)],
         [(0, 0), (1, 2), (2, 0), (3, 0)],
         [(0, 0), (1, 3), (2, 0), (3, 0)],
         [(0, 0), (1, 0), (2, 1), (3, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0)],
         [(0, 0), (1, 0), (2, 3), (3, 0)],
         [(0, 0), (1, 2), (2, 3), (3, 0)],
         [(0, 0), (1, 3), (2, 3), (3, 0)],
         [(0, 0), (1, 0), (2, 0), (3, 1)],
         [(0, 0), (1, 2), (2, 0), (3, 1)],
         [(0, 0), (1, 0), (2, 1), (3, 1)],
         [(0, 0), (1, 0), (2, 3), (3, 1)],
         [(0, 0), (1, 0), (2, 0), (3, 2)],
         [(0, 0), (1, 2), (2, 0), (3, 2)],
         [(0, 0), (1, 3), (2, 0), (3, 2)],
         [(0, 0), (1, 0), (2, 1), (3, 2)]]
        sage: sz=3; X=var_list('x',sz+1); Y=var_list('y',sz+1); Z=var_list('z',sz+1); Ta=HM(sz+1,sz+1,'a'); A=HM(sz+1,sz+1,'zero'); DgA=[Ta[i,i] for i in rg(sz+1)]
        sage: for i in rg(sz+1):
        ....:     for j in rg(sz+1):
        ....:         A[i,j]=Ta[i,j]*X[abs(j-i)]*prod(Y[u] for u in rg(min(i,j), max(i,j)))*prod(Z[u] for u in rg(1+min(i,j), 1+max(i,j)))
        ....:
        sage: sum(prod(A[t[0],t[1]] for t in tp) for tp in TupleTreeFunctionList(sz))
        a00*a10*a20*x0*x1*x2*y0^2*y1*z1^2*z2 + a00*a12*a20*x0*x1*x2*y0*y1^2*z1*z2^2 + a00*a10*a21*x0*x1^2*y0*y1*z1*z2
 

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[sz for i in range(sz-1)]; Lf=[]
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=[0]+entry
        # Initialization of the adjacency martrix of the functional graph
        A = sum([Id[:,j]*Id[f[j],:] for j in range(sz)])
        # Initialization of the dierected Laplacian
        lA= diagonal_matrix((A*ones_matrix(A.nrows(),1)).list())-A
        # Testing treeness
        if Matrix(ZZ,sz-1,sz-1,[lA[u,v] for u in rg(1,sz) for v in rg(1,sz)]).det()*A[0,0] == 1:
            #Lf.append(f)
            Lf.append([(j, f[j]) for j in rg(sz)])
    return Lf

def RootedTreeFunctionList(sz):
    """
    Goes through all the functions and determines which ones
    are associated with trees.

    EXAMPLES:
    ::
        sage: RootedTreeFunctionList(3)
        [[0, 0, 0],
         [1, 1, 0],
         [0, 2, 0],
         [0, 0, 1],
         [1, 1, 1],
         [2, 1, 1],
         [2, 0, 2],
         [1, 2, 2],
         [2, 2, 2]]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Initialization of the adjacency martrix of the functional graph
        A = sum([Id[:,j]*Id[f[j],:] for j in range(sz)])
        # Initialization of the dierected Laplacian
        lA= diagonal_matrix((A*ones_matrix(A.nrows(),1)).list())-A
        # Initialization of the list of sumbratrices
        Lmtr=[Matrix(ZZ,sz-1,sz-1,[lA[u,v] for u in range(sz) for v in range(sz) if u!=t if v!=t]) for t in range(sz)]
        # Testing treeness
        if sum(A[t,t]*Lmtr[t].det() for t in range(sz)) == 1:
            # Appending the function to the list
            Lf.append(f)
    return Lf

def RootedTupleGracefulTreeFunctionList(sz):
    """
    Goes through all the functions and determines which ones
    are associated with gracefully labeled trees.

    EXAMPLES:
    ::
        sage: RootedTupleGracefulTreeFunctionList(3)
        [[(0, 0), (1, 0), (2, 0)],
         [(0, 1), (1, 1), (2, 0)],
         [(0, 0), (1, 2), (2, 0)],
         [(0, 2), (1, 1), (2, 1)],
         [(0, 2), (1, 0), (2, 2)],
         [(0, 2), (1, 2), (2, 2)]]
        sage: sz=3; X=var_list('x',sz+1); Y=var_list('y',sz+1); Z=var_list('z',sz+1); Ta=HM(sz+1,sz+1,'a'); A=HM(sz+1,sz+1,'zero'); DgA=[Ta[i,i] for i in rg(sz+1)]
        sage: for i in rg(sz+1):
        ....:     for j in rg(sz+1):
        ....:         A[i,j]=Ta[i,j]*X[abs(j-i)]*prod(Y[u] for u in rg(min(i,j), max(i,j)))*prod(Z[u] for u in rg(1+min(i,j), 1+max(i,j)))
        ....:
        sage: sum(prod(A[t[0],t[1]] for t in tp) for tp in RootedTupleGracefulTreeFunctionList(sz))
        a00*a10*a20*x0*x1*x2*y0^2*y1*z1^2*z2 + a01*a11*a20*x0*x1*x2*y0^2*y1*z1^2*z2 + a02*a10*a22*x0*x1*x2*y0^2*y1*z1^2*z2 + a00*a12*a20*x0*x1*x2*y0*y1^2*z1*z2^2 + a02*a11*a21*x0*x1*x2*y0*y1^2*z1*z2^2 + a02*a12*a22*x0*x1*x2*y0*y1^2*z1*z2^2 + a00*a10*a21*x0*x1^2*y0*y1*z1*z2 + a01*a11*a21*x0*x1^2*y0*y1*z1*z2 + a01*a12*a22*x0*x1^2*y0*y1*z1*z2


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Initialization of the adjacency martrix of the functional graph
        A = sum([Id[:,j]*Id[f[j],:] for j in range(sz)])
        # Initialization of the dierected Laplacian
        lA= diagonal_matrix((A*ones_matrix(A.nrows(),1)).list())-A
        # Initialization of the list of sumbratrices
        Lmtr=[Matrix(ZZ,sz-1,sz-1,[lA[u,v] for u in range(sz) for v in range(sz) if u!=t if v!=t]) for t in range(sz)]
        # Testing treeness
        if sum(A[t,t]*Lmtr[t].det() for t in range(sz)) == 1 and Set([abs(f[i]-i) for i in range(sz)]).cardinality() == sz:
            # Appending the function to the list
            Lf.append([(i,f[i]) for i in range(sz)])
    return Lf

def RootedTupleInducedTreeFunctionList(sz, induced_edge_label_sequence):
    """
    Goes through all the functions and determines which ones
    are associated with a given edge labeled sequence sorted
    in non decreasing order

    EXAMPLES:
    ::
        sage: RootedTupleInducedTreeFunctionList(3,[0,1,2])
        [[(0, 0), (1, 0), (2, 0)],
         [(0, 1), (1, 1), (2, 0)],
         [(0, 0), (1, 2), (2, 0)],
         [(0, 2), (1, 1), (2, 1)],
         [(0, 2), (1, 0), (2, 2)],
         [(0, 2), (1, 2), (2, 2)]]
        sage: sz=3; X=var_list('x',sz+1); Y=var_list('y',sz+1); Z=var_list('z',sz+1); Ta=HM(sz+1,sz+1,'a'); A=HM(sz+1,sz+1,'zero'); DgA=[Ta[i,i] for i in rg(sz+1)]
        sage: for i in rg(sz+1):
        ....:     for j in rg(sz+1):
        ....:         A[i,j]=Ta[i,j]*X[abs(j-i)]*prod(Y[u] for u in rg(min(i,j), max(i,j)))*prod(Z[u] for u in rg(1+min(i,j), 1+max(i,j)))
        ....:
        sage: sum(prod(A[t[0],t[1]] for t in tp) for tp in RootedTupleInducedTreeFunctionList(sz,[0,1,2]))
        a00*a10*a20*x0*x1*x2*y0^2*y1*z1^2*z2 + a01*a11*a20*x0*x1*x2*y0^2*y1*z1^2*z2 + a02*a10*a22*x0*x1*x2*y0^2*y1*z1^2*z2 + a00*a12*a20*x0*x1*x2*y0*y1^2*z1*z2^2 + a02*a11*a21*x0*x1*x2*y0*y1^2*z1*z2^2 + a02*a12*a22*x0*x1*x2*y0*y1^2*z1*z2^2 + a00*a10*a21*x0*x1^2*y0*y1*z1*z2 + a01*a11*a21*x0*x1^2*y0*y1*z1*z2 + a01*a12*a22*x0*x1^2*y0*y1*z1*z2


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Initialization of the adjacency martrix of the functional graph
        A = sum([Id[:,j]*Id[f[j],:] for j in range(sz)])
        # Initialization of the dierected Laplacian
        lA= diagonal_matrix((A*ones_matrix(A.nrows(),1)).list())-A
        # Initialization of the list of sumbratrices
        Lmtr=[Matrix(ZZ,sz-1,sz-1,[lA[u,v] for u in range(sz) for v in range(sz) if u!=t if v!=t]) for t in range(sz)]
        # Testing treeness
        EdgLblSeq=[abs(f[i]-i) for i in range(sz)]; EdgLblSeq.sort()
        if sum(A[t,t]*Lmtr[t].det() for t in range(sz)) == 1 and  EdgLblSeq == induced_edge_label_sequence:
            # Appending the function to the list
            Lf.append([(i,f[i]) for i in range(sz)])
    return Lf

def RootedTupleGracefulFunctionList(sz):
    """
    Goes through all the functions and determines which ones
    are associated with gracefully labeled functional directed graphs.

    EXAMPLES:
    ::
        sage: RootedTupleGracefulFunctionList(3)
        [[(0, 0), (1, 0), (2, 0)],
         [(0, 1), (1, 1), (2, 0)],
         [(0, 0), (1, 2), (2, 0)],
         [(0, 0), (1, 0), (2, 1)],
         [(0, 1), (1, 1), (2, 1)],
         [(0, 2), (1, 1), (2, 1)],
         [(0, 2), (1, 0), (2, 2)],
         [(0, 1), (1, 2), (2, 2)],
         [(0, 2), (1, 2), (2, 2)]]
        sage: sz=3; X=var_list('x',sz+1); Y=var_list('y',sz+1); Z=var_list('z',sz+1); Ta=HM(sz+1,sz+1,'a'); A=HM(sz+1,sz+1,'zero'); DgA=[Ta[i,i] for i in rg(sz+1)]
        sage: for i in rg(sz+1):
        ....:     for j in rg(sz+1):
        ....:         A[i,j]=Ta[i,j]*X[abs(j-i)]*prod(Y[u] for u in rg(min(i,j), max(i,j)))*prod(Z[u] for u in rg(1+min(i,j), 1+max(i,j)))
        ....:
        sage: sum(prod(A[t[0],t[1]] for t in tp) for tp in RootedTupleGracefulFunctionList(sz))
        a00*a10*a20*x0*x1*x2*y0^2*y1*z1^2*z2 + a01*a11*a20*x0*x1*x2*y0^2*y1*z1^2*z2 + a02*a10*a22*x0*x1*x2*y0^2*y1*z1^2*z2 + a00*a12*a20*x0*x1*x2*y0*y1^2*z1*z2^2 + a02*a11*a21*x0*x1*x2*y0*y1^2*z1*z2^2 + a02*a12*a22*x0*x1*x2*y0*y1^2*z1*z2^2 + a00*a10*a21*x0*x1^2*y0*y1*z1*z2 + a01*a11*a21*x0*x1^2*y0*y1*z1*z2 + a01*a12*a22*x0*x1^2*y0*y1*z1*z2


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Testing Gracefulness
        if Set([abs(f[i]-i) for i in range(sz)]).cardinality() == sz:
            # Appending the function to the list
            Lf.append([(i,f[i]) for i in range(sz)])
    return Lf

def NonIsomorphicRootedTreeFunctionList(sz):
    """
    Obtains representative equivalence classes of tuple
    function rooted.


    EXAMPLES:

    ::


        sage: NonIsomorphicRootedTreeFunctionList(4)
        [[(0, 0), (1, 0), (2, 0), (3, 0)],
         [(0, 1), (1, 1), (2, 0), (3, 0)],
         [(0, 0), (1, 2), (2, 0), (3, 0)],
         [(0, 2), (1, 1), (2, 1), (3, 0)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=[T2SP(T) for T in RootedTupleTreeFunctionList(sz)]
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2DiGraphII([(j,j+s[j]) for j in range(sz)],sz).is_isomorphic(T2DiGraphII([(j,j+cL[i][j]) for j in range(sz)],sz)):
                nwT=False
                break
        if nwT==True:
            cL.append(s)
    return [SP2T(sg) for sg in cL]

def GracefulTupleFunctionList(sz):
    """
    Returns list of edge tuple desctiption for all 
    functional derected graphs.


    EXAMPLES:
    ::
        sage:GracefulTupleFunctionList(3)
        [[(0, 0), (1, 0), (2, 0), (3, 0)],
         [(0, 1), (1, 1), (2, 0), (3, 0)],
         [(0, 0), (1, 2), (2, 0), (3, 0)],
         [(0, 2), (1, 1), (2, 1), (3, 0)],
         [(0, 0), (1, 3), (2, 1), (3, 0)],
         [(0, 2), (1, 0), (2, 2), (3, 0)],
         [(0, 2), (1, 2), (2, 2), (3, 0)],
         [(0, 1), (1, 3), (2, 2), (3, 0)],
         [(0, 2), (1, 1), (2, 3), (3, 0)],
         [(0, 0), (1, 3), (2, 3), (3, 0)],
         [(0, 3), (1, 1), (2, 1), (3, 1)],
         [(0, 3), (1, 0), (2, 2), (3, 1)],
         [(0, 3), (1, 2), (2, 2), (3, 1)],
         [(0, 3), (1, 1), (2, 3), (3, 1)],
         [(0, 3), (1, 1), (2, 0), (3, 2)],
         [(0, 3), (1, 3), (2, 2), (3, 2)],
         [(0, 3), (1, 0), (2, 0), (3, 3)],
         [(0, 3), (1, 2), (2, 0), (3, 3)],
         [(0, 3), (1, 3), (2, 1), (3, 3)],
         [(0, 3), (1, 3), (2, 3), (3, 3)]]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Appending the function to the list
        if Set([abs(j-f[j]) for j in rg(sz)]) == Set(rg(sz)):
            Lf.append([(j,f[j]) for j in rg(sz)])
    return Lf

def PseudoTreeFunctionList(n):
    """
    Goes through all the functions and determines which ones
    are such that the image of 1,2,...,n-1 has size n-1

    EXAMPLES:
    ::
        sage: PseudoTreeFunctionList(4)
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
         [2, 3, 1],
         [0, 0, 2],
         [2, 0, 2],
         [3, 0, 2],
         [0, 1, 2],
         [3, 1, 2]]


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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # testing self image
        self_image=False
        for z in rg(len(f)):
            if f[z]==z+1:
                self_image=True
                break
        for z in rg(len(f)):
            if (min(f[z],z+1), max(f[z],z+1)) in [(min(f[y],y+1), max(f[y],y+1)) for y in rg(len(f)) if y!=z] :
                self_image=True
                break

        # Appending the function to the list
        #if Set(f).cardinality() == n-1 and self_image==False:
        if self_image==False:
            Lf.append(f)
    return Lf

def NonTreeFunctionList(n):
    """
    Goes through all the functions and determines which ones
    are associated with trees.

    EXAMPLES:
    ::
        sage: NonTreeFunctionList(3)
        [[1, 0], [1, 1], [2, 1], [0, 2], [1, 2], [2, 2]]

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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Appending the function to the list
        if not is_Tree(sum([Id[:,j]*Id[f[j-1],:]+Id[:,f[j-1]]*Id[j,:] for j in range(1,n)])):
            Lf.append(f)
    return Lf

def TreeHypermatrix(sz):
    """
    Goes through all the functions and determines which ones
    are associated with trees.

    EXAMPLES:
    ::
        sage: TreeHypermatrix(2).printHM()
        [:, :]=
        [1 0]
        [0 0]
        sage: TreeHypermatrix(3).printHM()
        [:, :, 0]=
        [1 0 1]
        [0 0 0]
        <BLANKLINE>
        [:, :, 1]=
        [1 0 0]
        [0 0 0]
        <BLANKLINE>
        [:, :, 2]=
        [0 0 0]
        [0 0 0]
        <BLANKLINE>


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[sz for i in range(sz-1)]; Lf=[]
    # Initialization of the hypermatrix
    Ha=apply(HM,[2]+l+['zero'])
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Appending the function to the list
        if is_Tree(sum([Id[:,j]*Id[f[j-1],:]+Id[:,f[j-1]]*Id[j,:] for j in range(1,sz)])):
            Lf.append(f)
            Ha[tuple([0]+f)]=1
    return Ha

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

def K_arryTreeSignedPermutationClasses(n, K):
    """
    Obtains representative equivalence classes of 
    signed function of tree with vertices of degree >=K
    defined by the isomorphism classes of associated trees.
    The tree isomorphism function is doing the heavy lifting
    here. The difference with the implementation above is 
    that we only store one tree per equivalence class.


    EXAMPLES:

    ::


        sage: K_arryTreeSignedPermutationClasses(4,2)
        [[0, 1, -2, -3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=TreeSignedPermutations(n)
    # Initialization of the list storing the equivalence class of trees.
    #cL=[ L[0] ]
    if Set(SignedPermutationListDegreeSequence(L[0])).issubset(Set(range(1,K+1))): 
        cL=[ L[0] ]
    else:
        cL=[ ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2Graph([(j,j+s[j]) for j in range(n)]).is_isomorphic(T2Graph([(j,j+cL[i][j]) for j in range(n)])):
                nwT=False
                break
        if nwT==True and Set(SignedPermutationListDegreeSequence(s)).issubset(Set(range(1,K+1))):
            cL.append(s)
    return cL

def Isomorphic_SignedPermutationList(L0,L1):
    """
    Compares the input list to determine if they
    are made up of isomorphic list of signed
    permutations.

    EXAMPLES:

    ::


        sage: L0=K_arryTreeSignedPermutationClasses(4,2) 
        sage: L1=K_arryNodeSplittingSignedPermutationClasses(4,2)
        sage: Isomorphic_SignedPermutationList(L0,L1)
        True
        sage: Isomorphic_SignedPermutationList(K_arryNodeSplittingSignedPermutationClasses(4,3),K_arryTreeSignedPermutationClasses(4,3))
        True


    AUTHORS:
    - Edinah K. Gnang
    """
    # Loop perfomring the isomorphism binning.
    for s0 in L0:
        nwT=True
        for s1 in L1:
            if T2Graph([(j,j+s0[j]) for j in range(len(s0))]).is_isomorphic(T2Graph([(j,j+s1[j]) for j in range(len(s1))])):
                nwT=False
                break
        if nwT==True:
            print 'A graph isomorphic to the signed permutation ', s0,' is missing.'
            return False 
    return True

def TupleSubgraphList(T,k):
    """
    Returns a list of tuple description of subgraph on
    k vertices. Ignores loop edges.
    

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
    # Initialization of the edge set
    Edge=Set(T[1:])
    # Initialization of the list
    L=[]
    for s in Edge.subsets(k):
        L.append(s.list())
    # Sorting the list
    for i in range(len(L)):
        L[i].sort()
    return L

def TupleSubgraphListII(T,k):
    """
    Returns a list of tuple description of subgraph on
    k vertices. The difference with the implementation
    above is that it includes the loop edge
    

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
    # Initialization of the edge set
    Edge=Set(T)
    # Initialization of the list
    L=[]
    for s in Edge.subsets(k):
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
        [1, 1, 1, 1, 0]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    return Graph(sum(Id[:,t[0]]*Id[t[1],:] for t in T)+sum(Id[:,t[0]]*Id[t[1],:] for t in T).transpose())

def T2DiGraphII(T,sz):
    """
    The method returns a directed graph object associated with 
    with the tuple list description of the directed graph on sz
    vertices.

    EXAMPLES:

    ::

        sage: T2DiGraphII([(1, 2), (3, 0)],5).degree_sequence()
        [1, 1, 1, 1, 0]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    return DiGraph(sum(Id[:,t[0]]*Id[t[1],:] for t in T))

def generate_script_fast_graceful_list(sz):
    """
    The produces an implementation an optimal sage script for
    displaying and listing gracefully labled graph on sz vertices.
    sz must be greater than 4. Make sure to include the
    Hypermatrix Package in the working directory

    EXAMPLES:

    ::

        sage: from subprocess import call
        sage: call("touch fast_listing_of_graceful_trees_on_10_vertices.sage", shell=True)
        0
        sage: call("rm fast_listing_of_graceful_trees_on_10_vertices.sage", shell=True)
        0
        sage: generate_script_fast_graceful_list(10)
        sage: call("ls fast_listing_of_graceful_trees_on_10_vertices.sage", shell=True)
        fast_listing_of_graceful_trees_on_10_vertices.sage
        0
        

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

        sage: from subprocess import call
        sage: call("touch graceful_labeling_of_a_trees_on_6_vertices.sage", shell=True)
        0
        sage: call("rm graceful_labeling_of_a_trees_on_6_vertices.sage", shell=True)
        0
        sage: generate_labeling_script([(0,0),(1,2),(2,3),(3,4),(4,5),(5,0)])
        sage: call("ls graceful_labeling_of_a_trees_on_6_vertices.sage", shell=True)
        graceful_labeling_of_a_trees_on_6_vertices.sage
        0
        

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
    the working directory. The difference with the code above
    is that it does not produce png files for each tree and has
    a properly wrapped up code. I also appear to stop as soon 
    as it finds one


    EXAMPLES:

    ::

        sage: from subprocess import call
        sage: call("touch graceful_labeling_of_a_trees_on_6_vertices.sage", shell=True)
        0
        sage: call("rm graceful_labeling_of_a_trees_on_6_vertices.sage", shell=True)
        0
        sage: generate_labeling_scriptII([(0,0),(1,2),(2,3),(3,4),(4,5),(5,0)])
        sage: call("ls graceful_labeling_of_a_trees_on_6_vertices.sage", shell=True)
        graceful_labeling_of_a_trees_on_6_vertices.sage
        0
        

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

        sage: from subprocess import call
        sage: call("touch fast_listing_of_partial_graceful_trees_on_10_vertices.sage", shell=True)
        0
        sage: call("rm fast_listing_of_partial_graceful_trees_on_10_vertices.sage", shell=True)
        0
        sage: generate_script_fast_partial_graceful_list(10)
        sage: call("ls fast_listing_of_partial_graceful_trees_on_10_vertices.sage", shell=True)
        fast_listing_of_partial_graceful_trees_on_10_vertices.sage
        0
        

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

        sage: from subprocess import call
        sage: call("touch fast_listing_of_partial_graceful_trees_on_10_verticesII.sage", shell=True)
        0
        sage: call("rm fast_listing_of_partial_graceful_trees_on_10_verticesII.sage", shell=True)
        0
        sage: generate_script_fast_partial_graceful_listII(10)
        sage: call("ls fast_listing_of_partial_graceful_trees_on_10_verticesII.sage", shell=True)
        fast_listing_of_partial_graceful_trees_on_10_verticesII.sage
        0
        

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

def generate_script_graceful_trees(sz):
    """
    The produces an implementation an optimal sage script for
    displaying and listing graceful functional directed graphs
    which never cross the axis y=x.


    EXAMPLES:

    ::

        sage: from subprocess import call
        sage: call("touch fast_listing_of_graceful_trees_on_10_vertices.sage", shell=True)
        0
        sage: call("rm fast_listing_of_graceful_trees_on_10_vertices.sage", shell=True)
        0
        sage: generate_script_graceful_trees(4)
        sage: call("ls fast_listing_of_graceful_trees_on_4_vertices.sage", shell=True)
        fast_listing_of_graceful_trees_on_4_vertices.sage
        0
        sage: load('fast_listing_of_graceful_trees_on_4_vertices.sage') 
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        sage: f
        12*a00*a10*a20*a30 + 4*a00*a10*a21*a30 + 4*a00*a10*a20*a31 + 12*a00*a10*a21*a31 + 4*a00*a10*a20*a32 + 4*a00*a10*a21*a32

 

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
    f.write('# Initialization of the symbolic variables\n')
    f.write("A=HM(sz,sz,'a'); Y=var_list('y',sz); K=var_list('k',sz)\n\n")
    f.write('# Initialization of the Permutations\n')
    f.write('P=Permutations(sz)\n\n')
    f.write('# Computing the sum over gracefully labeled functional directed graphs\n')
    f.write('f=0 \n\n')
    f.write('# Main loop\n')
    f.write('for t in NonIncreasingFunctionList(sz)[0]:\n')
    # variable storing the indentation spaces
    sp = '    '
    for i in range(sz):
        tmpString='for K['+str(i)+'] in rg(2):\n'
        f.write(sp+tmpString)
        sp=sp+'    '
    f.write(sp+'if prod(Y[j-t[j][1]]^K[j]*Y[sz-1-t[j][1]]^(1-K[j]) for j in rg(sz))==prod(Y):\n')
    sp=sp+'    '
    f.write(sp+'for p in P:\n')
    sp=sp+'    '
    f.write(sp+'if prod(p[k-t[k][1]]*K[k]+p[sz-1-t[k][1]]*(1-K[k]) >= p[k-t[k][1]]*(1-K[k])+p[sz-1-t[k][1]]*K[k] for k in rg(sz)):\n')
    sp=sp+'    '
    f.write(sp+'f=f+prod(A[p[sz-1-t[i][1]]-1,p[i-t[i][1]]-1]^(1-K[i])*A[p[i-t[i][1]]-1,p[sz-1-t[i][1]]-1]^K[i] for i in rg(sz))/2\n\n')
    f.write('# Initialization of the list of graphs\n')
    f.write('LtG = [Monomial2TII(mnm, A.list(), sz) for mnm in f.operands()]\n\n')
    f.write('# Plotting the trees \n')
    f.write('for g in LtG:\n')
    f.write('    T2DiGraphII(g,sz).plot().show()\n')
    # Closing the file
    f.close()

def generate_script_graceful_treesII(sz):
    """
    The produces an implementation an optimal sage script for
    displaying and listing graceful functional directed graphs
    which never cross the axis y=x.


    EXAMPLES:

    ::

        sage: from subprocess import call
        sage: call("touch fast_listing_of_graceful_trees_on_10_vertices.sage", shell=True)
        0
        sage: call("rm fast_listing_of_graceful_trees_on_10_vertices.sage", shell=True)
        0
        sage: generate_script_graceful_trees(4)
        sage: call("ls fast_listing_of_graceful_trees_on_4_vertices.sage", shell=True)
        fast_listing_of_graceful_trees_on_4_vertices.sage
        0
        sage: load('fast_listing_of_graceful_trees_on_4_vertices.sage') 
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        sage: f
        12*a00*a10*a20*a30 + 4*a00*a10*a21*a30 + 4*a00*a10*a20*a31 + 12*a00*a10*a21*a31 + 4*a00*a10*a20*a32 + 4*a00*a10*a21*a32

 

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
    f.write('# Initialization of the symbolic variables\n')
    f.write("A=HM(sz,sz,'a'); Y=var_list('y',sz); K=var_list('k',sz)\n\n")
    f.write('# Initialization of the Permutations\n')
    f.write('P=Permutations(sz)\n\n')
    f.write('# Computing the sum over gracefully labeled functional directed graphs\n')
    f.write('f=0 \n\n')
    f.write('# Main loop\n')
    f.write('for t in DecreasingFunctionList(sz):\n')
    # variable storing the indentation spaces
    sp = '    '
    for i in range(sz):
        tmpString='for K['+str(i)+'] in rg(2):\n'
        f.write(sp+tmpString)
        sp=sp+'    '
    f.write(sp+'if prod(Y[j-t[j][1]]^K[j]*Y[sz-1-t[j][1]]^(1-K[j]) for j in rg(sz))==prod(Y):\n')
    sp=sp+'    '
    f.write(sp+'for p in P:\n')
    sp=sp+'    '
    f.write(sp+'if prod(p[k-t[k][1]]*K[k]+p[sz-1-t[k][1]]*(1-K[k]) >= p[k-t[k][1]]*(1-K[k])+p[sz-1-t[k][1]]*K[k] for k in rg(sz)):\n')
    sp=sp+'    '
    f.write(sp+'f=f+prod(A[p[sz-1-t[i][1]]-1,p[i-t[i][1]]-1]^(1-K[i])*A[p[i-t[i][1]]-1,p[sz-1-t[i][1]]-1]^K[i] for i in rg(sz))/2\n\n')
    f.write('# Initialization of the list of graphs\n')
    f.write('LtG = [Monomial2TII(mnm, A.list(), sz) for mnm in f.operands()]\n\n')
    f.write('# Plotting the trees \n')
    f.write('for g in LtG:\n')
    f.write('    T2DiGraphII(g,sz).plot().show()\n')
    # Closing the file
    f.close()

def generate_script_gracefully_labeled_functional_digraphs(sz):
    """
    The produces an implementation an optimal sage script for
    displaying and listing gracefully labeled functional directed graphs
    on sz vertices.


    EXAMPLES:

    ::

        sage: from subprocess import call
        sage: call("touch fast_listing_of_graceful_graph_on_10_vertices.sage", shell=True)
        0
        sage: call("rm fast_listing_of_graceful_graph_on_10_vertices.sage", shell=True)
        0
        sage: generate_script_gracefully_labeled_functional_digraphs(4)
        sage: call("ls fast_listing_of_graceful_graphs_on_4_vertices.sage", shell=True)
        fast_listing_of_graceful_trees_on_4_vertices.sage
        0
        sage: load('fast_listing_of_graceful_graphs_on_4_vertices.sage')
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives 
        sage: f
        a00*a10*a20*a30 + a01*a11*a20*a30 + a00*a12*a20*a30 + a02*a11*a21*a30 + a00*a13*a21*a30 + a02*a10*a22*a30 + a02*a12*a22*a30 + a01*a13*a22*a30 + a02*a11*a23*a30 + a00*a13*a23*a30 + a03*a11*a21*a31 + a03*a10*a22*a31 + a03*a12*a22*a31 + a03*a11*a23*a31 + a03*a11*a20*a32 + a03*a13*a22*a32 + a03*a10*a20*a33 + a03*a12*a20*a33 + a03*a13*a21*a33 + a03*a13*a23*a33

 

    AUTHORS:
    - Edinah K. Gnang
    """
    # Creating the string corresponding to the file name
    filename = 'fast_listing_of_graceful_graphs_on_'+str(sz)+'_vertices.sage'
    # Opening the file
    f = open(filename,'w')
    f.write('# Loading the Hypermatrix Package\n')
    f.write("load('./Hypermatrix_Algebra_tst.sage')\n\n")
    f.write('# Loading the Graceful Graph Package\n')
    f.write("load('./graceful_graph_package.sage')\n\n")
    f.write('# Initialization of the size parameter\n')
    f.write('sz='+str(sz)+'\n\n')
    f.write('# Initialization of the symbolic variables\n')
    f.write("A=HM(sz,sz,'a'); Y=var_list('y',sz); K=var_list('k',sz)\n\n")
    f.write('# Initialization of the Permutations\n')
    f.write('P=Permutations(sz)\n\n')
    f.write('# Computing the sum over gracefully labeled functional directed graphs\n')
    f.write('f=0 \n\n')
    f.write('# Main loop\n')
    f.write('for t in NonIncreasingFunctionList(sz)[0]:\n')
    # variable storing the indentation spaces
    sp = '    '
    for i in range(sz):
        tmpString='for K['+str(i)+'] in rg(2):\n'
        f.write(sp+tmpString)
        sp=sp+'    '
    f.write(sp+'if prod(Y[j-t[j][1]]^K[j]*Y[sz-1-t[j][1]]^(1-K[j]) for j in rg(sz))==prod(Y):\n')
    sp=sp+'    '
    f.write(sp+'f=f+prod(A[sz-1-t[i][1], i-t[i][1]]^(1-K[i])*A[i-t[i][1], sz-1-t[i][1]]^K[i] for i in rg(sz))/2\n\n')
    f.write('# Initialization of the list of graphs\n')
    f.write('LtG = [Monomial2TII(mnm, A.list(), sz) for mnm in f.operands()]\n\n')
    f.write('# Plotting the graphs \n')
    f.write('for g in LtG:\n')
    f.write('    T2DiGraphII(g,sz).plot().show()\n')
    # Closing the file
    f.close()

def generate_script_gracefully_labeled_functional_digraphsII(sz):
    """
    The produces an implementation an optimal sage script for
    displaying and listing gracefully labeled functional directed graphs
    on sz vertices. The difference with the implementation above is that
    we start from trees.


    EXAMPLES:

    ::

        sage: from subprocess import call
        sage: call("touch fast_listing_of_graceful_graph_on_10_vertices.sage", shell=True)
        0
        sage: call("rm fast_listing_of_graceful_graph_on_10_vertices.sage", shell=True)
        0
        sage: generate_script_gracefully_labeled_functional_digraphsII(4)
        sage: call("ls fast_listing_of_graceful_graphs_on_4_vertices.sage", shell=True)
        fast_listing_of_graceful_trees_on_4_vertices.sage
        0
        sage: load('fast_listing_of_graceful_graphs_on_4_vertices.sage')
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        Launched png viewer for Graphics object consisting of 12 graphics primitives
        sage: f
        a00*a10*a20*a30 + a00*a12*a20*a30 + a00*a13*a21*a30 + a00*a13*a23*a30
 

    AUTHORS:
    - Edinah K. Gnang
    """
    # Creating the string corresponding to the file name
    filename = 'fast_listing_of_graceful_graphs_on_'+str(sz)+'_vertices.sage'
    # Opening the file
    f = open(filename,'w')
    f.write('# Loading the Hypermatrix Package\n')
    f.write("load('./Hypermatrix_Algebra_tst.sage')\n\n")
    f.write('# Loading the Graceful Graph Package\n')
    f.write("load('./graceful_graph_package.sage')\n\n")
    f.write('# Initialization of the size parameter\n')
    f.write('sz='+str(sz)+'\n\n')
    f.write('# Initialization of the symbolic variables\n')
    f.write("A=HM(sz,sz,'a'); Y=var_list('y',sz); K=var_list('k',sz)\n\n")
    f.write('# Initialization of the Permutations\n')
    f.write('P=Permutations(sz)\n\n')
    f.write('# Computing the sum over gracefully labeled functional directed graphs\n')
    f.write('f=0 \n\n')
    f.write('# Main loop\n')
    f.write('for t in IncreasingFunctionList(sz):\n')
    # variable storing the indentation spaces
    sp = '    '
    for i in range(sz-1):
        tmpString='for K['+str(i)+'] in rg(2):\n'
        f.write(sp+tmpString)
        sp=sp+'    '
    f.write(sp+'K[sz-1]=0\n')
    f.write(sp+'if prod(Y[j-t[j][1]]^K[j]*Y[sz-1-t[j][1]]^(1-K[j]) for j in rg(sz))==prod(Y):\n')
    sp=sp+'    '
    f.write(sp+'f=f+prod(A[sz-1-t[i][1], i-t[i][1]]^(1-K[i])*A[i-t[i][1], sz-1-t[i][1]]^K[i] for i in rg(sz))\n\n')
    f.write('# Initialization of the list of graphs\n')
    f.write('LtG = [Monomial2TII(mnm, A.list(), sz) for mnm in f.operands()]\n\n')
    f.write('# Plotting the graphs \n')
    f.write('for g in LtG:\n')
    f.write('    T2DiGraphII(g,sz).plot().show()\n')
    # Closing the file
    f.close()

def RandomGraphSignedPermutationClasses(n):
    """
    Returns a random graph selected uniformly among the
    unlabeled graphs which admits a graceful labelin.


    EXAMPLES:

    ::


        sage: RandomGraphSignedPermutationClasses(2)
        [0, -1]


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


        sage: RandomTreeSignedPermutationClasses(2)
        [0, -1]


    AUTHORS:
    - Edinah K. Gnang
    """
    L=TreeSignedPermutationClassesII(n)
    return L[randint(0,len(L)-1)]

def RandomTreeTupleClasses(n):
    """
    Returns a random tree selected uniformly among the
    unlabeled trees which admits a graceful labelin.


    EXAMPLES:

    ::


        sage: RandomTreeTupleClasses(2)
        [(0, 0), (1, 0)]


    AUTHORS:
    - Edinah K. Gnang
    """
    return SP2T(RandomTreeSignedPermutationClasses(n))

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
    Implements the counting formula derived via Cramer's rule
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

def TupleLabelInvolution(T):
    """
    Returns the tuple encoding of the involuted labeling


    EXAMPLES:

    ::

        sage: T=[(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)]
        sage: TupleLabelInvolution(T)
        [(0, 0), (1, 4), (2, 0), (3, 2), (4, 0)]


    AUTHORS:
    - Edinah K. Gnang
    """
    return [(0, 0)]+[(len(T)-1-T[i][0], len(T)-1-T[i][1]) for i in range(len(T)-2,0,-1)]+[(len(T)-1,0)]

def TupleLabelInvolutionII(T):
    """
    Returns the tuple encoding of the involuted labeling
    This implementation does not assume that the tuple is rooted at 0.

    EXAMPLES:

    ::

        sage: T=[(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)]
        sage: TupleLabelInvolutionII(T)
        [(0, 4), (1, 4), (2, 0), (3, 2), (4, 4)]


    AUTHORS:
    - Edinah K. Gnang
    """
    sz=len(T)
    tp=[(sz-1-t[0], sz-1-t[1]) for t in T]
    tp.sort()
    return tp

def SignedPermutationLabelInvolution(s):
    """
    Returns the signed permutations encoding of the involuted labeling


    EXAMPLES:

    ::

        sage: T=[(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)]
        sage: SignedPermutationLabelInvolution(T2SP(T))
        [0, 3, -2, -1, -4]


    AUTHORS:
    - Edinah K. Gnang
    """
    T=SP2T(s)
    return T2SP([(0, 0)]+[(len(T)-1-T[i][0], len(T)-1-T[i][1]) for i in range(len(T)-2,0,-1)]+[(len(T)-1,0)])

def incidenceHM(sigma):
    """
    Returns the transpose of the conventional incidence matrix associated with the input
    tree specified by a signed permutation. 


    EXAMPLES:

    ::

        sage: incidenceHM([0, 1, -2])
        [[0, 1, -1], [-1, 0, 1]]     


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
    Eq = [Le[i] == Y[i] for i in range(1,sz)]
    # Formating the constraints to extract the incidence
    # matrix and the edge vector on the right hand side.
    return ConstraintFormatorHM(Eq,X)[0]

def incidenceHMII(T):
    """
    Returns the transpose of the conventional incidence matrix associated with the input
    tree specified by a list of tuples associated with the function. 


    EXAMPLES:

    ::

        sage: incidenceHMII([(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)])
        [[1, 0, 0, 0, -1], [0, 1, 0, 0, -1], [0, 0, 1, 0, -1], [0, 0, 0, 1, -1], [0, 0, 0, 0, 0]]


    AUTHORS:
    - Edinah K. Gnang
    """
    sz = len(T)
    # Initialization of the list of variables
    X = var_list('x',sz); Y = var_list('y',sz)
    # Initializing the permutation performing the involution
    q=[sz-1-i for i in range(sz)]
    # Initialization of the list of tuples
    Tp=[(q[t[0]], q[t[1]]) for t in T[::-1]]
    # Initialization of the list of constraints
    Eq=[(Tuple2EdgeList(Tp,'x'))[i] == Y[i] for i in range(sz-1)]+[0*X[sz-1] == Y[sz-1]]
    return ConstraintFormatorHM(Eq,X)[0]

def incidenceHMIII(T):
    """
    Returns the transpose of the conventional incidence matrix associated with the input
    tree specified by a list of tuples associated with the function. 


    EXAMPLES:

    ::

        sage: incidenceHMIII([(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)])
        [[0, 0, 0, 0, 0], [1, -1, 0, 0, 0], [1, 0, -1, 0, 0], [1, 0, 0, -1, 0], [1, 0, 0, 0, -1]]


    AUTHORS:
    - Edinah K. Gnang
    """
    sz=len(T)
    # Initialization of the list of variables
    X=var_list('x',sz); Y=var_list('y',sz)
    # Initialization of the list of tuples
    Tp=[(t[0],t[1]) for t in T]
    # Initialization of the list of constraints
    Eq=[SR(Tuple2EdgeList(Tp,'x')[i]) == Y[i] for i in rg(sz)]
    return ConstraintFormatorHM(Eq,X)[0]

def adjacency_polynomial_constructionIa(X, tp):
    """
    Returns the polynomial adjacency polynomial construction over roots of
    unity which certifies the existence of a graceful labeling. The inputs
    to the functions are the list of variables associaated with vertices
    and a tuple edge list decription of the input graph.


    EXAMPLES:

    ::

        sage: sz=4; FnLst=[SP2T(l) for l in GraphSignedPermutationClassesII(sz)]
        sage: tp=FnLst[0]; X=var_list('x',sz); adjacency_polynomial_constructionIa(X, tp)
        [-x2/x1 + x0/x3,
         -x1/x2 + x0/x3,
         x0/x2 - x2/x1,
         1,
         -x2/x0 + x0/x3,
         (x0 - x1)*(x0 - x2)*(x0 - x3)*(x1 - x2)*(x1 - x3)*(x2 - x3),
         (w^3 - x0)*(w^3 - x1)*(w^3 - x2)*(w^3 - x3)*(w^2 - x0)*(w^2 - x1)*(w^2 - x2)*(w^2 - x3)*(w - x0)*(w - x1)*(w - x2)*(w - x3)]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter associated with the number of vertices
    sz=len(X)
    # Initialization of the list of roots of unity
    w=var('w')
    #W=[exp(I*2*pi*k/sz) for k in rg(sz)]
    # Initialization of the polynomial
    f0=prod( (X[tp[j][1]]/X[tp[j][0]] - X[tp[i][1]]/X[tp[i][0]]) for i in rg(1,sz) for j in rg(1,sz) if j>i and tp[j][1] != tp[i][1] and tp[i][1] != tp[j][0] and tp[j][1] != tp[i][0] )
    f1=prod( (X[tp[j][1]]/X[tp[j][0]] - X[tp[i][0]]/X[tp[i][1]]) for i in rg(1,sz) for j in rg(1,sz) if j>i and tp[j][1] != tp[i][1] and tp[i][1] != tp[j][0] and tp[j][1] != tp[i][0] )
    f2=prod( (X[tp[j][1]]/X[tp[j][0]] - X[tp[i][1]]/X[tp[i][0]]) for i in rg(1,sz) for j in rg(1,sz) if j>i and (tp[i][1] == tp[j][0] or tp[j][1] == tp[i][0]))
    f3=prod( (X[tp[j][1]]/X[tp[j][0]] - X[tp[i][0]]/X[tp[i][1]]) for i in rg(1,sz) for j in rg(1,sz) if j>i and tp[j][1] == tp[i][1] )
    f4=prod( (X[j] - X[i]) for i in rg(sz) for j in rg(sz) if j>i )
    f5=prod( (X[s] - w^t) for s in rg(sz) for t in rg(1,sz) )
    #f6=prod( (X[s] - 1/W[t]) for s in rg(sz) for t in rg(1,sz) )
    return [f0, f1, f2, f3, f4, f5]

def adjacency_polynomial_constructionIb(X, tp):
    """
    Returns the adjacency polynomial construction over roots of unity 
    which certifies the existence of a graceful labeling. The inputs
    to the functions are the list of variables associated with 
    vertices and a tuple edge list decription of the input graph.


    EXAMPLES:

    ::

        sage: sz=4; FnLst=[SP2T(l) for l in GraphSignedPermutationClassesII(sz)]
        sage: tp=FnLst[0]; X=var_list('x',sz); L=adjacency_polynomial_constructionIb(X, tp); L
        [-x1^6*x2 + x0*x3^6,
         -x1*x2^6 + x0*x3^6,
         -x1^6*x2 + x0*x2^6,
         1,
         -x0^6*x2 + x0*x3^6,
         (x0 - x1)*(x0 - x2)*(x0 - x3)*(x1 - x2)*(x1 - x3)*(x2 - x3),
         (w^3 - x0)*(w^3 - x1)*(w^3 - x2)*(w^3 - x3)*(w^2 - x0)*(w^2 - x1)*(w^2 - x2)*(w^2 - x3)*(w - x0)*(w - x1)*(w - x2)*(w - x3)]
        sage: f=prod(L[:-2]); f
        (x0^6*x2 - x0*x3^6)*(x1^6*x2 - x0*x2^6)*(x1^6*x2 - x0*x3^6)*(x1*x2^6 - x0*x3^6)
        sage: reduce_polynomial(f.expand(), X)
        x0^6*x1^6*x2^2 - x0^3*x1^6*x2*x3^4 + x0^4*x2^6*x3^4 - x0^3*x1^6*x3^5 + x0^2*x1^5*x2^2*x3^5 + x0*x1^6*x2^2*x3^5 - x0^3*x1*x2^5*x3^5 - x1^5*x2^3*x3^6 + x0^2*x2^6*x3^6 + x0*x1*x2^6*x3^6 - x2*x3^6 - 1


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter associated with the number of vertices
    sz=len(X)
    # Initialization of the variable associated with the primitive root of unity
    w=var('w')
    #W=[exp(I*2*pi*k/sz) for k in rg(sz)]
    d=2*(sz-1)
    # Initialization of the polynomial
    f0=prod( (X[tp[j][1]]*X[tp[j][0]]^d - X[tp[i][1]]*X[tp[i][0]]^d) for i in rg(1,sz) for j in rg(1,sz) if j>i and tp[j][1] != tp[i][1] and tp[i][1] != tp[j][0] and tp[j][1] != tp[i][0] )
    f1=prod( (X[tp[j][1]]*X[tp[j][0]]^d - X[tp[i][0]]*X[tp[i][1]]^d) for i in rg(1,sz) for j in rg(1,sz) if j>i and tp[j][1] != tp[i][1] and tp[i][1] != tp[j][0] and tp[j][1] != tp[i][0] )
    f2=prod( (X[tp[j][1]]*X[tp[j][0]]^d - X[tp[i][1]]*X[tp[i][0]]^d) for i in rg(1,sz) for j in rg(1,sz) if j>i and (tp[i][1] == tp[j][0] or tp[j][1] == tp[i][0]) )
    f3=prod( (X[tp[j][1]]*X[tp[j][0]]^d - X[tp[i][0]]*X[tp[i][1]]^d) for i in rg(1,sz) for j in rg(1,sz) if j>i and tp[j][1] == tp[i][1] )
    f4=prod( (X[j] - X[i]) for i in rg(sz) for j in rg(sz) if j>i )
    f5=prod( (X[s] - w^t) for s in rg(sz) for t in rg(1,sz) )
    #f6=prod( (X[s] - 1/W[t]) for s in rg(sz) for t in rg(1,sz) )
    return [f0, f1, f2, f3, f4, f5]

def adjacency_polynomial_constructionIc(X, tp):
    """
    Returns the adjacency polynomial construction over roots of unity 
    which certifies the existence of a graceful labeling. The inputs
    to the functions are the list of variables associaated with 
    vertices and a tuple edge list decription of the input graph.


    EXAMPLES:

    ::

        sage: sz=4; FnLst=[SP2T(l) for l in GraphSignedPermutationClassesII(sz)]
        sage: tp=FnLst[0]; X=var_list('x',sz); L=adjacency_polynomial_constructionIc(X, tp); L
        [x0*x1 - x2*x3,
         x0*x2 - x1*x3,
         x0*x1 - x2^2,
         1,
         x0^2 - x2*x3,
         (x0 - x1)*(x0 - x2)*(x0 - x3)*(x1 - x2)*(x1 - x3)*(x2 - x3),
         (w^3 - x0)*(w^3 - x1)*(w^3 - x2)*(w^3 - x3)*(w^2 - x0)*(w^2 - x1)*(w^2 - x2)*(w^2 - x3)*(w - x0)*(w - x1)*(w - x2)*(w - x3)]
        sage: f=prod(L[:-2]); f
        (x0^2 - x2*x3)*(x0*x1 - x2^2)*(x0*x1 - x2*x3)*(x0*x2 - x1*x3)
        sage: reduce_polynomial(f.expand(), X)
        x0^5*x1^2*x2 - x0^4*x1*x2^3 - x0^4*x1^3*x3 - x0^4*x1*x2^2*x3 + x0^3*x2^4*x3 + x0^2*x1*x2^4*x3 + x0^3*x1^2*x2*x3^2 + x0^2*x1^3*x2*x3^2 - x0*x1^2*x2^3*x3^2 - x0*x2^5*x3^2 - x0*x1^2*x2^2*x3^3 + x1*x2^4*x3^3


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter associated with the number of vertices
    sz=len(X)
    # Initialization of the variable associated with the primitive root of unity
    w=var('w')
    #W=[exp(I*2*pi*k/sz) for k in rg(sz)]
    d=2*(sz-1)
    # Initialization of the polynomial
    f0=prod( (X[tp[j][1]]*X[tp[i][0]] - X[tp[i][1]]*X[tp[j][0]]) for i in rg(1,sz) for j in rg(1,sz) if j>i and tp[j][1] != tp[i][1] and tp[i][1] != tp[j][0] and tp[j][1] != tp[i][0] )
    f1=prod( (X[tp[j][1]]*X[tp[i][1]] - X[tp[i][0]]*X[tp[j][0]]) for i in rg(1,sz) for j in rg(1,sz) if j>i and tp[j][1] != tp[i][1] and tp[i][1] != tp[j][0] and tp[j][1] != tp[i][0] )
    f2=prod( (X[tp[j][1]]*X[tp[i][0]] - X[tp[i][1]]*X[tp[j][0]]) for i in rg(1,sz) for j in rg(1,sz) if j>i and tp[i][1] == tp[j][0] )
    f3=prod( (X[tp[j][1]]*X[tp[i][0]] - X[tp[i][1]]*X[tp[j][0]]) for i in rg(1,sz) for j in rg(1,sz) if j>i and tp[j][1] == tp[i][0] )
    f4=prod( (X[tp[j][1]]*X[tp[i][1]] - X[tp[i][0]]*X[tp[j][0]]) for i in rg(1,sz) for j in rg(1,sz) if j>i and tp[j][1] == tp[i][1] )
    f5=prod( (X[j] - X[i]) for i in rg(sz) for j in rg(sz) if j>i )
    f6=prod( (X[s] - w^t) for s in rg(sz) for t in rg(1,sz) )
    #f6=prod( (X[s] - 1/W[t]) for s in rg(sz) for t in rg(1,sz) )
    return [f0, f1, f2, f3, f4, f5, f6]

def adjacency_polynomial_constructionII(X, tp):
    """
    Returns the polynomial adjacency polynomial construction over integers
    which certifies the existence of a graceful labeling.  The inputs to 
    the functions are the list of variables associaated with vertices and
    a tuple edge list decription of the input graph.


    EXAMPLES:

    ::

        sage: sz=4; FnLst=[SP2T(l) for l in GraphSignedPermutationClassesII(sz)]
        sage: tp=FnLst[0]; adjacency_polynomial_constructionII(var_list('x',sz), tp)
        [x0 + x1 - x2 - x3,
         x0 - x1 + x2 - x3,
         x0 + x1 - 2*x2,
         1,
         2*x0 - x2 - x3,
         (x0 - x1)*(x0 - x2)*(x0 - x3)*(x1 - x2)*(x1 - x3)*(x2 - x3),
         (x0 + 3)*(x0 + 2)*(x0 + 1)*(x1 + 3)*(x1 + 2)*(x1 + 1)*(x2 + 3)*(x2 + 2)*(x2 + 1)*(x3 + 3)*(x3 + 2)*(x3 + 1)]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter associated with the number of vertices
    sz=len(X)
    # Initializing the degree 
    d=2*(sz-1)
    # Initialization of the polynomial
    f0=prod( (X[tp[j][1]]-X[tp[j][0]] - (X[tp[i][1]]-X[tp[i][0]])) for i in rg(1,sz) for j in rg(1,sz) if j>i and tp[j][1] != tp[i][1] and tp[i][1] != tp[j][0] and tp[j][1] != tp[i][0] )
    f1=prod( (X[tp[j][1]]-X[tp[j][0]] - (X[tp[i][0]]-X[tp[i][1]])) for i in rg(1,sz) for j in rg(1,sz) if j>i and tp[j][1] != tp[i][1] and tp[i][1] != tp[j][0] and tp[j][1] != tp[i][0] )
    f2=prod( (X[tp[j][1]]-X[tp[j][0]] - (X[tp[i][1]]-X[tp[i][0]])) for i in rg(1,sz) for j in rg(1,sz) if j>i and (tp[i][1] == tp[j][0] or tp[j][1] == tp[i][0]))
    f3=prod( (X[tp[j][1]]-X[tp[j][0]] - (X[tp[i][0]]-X[tp[i][1]])) for i in rg(1,sz) for j in rg(1,sz) if j>i and tp[j][1] == tp[i][1] )
    f4=prod( (X[j] - X[i]) for i in rg(sz) for j in rg(sz) if j>i )
    f5=prod( (X[s] + t) for s in rg(sz) for t in rg(1,sz) )
    return [f0, f1, f2, f3, f4, f5]

def adjacency_polynomial_constructionIII(X, tp):
    """
    Returns the polynomial adjacency polynomial construction over integers
    which certifies the existence of a graceful labeling.  The inputs to 
    the functions are the list of variables associated with vertices and
    a tuple edge list decription of the input graph. This is the redundant
    version because it does not enforce i<j.


    EXAMPLES:

    ::

        sage: sz=4; FnLst=[SP2T(l) for l in GraphSignedPermutationClassesII(sz)]
        sage: tp=FnLst[0]; X=var_list('x',sz); adjacency_polynomial_constructionIII(X, tp)
        [x0 + x1 - x2 - x3,
         x0 - x1 + x2 - x3,
         x0 + x1 - 2*x2,
         1,
         2*x0 - x2 - x3,
         (x0 - x1)*(x0 - x2)*(x0 - x3)*(x1 - x2)*(x1 - x3)*(x2 - x3),
         (x0 + 3)*(x0 + 2)*(x0 + 1)*(x1 + 3)*(x1 + 2)*(x1 + 1)*(x2 + 3)*(x2 + 2)*(x2 + 1)*(x3 + 3)*(x3 + 2)*(x3 + 1)]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter associated with the number of vertices
    sz=len(X)
    # Initializing the degree 
    d=2*(sz-1)
    # Initialization of the polynomial
    f0=prod( (X[tp[j][1]]-X[tp[j][0]] - (X[tp[i][1]]-X[tp[i][0]])) for i in rg(1,sz) for j in rg(1,sz) if j != i and tp[j][1] != tp[i][1] and tp[i][1] != tp[j][0] and tp[j][1] != tp[i][0] )
    f1=prod( (X[tp[j][1]]-X[tp[j][0]] - (X[tp[i][0]]-X[tp[i][1]])) for i in rg(1,sz) for j in rg(1,sz) if j != i and tp[j][1] != tp[i][1] and tp[i][1] != tp[j][0] and tp[j][1] != tp[i][0] )
    f2=prod( (X[tp[j][1]]-X[tp[j][0]] - (X[tp[i][1]]-X[tp[i][0]])) for i in rg(1,sz) for j in rg(1,sz) if j != i and (tp[i][1] == tp[j][0] or tp[j][1] == tp[i][0]) )
    f3=prod( (X[tp[j][1]]-X[tp[j][0]] - (X[tp[i][0]]-X[tp[i][1]])) for i in rg(1,sz) for j in rg(1,sz) if j != i and tp[j][1] == tp[i][1] )
    return [f0, f1, f2, f3]

def reduce_polynomial(f, X):
    """
    Returns the remainder of the input polynomial modulo the algebraic
    relation in the variable entries of the input list X. The algebraic
    relation are associated with the roots of unity.


    EXAMPLES:

    ::

        sage: sz=4; FnLst=[SP2T(l) for l in GraphSignedPermutationClassesII(sz)]
        sage: tp=FnLst[0]; X=var_list('x',sz); L=adjacency_polynomial_constructionIb(X, tp)
        sage: f=expand(prod(L[:-2]))
        sage: reduce_polynomial(f, X)
        x0^6*x1^6*x2^2 - x0^3*x1^6*x2*x3^4 + x0^4*x2^6*x3^4 - x0^3*x1^6*x3^5 + x0^2*x1^5*x2^2*x3^5 + x0*x1^6*x2^2*x3^5 - x0^3*x1*x2^5*x3^5 - x1^5*x2^3*x3^6 + x0^2*x2^6*x3^6 + x0*x1*x2^6*x3^6 - x2*x3^6 - 1


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter
    sz=len(X)
    # Initialization of the temparary function
    ftmp = f
    # Obtaining the degree of each variables in the canonical monomial
    DgL=[ftmp.degree(X[i]) for i in rg(len(X))]
    # While loop performing the desired reduction
    while max(DgL) >= 2*sz-1:
        ftmp=fast_reduce(ftmp, [X[i]^DgL[i] for i in rg(sz) if DgL[i]>=2*sz-1],[X[i]^Integer(mod(DgL[i],2*sz-1)) for i in rg(sz) if DgL[i]>=2*sz-1])
        # Updating the list of degree list
        DgL=[ftmp.degree(X[i]) for i in rg(len(X))]
        #print DgL
    return ftmp

def matrix_tree_generating_polynomial(sz):
    """
    Uses the matrix tree theorem to derive the generating polynomial 
    which seprates the graphs edge weight patterns. The number of 
    graceful labelings is given the coefficient of 
    x^((sz^(sz-1)-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=4; matrix_tree_generating_polynomial(sz)
        sage: f.coefficient(x^((sz^(sz-1)-1)/(sz-1)))
        4


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the variables.
    x=var('x'); X=HM(sz,sz,'a')
    # Initialization of the hypermatrix.
    A=HM(sz,sz,[sqrt(x^(sz^(abs(i-j)-1))) for j in rg(sz) for i in rg(sz)])
    # Initialization of the symbolic incidence matrix.
    Hb=HM(sz,binomial(sz,2),'zero')
    # Initializing the column index
    clidx=0
    # Loop performing filling the incidence matrix.
    # Note that in an edge (i,j) we have i < j.
    for j in rg(1,sz):
        for i in rg(j):
            Hb[i,clidx]=A[i,j]; Hb[j,clidx]=-A[i,j]
            clidx=clidx+1
    # We ground the graph at the vertex labeled 0.
    B=HM(sz-1,binomial(sz,2),[Hb[i,j] for j in rg(binomial(sz,2)) for i in rg(1,sz)])
    # Computing the generating function.
    f=expand((B*B.transpose()).det())
    return f 

def matrix_tree_list_generating_polynomial(sz):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the generating polynomial which seprates the graphs edge weight
    patterns. The edge list of graceful labelings is given the 
    coefficient of x^((sz^(sz-1)-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=4; f=matrix_tree_list_generating_polynomial(sz)
        sage: f.coefficient(x^((sz^(sz-1)-1)/(sz-1)))
        a01*a02*a03 + a02*a03*a12 + a03*a12*a13 + a03*a13*a23


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the variables.
    x=var('x'); X=HM(sz,sz,'a')
    # Initialization of the hypermatrix.
    A=HM(sz,sz,[sqrt(X[i,j]*x^(sz^(abs(i-j)-1))) for j in rg(sz) for i in rg(sz)])
    # Initialization of the symbolic incidence matrix.
    Hb=HM(sz,binomial(sz,2),'zero')
    # Initializing the column index
    clidx=0
    # Loop performing filling the incidence matrix.
    # Note that in an edge (i,j) we have i < j.
    for j in rg(1,sz):
        for i in rg(j):
            Hb[i,clidx]=A[i,j]; Hb[j,clidx]=-A[i,j]
            clidx=clidx+1
    # We ground the graph at the vertex labeled 0.
    B=HM(sz-1,binomial(sz,2),[Hb[i,j] for j in rg(binomial(sz,2)) for i in rg(1,sz)])
    # Computing the generating function.
    f=expand((B*B.transpose()).det())
    return f

def matrix_tree_list_generating_polynomialII(sz):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the generating polynomial which seprates the graphs edge weight
    patterns. The edge list of graceful labelings is given the 
    coefficient of x^((sz^sz-1)/(sz-1)). The difference with the
    implementation above is that the degrees are higher.


    EXAMPLES:

    ::


        sage: sz=4; f=matrix_tree_list_generating_polynomialII(sz)
        sage: f.coefficient(x^((sz^sz-1)/(sz-1)))
        a01*a02*a03 + a02*a03*a12 + a03*a12*a13 + a03*a13*a23


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the variables.
    x=var('x'); X=HM(sz,sz,'a')
    # Initialization of the hypermatrix.
    A=HM(sz,sz,[sqrt(X[i,j]*x^(sz^abs(i-j))) for j in rg(sz) for i in rg(sz)])
    # Initialization of the symbolic incidence matrix.
    Hb=HM(sz,binomial(sz,2),'zero')
    # Initializing the column index
    clidx=0
    # Loop performing filling the incidence matrix.
    # Note that in an edge (i,j) we have i < j.
    for j in rg(1,sz):
        for i in rg(j):
            Hb[i,clidx]=A[i,j]; Hb[j,clidx]=-A[i,j]
            clidx=clidx+1
    # We ground the graph at the vertex labeled 0.
    B=HM(sz-1,binomial(sz,2),[Hb[i,j] for j in rg(binomial(sz,2)) for i in rg(1,sz)])
    # Computing the generating function.
    f=expand((B*B.transpose()).det())
    return f

def list_generating_polynomial(sz):
    """
    Uses a variant of the matrix tree theorem to construct the
    generating polynomial which seprates the graphs edge weight
    patterns. The edge list of graceful labelings is given the 
    coefficient of x^((sz^sz-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=3; f=list_generating_polynomial(sz)
        sage: factor(f.coefficient(x^((sz^sz-1)/(sz-1))))
        (a00 + a11 + a22)*(a01 + a12)*a02


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the variables.
    x=var('x'); X=HM(sz,sz,'a')
    # Initialization of the hypermatrix.
    A=HM(sz,sz,[(X[i,j]*x^(sz^(abs(i-j)))) for j in rg(sz) for i in rg(sz)])
    # Initialization of the symbolic incidence matrix.
    Hb=HM(sz,sz+binomial(sz,2),'one')
    # Initializing the column index
    clidx=0
    # Loop performing filling the incidence matrix.
    # Note that in an edge (i,j) we have i < j.
    for j in rg(sz):
        for i in rg(j+1):
            Hb[i,clidx]=A[i,j]
            clidx=clidx+1
    # Initializing the list of submatrices
    Lb = Hb.side_length_subhypermatrices(sz)
    # Computing and returning the generating function.
    return sum(prod(mtr.list()) for mtr in Lb)

def generating_polynomial(sz):
    """
    Uses a variant of the matrix tree theorem to construct the
    generating polynomial which seprates the graphs edge weight
    patterns. The number of graceful labelings is given the 
    coefficient of x^((sz^sz-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=3; f=generating_polynomial(sz)
        sage: f.coefficient(x^((sz^sz-1)/(sz-1)))
        3


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the variables.
    x=var('x'); X=HM(sz,sz,'a')
    # Initialization of the hypermatrix.
    A=HM(sz,sz,[x^(sz^(abs(i-j))) for j in rg(sz) for i in rg(sz)])
    # Initialization of the symbolic incidence matrix.
    Hb=HM(sz,sz+binomial(sz,2),'one')
    # Initializing the column index
    clidx=0
    # Loop performing filling the incidence matrix.
    # Note that in an edge (i,j) we have i < j.
    for j in rg(sz):
        for i in rg(j+1):
            Hb[i,clidx]=A[i,j]
            clidx=clidx+1
    # Initializing the list of submatrices
    Lb = Hb.side_length_subhypermatrices(sz)
    # Computing and returning the generating function.
    return sum(prod(mtr.list()) for mtr in Lb)

def list_generating_polynomialII(sz):
    """
    Uses a variant of the matrix tree theorem to construct the
    generating polynomial which seprates the graphs edge weight
    patterns. The edge list of graceful labelings is given the 
    coefficient of x^((sz^(sz-1)-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=4; f=list_generating_polynomialII(sz)
        sage: f.coefficient(x^((sz^(sz-1)-1)/(sz-1)))
        4*a01*a02*a03 + 4*a02*a03*a12 + 4*a01*a03*a13 + 4*a03*a12*a13 + 4*a02*a03*a23 + 4*a03*a13*a23


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the variables.
    x=var('x'); X=HM(sz,sz,'a')
    # Initialization of the hypermatrix.
    A=HM(sz,sz,[(X[i,j]*x^(sz^(abs(i-j)-1))) for j in rg(sz) for i in rg(sz)])
    # Initialization of the symbolic incidence matrix.
    Hb=HM(sz,binomial(sz,2),'one'); Hc=HM(sz,binomial(sz,2),'one')
    # Initializing the column index
    clidx=0
    # Loop performing filling the incidence matrix.
    # Note that in an edge (i,j) we have i < j.
    for j in rg(1,sz):
        for i in rg(j):
            Hb[i,clidx]=A[i,j]
            Hc[j,clidx]=A[i,j]
            clidx=clidx+1
    # Initializing the list of submatrices
    Lb = Hb.side_length_subhypermatrices(sz-1)
    Lc = Hc.side_length_subhypermatrices(sz-1)
    # Computing the generating function.
    f=sum(prod(mtr.list()) for mtr in Lb)
    g=sum(prod(mtr.list()) for mtr in Lc)
    return f+g

def generating_polynomialII(sz):
    """
    Uses a variant of the matrix tree theorem to construct the
    generating polynomial which seprates the graphs edge weight
    patterns. The number of graceful labelings is given the 
    coefficient of x^((sz^(sz-1)-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=4; f=generating_polynomialII(sz)
        sage: f.coefficient(x^((sz^(sz-1)-1)/(sz-1)))
        24


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the variables.
    x=var('x'); X=HM(sz,sz,'a')
    # Initialization of the hypermatrix.
    A=HM(sz,sz,[(x^(sz^(abs(i-j)-1))) for j in rg(sz) for i in rg(sz)])
    # Initialization of the symbolic incidence matrix.
    Hb=HM(sz,binomial(sz,2),'one'); Hc=HM(sz,binomial(sz,2),'one')
    # Initializing the column index
    clidx=0
    # Loop performing filling the incidence matrix.
    # Note that in an edge (i,j) we have i < j.
    for j in rg(1,sz):
        for i in rg(j):
            Hb[i,clidx]=A[i,j]
            Hc[j,clidx]=A[i,j]
            clidx=clidx+1
    # Initializing the list of submatrices
    Lb = Hb.side_length_subhypermatrices(sz-1)
    Lc = Hc.side_length_subhypermatrices(sz-1)
    # Computing the generating function.
    f=sum(prod(mtr.list()) for mtr in Lb)
    g=sum(prod(mtr.list()) for mtr in Lc)
    return f+g

def tree_sum_list_generating_polynomial(sz):
    """
    Computing the sum over non isomorphic trees and accounting
    for the action of the symetric group on the vertices to obtain
    generating polynomial which seprates the graphs edge weight
    patterns. The edge monomial list associated with graceful 
    labelings is given the coefficient of x^((sz^(sz-1)-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=4; f=tree_sum_list_generating_polynomial(sz)
        sage: f.coefficient(x^((sz^(sz-1)-1)/(sz-1)))
        a01*a02*a03 + a02*a03*a12 + a03*a12*a13 + a03*a13*a23


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the matrix
    A=HM(2,sz,'a','sym')
    # Initization of list of tree representative
    FnLst=[SP2T(l) for l in TreeSignedPermutationClassesII(sz)]
    # Initialization of the function and the Permutations
    f=0; P=Permutations(sz)
    # Looping through permutations
    for q in P:
        # fixing the permutation index
        p=[q[i]-1 for i in rg(sz)]
        # Looping through the tree classes
        for tp in FnLst:
            # Initialization of the order of the automorphism group.
            Ng = T2GraphII(tp[1:],sz).automorphism_group().order()
            f = f + prod(A[p[tp[i][0]],p[tp[i][1]]]*x^(sz^(abs(p[tp[i][1]]-p[tp[i][0]])-1)) for i in rg(1,sz))/Ng
    return f
 
def tree_sum_generating_polynomial(sz):
    """
    Computing the sum over non isomorphic trees and accounting
    for the action of the symetric group on the vertices to obtain
    generating polynomial which seprates the graphs edge weight
    patterns. The number of graceful labelings is given the 
    coefficient of x^((sz^(sz-1)-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=4; f=tree_sum_generating_polynomial(sz)
        sage: f.coefficient(x^((sz^(sz-1)-1)/(sz-1)))
        4


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initization of list of tree representative
    FnLst=[SP2T(l) for l in TreeSignedPermutationClassesII(sz)]
    # Initialization of the function and the Permutations
    f=0; P=Permutations(sz)
    # Looping through permutations
    for q in P:
        # fixing the permutation index
        p=[q[i]-1 for i in rg(sz)]
        # Looping through the tree classes
        for tp in FnLst:
            # Initialization of the order of the automorphism group.
            Ng = T2GraphII(tp[1:],sz).automorphism_group().order()
            f = f + prod(x^(sz^(abs(p[tp[i][1]]-p[tp[i][0]])-1)) for i in rg(1,sz))/Ng
    return f

def orbit_tree_sum_list_generating_polynomial(sz):
    """
    Computing the sum over non isomorphic trees and accounting
    for the action of the symetric group on the vertices to obtain
    generating polynomial which seprates the graphs edge weight
    patterns. The edge monomial list associated with graceful 
    labelings is given the coefficient of x^((sz^(sz-1)-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=4; f=orbit_tree_sum_list_generating_polynomial(sz)
        sage: f.coefficient(x^((sz^(sz-1)-1)/(sz-1)))
        24*a01*a02*a03 + 24*a02*a03*a12 + 24*a03*a12*a13 + 24*a03*a13*a23


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the matrix
    A=HM(2,sz,'a','sym')
    # Initization of list of all trees
    FnLst=[[(0,0)]+[(i+1,l[i]) for i in rg(len(l))] for l in TreeFunctionList(sz)]
    # Initialization of the function and the Permutations
    f=0; P=Permutations(sz)
    # Looping through permutations
    for q in P:
        # fixing the permutation index
        p=[q[i]-1 for i in rg(sz)]
        # Looping through the tree classes
        for tp in FnLst:
            # Updating the polynomial f
            f = f + prod(A[p[tp[i][0]],p[tp[i][1]]]*x^(sz^(abs(p[tp[i][1]]-p[tp[i][0]])-1)) for i in rg(1,sz))
    return f
 
def orbit_tree_sum_generating_polynomial(sz):
    """
    Computing the sum over non isomorphic trees and accounting
    for the action of the symetric group on the vertices to obtain
    generating polynomial which seprates the graphs edge weight
    patterns. The number of graceful labelings is given the 
    coefficient of x^((sz^(sz-1)-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=4; f=orbit_tree_sum_generating_polynomial(sz)
        sage: f.coefficient(x^((sz^(sz-1)-1)/(sz-1)))
        96


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initization of list of all trees
    FnLst=[[(0,0)]+[(i+1,l[i]) for i in rg(len(l))] for l in TreeFunctionList(sz)]
    # Initialization of the function and the Permutations
    f=0; P=Permutations(sz)
    # Looping through permutations
    for q in P:
        # fixing the permutation index
        p=[q[i]-1 for i in rg(sz)]
        # Looping through the tree classes
        for tp in FnLst:
            # Updating the polynomial f
            f = f + prod(x^(sz^(abs(p[tp[i][1]]-p[tp[i][0]])-1)) for i in rg(1,sz))
    return f

def orbit_non_tree_sum_list_generating_polynomial(sz):
    """
    Computing the sum over non isomorphic trees and accounting
    for the action of the symetric group on the vertices to obtain
    generating polynomial which seprates the graphs edge weight
    patterns. The edge monomial list associated with graceful 
    labelings is given the coefficient of x^((sz^(sz-1)-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=4; f=orbit_non_tree_sum_list_generating_polynomial(sz)
        sage: f.coefficient(x^((sz^sz-1)/(sz-1)))
        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the matrix
    A=HM(2,sz,'a','sym')
    # Initization of list of all trees
    FnLst=[[(0,0)]+[(i+1,l[i]) for i in rg(len(l))] for l in NonTreeFunctionList(sz)]
    # Initialization of the function and the Permutations
    f=0; P=Permutations(sz)
    # Looping through permutations
    for q in P:
        # fixing the permutation index
        p=[q[i]-1 for i in rg(sz)]
        # Looping through the tree classes
        for tp in FnLst:
            # Updating the polynomial f
            f = f + prod(A[p[tp[i][0]],p[tp[i][1]]]*x^(sz^abs(p[tp[i][1]]-p[tp[i][0]])) for i in rg(sz))
    return f
 
def orbit_non_tree_sum_generating_polynomial(sz):
    """
    Computing the sum over non isomorphic trees and accounting
    for the action of the symetric group on the vertices to obtain
    generating polynomial which seprates the graphs edge weight
    patterns. The number of graceful labelings is given the 
    coefficient of x^((sz^(sz-1)-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=4; f=orbit_non_tree_sum_generating_polynomial(sz)
        sage: f.coefficient(x^((sz^(sz-1)-1)/(sz-1)))
        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initization of list of all trees
    FnLst=[[(0,0)]+[(i+1,l[i]) for i in rg(len(l))] for l in NonTreeFunctionList(sz)]
    # Initialization of the function and the Permutations
    f=0; P=Permutations(sz)
    # Looping through permutations
    for q in P:
        # fixing the permutation index
        p=[q[i]-1 for i in rg(sz)]
        # Looping through the tree classes
        for tp in FnLst:
            # Updating the polynomial f
            f = f + prod(x^(sz^(abs(p[tp[i][1]]-p[tp[i][0]])-1)) for i in rg(1,sz))
    return f

def mDetII(A, tp):
    """
    Computes symbolically the determinant of a 
    square matrix using the sum over permutation
    formula.

    EXAMPLES:

    ::

        sage: sz=4; tp=[(0,0)]+[(0,i) for i in rg(1,sz)]; X=HM(sz,sz,[x^(sz^abs(i-j)) for j in rg(sz) for i in rg(sz)])
        sage: mDetII(X,tp).coefficient(x^((sz^sz-1)/(sz-1)))
        0
        sage: sz=4; mDetII(HM(sz,sz,'a'), [(i,0) for i in rg(sz)])
        0


    AUTHORS:
    - Edinah K. Gnang
    """
    sz=A.n(0)
    # Initializing the permutations
    P = Permutations(sz)
    # Initialization of the function
    f=0
    for q in P:
        # Initialization of the inverse
        qi=q.inverse()
        # fixing the permutation index
        p =[q[i]-1 for i in rg(sz)]
        piv=[qi[i]-1 for i in rg(sz)]
        f=f+q.signature()*prod([A[p[tp[k][0]],p[tp[k][1]]] for k in range(A.nrows())])
    return f

def mPerII(A,tp):
    """
    Computes symbolically the determinant of a 
    square matrix using the sum over permutation
    formula.

    EXAMPLES:

    ::

        sage: sz=4; L=[[(0,0)]+[(i+1,l[i]) for i in rg(len(l))] for l in TreeFunctionList(sz)]; tp=L[1]
        sage: A=HM(sz,sz,'a').elementwise_product(HM(sz,sz,[x^(sz^abs(j-i)) for j in rg(sz) for i in rg(sz)]))
        sage: mPerII(A,tp).coefficient(x^((sz^sz-1)/(sz-1)))
        a00*a12*a20*a30 + a02*a12*a22*a30 + a03*a11*a21*a31 + a03*a13*a21*a33
        sage: sz=4; mPerII(HM(sz,sz,'a'), [(i,0) for i in rg(sz)])
        6*a00*a10*a20*a30 + 6*a01*a11*a21*a31 + 6*a02*a12*a22*a32 + 6*a03*a13*a23*a33
       
 

    AUTHORS:
    - Edinah K. Gnang
    """
    sz=A.n(0)
    # Initializing the permutations
    P = Permutations(sz)
    # Initialization of the function
    f=0
    for q in P:
        # fixing the permutation index
        p=[q[i]-1 for i in rg(sz)]
        f=f+prod([A[p[tp[k][0]],p[tp[k][1]]] for k in range(A.nrows())])
    return f

def mPerIII(A,tp):
    """
    Computes symbolically the modified permanent
    by summing over all permutations and where the
    vertex with the largest label mapping to the image
    of its orginal image


    EXAMPLES:

    ::

        sage: sz=4; L=[[(0,0)]+[(i+1,l[i]) for i in rg(len(l))] for l in TreeFunctionList(sz)]; tp=L[1]
        sage: A=HM(sz,sz,'a').elementwise_product(HM(sz,sz,[x^(sz^abs(j-i)) for j in rg(sz) for i in rg(sz)]))
        sage: mPerII(A,tp).coefficient(x^((sz^sz-1)/(sz-1)))
        a00*a12*a20*a30 + a02*a12*a22*a30 + a03*a11*a21*a31 + a03*a13*a21*a33
        sage: sz=4; mPerIII(HM(sz,sz,'a'), [(i,0) for i in rg(sz)])
        6*a00*a10*a20*a30 + 6*a01*a11*a21*a31 + 6*a02*a12*a22*a32 + 6*a03*a13*a23*a33
       
 

    AUTHORS:
    - Edinah K. Gnang
    """
    sz=A.n(0)
    # Initializing the permutations
    P = Permutations(sz)
    # Initialization of the function
    f=0
    for q in P:
        # fixing the permutation index
        p=[q[i]-1 for i in rg(sz)]
        f=f+prod(A[p[tp[k][0]],p[tp[k][1]]] for k in rg(A.nrows()-1))*A[p[tp[sz-1][0]],p[tp[tp[sz-1][1]][1]]]
    return f

def per_generating_polynomial(sz):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the generating polynomial which seprates the graphs edge weight
    patterns. The number of graceful labelings is given the 
    coefficient of x^((sz^(sz-1)-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=4; f=per_generating_polynomial(sz)

        sage: f.coefficient(x^((sz^(sz-1)-1)/(sz-1)))
        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the variables.
    x=var('x'); X=HM(sz,sz,'a')
    # Initialization of the hypermatrix.
    A=HM(sz,sz,[(x^(sz^(abs(i-j)-1))) for j in rg(sz) for i in rg(sz)])
    # Initialization of the symbolic incidence matrix.
    Hb=HM(sz,binomial(sz,2),'one')
    # Initializing the column index
    clidx=0
    # Loop performing filling the incidence matrix.
    # Note that in an edge (i,j) we have i < j.
    for j in rg(1,sz):
        for i in rg(j):
            Hb[i,clidx]=A[i,j]; Hb[j,clidx]=A[i,j]
            clidx=clidx+1
    # We ground the graph at the vertex labeled 0.
    B=HM(sz-1,binomial(sz,2),[Hb[i,j] for j in rg(binomial(sz,2)) for i in rg(1,sz)])
    # Initializing the list of submatrices
    Lb = B.side_length_subhypermatrices(sz-1)
    # Computing the generating function.
    f=sum(m.per() for m in Lb)
    return f

def per_list_generating_polynomial(sz):
    """
    Goes through all the permutation of n > 1 elements and outputs
    the generating polynomial which seprates the graphs edge weight
    patterns. The edge list of graceful labelings is given the 
    coefficient of x^((sz^(sz-1)-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=4; f=per_list_generating_polynomial(sz)

        sage: f.coefficient(x^((sz^(sz-1)-1)/(sz-1)))
        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the variables.
    x=var('x'); X=HM(sz,sz,'a')
    # Initialization of the hypermatrix.
    A=HM(sz,sz,[(X[i,j]*x^(sz^(abs(i-j)-1))) for j in rg(sz) for i in rg(sz)])
    # Initialization of the symbolic incidence matrix.
    Hb=HM(sz,binomial(sz,2),'one')
    # Initializing the column index
    clidx=0
    # Loop performing filling the incidence matrix.
    # Note that in an edge (i,j) we have i < j.
    for j in rg(1,sz):
        for i in rg(j):
            Hb[i,clidx]=A[i,j]; Hb[j,clidx]=A[i,j]
            clidx=clidx+1
    # We ground the graph at the vertex labeled 0.
    B=HM(sz-1,binomial(sz,2),[Hb[i,j] for j in rg(binomial(sz,2)) for i in rg(1,sz)])
    # Initializing the list of submatrices
    Lb = B.side_length_subhypermatrices(sz-1)
    # Computing the generating function.
    f=sum(m.per() for m in Lb)
    return f

def tree_junction_polynomial(sz):
    """
    The method performs the tree junction starting from Isolated vertices.
    The edge list of graceful labelings is given the coefficient of 
    x^((sz^(sz-1)-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=4; f=tree_junction_polynomial(sz)

        sage: f.coefficient(x^((sz^(sz-1)-1)/(sz-1)))
        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of alphabet
    aL=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',\
    'p', 'q', 'r', 's', 't', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
    # Initializing the list of symbolic matrices
    L=[HM(2,sz,c,'sym') for c in aL[:sz]]
    # Initialization of the list of variables
    Lv=[]
    for t in rg(sz):
        for i in rg(sz):
            for j in rg(i+1):
                Lv.append(L[t][i,j])
    # Initialization of the incidence matrix
    Hb=HM(sz,binomial(sz,2),'zero'); clidx=0
    for i in rg(sz):
        for j in rg(sz):
            if i < j:
                Hb[i,clidx]= sqrt( L[i][i,j]*x^(sz^abs(j-i)) )
                Hb[j,clidx]=-sqrt( L[j][i,j]*x^(sz^abs(i-j)) )
                clidx=clidx+1
    # Grounding at the vertex A
    M=HM(sz-1,binomial(sz,2),[Hb[i,j] for j in rg(binomial(sz,2)) for i in rg(1,sz)])
    # Initialization of the function
    F=expand(L[0][0,0]*x*sum(m.det()^2 for m in M.side_length_subhypermatrices(sz-1)[:binomial(binomial(sz,2),sz-1)-1]))
    return F

def tree_junction_polynomialII(sz):
    """
    The method performs the tree junction starting from Isolated vertices.
    The edge list of graceful labelings is given the coefficient of 
    x^((sz^(sz-1)-1)/(sz-1)).


    EXAMPLES:

    ::


        sage: sz=4; f=tree_junction_polynomialII(sz)

        sage: f.coefficient(x^((sz^(sz-1)-1)/(sz-1)))
        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of alphabet
    aL=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',\
    'p', 'q', 'r', 's', 't', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
    # Initializing the list of symbolic matrices
    L=[HM(2,sz,c,'sym') for c in aL[:sz]]
    # Initialization of the list of variables
    Lv=[]
    for t in rg(sz):
        for i in rg(sz):
            for j in rg(i+1):
                Lv.append(L[t][i,j])
    # Initialization of the incidence matrix
    Hb=HM(sz,binomial(sz,2),'zero'); clidx=0
    for i in rg(sz):
        for j in rg(sz):
            if i < j:
                Hb[i,clidx]= sqrt(sum(L[i][t,i]*x^(sz^abs(t-i)) for t in rg(sz) if t != i))
                Hb[j,clidx]=-sqrt(sum(L[j][t,j]*x^(sz^abs(t-j)) for t in rg(sz) if t != j))
                clidx=clidx+1
    # Grounding at the vertex A
    M=HM(sz-1,binomial(sz,2),[Hb[i,j] for j in rg(binomial(sz,2)) for i in rg(1,sz)])
    # Initialization of the function
    F=expand(L[0][0,0]*x*sum(m.det()^2 for m in M.side_length_subhypermatrices(sz-1)[:binomial(binomial(sz,2),sz-1)-1]))
    return F

def tree_edge_weight_patterns(sz):
    """
    The method returns the list description of the edges weight patterns
    associated with non-isomorphic trees on sz vertices.


    EXAMPLES:

    ::


        sage: sz=4; tree_edge_weight_patterns(sz)
        [[[1, [1, 0, 2, 1], [(0, 0), (1, 2), (2, 0), (3, 0)]],
          [2, [1, 1, 1, 1], [(0, 0), (1, 2), (2, 0), (3, 0)]],
          [3, [1, 2, 0, 1], [(0, 0), (1, 2), (2, 0), (3, 0)]],
          [3, [1, 1, 2, 0], [(0, 0), (1, 2), (2, 0), (3, 0)]],
          [2, [1, 2, 1, 0], [(0, 0), (1, 2), (2, 0), (3, 0)]],
          [1, [1, 3, 0, 0], [(0, 0), (1, 2), (2, 0), (3, 0)]]],
         [[2, [1, 1, 1, 1], [(0, 0), (1, 0), (2, 0), (3, 0)]],
          [2, [1, 2, 1, 0], [(0, 0), (1, 0), (2, 0), (3, 0)]]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    if sz == 2:
        return [[1, [1, 1], [(0, 0), (1, 0)]]]
    else:
        # Initialization of the symbolic matrix of edge weights associated with the complete graph
        X=HM(sz,sz,[x^(sz^abs(j-i)) for j in rg(sz) for i in rg(sz)])
        # Obtaining representatives for equivalence classes of non-isomorphic trees on sz vertices
        FnL=[SP2T(s) for s in TreeSignedPermutationClassesIII(sz)]
        #return [[[f.subs(x==1), Integer((f/f.subs(x==1)).operands()[1]).digits(base=sz)+[0 for i in rg(sz-len(Integer((f/f.subs(x==1)).operands()[1]).digits(base=sz)))], tp] for f in mPerII(X,tp).operands()] for tp in FnL]
        return [[[f.subs(x==1)/T2GraphII(tp[1:],sz).automorphism_group().order(), Integer((f/f.subs(x==1)).operands()[1]).digits(base=sz)+[0 for i in rg(sz-len(Integer((f/f.subs(x==1)).operands()[1]).digits(base=sz)))], tp] for f in mPerII(X,tp).operands()] for tp in FnL]

def tree_edge_weight_polynomial_list(sz):
    """
    The method returns the list description of the edges weight patterns
    associated with non-isomorphic trees on sz vertices.


    EXAMPLES:

    ::


        sage: sz=4; tree_edge_weight_polynomial_list(sz)



    AUTHORS:
    - Edinah K. Gnang
    """
    if sz == 2:
        return [[1, [1, 1], [(0, 0), (1, 0)]]]
    else:
        # Initialization of the symbolic matrix of edge weights associated with the complete graph
        X=HM(sz,sz,[x^(sz^abs(j-i)) for j in rg(sz) for i in rg(sz)])
        # Obtaining representatives for equivalence classes of non-isomorphic trees on sz vertices
        FnL=[SP2T(s) for s in TreeSignedPermutationClassesIII(sz)]
        # Initialization of the total list of mononomials
        S=Set([f/f.subs(x==1) for f in sum(mPerII(X,tp) for tp in FnL).operands()])
        # Initialization of the list
        L=[]; y=var('y')
        for tp in FnL:
            Sf=Set([f/f.subs(x==1) for f in mPerII(X,tp).operands()])
            L.append(sum(prod((y-b)/(a-b) for b in S.difference(Set([a]))) for a in Sf))
        return L

def Monomial2T(mnm, VrbL, sz):
    """
    Outputs the tuple edge list description of the input monomial mnm.
    The variables are obtain by listing the entries of a hypermatrix.
    This function returns pairs of vertex bales corresponding to the
    edges where each pair is in non-decreasing order.


    EXAMPLES:
 
    ::

        sage: sz=4; A=HM(sz,sz,'a'); Tp=Monomial2T(prod(A[0,i] for i in rg(sz)), A.list(), sz)
        sage: Tp
        [(0, 0), (0, 1), (0, 2), (0, 3)]
        


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Getting rid of the coefficient
    f=mnm/(mnm.subs([v==1 for v in VrbL]))
    # Initialization of the dictionary
    EdgeDct=dict([(VrbL[j*sz+i],(min(i,j),max(i,j))) for j in rg(sz) for i in rg(sz)])
    if f in VrbL:
        return [EdgeDct[f]] 
    else:
        Tp=[EdgeDct[g] for g in f.operands()]; Tp.sort()
    return Tp

def Monomial2TII(mnm, VrbL, sz):
    """
    Outputs the tuple edge list description of the input monomial mnm.
    The variables are obtain by listing the entries of a hypermatrix.
    directed graph case.


    EXAMPLES:
 
    ::

        sage: sz=4; A=HM(sz,sz,'a'); Tp=Monomial2TII(prod(A[0,i] for i in rg(sz)), A.list(), sz)
        sage: Tp
        [(0, 0), (0, 1), (0, 2), (0, 3)]
        


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Getting rid of the coefficient
    f=mnm/(mnm.subs([v==1 for v in VrbL]))
    # Initialization of the dictionary
    EdgeDct=dict([(VrbL[j*sz+i],(i,j)) for j in rg(sz) for i in rg(sz)])
    if f in VrbL:
        return [EdgeDct[f]] 
    else:
        Tp=[EdgeDct[g] for g in f.operands()]; Tp.sort()
    return Tp

def is_subgraph_isomorphic(smlTp, Tp, sz):
    """
    Determins if the smlTp is isomorphic to a subgraph in Tp
    both input graph smlTp and Tp are specified as edge list



    EXAMPLES:
 
    ::

        sage: sz=3; smlTp=[(1,1)]; Tp=[(0,0),(0,1),(0,2)]
        sage: is_subgraph_isomorphic(smlTp, Tp, sz)
        True
        


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the list of subgraph
    TpL=TupleSubgraphListII(Tp, len(smlTp))
    # Checking the isomorphism sage is really doing the heavy lifting here
    subGraph=False
    for t in TpL:
        if T2GraphII(t,sz).is_isomorphic(T2GraphII(smlTp,sz)):
            subGraph=True; break
    return subGraph

def is_subgraph_isomorphicII(smlTp, Tp, sz):
    """
    Determins if the smlTp is isomorphic to a subgraph in Tp
    both input graph smlTp and Tp are specified as edge list
    Directed graph case


    EXAMPLES:
 
    ::

        sage: sz=3; smlTp=[(1,1)]; Tp=[(0,0),(0,1),(0,2)]
        sage: is_subgraph_isomorphicII(smlTp, Tp, sz)
        True
        


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the list of subgraph
    TpL=TupleSubgraphListII(Tp, len(smlTp))
    # Checking the isomorphism sage is really doing the heavy lifting here
    subGraph=False
    for t in TpL:
        if T2DiGraphII(t,sz).is_isomorphic(T2DiGraphII(smlTp,sz)):
            subGraph=True; break
    return subGraph

def mper_gaussian_elimination_ReductionHM(Cf, VrbL, Rlts, Tp):
    """
    Outputs the row echelon form of the input second order hypermatrix and the right hand side.
    Does not normalize the rows to ensure that the first non zero entry of a non zero row = 1
    The algorithm perform the reduction assuming that the leading term in each relations
    is a monic powers in a distinct variable as illustrated in the example bellow.
    The algorithm also only retains trees which are subtrees of the input input undirected graph
    Tp specified as a list of ordered edge pairs. If Tp contains all possible edges then all
    trees will be included. It must also be emphasized that the sage graph implementaiton is
    doing the heavy lifiting by solving subgraph isomorphism instance.
    The computation is performed in such a way that the last diagonal entry holds the determinant
    of the whole matrix it is also true that each diagonal entry corresponds to the determinant
    of the corresponding top diagonal block of the matrix. Note that the relation in Rlts are
    assumed to be univariate leading term. This implementaion also gets the sign right for the
    determinant. This implementation is considerably more efficient then the previous one above
    because the extra multiplicative factor is kept minimal.


    EXAMPLES:
 
    ::

        sage: od=2; sz=4 # Initialization of the order and size parameter
        sage: A=HM(sz,sz,'a') # Initialization of the matrix
        sage: v=HM(sz,1,'one') # Initialization of the all one vector
        sage: Hb=HM(od,(A*v).list(),'diag')-A # Initialization of the Laplacian
        sage: t=0; B=HM(sz-1,sz-1,[Hb[i,j] for j in rg(Hb.n(1)) for i in rg(Hb.n(0)) if j!=t if i!=t]) # Grounding at the vertex t.
        sage: M=(HM(od,[A[t,t]]+[1 for i in rg(sz-2)],'diag')*B).expand() # Initialization of the fundamental matrix
        sage: Rh=mper_gaussian_elimination_ReductionHM(M, A.list(), [v^2-v for v in A.list()], [(i,0) for i in rg(sz)]) # Gaussian Elimination
        sage: Rh.printHM()
        [:, :]=
        [a00*a10 + a00*a12 + a00*a13                    -a00*a12                    -a00*a13]
        [                          0                 a00*a10*a20                           0]
        [                          0                           0             a00*a10*a20*a30]       
        sage: Rh=mper_gaussian_elimination_ReductionHM(M, A.list(), [v^2-v for v in A.list()], [(i,i+1) for i in rg(sz-1)]+[(sz-1,sz-1)]) # Gaussian Elimination
        sage: Rh[sz-2,sz-2]
        a00*a13*a21*a30 + a00*a12*a23*a30 + a00*a12*a20*a31 + a00*a10*a23*a31 + a00*a13*a20*a32 + a00*a10*a21*a32
        sage: f=0 # Initialization of the polynomial which will store the modified permanent.
        sage: for t in rg(sz):
        ....:     B=HM(sz-1,sz-1,[Hb[i,j] for j in rg(Hb.n(1)) for i in rg(Hb.n(0)) if j!=t if i!=t]) # Grounding at the vertex t.
        ....:     M=(HM(od,[A[t,t]]+[1 for i in rg(sz-2)],'diag')*B).expand() # Initialization of the fundamental matrix
        ....:     Rh=mper_gaussian_elimination_ReductionHM(M, A.list(), [v^2-v for v in A.list()], [(i,i+1) for i in rg(sz-1)]+[(sz-1,sz-1)]) # Gaussian Elimination
        ....:     f=f+Rh[sz-2,sz-2]
        ....:
        sage: f
        a02*a11*a21*a30 + a00*a13*a21*a30 + a01*a12*a22*a30 + a02*a13*a22*a30 + a01*a11*a23*a30 + a00*a12*a23*a30 + a03*a11*a20*a31 + a00*a12*a20*a31 + a02*a10*a22*a31 + a03*a12*a22*a31 + a00*a10*a23*a31 + a02*a11*a23*a31 + a01*a11*a20*a32 + a00*a13*a20*a32 + a00*a10*a21*a32 + a03*a11*a21*a32 + a03*a10*a22*a32 + a01*a13*a22*a32 + a03*a12*a20*a33 + a01*a13*a20*a33 + a03*a10*a21*a33 + a02*a13*a21*a33 + a02*a10*a23*a33 + a01*a12*a23*a33
        sage: sz=4 # Initialization of the size parametery
        sage: A=HM(sz,sz,'a') # Initialization of the matrix
        sage: v=HM(sz,1,'one') # Initialization of the all one vector
        sage: Hb=HM(od,(A*v).list(),'diag')-A # Initialization of the Laplacian
        sage: g=0 # Initialization of the polynomial which will store the modified permanent.
        sage: for t in rg(sz):
        ....:     B=HM(sz-1,sz-1,[Hb[i,j] for j in rg(Hb.n(1)) for i in rg(Hb.n(0)) if j!=t if i!=t]) # Grounding at the vertex t.
        ....:     M=(HM(od,[A[t,t]]+[1 for i in rg(sz-2)],'diag')*B).expand() # Initialization of the fundamental matrix
        ....:     Rh=mper_gaussian_elimination_ReductionHM(M, A.list(), [v^2-v for v in A.list()], [(i,0) for i in rg(sz)]) # Gaussian Elimination
        ....:     g=g+Rh[sz-2,sz-2]
        ....:
        sage: g
        a00*a10*a20*a30 + a01*a11*a21*a31 + a02*a12*a22*a32 + a03*a13*a23*a33
        sage: sz=5 # Initialization of the size parametery
        sage: A=HM(sz,sz,'a') # Initialization of the matrix
        sage: v=HM(sz,1,'one') # Initialization of the all one vector
        sage: Hb=HM(od,(A*v).list(),'diag')-A # Initialization of the Laplacian
        sage: g=0 # Initialization of the polynomial which will store the modified permanent.
        sage: for t in rg(sz):
        ....:     B=HM(sz-1,sz-1,[Hb[i,j] for j in rg(Hb.n(1)) for i in rg(Hb.n(0)) if j!=t if i!=t]) # Grounding at the vertex t.
        ....:     M=(HM(od,[A[t,t]]+[1 for i in rg(sz-2)],'diag')*B).expand() # Initialization of the fundamental matrix
        ....:     Rh=mper_gaussian_elimination_ReductionHM(M, A.list(), [v^2-v for v in A.list()], [(i,0) for i in rg(sz)]) # Gaussian Elimination
        ....:     g=g+Rh[sz-2,sz-2]
        ....:
        sage: g
        a00*a10*a20*a30*a40 + a01*a11*a21*a31*a41 + a02*a12*a22*a32*a42 + a03*a13*a23*a33*a43 + a04*a14*a24*a34*a44



    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing a copy of the input second order hypermatrices.
    A=Cf.copy()
    # Initialization of the row and column index
    i=0; j=0
    while i < A.n(0) and j < A.n(1):
        while HM(A.n(0)-i, 1, [A[i0,j] for i0 in range(i,A.n(0))]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        if HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))]).is_zero()==False:
            while A[i,j].is_zero(): 
                Ta=HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))])
                # Initializing the cyclic shift permutation matrix
                Id=HM(2, Ta.n(0), 'kronecker')
                P=sum([HM(Ta.n(0),1,[Id[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Id[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                Ta=P*Ta
                for i0 in range(Ta.n(0)):
                    for j0 in range(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
            # Performing the row operations.
            cf1=A[i,j]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if HM(1,A.n(1),[A[r,j0] for j0 in range(A.n(1))]).is_zero():
                    r=r+1
                else:
                    # Initialization of the coefficient
                    cf2=A[r,j]
                    if r==j+1:
                        for j0 in range(j,A.n(1)):
                            # Performing the reduction
                            f=(-cf2*A[i,j0]+cf1*A[r,j0]).numerator().expand()
                            for v in range(len(VrbL)):
                                #f=f.maxima_methods().divide(Rlts[v])[1]
                                for d in range(f.degree(VrbL[v])-Rlts[v].degree(VrbL[v]),-1,-1):
                                    f=expand(fast_reduce(f,[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))],[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))-expand(Rlts[v]*VrbL[v]^d)]))
                            f=SR(sum([trm for trm in f.operands() if is_subgraph_isomorphicII(Monomial2TII(trm.subs(x==1),HM(sz,sz,'a').list(),sz), Tp, sz)]))
                            A[r,j0]=f
                    else:
                        for j0 in range(j,A.n(1)):
                            # Performing the reduction
                            g=expand((-cf2*A[i,j0]+cf1*A[r,j0])/cf1)
                            A[r,j0]=g
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return A

def tpl_vertex_bipartition(tp):
    """
    returns the vertex bipartition associated with a tree
    

    EXAMPLES:

    ::

        sage: tpl_vertex_bipartition([(0, 1), (1, 2), (2, 3), (4, 4)])
        [[4, 2, 0], [3, 1]]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    sz=len(tp)
    # Initialization of the root node
    i=find_sink(tp)
    # Initialization of the bipartite partition
    lL=[i]; rL=[j for j in tpl_pre_image_set(tp,i) if j !=i]
    TplL=copy(lL); TprL=copy(rL)
    while len(lL+rL) < sz:
        TtplL=[]
        for v in TprL:
            TtplL=TtplL+tpl_pre_image_set(tp,v)
        lL=lL+TtplL; TplL=copy(TtplL)
        #print 'lL=',lL
        TtprL=[]
        for v in TplL:
            TtprL=TtprL+tpl_pre_image_set(tp,v)
        rL=rL+TtprL; TprL=copy(TtprL)
        #print 'rL=',rL
    return [lL,rL]

def FindTreeTupleComponents(T):
    """
    Returns a tuple list each of which corresponds to a pair
    made up of a vertex and the associated sink vertex.
    This implementation assume the input is a tree.
    This implementation assume that the given tree is
    rooted at 0


    EXAMPLES:

    ::


        sage: FindTreeTupleComponents([(0, 0), (1, 2), (2, 4), (3, 7), (4, 4), (5, 0), (6, 0), (7, 0)])
        [(0, 0), (1, 4), (2, 4), (3, 0), (4, 4), (5, 0), (6, 0), (7, 0)]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Classifying the vertices
    C=[(0,0)]
    for j in range(1,len(T)):
        i=j; cnt=0
        while T[i][1] != T[i][0] and cnt <= len(T)+1:
            i=T[i][1]; cnt=cnt+1
        if T[i][1] == T[i][0]:
            C.append((j,T[i][0]))
        else:
            raise ValueError, "Expected a tree"
    return C

def graph_function_orbit(T):
    """
    Obtain the set of non isomorphic functional
    directed graphs resulting from switching the
    sink node.


    EXAMPLES:

    ::


        sage: graph_function_orbit([(0, 0), (1, 0), (2, 1), (3, 2), (4, 3)])
        [[(0, 0), (1, 0), (2, 1), (3, 2), (4, 3)],
         [(0, 1), (1, 1), (2, 1), (3, 2), (4, 3)],
         [(0, 1), (1, 2), (2, 2), (3, 2), (4, 3)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=[T]+[switch_sink(T, i) for i in rg(1,len(T))]
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2GraphII(s,len(s)).is_isomorphic(T2GraphII(cL[i],len(T))):
                nwT=False; break
        if nwT==True:
            cL.append(s)
    return cL

def PartialmPerII(A,tq):
    """
    Computes symbolically the partial modified permanent by summing only over
    representatives of the coeset.


    EXAMPLES:

    ::

        sage: sz=4; L=[[(0,0)]+[(i+1,l[i]) for i in rg(len(l))] for l in TreeFunctionList(sz)]; tp=L[1]
        sage: A=HM(sz,sz,'a').elementwise_product(HM(sz,sz,[x^(sz^abs(j-i)) for j in rg(sz) for i in rg(sz)]))
        sage: PartialmPerII(A,tp).coefficient(x^((sz^sz-1)/(sz-1)))
        a00*a12*a20*a30 + a00*a03*a21*a31
        sage: sz=4; PartialmPerII(HM(sz,sz,'a'), [(i,0) for i in rg(sz)])
        a00*a03*a13*a23 + a00*a10*a20*a30 + a00*a01*a21*a31 + a00*a02*a12*a32
       
 

    AUTHORS:
    - Edinah K. Gnang
    """
    sz=A.n(0)
    tp=[(1+tq[i][0], 1+tq[i][1]) for i in rg(len(tq))]
    # Initializing the permutations
    P = Permutations(sz); S=SymmetricGroup(sz)
    # Initializing the graph
    grph=T2GraphII(tp[1:],sz+1)
    # Initializing the automorphism group
    AutGrp=grph.automorphism_group()
    # Initializing representatives of Left coset as strings
    Lcstr=[CstL[0].cycle_string() for CstL in S.cosets(AutGrp)]
    # Loop enumerating the number of graceful labelings
    fctr=0
    # Initialization of the function
    f=0
    for p in P:
        if p.cycle_string() in Lcstr:
            # fixing the permutation index
            q=[p[i]-1 for i in rg(sz)]
            f=f+A[0,0]*prod([A[q[i],q[tp[i][1]-1]] for i in range(1,A.nrows())])
    return f

def Gpoly(sz,sf,c='a'):
    """
    returns the polynomial construction
    The input size correspond to the number of vertices
    The input sf is the index shifting function

    EXAMPLES:

    ::

        sage: Gpoly(3,0)
        a00*a10*a21*x0*x1^2 + a00*a10*a20*x0*x1*x2 + a00*a12*a20*x0*x1*x2
        

    AUTHORS:
    - Edinah K. Gnang
    """
    if sz==0:
        return 0
    elif sz==1:
        # Initializing the symbolic matrix
        A=HM(sz+sf,sz+sf,c)
        #A=HM(sz+sf,sz+sf,[var(c+str(i)+'_'+str(j)) for j in rg(sz+sf) for i in rg(sz+sf)])
        return A[sf,sf]*x0
    elif sz > 1:
        # Initializing the symbolic matrix
        A=HM(sz+sf,sz+sf,c)
        #A=HM(sz+sf,sz+sf,[var(c+str(i)+'_'+str(j)) for j in rg(sz+sf) for i in rg(sz+sf)])
        # Initialization of thelist of variables
        Lx=var_list('x',sz)
        # Initialization of the list of rooted tree tuple functions.
        #FnLst=[[(sf+tp[i][0], sf+tp[i][1]) for i in rg(sz)] for tp in RootedTupleTreeFunctionList(sz)]
        # Initialization of the list of tree tuple functions rooted at 0.
        FnLst=[[(sf+0,sf+0)]+[(sf+i+1,sf+l[i]) for i in rg(len(l))] for l in TreeFunctionList(sz)]
        # Initialization of the list of monomials
        Lm=[prod(Lx[abs(k-j)] for j in rg(sz)) for k in rg(1+ceil((sz-1)/2))]
        # Initiatization of the polynomial
        g=0
        for tp in FnLst:
            if prod(Lx[abs(tp[i][1]-tp[i][0])] for i in rg(sz))==prod(Lx):
                g=g+prod(A[tp[i][0],tp[i][1]]*Lx[abs(tp[i][1]-tp[i][0])] for i in rg(sz))
        return g

def GpolyII(sz,sf,c='a'):
    """
    returns the polynomial construction
    The input size correspond to the number of vertices
    The input sf is the index shifting function

    EXAMPLES:

    ::

        sage: GpolyII(3,0)
        a00*a10*a21*x0*x1^2 + a00*a10*a20*x0*x1*x2 + a00*a12*a20*x0*x1*x2
        

    AUTHORS:
    - Edinah K. Gnang
    """
    if sz==0:
        return 0
    elif sz==1:
        # Initializing the symbolic matrix
        A=HM(sz+sf,sz+sf,c)
        #A=HM(sz+sf,sz+sf,[var(c+str(i)+'_'+str(j)) for j in rg(sz+sf) for i in rg(sz+sf)])
        return A[sf,sf]*x0
    elif sz > 1:
        # Initializing the symbolic matrix
        A=HM(sz+sf,sz+sf,c)
        #A=HM(sz+sf,sz+sf,[var(c+str(i)+'_'+str(j)) for j in rg(sz+sf) for i in rg(sz+sf)])
        # Initialization of thelist of variables
        Lx=var_list('x',sz)
        # Initialization of the list of tree tuple functions rooted at 0.
        FnLst=[[(sf+0,sf+0)]+[(sf+i+1,sf+l[i]) for i in rg(len(l))] for l in TreeFunctionList(sz)]
        #print FnLst
        #print Lm
        # Initiatization of the polynomial
        g=0
        for tp in FnLst:
            Ldg=[prod(Lx[abs(tp[i][1]-tp[i][0])] for i in rg(sz)).degree(v) for v in Lx[1:]]; Ldg.sort(); Ldg.reverse()
            #print 'Ldg = ', Ldg
            #print '[prod(Lx[abs(tp[i][1]-tp[i][0])] for i in rg(sz)).degree(v) for v in Lx[:1]] = ', [prod(Lx[abs(tp[i][1]-tp[i][0])] for i in rg(sz)).degree(v) for v in Lx[:1]]
            if Ldg == [prod(Lx[abs(tp[i][1]-tp[i][0])] for i in rg(sz)).degree(v) for v in Lx[1:]]:
                g=g+prod(A[tp[i][0],tp[i][1]]*Lx[abs(tp[i][1]-tp[i][0])] for i in rg(sz))
        return g

def CountInversions(tp):
    """
    returns the number of inversion occuring in the input
    function specified as a tuple list

    EXAMPLES:

    ::

        sage: CountInversions([(0, 1), (1, 2), (2, 3), (3, 0)])
        3
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the count and number of vertices
    cnt=0; sz=len(tp)
    for j in rg(1,sz):
        for i in rg(j):
            if tp[i][1]>tp[j][1]:
                cnt=cnt+1
    return cnt

def SgnTp(tp):
    """
    returns the sign of the function specified as a tuple list.

    EXAMPLES:

    ::

        sage: SgnTp([(0, 1), (1, 2), (2, 3), (3, 0)])
        -1
        

    AUTHORS:
    - Edinah K. Gnang
    """
    return (-1)^CountInversions(tp)

def GrLmPer(A,tq):
    """
    Computes symbolically the partial modified permanent by summing only over
    representatives of the coeset and graceful labelings.


    EXAMPLES:

    ::

        sage: sz=4; L=[[(0,0)]+[(i+1,l[i]) for i in rg(len(l))] for l in TreeFunctionList(sz)]; tp=L[1]
        sage: A=HM(sz,sz,'a').elementwise_product(HM(sz,sz,[x^(sz^abs(j-i)) for j in rg(sz) for i in rg(sz)]))
        sage: GrLmPer(A,tp)
        a00*a12*a20*a30*x^85 + a02*a12*a22*a30*x^85 + a03*a11*a21*a31*x^85 + a03*a13*a21*a33*x^85
        sage: sz=4; GrLmPer(HM(sz,sz,'a'), [(i,0) for i in rg(sz)])
        a00*a10*a20*a30 + a03*a13*a23*a33
       

    AUTHORS:
    - Edinah K. Gnang
    """
    sz=A.n(0)
    tp=[(1+tq[i][0], 1+tq[i][1]) for i in rg(len(tq))]
    # Initializing the permutations
    P = Permutations(sz); S=SymmetricGroup(sz)
    # Initializing the graph
    grph=T2DiGraphII(tp, sz+1)
    # Initializing the automorphism group
    AutGrp=grph.automorphism_group()
    # Initializing representatives of Left coset as strings
    Lcstr=[CstL[0].cycle_string() for CstL in S.cosets(AutGrp)]
    # Loop enumerating the number of graceful labelings
    fctr=0
    # Initialization of the function
    f=0
    for p in P:
        if p.cycle_string() in Lcstr:
            # Initializing the inverse
            pinv=p.inverse()
            # fixing the permutations index
            q =[p[i]-1 for i in rg(sz)]
            qi=[pinv[i]-1 for i in rg(sz)]
            if len(Set([abs(j-q[tp[qi[j]][1]-1]) for j in rg(sz)]).list())==sz:
                f=f+prod([A[j, q[tp[qi[j]][1]-1]] for j in range(sz)])
    return f

def mDet(A,tq):
    """
    Computes symbolically the partial modified permanent by summing only over
    representatives of the coeset.


    EXAMPLES:

    ::

        sage: sz=4; L=[[(0,0)]+[(i+1,l[i]) for i in rg(len(l))] for l in TreeFunctionList(sz)]; tp=L[1]
        sage: A=HM(sz,sz,'a').elementwise_product(HM(sz,sz,[x^(sz^abs(j-i)) for j in rg(sz) for i in rg(sz)]))
        sage: mDet(A,tp).coefficient(x^((sz^sz-1)/(sz-1)))
        a00*a12*a20*a30 - a03*a11*a21*a31
        sage: sz=4; mDet(HM(sz,sz,'a'), [(i,0) for i in rg(sz)])
        a00*a10*a20*a30 + a01*a11*a21*a31 + a02*a12*a22*a32 + a03*a13*a23*a33
       
 

    AUTHORS:
    - Edinah K. Gnang
    """
    sz=A.n(0)
    tp=[(1+tq[i][0], 1+tq[i][1]) for i in rg(len(tq))]
    # Initializing the permutations
    P = Permutations(sz); S=SymmetricGroup(sz)
    # Initializing the graph
    grph=T2DiGraphII(tp, sz+1)
    # Initializing the automorphism group
    AutGrp=grph.automorphism_group()
    # Initializing representatives of Left coset as strings
    Lcstr=[CstL[0].cycle_string() for CstL in S.cosets(AutGrp)]
    # Loop enumerating the number of graceful labelings
    fctr=0
    # Initialization of the function
    f=0
    for p in P:
        if p.cycle_string() in Lcstr:
            # Initializing the inverse
            pinv=p.inverse()
            # fixing the permutations index
            q =[p[i]-1 for i in rg(sz)]
            qi=[pinv[i]-1 for i in rg(sz)]
            f=f+SgnTp([(i, q[tp[qi[i]][1]-1]) for i in rg(sz)])*prod([A[j, q[tp[qi[j]][1]-1]] for j in range(sz)])
    return f

def Permutation2pseudoTuple(M):
    """
    Returns list of edge tuple desctiption associated with the 
    input permutation. The encoding is based on the diagonals


    EXAMPLES:
    ::
        sage: Permutation2pseudoTuple(HM([[1,0,0],[0,1,0],[0,0,1]]))
        [[(0, 0)], [(0, 1), (1, 0)], [(2, 0), (0, 2)]]
        sage: Permutation2pseudoTuple(HM([[0,0,1],[0,1,0],[1,0,0]]))
        [[(2, 2)], [(1, 2), (2, 1)], [(2, 0), (0, 2)]]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of running edge weight parameter and the Tuple list
    edgw=0; tp=[]
    while M.n(0) > 1:
        for i in rg(M.n(0)):
            if M[0,i] == 1:
                M=SecondOrderSlicer(SecondOrderSlicer(M, [j  for j in rg(1,M.n(0))], 'row'), [k for k in rg(M.n(1)) if k != i], 'col')
                if edgw == 0:
                    tp.append([(i, i+edgw)])
                else:
                    tp.append([(i, i+edgw), (i+edgw, i)])
                edgw=edgw+1
                break
    tp.append([(edgw,0), (0,edgw)])
    return tp

def Tuple2pseudoTuple(tp):
    """
    Returns list of unidrected edge tuple desctiption associated with the 
    input tuple.


    EXAMPLES:
    ::
        sage: Tuple2pseudoTuple([(0, 0), (1, 0), (2, 0)])
        [[(0, 0)], [(1, 0), (0, 1)], [(2, 0), (0, 2)]]
        sage: Tuple2pseudoTuple([(0, 2), (1, 2), (2, 2)])
        [[(2, 2)], [(1, 2), (2, 1)], [(0, 2), (2, 0)]]



    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the pseudo tuple list
    psdT = []
    for i in rg(len(tp)):
        if tp[i][0]==tp[i][1]:
            psdT.append([(tp[i][0],tp[i][1])])
        else:
            psdT.append([(tp[i][0], tp[i][1]), (tp[i][1], tp[i][0])])
    T=[]
    for i in rg(len(tp)):
        for j in rg(len(tp)):
            if abs(psdT[j][0][0]-psdT[j][0][1]) == i:
                T.append(psdT[j])
    return T

def PseudoTuple2Permutation(psdT):
    """
    Returns a permutation matrix associated with list of edge specified as a pseudo tuple edge list.


    EXAMPLES:
    ::
        sage: PseudoTuple2Permutation([[(0, 0)], [(0, 1), (1, 0)], [(2, 0), (0, 2)]])
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        sage: PseudoTuple2Permutation([[(2, 2)], [(1, 2), (2, 1)], [(2, 0), (0, 2)]])
        [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
        sage: PseudoTuple2Permutation([[(1, 1)], [(1, 2), (2, 1)], [(2, 0), (0, 2)]])
        [[0, 1, 0], [0, 0, 1], [1, 0, 0]]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the input matrix
    M=HM(len(psdT),len(psdT),'zero'); M[0,psdT[0][0][0]]=1
    # Initialization of the forbiden index
    allowed_row_index=rg(1,len(psdT))
    #print 'allowed_row_index =', allowed_row_index
    allowed_col_index=[i for i in rg(len(psdT)) if i != psdT[0][0][0]]
    #print 'allowed_col_index =', allowed_col_index
    # Initialization fo the counter
    for indx in rg(1,len(psdT)):
        #print 'indx = ', indx
        M[indx, allowed_col_index[min(psdT[indx][0])]]=1
        # Updating the allowable indices
        allowed_row_index.remove(indx)
        #print 'allowed_row_index =', allowed_row_index
        allowed_col_index.remove(allowed_col_index[min(psdT[indx][0])])
        #print 'allowed_col_index =', allowed_col_index
    return M

def pseudoTuple2Graph(T):
    """
    The method returns an undirected graph object associated with 
    with the pseudo tuple list description of the directed graph

    EXAMPLES:

    ::

        sage: pseudoTuple2Graph([[(0, 0)], [(0, 1), (1, 0)], [(2, 0), (0, 2)]]).degree_sequence()
        [4, 1, 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(len(T))
    return Graph(sum(sum(Id[:,t[0]] * Id[t[1],:] for t in Lt) for Lt in T))

def ConjugateGraphTerm(tp, p, A):
    """
    Performs the conjugation of the graph via the
    permutation p in the graph specified by the
    input tuple edge list tp. Note that p is a list
    of index starting from 1 and not from zero.
    objects returned by Permutations(len(tp)) as a result
    the index are shifted by one the code fixes this.



    EXAMPLES:
    ::


        sage: sz=4; p=Permutations(sz)[3]; A=HM(sz,sz,'a')
        sage: ConjugateGraphTerm([(0,1), (1,2), (2,3), (3,3)], p, A)
        a02*a11*a23*a31


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of primes
    sz=len(tp)
    # Initialization of the inverse
    return prod(A[j,p[tp[p.inverse()[j]-1][1]]-1] for j in rg(sz))

def ConjugateGraphTermII(tp, p, A):
    """
    Performs the conjugation of the graph via the
    permutation p in the graph specified by the
    input tuple edge list tp. Note that p is a list
    of index starting from 1 and not from zero.
    objects returned by Permutations(len(tp)) as a result
    the index are shifted by one the code fixes this.



    EXAMPLES:
    ::


        sage: sz=4; p=[0, 2, 3, 1]; A=HM(sz,sz,'a')
        sage: ConjugateGraphTermII([(0,1), (1,2), (2,3), (3,3)], p, A)
        a02*a11*a23*a31


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of primes
    sz=len(tp)
    return prod(A[p[tp[j][0]],p[tp[j][1]]] for j in rg(sz))

def Polynomial_Construction(tq, A, B):
    """
    Implements the polynomial construction associated
    with the strong form of the graceful labeling conjecture
    it also enables us to check the symmetry property. The 
    input tq is rooted tree tuple.


    EXAMPLES:
    ::


        sage: sz=3; tq=[(i, 0) for i in rg(sz)]; A=HM(sz,sz,'a')
        sage: Polynomial_Construction(tq, A, A)
        (a00*a10*a20 - a01*a11*a20)*(a00*a10*a20 - a00*a12*a20)*(a00*a10*a20 - a00*a10*a21)*(a00*a10*a20 - a02*a11*a21)*(a00*a10*a20 - a02*a10*a22)*(a00*a10*a20 - a01*a12*a22) - (a01*a11*a20 - a01*a11*a21)*(a00*a12*a20 - a01*a11*a21)*(a00*a10*a21 - a01*a11*a21)*(a01*a11*a21 - a02*a11*a21)*(a01*a11*a21 - a02*a10*a22)*(a01*a11*a21 - a01*a12*a22) + (a01*a11*a20 - a02*a12*a22)*(a00*a12*a20 - a02*a12*a22)*(a00*a10*a21 - a02*a12*a22)*(a02*a11*a21 - a02*a12*a22)*(a02*a10*a22 - a02*a12*a22)*(a01*a12*a22 - a02*a12*a22)


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of all functional directed graphs and the corresponding monomials
    sz=len(tq); FnL=NonIsomorphicRootedTreeFunctionList(sz)

    # Initialization of the list of permutations
    P=Permutations(sz); S=SymmetricGroup(sz)
   
    # Initializing the graph
    grph1=T2DiGraphII([(tq[i][0]+1, tq[i][1]+1) for i in rg(sz)], sz+1)
    # Initializing the automorphism group
    AutGrp1=grph1.automorphism_group()
    # Initializing representatives of Left coset as strings
    Lcstr1=[CstL[0].cycle_string() for CstL in S.cosets(AutGrp1)]

    # Initialization of the function which will store the polynomial
    f=0
    for q in P:
        if q.cycle_string() in Lcstr1:
            for tp in FnL:
                # Initializing the graph
                grph2=T2DiGraphII([(tp[i][0]+1, tp[i][1]+1) for i in rg(sz)], sz+1)
                # Initializing the automorphism group
                AutGrp2=grph2.automorphism_group()
                # Initializing representatives of Left coset as strings
                Lcstr2=[CstL[0].cycle_string() for CstL in S.cosets(AutGrp2)]
                # Initialization of the multiplicative factor polynomial
                g=1
                for p in P:
                    if p.cycle_string() in Lcstr2:
                        g=g*( ConjugateGraphTerm(tq, q, A) - ConjugateGraphTerm(tp, p, B) )
                f=f+g 
    return f

def Polynomial_ConstructionII(tq, A, B):
    """
    Implements the polynomial construction
    to check the symmetry property. tq is
    rooted tree tuple. The difference with
    the previous function is that the functional
    star trees are removed from the sum
    in order to investigate subtle differences
    betwen the functional star trees and the 
    functional pseudo star trees.


    EXAMPLES:
    ::


        sage: sz=5; tq=[(i,0) for i in rg(sz-1)]+[(sz-1, sz-2)]; tp=[(sz-1-i,sz-1) for i in rg(sz-1)]+[(0,1)]; A=HM(sz,sz,'a')
        sage: f=Polynomial_Construction(tq, A, A); g=Polynomial_ConstructionII(tp, A, A)
        sage: (f-g).is_zero()
        False



    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of all functional directed graphs and the corresponding monomials
    sz=len(tq); FnL=[T for T in NonIsomorphicRootedTreeFunctionList(sz) \
if not (mPerII(HM(sz,sz,[x^(sz^abs(i-j)) for j in rg(sz) for i in rg(sz)]),T) - \
mPerII(HM(sz,sz,[x^(sz^abs(i-j)) for j in rg(sz) for i in rg(sz)]),[(i,0) for i in rg(sz)])).is_zero()]
    
    # Initialization of the list of permutations
    P=Permutations(sz); S=SymmetricGroup(sz)
   
    # Initializing the graph
    grph1=T2DiGraphII([(tq[i][0]+1, tq[i][1]+1) for i in rg(sz)], sz+1)
    # Initializing the automorphism group
    AutGrp1=grph1.automorphism_group()
    # Initializing representatives of Left coset as strings
    Lcstr1=[CstL[0].cycle_string() for CstL in S.cosets(AutGrp1)]

    # Initialization of the function which will store the polynomial
    f=0
    for q in P:
        if q.cycle_string() in Lcstr1:
            for tp in FnL:
                # Initializing the graph
                grph2=T2DiGraphII([(tp[i][0]+1, tp[i][1]+1) for i in rg(sz)], sz+1)
                # Initializing the automorphism group
                AutGrp2=grph2.automorphism_group()
                # Initializing representatives of Left coset as strings
                Lcstr2=[CstL[0].cycle_string() for CstL in S.cosets(AutGrp2)]
                # Initialization of the multiplicative factor polynomial
                g=1
                for p in P:
                    if p.cycle_string() in Lcstr2:
                        g=g*( ConjugateGraphTerm(tq, q, A) - ConjugateGraphTerm(tp, p, B) )
                f=f+g 
    return f

def comp(tp, tq):
    """
    Performs the composition of functions. The inputs to this method
    are two tuple lists of the same size. The current version does
    not check that the list are in fact of the same size.


    EXAMPLES:

    ::

        sage: tp=[(0,1), (1,2), (2,3), (3,0)]; tq=[(0,0), (1,0), (2,0), (3,0)]
        sage: comp(tp, tq)
        [(0, 1), (1, 1), (2, 1), (3, 1)]
       

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter
    return [(i,tp[tq[i][1]][1]) for i in rg(len(tp))]


def tree_conjugation(sz):
    """
    Runs the conjugation algorithm to produce the orbit
    of trees which admit a graceful labeling. The algorithm
    returns the list of monomial obtained in the orbit
    and the polynomial on which the conjugation is being
    performed


    EXAMPLES:

    ::

        sage: sz=4; [LtRm, Fa] = tree_conjugation(sz)
        sage: len(LtRm)
        64
        sage: Fa
        a00*a10*a20*a30*x^85 + a01*a11*a20*a30*x^85 + a00*a12*a20*a30*x^85 + a02*a11*a21*a30*x^85 + a00*a13*a21*a30*x^85 + a02*a10*a22*a30*x^85 + a02*a12*a22*a30*x^85 + a00*a13*a23*a30*x^85 + a03*a11*a21*a31*x^85 + a03*a12*a22*a31*x^85 + a03*a11*a23*a31*x^85 + a03*a13*a22*a32*x^85 + a03*a10*a20*a33*x^85 + a03*a12*a20*a33*x^85 + a03*a13*a21*a33*x^85 + a03*a13*a23*a33*x^85 + a00*a10*a21*a31*x^25 + a01*a11*a21*a31*x^25 + a01*a12*a22*a31*x^25 + a00*a10*a23*a31*x^25 + a01*a11*a23*a31*x^25 + a00*a10*a20*a32*x^25 + a01*a11*a20*a32*x^25 + a00*a12*a20*a32*x^25 + a02*a11*a21*a32*x^25 + a02*a10*a22*a32*x^25 + a02*a12*a22*a32*x^25 + a01*a13*a22*a32*x^25 + a01*a13*a21*a33*x^25 + a02*a10*a23*a33*x^25 + a02*a12*a23*a33*x^25 + a01*a13*a23*a33*x^25
       

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the order parameter
    od=2
    # Initialization of the variables
    x=var('x')
    # Initialization of the second order hypermatrices
    TpA=HM(sz,sz,'a'); X=HM(sz,sz,[x^(sz^abs(j-i)) for j in rg(sz) for i in rg(sz)])
    A=TpA.elementwise_product(X)
    # Initialization of the list of monomials associated with induced edge labelings of the star
    Lm=[prod(x^(sz^abs(t-i)) for i in rg(sz)) for t in rg(ceil(sz/2))]
    #print 'p = ', rg(1,sz+1)
    # Initialization of the polynomials and the list which will store the trees
    Fa=0; La=[]; LtRm=[]
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]
    # Initialization of the identity matrix
    Id=HM(od,sz,'kronecker')
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Initialization of the adjacency martrix of the functional graph
        B = Matrix2HM(sum([Id.matrix()[:,j]*Id.matrix()[f[j],:] for j in rg(sz)]))
        # Initialization of the dierected Laplacian
        lB= HM(od,(B*HM(B.nrows(),1,'one')).list(),'diag')-B
        # Initialization of the list of sumbratrices
        Lmtr=[HM(sz-1,sz-1,[lB[u,v] for u in range(sz) for v in range(sz) if u!=t if v!=t]) for t in range(sz)]
        # Testing treeness
        if sum(B[t,t]*Lmtr[t].det() for t in range(sz)) == 1:
            # Initialization of the tree truple description
            tq=[(i, f[i]) for i in range(sz)]
            if prod(x^(sz^abs(tq[i][1]-tq[i][0])) for i in rg(sz)) in Lm:
                La.append(tq)
                LtRm.append(prod(A[tq[i][0],tq[i][1]] for i in rg(sz)))
                Fa=Fa+prod(A[tq[i][0],tq[i][1]] for i in rg(sz))
    #print 'New monomials = ', LtRm
    #print 'Number of terms thus far = ', len(Fa.operands())
    # Initialization of the list of permutations
    P=Permutations(sz)
    # Summing over representative of the quotient by the automorphism group of the Fa
    # which have the same automorphism group made up only of the identity permutation.
    # We exclude the identity from the sum.
    F=Fa
    for p in P[1:]:
        #print 'p = ', p
        Ga=sum(ConjugateGraphTerm(tq,p,A) for tq in La).subs([trm==0 for trm in LtRm])
        if Ga.is_zero() and len(LtRm)==sz^(sz-1):
            break
        elif Ga.subs([v==1 for v in [x]+TpA.list()]) == 1:
            # Updading LtRm with the new monomials
            LtRm=LtRm+[Ga]
            #print 'New monomial = ', [Ga.subs(x==1)]
    	# Updating the sum
            F=F+Ga
            #print 'Number of terms thus far = ', len(F.operands())
        else:
            # Updading LtRm with the new monomials
            LtRm=Set(LtRm+[mn for mn in Ga.operands()]).list()
            #print 'New monomials = ', Set([mn for mn in Ga.operands()]).list()
            # Updating the generator set with the new tree monomial summands.
            F=F+Ga
            #print 'Number of terms thus far = ', len(F.operands())
    return [LtRm, Fa]

def tree_conjugationII(sz):
    """
    Runs the conjugation algorithm to produce the orbit
    of trees which admit a graceful labeling. The algorithm
    returns the list of monomial obtained in the orbit
    and the polynomial on which the conjugation is being
    performed


    EXAMPLES:

    ::

        sage: sz=3; [LtRm, Fa] = tree_conjugationII(sz)
        sage: len(LtRm)
        9
        sage: sum(LtRm)
        a00*a10*a20*x^13 + a01*a11*a20*x^13 + a00*a12*a20*x^13 + a02*a11*a21*x^13 + a02*a10*a22*x^13 + a02*a12*a22*x^13 + a00*a10*a21*x^7 + a01*a11*a21*x^7 + a01*a12*a22*x^7
        sage: Fa
        a00*a10*a20*x^13 + a01*a11*a20*x^13 + a00*a12*a20*x^13 + a02*a11*a21*x^13 + a02*a10*a22*x^13 + a02*a12*a22*x^13

       

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the order parameter
    od=2
    # Initialization of the variables
    x=var('x')
    # Initialization of the second order hypermatrices
    TpA=HM(sz,sz,'a'); X=HM(sz,sz,[x^(sz^abs(j-i)) for j in rg(sz) for i in rg(sz)])
    A=TpA.elementwise_product(X)
    #print 'p = ', rg(1,sz+1)
    # Initialization of the polynomials and the list which will store the trees
    Fa=0; La=[]; LtRm=[]
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Initialization of the adjacency martrix of the functional graph
        B=Matrix2HM(sum([Id.matrix()[:,j]*Id.matrix()[f[j],:] for j in rg(sz)]))
        # Initialization of the dierected Laplacian
        lB=HM(od,(B*HM(B.nrows(),1,'one')).list(),'diag')-B
        # Initialization of the list of sumbratrices
        Lmtr=[HM(sz-1,sz-1,[lB[u,v] for u in rg(sz) for v in rg(sz) if u!=t if v!=t]) for t in rg(sz)]
        # Testing treeness
        if sum(B[t,t]*Lmtr[t].det() for t in rg(sz)) == 1:
            # Initialization of the tree truple description
            tq=[(i,f[i]) for i in range(sz)]
            if prod(x^(sz^abs(tq[i][1]-tq[i][0])) for i in rg(sz)) == prod(x^(sz^abs(0-i)) for i in rg(sz)):
                La.append(tq)
                LtRm.append(prod(A[tq[i][0],tq[i][1]] for i in rg(sz)))
                Fa=Fa+prod(A[tq[i][0],tq[i][1]] for i in rg(sz))
    #print 'New monomials = ', LtRm
    #print 'Number of terms thus far = ', len(Fa.operands())
    # Initialization of the list of permutations
    P=Permutations(sz)
    # Summing over representative of the quotient by the automorphism group of the Fa
    # which have the same automorphism group made up only of the identity permutation.
    # We exclude the identity from the sum.
    F=Fa
    for p in P[1:]:
        #print 'p = ', p
        Ga=sum(ConjugateGraphTerm(tq,p,A) for tq in La).subs([trm==0 for trm in LtRm])
        if Ga.is_zero() and len(LtRm)==sz^sz:
            break
        elif Ga.subs([v==1 for v in [x]+TpA.list()]) == 1:
            # Updading LtRm with the new monomials
            LtRm=LtRm+[Ga]
            #print 'New monomial = ', [Ga.subs(x==1)]
    	# Updating the sum
            F=F+Ga
            #print 'Number of terms thus far = ', len(F.operands())
        else:
            # Updading LtRm with the new monomials
            LtRm=Set(LtRm+[mn for mn in Ga.operands()]).list()
            #print 'New monomials = ', Set([mn for mn in Ga.operands()]).list()
            # Updating the generator set with the new tree monomial summands.
            F=F+Ga
            #print 'Number of terms thus far = ', len(F.operands())
    return [LtRm, Fa]

def graph_conjugation(sz):
    """
    Runs the conjugation algorithm to produce the orbit
    of trees which admit a graceful labeling. The algorithm
    returns the list of monomial obtained in the orbit
    and the polynomial on which the conjugation is being
    performed


    EXAMPLES:

    ::

        sage: sz=3; [LtRm, Fa] = graph_conjugation(sz)
        sage: len(LtRm)
        9
        sage: sum(LtRm)
        a00*a10*a20*x^13 + a01*a11*a20*x^13 + a00*a12*a20*x^13 + a02*a11*a21*x^13 + a02*a10*a22*x^13 + a02*a12*a22*x^13 + a00*a10*a21*x^7 + a01*a11*a21*x^7 + a01*a12*a22*x^7
        sage: Fa
        a00*a10*a20*x^13 + a01*a11*a20*x^13 + a00*a12*a20*x^13 + a02*a11*a21*x^13 + a02*a10*a22*x^13 + a02*a12*a22*x^13
        sage: sz=4; [LtRm1, Fa] = tree_conjugation(sz)
        sage: L=Set(LtRm).difference(Set(LtRm1)).list(); L
        [a00*a12*a23*a31*x^25,
         a02*a10*a21*a33*x^25,
         a00*a13*a21*a32*x^25,
         a03*a11*a20*a32*x^85,
         a02*a11*a23*a30*x^85,
         a01*a13*a22*a30*x^85,
         a03*a10*a22*a31*x^85,
         a01*a12*a20*a33*x^25]
       

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the order parameter
    od=2
    # Initialization of the variables
    x=var('x')
    # Initialization of the second order hypermatrices
    TpA=HM(sz,sz,'a'); X=HM(sz,sz,[x^(sz^abs(j-i)) for j in rg(sz) for i in rg(sz)])
    A=TpA.elementwise_product(X)
    #print 'p = ', rg(1,sz+1)
    # Initialization of the polynomials and the list which will store the trees
    Fa=0; La=[]; LtRm=[]
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Initialization of the tree truple description
        tq=[(i,f[i]) for i in range(sz)]
        if prod(x^(sz^abs(tq[i][1]-tq[i][0])) for i in rg(sz)) == prod(x^(sz^abs(0-i)) for i in rg(sz)):
            La.append(tq)
            LtRm.append(prod(A[tq[i][0],tq[i][1]] for i in rg(sz)))
            Fa=Fa+prod(A[tq[i][0],tq[i][1]] for i in rg(sz))
    #print 'New monomials = ', LtRm
    #print 'Number of terms thus far = ', len(Fa.operands())
    # Initialization of the list of permutations
    P=Permutations(sz)
    # Summing over representative of the quotient by the automorphism group of the Fa
    # which have the same automorphism group made up only of the identity permutation.
    # We exclude the identity from the sum.
    F=Fa
    for p in P[1:]:
        #print 'p = ', p
        Ga=sum(ConjugateGraphTerm(tq,p,A) for tq in La).subs([trm==0 for trm in LtRm])
        if Ga.is_zero() and len(LtRm)==sz^sz:
            break
        elif Ga.subs([v==1 for v in [x]+TpA.list()]) == 1:
            # Updading LtRm with the new monomials
            LtRm=LtRm+[Ga]
            #print 'New monomial = ', [Ga.subs(x==1)]
    	# Updating the sum
            F=F+Ga
            #print 'Number of terms thus far = ', len(F.operands())
        else:
            # Updading LtRm with the new monomials
            LtRm=Set(LtRm+[mn for mn in Ga.operands()]).list()
            #print 'New monomials = ', Set([mn for mn in Ga.operands()]).list()
            # Updating the generator set with the new tree monomial summands.
            F=F+Ga
            #print 'Number of terms thus far = ', len(F.operands())
    return [LtRm, Fa]

def creative_stabilizing(sz):
    """
    This function implements the main argument of the proof.
    The function starts by collecting gracefully labeled trees
    on sz vertices and proceed to generate by conjugation 
    the orbit of all graceful trees on sz+1 vertices.


    EXAMPLES:

    ::

        sage: sz=2; [LtRm, Fa] = creative_stabilizing(sz)
        sage: len(LtRm)
        9
        sage: sum(LtRm)
        a00*a10*a20*x^13 + a01*a11*a20*x^13 + a00*a12*a20*x^13 + a02*a11*a21*x^13 + a02*a10*a22*x^13 + a02*a12*a22*x^13 + a00*a10*a21*x^7 + a01*a11*a21*x^7 + a01*a12*a22*x^7
        sage: Fa
        a00*a10*a20*x^13 + a01*a11*a20*x^13 + a00*a12*a20*x^13 + a02*a11*a21*x^13 + a02*a10*a22*x^13 + a02*a12*a22*x^13
        sage: sz=3; [LtRm, Fa] = creative_stabilizing(sz)
        sage: TpA=HM(sz+1,sz+1,'a'); X=HM(sz+1,sz+1,[x^((sz+1)^abs(j-i)) for j in rg(sz+1) for i in rg(sz+1)]); A=TpA.elementwise_product(X)
        sage: sum(LtRm)-sum(prod(A[tq[j][0],tq[j][1]] for j in rg(sz+1)) for tq in RootedTupleTreeFunctionList(sz+1))
        0
       

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the order parameter
    od=2
    # Initialization of the variables
    x=var('x')
    # Initialization of the second order hypermatrices
    TpA=HM(sz+1,sz+1,'a'); X=HM(sz+1,sz+1,[x^((sz+1)^abs(j-i)) for j in rg(sz+1) for i in rg(sz+1)])
    A=TpA.elementwise_product(X)
    #print 'p = ', rg(1,sz+1)
    # Initialization of the polynomials and the list which will store the trees
    Fa=0; La=[]; LtRm=[]
    # Initialization of the lists
    l=[sz for i in rg(sz)]; Lf=[]
    # Initialization of the identity matrix
    Id=HM(od,sz,'kronecker')
    # Main loop performing the collecting the functions.
    for i in rg(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in rg(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Initialization of the adjacency martrix of the functional graph
        B=Matrix2HM(sum([Id.matrix()[:,j]*Id.matrix()[f[j],:] for j in rg(sz)]))
        # Initialization of the dierected Laplacian
        lB=HM(od,(B*HM(B.nrows(),1,'one')).list(),'diag')-B
        # Initialization of the list of sumbratrices
        Lmtr=[HM(sz-1,sz-1,[lB[u,v] for u in rg(sz) for v in rg(sz) if u!=t if v!=t]) for t in rg(sz)]
        # Testing treeness
        if sum(B[t,t]*Lmtr[t].det() for t in rg(sz)) == 1:
            # Initialization of the tree truple description
            tq=[(j,f[j]) for j in rg(sz)]+[(sz,0)]
            #if prod(x^((sz+1)^abs(tq[j][1]-tq[j][0])) for j in rg(sz+1)) == prod(x^((sz+1)^abs(0-j)) for j in rg(sz+1)):
            if prod(x^((sz+1)^abs(tq[j][1]-tq[j][0])) for j in rg(sz+1)) == prod(x^((sz+1)^abs(0-j)) for j in rg(sz+1)):
                # Initializing the corresponding tree rooted at sz.
                tp=switch_sink(tq,sz)
                # Appending the extended gracefully labeled trees
                La.append(tq) 
                LtRm.append(prod(A[tq[j][0],tq[j][1]] for j in rg(sz+1)))
                if ( prod(A[tp[j][0],tp[j][1]] for j in rg(sz+1)) in LtRm ) == False:
                    La.append(tp)
                    LtRm.append(prod(A[tp[j][0],tp[j][1]] for j in rg(sz+1)))
                    Fa=Fa + prod(A[tq[j][0],tq[j][1]] for j in rg(sz+1)) + prod(A[tp[j][0],tp[j][1]] for j in rg(sz+1))
                else:
                    Fa=Fa+prod(A[tq[j][0],tq[j][1]] for j in rg(sz+1))
    #print 'New monomials = ', LtRm
    #print 'Number of terms thus far = ', len(Fa.operands())
    # Initialization of the list of permutations
    P=Permutations(sz+1)
    # Summing over representative of the quotient by the automorphism group of the Fa
    # which have the same automorphism group made up only of the identity permutation.
    # We exclude the identity from the sum.
    F=Fa
    for p in P[1:]:
        #print 'p = ', p
        Ga=sum(ConjugateGraphTerm(tq,p,A) for tq in La).subs([trm==0 for trm in LtRm])
        if Ga.is_zero() and len(LtRm)==(sz+1)^sz:
            break
        elif Ga.subs([v==1 for v in [x]+TpA.list()]) == 1:
            # Updading LtRm with the new monomials
            LtRm=LtRm+[Ga]
            #print 'New monomial = ', [Ga.subs(x==1)]
    	# Updating the sum
            F=F+Ga
            #print 'Number of terms thus far = ', len(F.operands())
        else:
            # Updading LtRm with the new monomials
            LtRm=Set(LtRm+[mn for mn in Ga.operands()]).list()
            #print 'New monomials = ', Set([mn for mn in Ga.operands()]).list()
            # Updating the generator set with the new tree monomial summands.
            F=F+Ga
            #print 'Number of terms thus far = ', len(F.operands())
    return [LtRm, Fa]

def undirected_graceful_graph_polynomial(sz):
    """
    This function implements the polynomial construction
    which list all undirected graphs having sz-1 non loop
    edges which admit a graceful labeling. The coefficient
    of each term counts ???


    EXAMPLES:

    ::

        sage: sz=3; factor(undirected_graceful_graph_polynomial(sz))
        4*(a01*a02 + a01*a12 + a02*a12)*(a00 + a11 + a22)
        sage: sz=4; factor(undirected_graceful_graph_polynomial(sz))
        4*(3*a01*a02*a03 + 3*a01*a02*a12 + a01*a03*a12 + a02*a03*a12 + a01*a02*a13 + 3*a01*a03*a13 + a02*a03*a13 + 3*a01*a12*a13 + a02*a12*a13 + a03*a12*a13 + a01*a02*a23 + a01*a03*a23 + 3*a02*a03*a23 + a01*a12*a23 + 3*a02*a12*a23 + a03*a12*a23 + a01*a13*a23 + a02*a13*a23 + 3*a03*a13*a23 + 3*a12*a13*a23)*(a00 + a11 + a22 + a33)


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the order parameter
    od=2
    # Initialization of the matrix
    A=HM(2*sz-1,2*sz-1,'zero').fill_with(HM(od,sz,'a','sym'))
    # Going through all the permutations
    h=0
    for p in Permutations(sz):
        # Shifting back the index
        q=[p[i]-1 for i in rg(sz)]+[i for i in rg(sz,2*sz-1)]
        h=h+sum(prod(A[q[sz-1-t[0]+t[1]],q[t[1]]] for t in tp) for tp in TupleFunctionList(sz))
    return h

def undirected_subgraph_polynomial(sz):
    """
    This function implements the polynomial construction
    which list all undirected graphs having sz-1 non loop
    edges and exactly one loop edge. 


    EXAMPLES:

    ::

        sage: sz=3; undirected_subgraph_polynomial(sz)
        (a01*a02 + a01*a12 + a02*a12)*(a00 + a11 + a22)
        sage: sz=4; undirected_subgraph_polynomial(sz)
        (a01*a02*a03 + a01*a02*a12 + a01*a03*a12 + a02*a03*a12 + a01*a02*a13 + a01*a03*a13 + a02*a03*a13 + a01*a12*a13 + a02*a12*a13 + a03*a12*a13 + a01*a02*a23 + a01*a03*a23 + a02*a03*a23 + a01*a12*a23 + a02*a12*a23 + a03*a12*a23 + a01*a13*a23 + a02*a13*a23 + a03*a13*a23 + a12*a13*a23)*(a00 + a11 + a22 + a33)
       

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the order parameter
    od=2
    # Initialization of the matrix
    A=HM(od,sz,'a','sym')
    return sum(A[i,i] for i in rg(sz))*sum(prod(st) for st in Set([A[i,j] for j in rg(sz) for i in rg(sz) if i<j]).subsets(sz-1))

def undirected_ungraceful_graph_polynomial(sz):
    """
    This function implements the polynomial construction
    which list all undirected graphs having sz-1 non loop
    edges which do not admit a graceful labeling.


    EXAMPLES:

    ::

        sage: sz=4; undirected_ungraceful_graph_polynomial(sz)
        0
        sage: sz=5; factor(undirected_ungraceful_graph_polynomial(sz))
        a00*a04*a12*a13*a23 + a04*a11*a12*a13*a23 + a00*a02*a03*a14*a23 + a00*a01*a04*a14*a23 + a02*a03*a11*a14*a23 + a01*a04*a11*a14*a23 + a04*a12*a13*a22*a23 + a02*a03*a14*a22*a23 + a01*a04*a14*a22*a23 + a00*a01*a03*a13*a24 + a00*a02*a04*a13*a24 + a01*a03*a11*a13*a24 + a02*a04*a11*a13*a24 + a00*a03*a12*a14*a24 + a03*a11*a12*a14*a24 + a01*a03*a13*a22*a24 + a02*a04*a13*a22*a24 + a03*a12*a14*a22*a24 + a04*a12*a13*a23*a33 + a02*a03*a14*a23*a33 + a01*a04*a14*a23*a33 + a01*a03*a13*a24*a33 + a02*a04*a13*a24*a33 + a03*a12*a14*a24*a33 + a00*a01*a02*a12*a34 + a00*a03*a04*a12*a34 + a01*a02*a11*a12*a34 + a03*a04*a11*a12*a34 + a00*a02*a13*a14*a34 + a02*a11*a13*a14*a34 + a01*a02*a12*a22*a34 + a03*a04*a12*a22*a34 + a02*a13*a14*a22*a34 + a00*a01*a23*a24*a34 + a01*a11*a23*a24*a34 + a01*a22*a23*a24*a34 + a01*a02*a12*a33*a34 + a03*a04*a12*a33*a34 + a02*a13*a14*a33*a34 + a01*a23*a24*a33*a34 + a04*a12*a13*a23*a44 + a02*a03*a14*a23*a44 + a01*a04*a14*a23*a44 + a01*a03*a13*a24*a44 + a02*a04*a13*a24*a44 + a03*a12*a14*a24*a44 + a01*a02*a12*a34*a44 + a03*a04*a12*a34*a44 + a02*a13*a14*a34*a44 + a01*a23*a24*a34*a44
       

    AUTHORS:
    - Edinah K. Gnang
    """
    A=HM(sz,sz,'a'); La=[A[i,j] for j in rg(sz) for i in rg(sz) if i<=j]
    # Initialization of the polynomial construction of undirected graceful polynomial
    f=undirected_graceful_graph_polynomial(sz)
    # Scalling the non vanishing coefficient to 1
    F=sum(trm/trm.subs([v==1 for v in La]) for trm in f.operands())
    return expand(undirected_subgraph_polynomial(sz))-F

def composition_generator_trees(sz):
    """
    This function goes through all rooted trees and determines which ones
    can not be generated by a self composition.


    EXAMPLES:
    ::


        sage: sz=4; composition_generator_trees(sz)
        [[(0, 1), (1, 3), (2, 1), (3, 3)],
         [(0, 0), (1, 2), (2, 0), (3, 2)],
         [(0, 3), (1, 1), (2, 1), (3, 2)],
         [(0, 3), (1, 0), (2, 0), (3, 3)],
         [(0, 1), (1, 1), (2, 3), (3, 0)],
         [(0, 2), (1, 3), (2, 1), (3, 3)],
         [(0, 3), (1, 3), (2, 2), (3, 2)],
         [(0, 2), (1, 0), (2, 3), (3, 3)],
         [(0, 1), (1, 1), (2, 0), (3, 2)],
         [(0, 0), (1, 2), (2, 3), (3, 0)],
         [(0, 0), (1, 2), (2, 0), (3, 1)],
         [(0, 0), (1, 3), (2, 1), (3, 0)],
         [(0, 0), (1, 0), (2, 3), (3, 1)],
         [(0, 3), (1, 2), (2, 2), (3, 1)],
         [(0, 2), (1, 1), (2, 1), (3, 0)],
         [(0, 1), (1, 3), (2, 2), (3, 2)],
         [(0, 1), (1, 2), (2, 2), (3, 0)],
         [(0, 3), (1, 0), (2, 1), (3, 3)],
         [(0, 3), (1, 0), (2, 2), (3, 2)],
         [(0, 0), (1, 0), (2, 1), (3, 2)],
         [(0, 1), (1, 1), (2, 0), (3, 0)],
         [(0, 2), (1, 2), (2, 3), (3, 3)],
         [(0, 2), (1, 0), (2, 2), (3, 0)],
         [(0, 3), (1, 1), (2, 3), (3, 1)],
         [(0, 2), (1, 3), (2, 2), (3, 0)],
         [(0, 3), (1, 2), (2, 0), (3, 3)],
         [(0, 3), (1, 1), (2, 0), (3, 1)],
         [(0, 0), (1, 3), (2, 3), (3, 0)],
         [(0, 2), (1, 1), (2, 3), (3, 1)],
         [(0, 1), (1, 2), (2, 2), (3, 1)],
         [(0, 0), (1, 3), (2, 0), (3, 2)],
         [(0, 2), (1, 0), (2, 2), (3, 1)],
         [(0, 1), (1, 2), (2, 3), (3, 3)],
         [(0, 0), (1, 0), (2, 1), (3, 1)],
         [(0, 2), (1, 1), (2, 1), (3, 2)],
         [(0, 1), (1, 3), (2, 0), (3, 3)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the variable
    z=var('z')
    # Initialization of the set of functions
    Sf=Set([expand(BasicLagrangeInterpolation(tp,z)) for tp in RootedTupleTreeFunctionList(sz)])
    # Initialization of Set of functions obtainaible by composition
    Sff=Set([expand(composition(f,z,2,sz)) for f in Sf])
    # Initialization of the generator set
    Sg=Sf.difference(Sff)
    # Converting SG back into tuples
    return [[(i,f.subs(z==i)) for i in rg(sz)] for f in Sg] 

def composition_generator_tree_classes(sz):
    """
    This function goes through all rooted trees and determines which ones
    can not be generated by a self composition up to isomorphism classes.


    EXAMPLES:
    ::


        sage: sz=5; composition_generator_tree_classes(sz)
        [[(0, 0), (1, 3), (2, 0), (3, 2), (4, 0)],
         [(0, 2), (1, 1), (2, 1), (3, 0), (4, 0)],
         [(0, 1), (1, 1), (2, 0), (3, 4), (4, 0)],
         [(0, 0), (1, 2), (2, 0), (3, 2), (4, 2)],
         [(0, 1), (1, 1), (2, 4), (3, 2), (4, 0)]]
        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations
    L=[T2SP(T) for T in composition_generator_trees(sz)]
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for s in L:
        nwT=True
        for i in range(len(cL)):
            if T2DiGraphII([(j,j+s[j]) for j in range(sz)],sz).is_isomorphic(T2DiGraphII([(j,j+cL[i][j]) for j in range(sz)],sz)):
                nwT=False
                break
        if nwT==True:
            cL.append(s)
    return [SP2T(sg) for sg in cL]

def composition_primitive_trees(sz):
    """
    This function goes through all rooted trees and determines which ones
    can not be generated by a combination of self composition and sink switching.


    EXAMPLES:
    ::


        sage: sz=4; composition_primitive_trees(sz)
        []
        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the variable
    z=var('z')
    # Initialization of the set of functions
    Sf=Set([expand(BasicLagrangeInterpolation(tp,z)) for tp in RootedTupleTreeFunctionList(sz)])
    # Initialization of Set of functions obtainaible by composition
    Sff=Set([])
    for f in Sf:
        for sink in rg(sz):
            # Initialization of the composition
            ff=expand(composition(f,z,2,sz))
            # Initialization of the tuple encoding
            Tf=[(i,ff.subs(z==i)) for i in rg(sz)]
            # Cheking the primitivity criteria
            if not expand(BasicLagrangeInterpolation(switch_sink(Tf,sink),z)) in Sff:
                Sff=Sff.union(Set([expand(BasicLagrangeInterpolation(switch_sink(Tf,sink),z))]))
                print 'Sff.cardinality()=', Sff.cardinality()
    # Initialization of the generator set
    if Sff.cardinality()==sz^(sz-2):
        return []
    else:
        # Converting SG back into tuples
        return [[(i,f.subs(z==i)) for i in rg(sz)] for f in Sf.difference(Sff)] 

def GeneratePartition(sz):
    """
    Creates Partition to be used for computing the permanent.
    This function follows a maple implementation suggested
    Harry Crane


    EXAMPLES:
    ::


        sage: GeneratePartition(5)
        [[1, 1, 1], [1, 1, 2], [1, 2, 1], [1, 2, 2], [1, 2, 3]]
        


    AUTHORS:
    - Edinah K. Gnang, Harry Crane
    """
    N = [[1]]
    if sz > 1:
        for i in range(1,sz):
            M = N
            r = len(M)
            N = []
            for j in range(r):
                mx = max(M[j])
                for k in range(1,mx+2):
                    N = N + [M[j]+[k]]
    return N

def Partition2Matrix(part):
    """
    Converts a partition into matrices. This function follows a maple implementation suggested
    Harry Crane


    EXAMPLES:
    ::


        sage: Partition2Matrix([1, 2, 1])
        [1 0 1]
        [0 1 0]
        [1 0 1]



    AUTHORS:
    - Edinah K. Gnang, Harry Crane
    """
    M = max(part)
    N = len(part)
    B = zero_matrix(N,N)
    for i in range(1,M+1):
        d  = zero_matrix(N,1)
        for j in range(N):
            if part[j] == i:
                d[j,0] = 1
        B = B + d*d.transpose()
    return B

def PermComp(A):
    """
    Partition based formula for computing the permanent.
    This function follows a maple implementation suggested
    Harry Crane


    EXAMPLES:
    ::


        sage: PermComp3(HM(3,3,'a').matrix())
        a02*a11*a20 + a01*a12*a20 + a02*a10*a21 + a00*a12*a21 + a01*a10*a22 + a00*a11*a22
        


    AUTHORS:
    - Edinah K. Gnang, Harry Crane
    """
    # Initializing the symbolic matrix
    sz = min(A.nrows(),A.ncols())
    P = GeneratePartition(sz)
    Pml = []
    for part in P:
        Pml.append((-1)^(max(part))*factorial(max(part))*((Partition2Matrix(part)).elementwise_product(A)).det())
    return expand(sum(Pml))*(-1)^A.nrows()

def PhiMatrix(n):
    """
    First of two matrices associated with a construction of Robin W. Whitty
    for listing gracefully labeled trees and graphs


    EXAMPLES:
    ::


        sage: PhiMatrix(3)
        [  0 x02 x01]
        [  0   0 x12]
        [  0   0   0]
        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the matrix
    M = Matrix(SR,n,n,zero_matrix(n,n))
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            if i < j:
                M[i,j] = var('x'+str(i)+str((n-1)-j+i+1))
    return M

def PhiHM(A):
    """
    First of two matrices associated with a construction of Robin W. Whitty
    for listing gracefully labeled trees and graphs


    EXAMPLES:
    ::


        sage: PhiHM(HM(3,3,'x'))
        [:, :]=
        [  0 x02 x01]
        [  0   0 x12]
        [  0   0   0]
        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the matrix
    M = HM(A.n(0),A.n(1),'zero')
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            if i < j:
                M[i,j] = A[i, (A.n(0)-1)-j+i+1]
    return M

def LambdaMatrix(n):
    """
    Second of two matrices associated with a construction of Robin W. Whitty
    for listing gracefully labeled trees and graphs


    EXAMPLES:
    ::


        sage: LambdaMatrix(3)
        [  0   0   0]
        [  0   0 x01]
        [  0 x02 x12]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the matrix
    M = Matrix(SR,n,n,zero_matrix(n,n))
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            if i > (n-1)-j:
                M[i,j] = var('x'+str(j-(n-1)+i-1)+str(i))
    return M

def LambdaHM(A):
    """
    Second of two matrices associated with a construction of Robin W. Whitty
    for listing gracefully labeled trees and graphs


    EXAMPLES:
    ::


        sage: LambdaHM(HM(3,3,'x'))
        [:, :]=
        [  0   0   0]
        [  0   0 x01]
        [  0 x02 x12]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the matrix
    M = HM(A.n(0),A.n(1),'zero')
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            if i > (A.n(0)-1)-j:
                M[i,j] = A[j-(A.n(0)-1)+i-1, i]
    return M

def ListingGracefulTrees(sz):
    """
    Listing gracefully labeled trees via the consruction of Robin W. Whitty.


    EXAMPLES:
    ::


        sage: ListingGracefulTrees(3)
        -x01*x02 + x02*x12


    AUTHORS:
    - Edinah K. Gnang
    """
    return expand(det((PhiMatrix(sz)-LambdaMatrix(sz))[1:sz,1:sz]))

def ListingGracefulGraphs(n):
    """
    Listing gracefully labeled graphs via the consruction of Robin W. Whitty.


    EXAMPLES:
    ::


        sage: ListingGracefulGraphs(3)
        x01*x02 + x02*x12


    AUTHORS:
    - Edinah K. Gnang, Harry Crane
    """
    return PermComp((PhiMatrix(n)+LambdaMatrix(n))[1:n,1:n])

def GrL(T):
    """
    Outputs the list of graceful relabling of the input functional
    directed graph T specified in tuple description.


    EXAMPLES:

    ::

        sage: sz=5; GrL([(i,0) for i in rg(sz)])
        [[(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)],  [(0, 4), (1, 4), (2, 4), (3, 4), (4, 4)]]
 

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the number of vertices
    sz = len(T)
    # Shifting the vertex indices to make the sage graph theory package happy
    tp = [(1+T[i][0], 1+T[i][1]) for i in rg(sz)]
    # Initializing the permutations
    P = Permutations(sz); S = SymmetricGroup(sz)
    # Initializing the graph
    grph = T2DiGraphII(tp,sz+1)
    # Initializing the automorphism group
    AutGrp = grph.automorphism_group()
    # Initializing representatives of Left coset as strings
    Lcst = [CstL[0] for CstL in S.cosets(AutGrp)]
    # Initializing the output list
    Lf = []
    # Looping through the coset representatives
    for p in Lcst:
        if Set([abs(p.dict()[tp[i][1]] - p.dict()[tp[i][0]]) for i in rg(sz)]) == Set(rg(sz)):
            Lf.append([(p.dict()[tp[i][0]]-1, p.dict()[tp[i][1]]-1) for i in rg(sz)])
    return Lf

def GrSp(T):
    """
    Outputs the list of triples sigma, sigma inverse, gamma, sign.
    directed graph T specified in tuple description.
    The function takes only one permutation per cosets
    of the automorphism.


    EXAMPLES:

    ::

        sage: sz=5; GrSp([(i,1) for i in rg(sz)])
        [[[(0, 0), (1, 4), (2, 1), (3, 2), (4, 3)],
          [(0, 0), (1, 2), (2, 3), (3, 4), (4, 1)],
          [(0, 4), (1, 0), (2, 3), (3, 2), (4, 1)],
          [(0, -1), (1, 0), (2, -1), (3, -1), (4, -1)]],
         [[(0, 1), (1, 0), (2, 2), (3, 3), (4, 4)],
          [(0, 1), (1, 0), (2, 2), (3, 3), (4, 4)],
          [(0, 1), (1, 0), (2, 2), (3, 3), (4, 4)],
          [(0, 1), (1, 0), (2, 1), (3, 1), (4, 1)]]]
 

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the number of vertices
    sz = len(T)
    # Shifting the vertex indices to make the sage graph theory package happy
    tp = [(1+T[i][0], 1+T[i][1]) for i in rg(sz)]
    # Initializing the permutations
    P = Permutations(sz); S = SymmetricGroup(sz)
    # Initializing the graph
    grph = T2DiGraphII(tp, sz+1)
    # Initializing the automorphism group
    AutGrp = grph.automorphism_group()
    # Initializing representatives of Left coset as strings
    Lcst = [CstL[0] for CstL in S.cosets(AutGrp)]
    # Initializing the output list
    Lf = []
    # Looping through the coset representatives
    for p in Lcst:
        if Set([abs(p.dict()[tp[i][1]] - p.dict()[tp[i][0]]) for i in rg(sz)]) == Set(rg(sz)):
            # Initialization of the inverse
            q = p.inverse()
            Lf.append([[(i, p.dict()[i+1]-1) for i in rg(sz)], [(i, q.dict()[i+1]-1) for i in rg(sz)], [(i, abs(p.dict()[tp[i][0]] - p.dict()[tp[i][1]])) for i in rg(sz)], [(i, sign(p.dict()[tp[i][0]] - p.dict()[tp[i][1]])) for i in rg(sz)]])
    return Lf

