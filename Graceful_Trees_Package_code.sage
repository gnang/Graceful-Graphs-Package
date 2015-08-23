@cached_function
def TstPerm(n):
    """
    Goes through all permutation of n elements and outputs
    the list of permutations from which a graceful graph
    can be deduced.

    EXAMPLES:
    ::
        sage: TstPerm(5)
        [[1, 2, 3, 4], [2, 1, 3, 4], [2, 3, 1, 4], [3, 2, 1, 4]]

    AUTHORS:
    - Edinah K. Gnang
    """
    P  = Permutations(range(1,n))
    L = []
    for p in P:
        b = 1
        for j in range(len(p)):
            if ((j+1)>p[j]) and (p[j]>(n-1)-(j+1)):
                b = 0
        if b == 1:
            L.append(p)
    return L

@cached_function
def CntFnct(n):
    """
    Provides a list of permutations followed by the corresponding count
    for the number of directed graceful graphs associated with the 
    graceful graph permutation construction.

    EXAMPLES:
    ::
        sage: CntFnct(5)
        [[[1, 2, 3, 4], 4], [[2, 1, 3, 4], 2], [[2, 3, 1, 4], 2], [[3, 2, 1, 4], 4]]

    AUTHORS:
    - Edinah K. Gnang
    """
    P  = Permutations(range(1,n))
    L = []
    for p in P:
        b = 1
        for j in range(len(p)):
            if ((j+1)>p[j]) and (p[j]>(n-1)-(j+1)):
                b = 0
        if b==1:
            c = 1
            for j in range(len(p)):
                if ((j+1) <= p[j]) and (p[j] <= (n-1)-(j+1)):
                    c = c*2
            L.append([p,c])
    return L

@cached_function
def GnrtFnct(n):
    """
    Provides a list of permutations followed by the list of second indices
    for directed graceful graphs associated with the graceful graph
    permutation construction.

    EXAMPLES:
    ::
        sage: GnrtFnct(5)
        [[[1, 2, 3, 4], [[0, 0, 0, 0], [2, 0, 0, 0], [0, 4, 0, 0], [2, 4, 0, 0]]],
        [[2, 1, 3, 4], [[1, 3, 0, 0], [3, 3, 0, 0]]],
        [[2, 3, 1, 4], [[1, 1, 4, 0], [3, 1, 4, 0]]],
        [[3, 2, 1, 4], [[2, 0, 4, 0], [4, 0, 4, 0], [2, 4, 4, 0], [4, 4, 4, 0]]]]

    AUTHORS:
    - Edinah K. Gnang
    """
    P  = Permutations(range(1,n))
    L = []
    for p in P:
        b = 1
        for j in range(len(p)):
            if ((j+1)>p[j]) and (p[j]>(n-1)-(j+1)):
                b = 0
        if b==1:
            c = []
            for j in range(len(p)):
                if ((j+1)<=p[j]) and (p[j]>(n-1)-(j+1)):
                    if len(c) == 0:
                        c.append([p[j]-(j+1)])
                    else:
                        for i in range(len(c)):
                            c[i] = c[i]+[p[j]-(j+1)]

                elif ((j+1) > p[j]) and (p[j] <= (n-1)-(j+1)):
                    if len(c)==0:
                        c.append([p[j]+(j+1)])
                    else:
                        for i in range(len(c)):
                            c[i] = c[i]+[p[j]+(j+1)]

                elif ((j+1) <= p[j]) and (p[j] <= (n-1)-(j+1)):
                    if len(c)==0:
                        c.append([p[j]-(j+1)])
                        c.append([p[j]+(j+1)])
                    else:
                        d = copy(c)
                        for i in range(len(c)):
                            c[i] = c[i]+[p[j]-(j+1)]
                        for i in range(len(d)):
                            d[i] = d[i]+[p[j]+(j+1)]
                        c = c+d
            L.append([p,c])
    return L

def get_permutation(la,lb):
    """
    Obtains a permutation list from two lists of the same size.
    No check is performed here the user must be very carefull
    to input lists of the same size
    
 
    EXAMPLES:
    ::
        sage: get_permutation([1, 2, 3, 4, 6, 8],[1, 2, 4, 8, 3, 6])
        [0, 1, 4, 2, 5, 3]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing the output
    L = list()

    # Loop performing the evaluation.
    for i1 in range(len(la)):
        for i2 in range(len(lb)):
            if la[i1] == lb[i2]:
                L.append(i2)
                break
   
    return L

@cached_function
def GnrtFnctII(n):
    """
    Provides a list of permutations followed by the list of second indices
    for directed graceful graphs associated with the graceful graph
    permutation construction.
    The difference between this and the implementation in GnrtFnct is the fact
    that we go through two list of permutation as opposed to a single list.

    EXAMPLES:
    ::
        sage: GnrtFnctII(5)
        [[[1, 2, 3, 4], [[0, 0, 0, 0], [2, 0, 0, 0], [0, 4, 0, 0], [2, 4, 0, 0]]],
        [[2, 1, 3, 4], [[2, 3, 0, 0], [3, 3, 0, 0]]],
        [[2, 3, 1, 4], [[3, 3, 4, 0], [2, 3, 4, 0]]],
        [[3, 2, 1, 4], [[2, 0, 4, 0], [4, 0, 4, 0], [2, 4, 4, 0], [4, 4, 4, 0]]]]

    AUTHORS:
    - Edinah K. Gnang
    """
    P  = Permutations(range(1,n))
    L = []
    for p in P:
        b = 1
        for j in range(len(p)):
            if ((j+1)>p[j]) and (p[j]>(n-1)-(j+1)):
                b = 0
        if b==1:
            # Initializing the list of functions
            c = []
            # Initializing the inverse permutations
            q = get_permutation(range(1,n),p)
            for j in range(len(p)):
                if ((j+1)<=p[j]) and (p[j]>(n-1)-(j+1)):
                    if len(c) == 0:
                        if (p[j]-(j+1)) == 0:
                            c.append([0])
                        else:
                            c.append([q[p[j]-(j+1)-1]+1])
                    else:
                        if (p[j]-(j+1)) == 0:
                            for i in range(len(c)):
                                c[i] = c[i]+[0]
                        else:
                            for i in range(len(c)):
                                c[i] = c[i]+[q[p[j]-(j+1)-1]+1]
                elif ((j+1) > p[j]) and (p[j] <= (n-1)-(j+1)):
                    if len(c)==0:
                        c.append([q[p[j]+(j+1)-1]+1])
                    else:
                        for i in range(len(c)):
                            c[i] = c[i]+[q[p[j]+(j+1)-1]+1]
                elif ((j+1) <= p[j]) and (p[j] <= n-(j+1)):
                    if len(c)==0:
                        if p[j]-(j+1) == 0:
                            c.append([0])
                        else:
                            c.append([q[p[j]-(j+1)-1]+1])
                        c.append([q[p[j]+(j+1)-1]+1])
                    else:
                        d = copy(c)
                        if p[j]-(j+1) == 0:
                            for i in range(len(c)):
                                c[i] = c[i]+[0]
                        else:
                            for i in range(len(c)):
                                c[i] = c[i]+[q[p[j]-(j+1)-1]+1]
                        for i in range(len(d)):
                            d[i] = d[i]+[q[p[j]+(j+1)-1]+1]
                        c = c+d
            L.append([p,c])
    return L

@cached_function
def GnrtFnctIII(n):
    """
    Provides a list of permutations followed by the list of second indices
    for directed graceful graphs associated with the graceful graph
    permutation construction.
    The difference between this and the implementation in GnrtFnctII is
    the fact that we are working with permutation of (n-1) elements.

    EXAMPLES:
    ::
        sage: GnrtFnctIII(5)
        [[[[1, 2, 3], [1, 2, 3]], [[0, 0, 0], [2, 0, 0]]],
        [[[1, 2, 3], [2, 1, 3]], [[3, 1, 0], [3, 3, 0]]],
        [[[1, 3, 2], [1, 3, 2]], [[0, 0, 0], [3, 0, 0]]],
        [[[1, 3, 2], [2, 3, 1]], [[2, 0, 1], [2, 0, 2]]],
        [[[2, 1, 3], [1, 2, 3]], [[2, 3, 0], [3, 3, 0]]],
        [[[2, 1, 3], [2, 1, 3]], [[0, 0, 0], [0, 1, 0]]],
        [[[2, 3, 1], [1, 3, 2]], [[3, 0, 2], [2, 0, 2]]],
        [[[2, 3, 1], [2, 3, 1]], [[0, 0, 0], [0, 0, 1]]],
        [[[3, 1, 2], [3, 1, 2]], [[0, 0, 0], [0, 3, 0]]],
        [[[3, 1, 2], [3, 2, 1]], [[0, 1, 2], [0, 1, 1]]],
        [[[3, 2, 1], [3, 1, 2]], [[0, 3, 1], [0, 1, 1]]],
        [[[3, 2, 1], [3, 2, 1]], [[0, 0, 0], [0, 0, 2]]]]

    AUTHORS:
    - Edinah K. Gnang
    """
    P  = Permutations(range(1,n))
    Q  = Permutations(range(1,n))
    L = []
    for p in P:
        for q in Q:
            b = 1
            for j in range(len(p)):
                if (q[j]>p[j]) and (p[j]>(n-1)-q[j]):
                    b = 0
            if b==1:
                # Initializing the list of functions
                c = []
    
                # Initializing the inverse permutations
                pinv = get_permutation(range(1,n),p)
    
                for j in range(len(p)):
                    if (q[j] <= p[j]) and (p[j] > (n-1)-q[j]):
                        if len(c) == 0:
                            if (p[j]-q[j]) == 0:
                                c.append([0])
                            else:
                                c.append([pinv[p[j]-q[j]-1]+1])
                        else:
                            if (p[j]-q[j]) == 0:
                                for i in range(len(c)):
                                    c[i] = c[i]+[0]
                            else:
                                for i in range(len(c)):
                                    c[i] = c[i]+[pinv[p[j]-q[j]-1]+1]
    
                    elif (q[j] > p[j]) and (p[j] <= (n-1)-q[j]):
                        if len(c)==0:
                            c.append([pinv[p[j]+q[j]-1]+1])
                        else:
                            for i in range(len(c)):
                                c[i] = c[i]+[pinv[p[j]+q[j]-1]+1]
    
                    elif (q[j] <= p[j]) and (p[j] <= n-q[j]):
                        if len(c)==0:
                            if p[j]-q[j] == 0:
                                c.append([0])
                            else:
                                c.append([pinv[p[j]-q[j]-1]+1])
                            c.append([pinv[p[j]+q[j]-1]+1])
                        else:
                            d = copy(c)
                            if p[j]-q[j] == 0:
                                for i in range(len(c)):
                                    c[i] = c[i]+[0]
                            else:
                                for i in range(len(c)):
                                    c[i] = c[i]+[pinv[p[j]-q[j]-1]+1]
                            for i in range(len(d)):
                                d[i] = d[i]+[pinv[p[j]+q[j]-1]+1]
                            c = c+d
                L.append([[p,q],c])
    return L

def List2Adj(L):
    """
    Starting from a list of second indices L and returns
    an adjacency matrices.

    EXAMPLES:
    ::
        sage: List2Adj([0,0,0,0])
        [0 1 1 1 1]
        [1 0 0 0 0]
        [1 0 0 0 0]
        [1 0 0 0 0]
        [1 0 0 0 0]

    AUTHORS:
    - Edinah K. Gnang
    """
    Id = identity_matrix(1+len(L))
    M  = zero_matrix(1+len(L),1+len(L)) 
    for i in range(len(L)): 
        M = M + (Id[i+1]).column()*(Id[L[i]]).row() +  (Id[L[i]]).column()*(Id[i+1]).row()
    return M

def is_Tree(A):
    """
    Returns an boolean value determining if the input unweighted adjacency matrix
    is associated with a tree. The implementation is based on a direct implementation
    of the matrix tree theorem.

    EXAMPLES:
    ::
        sage: is_Tree(List2Adj([0,0,0,0]))
        True

    AUTHORS:
    - Edinah K. Gnang
    """
    if ((diagonal_matrix((A*ones_matrix(A.nrows(),1)).list())-A)[0:A.nrows()-1,0:A.ncols()-1]).det()==1:
        return True
    else:
        return False
 
@cached_function
def CntPerm(n):
    """
    Provides a count for the number of permutation associated with the 
    graceful graph permutation construction.

    EXAMPLES:
    ::
        sage: CntPerm(5)
        96

    AUTHORS:
    - Edinah K. Gnang
    """
    P  = Permutations(range(1,n))
    cnt = 0
    for q in P:
        for p in P:
            b = 1
            for j in range(len(p)):
                if (q[j]>p[j]) and (p[j]>(n-1)-q[j]):
                    b = 0
            if b == 1:
                cnt = cnt+1
    return cnt

@cached_function
def CntFnctII(n):
    """
    Provides a count for the number of second indices list associated with the 
    graceful graph permutation construction.

    EXAMPLES:
    ::
        sage: CntFnctII(5)
        12

    AUTHORS:
    - Edinah K. Gnang
    """
    P  = Permutations(range(1,n))
    cnt = 0
    q = P[0]
    for p in P:
        b = 1
        for j in range(len(p)):
            if (q[j]>p[j]) and (p[j]>(n-1)-q[j]):
                b = 0
        if b==1:
            c = 1
            for j in range(len(p)):
                if (q[j] <= p[j]) and (p[j] <= (n-1)-q[j]):
                    c = c*2
            cnt = cnt+c
    return cnt

@cached_function
def CntFnctIII(n):
    """
    Provides a count for the number of second indices list associated with the 
    graceful graph permutation construction. The difference with CntFnctII
    is that the loop is here runs over two list of permutations.

    EXAMPLES:
    ::
        sage: CntFnctIII(5)
        288

    AUTHORS:
    - Edinah K. Gnang
    """
    P  = Permutations(range(1,n))
    cnt = 0
    for q in P:
        for p in P:
            b = 1
            for j in range(len(p)):
                if (q[j]>p[j]) and (p[j]>(n-1)-q[j]):
                    b = 0
            if b==1:
                c = 1
                for j in range(len(p)):
                    if (q[j] <= p[j]) and (p[j] <= (n-1)-q[j]):
                        c = c*2
                cnt = cnt+c
    return cnt

def GnrtAdj(n):
    """
    Returns a list of adjacency matrices each of which is followed by 
    a bollean value indicating if it is associated with a tree.

    EXAMPLES:
    ::
        sage: GnrtAdj(3)
        [[
        [0 1 1]      
        [1 0 0]      
        [1 0 0], True
        ],
         [
        [0 0 1]      
        [0 0 1]      
        [1 1 0], True
        ]]

    AUTHORS:
    - Edinah K. Gnang
    """
    Lg = GnrtFnctII(n) 
    return [[List2Adj(Lg[i][1][k]), is_Tree(List2Adj(Lg[i][1][k]))] for i in range(len(Lg)) for k in range(len(Lg[i]))]

def PlotGracefulGraphs(n):
    """
    Stores in the working directory png files associated with graceful graphs
    derived from the permutation construction. The names of the files yield
    permutation of the vertices which render the graph graceful.

    EXAMPLES:
    ::
        sage: PlotGracefulGraphs(3)
        

    AUTHORS:
    - Edinah K. Gnang
    """
    for p in Permutations(n-1):
        b = 1
        for j in range(len(p)):
            if ((j+1)>p[j]) and (p[j]>(n-1)-(j+1)):
                b = 0
        if b==1:
            # Initializing the list of functions
            c = []

            # Initializing the inverse permutations
            q = get_permutation(range(1,n),p)

            for j in range(len(p)):
                if ((j+1) <= p[j]) and (p[j] > (n-1)-(j+1)):
                    if len(c) == 0:
                        if (p[j]-(j+1)) == 0:
                            c.append([0])
                        else:
                            c.append([q[p[j]-(j+1)-1]+1])
                    else:
                        if (p[j]-(j+1)) == 0:
                            for i in range(len(c)):
                                c[i] = c[i]+[0]
                        else:
                            for i in range(len(c)):
                                c[i] = c[i]+[q[p[j]-(j+1)-1]+1]

                elif ((j+1) > p[j]) and (p[j] <= (n-1)-(j+1)):
                    if len(c)==0:
                        c.append([q[p[j]+(j+1)-1]+1])
                    else:
                        for i in range(len(c)):
                            c[i] = c[i]+[q[p[j]+(j+1)-1]+1]

                elif ((j+1) <= p[j]) and (p[j] <= n-(j+1)):
                    if len(c)==0:
                        if p[j]-(j+1) == 0:
                            c.append([0])
                        else:
                            c.append([q[p[j]-(j+1)-1]+1])
                        c.append([q[p[j]+(j+1)-1]+1])
                    else:
                        d = copy(c)
                        if p[j]-(j+1) == 0:
                            for i in range(len(c)):
                                c[i] = c[i]+[0]
                        else:
                            for i in range(len(c)):
                                c[i] = c[i]+[q[p[j]-(j+1)-1]+1]
                        for i in range(len(d)):
                            d[i] = d[i]+[q[p[j]+(j+1)-1]+1]
                        c = c+d
            for t in range(len(c)):
                Graph(List2Adj(c[t])).plot().save(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')

def PlotGracefulTrees(n):
    """
    Stores in the working directory png files associated with graceful trees
    derived from the permutation construction. The names of the files yield
    permutation of the vertices which render the graph graceful.

    EXAMPLES:
    ::
        sage: PlotGracefulTrees(3)
        

    AUTHORS:
    - Edinah K. Gnang
    """
    for p in Permutations(n-1):
        b = 1
        for j in range(len(p)):
            if ((j+1)>p[j]) and (p[j]>(n-1)-(j+1)):
                b = 0
        if b==1:
            # Initializing the list of functions
            c = []
            # Initializing the inverse permutations
            q = get_permutation(range(1,n),p)
            for j in range(len(p)):
                if ((j+1) <= p[j]) and (p[j] > (n-1)-(j+1)):
                    if len(c) == 0:
                        if (p[j]-(j+1)) == 0:
                            c.append([0])
                        else:
                            c.append([q[p[j]-(j+1)-1]+1])
                    else:
                        if (p[j]-(j+1)) == 0:
                            for i in range(len(c)):
                                c[i] = c[i]+[0]
                        else:
                            for i in range(len(c)):
                                c[i] = c[i]+[q[p[j]-(j+1)-1]+1]
                elif ((j+1) > p[j]) and (p[j] <= (n-1)-(j+1)):
                    if len(c)==0:
                        c.append([q[p[j]+(j+1)-1]+1])
                    else:
                        for i in range(len(c)):
                            c[i] = c[i]+[q[p[j]+(j+1)-1]+1]
                elif ((j+1) <= p[j]) and (p[j] <= n-(j+1)):
                    if len(c)==0:
                        if p[j]-(j+1) == 0:
                            c.append([0])
                        else:
                            c.append([q[p[j]-(j+1)-1]+1])
                        c.append([q[p[j]+(j+1)-1]+1])
                    else:
                        d = copy(c)
                        if p[j]-(j+1) == 0:
                            for i in range(len(c)):
                                c[i] = c[i]+[0]
                        else:
                            for i in range(len(c)):
                                c[i] = c[i]+[q[p[j]-(j+1)-1]+1]
                        for i in range(len(d)):
                            d[i] = d[i]+[q[p[j]+(j+1)-1]+1]
                        c = c+d
            for t in range(len(c)):
                # Testing Treenes
                if is_Tree(List2Adj(c[t])): 
                    Graph(List2Adj(c[t])).plot().save(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')


def PlotGracefulGraphsII(n):
    """
    Stores in the working directory png files associated with graceful graphs
    derived from the permutation construction. The names of the files yield
    permutation of the vertices which render the graph graceful. The difference
    with PlotGracefulGraphs is that we loop over two list of permutations as opposed
    to a single list of permutations. 

    EXAMPLES:
    ::
        sage: PlotGracefulGraphsII(3)
        

    AUTHORS:
    - Edinah K. Gnang
    """
    for p in Permutations(n-1):
        for q in Permutations(n-1):
            b = 1
            for j in range(len(p)):
                if (q[j]>p[j]) and (p[j]>(n-1)-q[j]):
                    b = 0
            if b==1:
                # Initializing the list of functions
                c = []
                # Initializing the inverse permutations
                pinv = get_permutation(range(1,n),p)
                for j in range(len(p)):
                    if (q[j] <= p[j]) and (p[j] > (n-1)-q[j]):
                        if len(c) == 0:
                            if (p[j]-q[j]) == 0:
                                c.append([0])
                            else:
                                c.append([pinv[p[j]-q[j]-1]+1])
                        else:
                            if (p[j]-q[j]) == 0:
                                for i in range(len(c)):
                                    c[i] = c[i]+[0]
                            else:
                                for i in range(len(c)):
                                    c[i] = c[i]+[pinv[p[j]-q[j]-1]+1]
                    elif (q[j] > p[j]) and (p[j] <= (n-1)-q[j]):
                        if len(c)==0:
                            c.append([pinv[p[j]+q[j]-1]+1])
                        else:
                            for i in range(len(c)):
                                c[i] = c[i]+[pinv[p[j]+q[j]-1]+1]
                    elif (q[j] <= p[j]) and (p[j] <= n-q[j]):
                        if len(c)==0:
                            if p[j]-q[j] == 0:
                                c.append([0])
                            else:
                                c.append([pinv[p[j]-q[j]-1]+1])
                            c.append([pinv[p[j]+q[j]-1]+1])
                        else:
                            d = copy(c)
                            if p[j]-q[j] == 0:
                                for i in range(len(c)):
                                    c[i] = c[i]+[0]
                            else:
                                for i in range(len(c)):
                                    c[i] = c[i]+[pinv[p[j]-q[j]-1]+1]
                            for i in range(len(d)):
                                d[i] = d[i]+[pinv[p[j]+q[j]-1]+1]
                            c = c+d
                for t in range(len(c)):
                    Graph(List2Adj(c[t])).plot().save(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')

def GracefulAdjacencyMatrices(n):
    """
    Returns a list of adjacency matrices associated with graceful graphs.
    Each of which are followed by the list of second indices and the permutation
    from which the list of second indices are deduced.

    EXAMPLES:
    ::
        sage: GracefulAdjacencyMatrices(3)
        [[
        [0 1 1]                
        [1 0 0]                
        [1 0 0], [1, 2], [0, 0]
        ],
        [
        [0 0 1]                
        [0 0 1]                
        [1 1 0], [1, 2], [2, 0]
        ],
        [
        [0 1 1]                
        [1 0 0]                
        [1 0 0], [2, 1], [0, 0]
        ],
        [
        [0 1 0]                
        [1 0 1]                
        [0 1 0], [2, 1], [0, 1]
        ]]

    AUTHORS:
    - Edinah K. Gnang
    """
    L = []
    for p in Permutations(n-1):
        for q in Permutations(n-1):
            b = 1
            for j in range(len(p)):
                if (q[j]>p[j]) and (p[j]>(n-1)-q[j]):
                    b = 0
            if b==1:
                # Initializing the list of functions
                c = []
                # Initializing the inverse permutations
                pinv = get_permutation(range(1,n),p)
                for j in range(len(p)):
                    if (q[j] <= p[j]) and (p[j] > (n-1)-q[j]):
                        if len(c) == 0:
                            if (p[j]-q[j]) == 0:
                                c.append([0])
                            else:
                                c.append([pinv[p[j]-q[j]-1]+1])
                        else:
                            if (p[j]-q[j]) == 0:
                                for i in range(len(c)):
                                    c[i] = c[i]+[0]
                            else:
                                for i in range(len(c)):
                                    c[i] = c[i]+[pinv[p[j]-q[j]-1]+1]
                    elif (q[j] > p[j]) and (p[j] <= (n-1)-q[j]):
                        if len(c)==0:
                            c.append([pinv[p[j]+q[j]-1]+1])
                        else:
                            for i in range(len(c)):
                                c[i] = c[i]+[pinv[p[j]+q[j]-1]+1]
                    elif (q[j] <= p[j]) and (p[j] <= n-q[j]):
                        if len(c)==0:
                            if p[j]-q[j] == 0:
                                c.append([0])
                            else:
                                c.append([pinv[p[j]-q[j]-1]+1])
                            c.append([pinv[p[j]+q[j]-1]+1])
                        else:
                            d = copy(c)
                            if p[j]-q[j] == 0:
                                for i in range(len(c)):
                                    c[i] = c[i]+[0]
                            else:
                                for i in range(len(c)):
                                    c[i] = c[i]+[pinv[p[j]-q[j]-1]+1]
                            for i in range(len(d)):
                                d[i] = d[i]+[pinv[p[j]+q[j]-1]+1]
                            c = c+d
                for t in range(len(c)):
                    L.append([List2Adj(c[t]), p, c[t]])
    return L 

def GracefulTreeAdjacencyMatrices(n):
    """
    Returns a list of adjacency matrices associated with graceful trees 

    EXAMPLES:
    ::
        sage: GracefulTreeAdjacencyMatrices(3)
        [[
        [0 1 1]                
        [1 0 0]                
        [1 0 0], [1, 2], [0, 0]
        ],
        [
        [0 0 1]                
        [0 0 1]                
        [1 1 0], [1, 2], [2, 0]
        ]]
    AUTHORS:
    - Edinah K. Gnang
    """
    L = []
    for p in Permutations(n-1):
        b = 1
        for j in range(len(p)):
            if ((j+1)>p[j]) and (p[j]>(n-1)-(j+1)):
                b = 0
        if b==1:
            # Initializing the list of functions
            c = []
            # Initializing the inverse permutations
            q = get_permutation(range(1,n),p)
            for j in range(len(p)):
                if ((j+1) <= p[j]) and (p[j] > (n-1)-(j+1)):
                    if len(c) == 0:
                        if (p[j]-(j+1)) == 0:
                            c.append([0])
                        else:
                            c.append([q[p[j]-(j+1)-1]+1])
                    else:
                        if (p[j]-(j+1)) == 0:
                            for i in range(len(c)):
                                c[i] = c[i]+[0]
                        else:
                            for i in range(len(c)):
                                c[i] = c[i]+[q[p[j]-(j+1)-1]+1]
                elif ((j+1) > p[j]) and (p[j] <= (n-1)-(j+1)):
                    if len(c)==0:
                        c.append([q[p[j]+(j+1)-1]+1])
                    else:
                        for i in range(len(c)):
                            c[i] = c[i]+[q[p[j]+(j+1)-1]+1]
                elif ((j+1) <= p[j]) and (p[j] <= n-(j+1)):
                    if len(c)==0:
                        if p[j]-(j+1) == 0:
                            c.append([0])
                        else:
                            c.append([q[p[j]-(j+1)-1]+1])
                        c.append([q[p[j]+(j+1)-1]+1])
                    else:
                        d = copy(c)
                        if p[j]-(j+1) == 0:
                            for i in range(len(c)):
                                c[i] = c[i]+[0]
                        else:
                            for i in range(len(c)):
                                c[i] = c[i]+[q[p[j]-(j+1)-1]+1]
                        for i in range(len(d)):
                            d[i] = d[i]+[q[p[j]+(j+1)-1]+1]
                        c = c+d
            for t in range(len(c)):
                # Testing Treenes
                if is_Tree(List2Adj(c[t])): 
                    L.append([List2Adj(c[t]), p, c[t]])
    return L 
