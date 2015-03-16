@cached_function
def TstPerm(n):
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
                if ((j+1) <= p[j]) and (p[j]>(n-1)-(j+1)):
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
            L.append([p,c])
    return L

@cached_function
def GnrtFnctIII(n):
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
    Id = identity_matrix(1+len(L))
    M  = zero_matrix(1+len(L),1+len(L)) 
    for i in range(len(L)): 
        M = M + (Id[i+1]).column()*(Id[L[i]]).row() +  (Id[L[i]]).column()*(Id[i+1]).row()
    return M

def is_Tree(A):
    return ((diagonal_matrix((A*ones_matrix(A.nrows(),1)).list())-A)[0:A.nrows()-1,0:A.ncols()-1]).det()==1

@cached_function
def CntPerm(n):
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
    Lg = GnrtFnctII(n) 
    return [[List2Adj(Lg[i][1][k]), is_Tree(LisLg[i][1][k])] for i in range(len(Lg)) for k in range(len(Lg[i]))]

def PlotGracefulGraphs(n):
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


def GracefulGraphsAdjacencyMatrices(n):
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
                    #Graph(List2Adj(c[t])).plot().save(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
                    L.append(List2Adj(c[t]))
    return L 


