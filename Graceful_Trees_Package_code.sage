#*************************************************************************#
#       Copyright (C) 2015 Edinah K. Gnang <kgnang@gmail.com>,            #
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
    graceful graphs with all the vertices have degree>0.

    EXAMPLES:
    ::
        sage: GracefulPermutations(4)
        [[1, 2, 3], [2, 1, 3]]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Permutations of elements from 1 to (n-1)
    P=Permutations(n-1)
    # Initialization of the sets
    S1=Set(range(n)); S2=Set(range(2,n))
    # Intialization of the list collecting the graceful permutations
    L=[]
    # Loop collecting the GracefulPermutations
    for q in P:
        # Appending 0 at end of the permutation
        #p=[q[i] for i in range(n-1)]+[0]
        p=[q[i] for i in range(n-1)]
        # Initialization of our boolean variable
        bl=True
        for j in range(n-1):
            # Verifying that the graceful criteria is met
            if ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                bl=False
                break
        if bl:
            L.append(p)
    return L

@cached_function
def CountGracefulFunctions(n):
    """
    Goes through all the permutation of n>1 elements and outputs
    the list of of functions derived from graceful permutations.

    EXAMPLES:
    ::
        sage: CountGracefulFunctions(4)
        [[2, [1, 2, 3]], [2, [2, 1, 3]]]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Graceful Permutations.
    P=GracefulPermutations(n)
    # Initialization of the sets.
    S1=Set(range(n)); S2=Set(range(2,n))
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Loop collecting the GracefulPermutations.
    for p in P:
        c=1
        for j in range(n-1):
            if ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) in S2 ):
                c=2*c
        L.append([c,p])
    return L

@cached_function
def GracefulFunctions(n):
    """
    Goes through all the permutation of n>1 elements and outputs
    the list of of functions on the vertices derived from 
    graceful permutations.

    EXAMPLES:
    ::
        sage: GracefulFunctions(3)
        [[[[0, 0], [2, 0]], [1, 2]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Graceful Permutations.
    P=GracefulPermutations(n)
    # Initialization of the sets.
    S1=Set(range(n)); S2=Set(range(2,n))
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Loop going through the GracefulPermutations.
    for p in P:
        # Initialization of the list of functions.
        c=[]
        for j in range(n-1):
            # testing that only the first criteria is met
            if ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                if len(c)==0:
                    c.append([(j+1)-p[j]])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[abs((j+1)-p[j])]
            # testing that only the second criteria is met
            elif ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) in S2 ):
                if len(c)==0:
                    c.append([(j+1)+p[j]])
                else:
                    for i in range(len(c)):
                        c[i]=c[i]+[(j+1)+p[j]]
            # testing that all two criterias are met
            elif ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) in S2 ):
                if len(c)==0:
                    c.append([(j+1)-p[j]])
                    c.append([(j+1)+p[j]])
                else:
                    d=copy(c)
                    for i in range(len(c)):
                        c[i]=c[i]+[(j+1)-p[j]]
                    for i in range(len(d)):
                        d[i]=d[i]+[(j+1)+p[j]]
                    c=c+d
        L.append([c,p])
    return L

def CountGracefulTree(n):
    """
    Goes through all the permutation of n>1 elements 
    and enumerates the graceful trees on n vertices.

    EXAMPLES:
    ::
        sage: CountGracefulTree(3)
        2
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the list of Graceful Permutations.
    P=Permutations(n-1)
    # Initialization of the sets.
    S1=Set(range(n)); S2=Set(range(2,n))
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Initialization of the counter
    cnt=0
    # Loop going through the GracefulPermutations.
    for q in P:
        # Appending 0 at end of the permutation
        p=[q[i] for i in range(n-1)]+[0]
        # Initialization of our boolean variable
        bl=True
        for j in range(n-1):
            # Verifying that the graceful criteria is met
            if ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                bl=False
                break
        if bl:
            L.append(p)
            # Initialization of the list of functions.
            c=[]
            for j in range(n-1):
                # testing that only the first criteria is met
                if ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[abs((j+1)-p[j])]
                # testing that only the second criteria is met
                elif ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)+p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)+p[j]]
                # testing that all two criterias are met
                elif ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                        c.append([(j+1)+p[j]])
                    else:
                        d=copy(c)
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)-p[j]]
                        for i in range(len(d)):
                            d[i]=d[i]+[(j+1)+p[j]]
                        c=c+d
            # Put here the stuff that is needed.
            for t in range(len(c)):
                # Initialization of the adjacency matrix.
                A=sum([Id[:,k]*Id[c[t][k-1],:]+Id[:,c[t][k-1]]*Id[k,:]  for k in range(1,n)])
                if is_Tree(A):
                    cnt=cnt+1
    return cnt

def GracefulGraphAdjacencyMatrixList(n):
    """
    Goes through all the permutation of graceful second index functions
    and outputs the list of undirected adjacency matrices of graceful 
    graphs.

    EXAMPLES:
    ::
        sage: GracefulGraphAdjacencyMatrixList(3)
        [
        [0 1 1]  [0 0 1]
        [1 0 0]  [0 0 1]
        [1 0 0], [1 1 0]
        ]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GracefulFunctions(n)
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        for f in l[0]:
            Rslt.append(sum([Id[:,j]*Id[f[j-1],:]+Id[:,f[j-1]]*Id[j,:] for j in range(1,n)]))
    return Rslt

def is_Tree(A):
    """
    Returns an boolean value determining if the input unweighted adjacency matrix
    is associated with a tree. The implementation is based on a direct implementation
    of the matrix tree theorem.

    EXAMPLES:
    ::
        sage: is_Tree(GracefulTreeAdjacencyMatrixList(3)[0])
        True

    AUTHORS:
    - Edinah K. Gnang
    """
    if ((diagonal_matrix((A*ones_matrix(A.nrows(),1)).list())-A)[0:A.nrows()-1,0:A.ncols()-1]).det()==1:
        return True
    else:
        return False

def GracefulTreeAdjacencyMatrixList(n):
    """
    Goes through all the permutation of graceful second index functions
    and outputs the list adjacency matrices undirected graceful tree.

    EXAMPLES:
    ::
        sage: GracefulTreeAdjacencyMatrixList(3)
        [
        [0 1 1]  [0 0 1]
        [1 0 0]  [0 0 1]
        [1 0 0], [1 1 0]
        ]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GracefulFunctions(n)
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        for f in l[0]:
            # Testing the treeness
            if is_Tree(sum([Id[:,j]*Id[f[j-1],:]+Id[:,f[j-1]]*Id[j,:]  for j in range(1,n)])):
                Rslt.append(sum([Id[:,j]*Id[f[j-1],:]+Id[:,f[j-1]]*Id[j,:] for j in range(1,n)]))
    return Rslt

def GracefulTreeCharPolyAdjacencyMatrixList(n, x0=0, x1=1):
    """
    Goes through all the permutation of graceful second index functions
    and outputs the list adjacency matrices undirected graceful tree.

    EXAMPLES:
    ::
        sage: GracefulTreeCharPolyAdjacencyMatrixList(3)
        [[
        [0 1 1]                      
        [1 0 0]                      
        [1 0 0], x^3 - 2*x, [1, 2]
        ]]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GracefulFunctions(n)
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        for f in l[0]:
            # Initialization of the adjacency matrix
            TmpA = sum([Id[:,j]*Id[f[j-1],:]+Id[:,f[j-1]]*Id[j,:]  for j in range(1,n)])
            A = TmpA*x1+(ones_matrix(n,n)-TmpA)*x0
            # Testing the treeness
            if is_Tree(TmpA):
                rpt=False
                for M in Rslt:
                    if M[1]==A.charpoly():
                        rpt=True
                        break
                if rpt==False:
                    Rslt.append([sum([Id[:,j]*Id[f[j-1],:]+Id[:,f[j-1]]*Id[j,:]  for j in range(1,n)]), A.charpoly(), l[1]])
    return Rslt

def GracefulTreeDrawingsI(n):
    """
    Goes through all the permutation of n>1 elements constructs
    the functions for constructing graceful graphs and saves 
    in the working directory the png files associated with
    drawings of gracefully labeled trees. This implemntation 
    never stores a list it just saves the images on the go

    EXAMPLES:
    ::
        sage: GracefulTreeDrawingsI(3)
        ['1_2_0__0.png', '1_2_0__1.png']
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of names
    List_of_Names=[]
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the list of Graceful Permutations.
    P=Permutations(n-1)
    # Initialization of the sets.
    S1=Set(range(n)); S2=Set(range(2,n))
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Loop going through the GracefulPermutations.
    for q in P:
        # Appending 0 at end of the permutation
        p=[q[i] for i in range(n-1)]+[0]
        # Initialization of our boolean variable
        bl=True
        for j in range(n-1):
            # Verifying that the graceful criteria is met
            if ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                bl=False
                break
        if bl:
            L.append(p)
            # Initialization of the list of functions.
            c=[]
            for j in range(n-1):
                # testing that only the first criteria is met
                if ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[abs((j+1)-p[j])]
                # testing that only the second criteria is met
                elif ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)+p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)+p[j]]
                # testing that all two criterias are met
                elif ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                        c.append([(j+1)+p[j]])
                    else:
                        d=copy(c)
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)-p[j]]
                        for i in range(len(d)):
                            d[i]=d[i]+[(j+1)+p[j]]
                        c=c+d
            # Put here the stuff that is needed.
            for t in range(len(c)):
                # Initialization of the adjacency matrix.
                A=sum([Id[:,k]*Id[c[t][k-1],:]+Id[:,c[t][k-1]]*Id[k,:]  for k in range(1,n)])
                if is_Tree(A):
                    Graph(A).plot().save(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
                    List_of_Names.append(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
    return List_of_Names

def GracefulTreeCharPolyRepresentativeDrawings(n, x0=0, x1=1):
    """
    Goes through all the permutation of n>1 elements constructs
    the functions for constructing graceful graphs and saves 
    in the working directory the png files associated with
    drawings of gracefully labeled trees. This implemntation 
    never stores a list it just saves the images on the go

    EXAMPLES:
    ::
        sage: L=GracefulTreeCharPolyRepresentativeDrawings(3)
        sage: L
        [[x^3 - 2*x, [1, 2, 0]]]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the list of Graceful Permutations.
    P=Permutations(n-1)
    # Initialization of the sets.
    S1=Set(range(n)); S2=Set(range(2,n))
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Initializaing the list storing the matrices
    Lm=[]
    # Loop going through the GracefulPermutations.
    for q in P:
        # Appending 0 at end of the permutation
        p=[q[i] for i in range(n-1)]+[0]
        # Initialization of our boolean variable
        bl=True
        for j in range(n-1):
            # Verifying that the graceful criteria is met
            if ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                bl=False
                break
        if bl:
            L.append(p)
            # Initialization of the list of functions.
            c=[]
            for j in range(n-1):
                # testing that only the first criteria is met
                if ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[abs((j+1)-p[j])]
                # testing that only the second criteria is met
                elif ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)+p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)+p[j]]
                # testing that all two criterias are met
                elif ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                        c.append([(j+1)+p[j]])
                    else:
                        d=copy(c)
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)-p[j]]
                        for i in range(len(d)):
                            d[i]=d[i]+[(j+1)+p[j]]
                        c=c+d
            # Put here the stuff that is needed.
            for t in range(len(c)):
                # Initialization of the adjacency matrix.
                TmpA=sum([Id[:,k]*Id[c[t][k-1],:]+Id[:,c[t][k-1]]*Id[k,:]  for k in range(1,n)])
                A=TmpA*x1+(ones_matrix(n,n)-TmpA)*x0
                if is_Tree(TmpA):
                    rpt=False
                    for M in Lm:
                        if M[0]==A.charpoly():
                            rpt=True
                            break
                    if rpt==False:
                        Lm.append([A.charpoly(), p])
                        Graph(TmpA).plot().save(str(len(Lm))+'__'+str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
    return Lm

def GracefulNonTreeAdjacencyMatrix(n):
    """
    Goes through all the permutation of graceful second index functions
    and outputs the list adjacency matrices undirected graceful tree.

    EXAMPLES:
    ::
        sage: GracefulNonTreeAdjacencyMatrix(7)
        [
        [0 0 0 0 0 1 1]  [0 0 0 0 0 1 1]  [0 0 0 0 1 0 1]  [0 0 0 0 1 0 1]
        [0 0 0 1 1 0 0]  [0 0 0 1 1 0 0]  [0 0 0 0 0 0 1]  [0 0 0 0 0 0 1]
        [0 0 0 0 0 0 1]  [0 0 0 0 0 0 1]  [0 0 0 1 0 1 0]  [0 0 0 1 0 1 0]
        [0 1 0 0 1 0 0]  [0 1 0 0 1 0 0]  [0 0 1 0 0 1 0]  [0 0 1 0 0 1 0]
        [0 1 0 1 0 0 0]  [0 1 0 1 0 0 0]  [1 0 0 0 0 0 0]  [1 0 0 0 0 0 0]
        [1 0 0 0 0 0 0]  [1 0 0 0 0 0 0]  [0 0 1 1 0 0 0]  [0 0 1 1 0 0 0]
        [1 0 1 0 0 0 0], [1 0 1 0 0 0 0], [1 1 0 0 0 0 0], [1 1 0 0 0 0 0]
        ]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GracefulFunctions(n)
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        for f in l[0]:
            # Testing the treeness
            if not is_Tree(sum([Id[:,j]*Id[f[j-1],:]+Id[:,f[j-1]]*Id[j,:]  for j in range(1,n)])):
                Rslt.append(sum([Id[:,j]*Id[f[j-1],:]+Id[:,f[j-1]]*Id[j,:] for j in range(1,n)]))
    return Rslt

def GracefulUndirectedBinaryTreeAdjacencyMatrix(n):
    """
    Goes through all the list of functions 
    derived from graceful permutations and
    construct directed adjacency matrices 
    of binary trees.

    EXAMPLES:
    ::
        sage: GracefulUndirectedBinaryTreeAdjacencyMatrix(3)
        [[[
        [0 1 1]  [0 0 1]
        [1 0 0]  [0 0 1]
        [1 0 0], [1 1 0]
        ], [1, 2]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GracefulFunctions(n)
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        Rslt.append([[]])
        for f in l[0]:
            # Initializing the adjacency matrix
            A = sum([Id[:,j]*Id[f[j-1],:]+Id[:,f[j-1]]*Id[j,:]  for j in range(1,n)])
            # Testing binary treeness
            if is_Tree(A) and (Set((A*ones_matrix(A.nrows(),1)).list())==Set([1,2]) or Set((A*ones_matrix(A.nrows(),1)).list())==Set([1,3]) or Set((A*ones_matrix(A.nrows(),1)).list())==Set([1,2,3])):
                Rslt[len(Rslt)-1][0].append(A)
        Rslt[len(Rslt)-1].append(l[len(l)-1])
    return Rslt


def GracefulBinaryTreeGenerator(n):
    """
    Goes through all the list of functions 
    derived from graceful permutations and
    construct directed adjacency matrices 
    of binary trees.

    EXAMPLES:
    ::
        sage: GracefulBinaryTreeGenerator(3)
        [[[], [1, 2]]]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GracefulFunctions(n)
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        Rslt.append([[]])
        for f in l[0]:
            # Initializing the adjacency matrix
            A = sum([Id[:,j]*Id[f[j-1],:]+Id[:,f[j-1]]*Id[j,:]  for j in range(1,n)])
            # Testing binary treeness
            if is_Tree(A) and (Set((A*ones_matrix(A.nrows(),1)).list())==Set([1,2]) or Set((A*ones_matrix(A.nrows(),1)).list())==Set([1,3]) or Set((A*ones_matrix(A.nrows(),1)).list())==Set([1,2,3])) and (A*ones_matrix(A.ncols(),1))[0,0]!=1 and (A*ones_matrix(A.ncols(),1))[A.nrows()-1,0]!=1:
                Rslt[len(Rslt)-1][0].append(A)
        Rslt[len(Rslt)-1].append(l[len(l)-1])
    return Rslt


def HalfPermutationUndirectedAdjacencyMatrixDrawings(n):
    """
    Returns adjacency matrices associated with the
    half permutation construction. The input must be 
    an odd positive integer.

    EXAMPLES:
    ::
        sage: HalfPermutationUndirectedAdjacencyMatrixDrawings(3)
        [[
        [0 0 1]        
        [0 0 1]        
        [1 1 0], [1, 2, 3]
        ]]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Permutations of elements from 1 to (n-1)
    P=Permutations(floor(n/2))
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initializing the list of adjacency matrices.
    L=[]
    # Loop constructing the adjacency matrices.
    for q in P:
        p=[q[i] for i in range(floor(n/2))]
        L.append([sum([Id[:,j]*Id[j+p[j-1],:]+Id[:,j+p[j-1]]*Id[j,:] for j in range(1,ceil(n/2))])+sum([Id[:,j]*Id[0,:]+Id[:,0]*Id[j,:] for j in range(ceil(n/2),n)]), [p[j-1] for j in range(1,ceil(n/2))]+range(ceil(n/2),n+1)])
        Graph(sum([Id[:,j]*Id[j+p[j-1],:]+Id[:,j+p[j-1]]*Id[j,:] for j in range(1,ceil(n/2))])+sum([Id[:,j]*Id[0,:]+Id[:,0]*Id[j,:] for j in range(ceil(n/2),n)])).plot().save(str([p[j-1] for j in range(1,ceil(n/2))]).replace(', ','_').replace('[','').replace(']','')+'.png')
    return L

def HalfPermutationGeneratorUndirectedAdjacencyMatrixDrawings(n):
    """
    Returns adjacency matrices associated with the
    half permutation construction. The input must be 
    an odd positive integer.

    EXAMPLES:
    ::
        sage: HalfPermutationGeneratorUndirectedAdjacencyMatrixDrawings(3)
        []

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Permutations of elements from 1 to (n-1)
    P=Permutations(floor(n/2))
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initializing the list of adjacency matrices.
    L=[]
    # Loop constructing the adjacency matrices.
    for q in P:
        p=[q[i] for i in range(floor(n/2))]
        # Initialization of the adjacency matrix
        TmpA=sum([Id[:,j]*Id[j+p[j-1],:]+Id[:,j+p[j-1]]*Id[j,:] for j in range(1,ceil(n/2))])+sum([Id[:,j]*Id[0,:]+Id[:,0]*Id[j,:] for j in range(ceil(n/2),n)])
        if (TmpA*ones_matrix(TmpA.ncols(),1))[0,0]!=1 and (TmpA*ones_matrix(TmpA.ncols(),1))[TmpA.nrows()-1,0]!=1:
            L.append([TmpA, [p[j-1] for j in range(1,ceil(n/2))]+range(ceil(n/2),n+1)])
            Graph(sum([Id[:,j]*Id[j+p[j-1],:]+Id[:,j+p[j-1]]*Id[j,:] for j in range(1,ceil(n/2))])+sum([Id[:,j]*Id[0,:]+Id[:,0]*Id[j,:] for j in range(ceil(n/2),n)])).plot().save(str([p[j-1] for j in range(1,ceil(n/2))]).replace(', ','_').replace('[','').replace(']','')+'.png')
    return L

def HalfPermutationUndirectedAdjacencyMatrixDrawingsII(n):
    """
    Returns adjacency matrices associated with the
    half permutation construction. The input must be 
    an odd positive integer.

    EXAMPLES:
    ::
        sage: HalfPermutationUndirectedAdjacencyMatrixDrawingsII(3)
        []
        [[
        [0 0 1]        
        [0 0 1]        
        [1 1 0], [1, 2]
        ]]


    AUTHORS:
    - Edinah K. Gnang
    """
    if n%2 == 1:
        # Initialization of the list of Permutations of elements from 1 to ((n-1)-floor(n/2)-1)
        P=Permutations((n-1)-floor(n/2)-1)
        # Initialization of the identity matrix
        Id=identity_matrix(n)
        # Initializing the list of adjacency matrices.
        L=[]
        # Loop constructing the adjacency matrices.
        for q in P:
            p=[q[i] for i in range((n-1)-floor(n/2)-1)]
            print p
            L.append([sum([Id[:,j]*Id[(n-1),:]+Id[:,(n-1)]*Id[j,:] for j in range(ceil(n/2))])+sum([Id[:,j]*Id[j-p[j-ceil(n/2)],:]+Id[:,j-p[j-ceil(n/2)]]*Id[j,:] for j in range(ceil(n/2),n-1)]),[(n-1)-j for j in range(1,ceil(n/2))]+[p[j] for j in range((n-1)-floor(n/2)-1)]+[n-1]])
            Graph(sum([Id[:,j]*Id[(n-1),:]+Id[:,(n-1)]*Id[j,:] for j in range(ceil(n/2))])+sum([Id[:,j]*Id[j-p[j-ceil(n/2)],:]+Id[:,j-p[j-ceil(n/2)]]*Id[j,:] for j in range(ceil(n/2),n-1)])).plot().save(str([(n-1)-j for j in range(1,ceil(n/2))]+[p[j] for j in range((n-1)-floor(n/2)-1)]+[n-1]).replace(', ','_').replace('[','').replace(']','')+'.png')
        return L
    else:
        # Initialization of the list of Permutations of elements from 1 to ((n-1)-floor(n/2)-1)
        P=Permutations(Integer(n/2)-1)
        # Initialization of the identity matrix
        Id=identity_matrix(n)
        # Initializing the list of adjacency matrices.
        L=[]
        # Loop constructing the adjacency matrices.
        for q in P:
            p=[q[i] for i in range(Integer(n/2)-1)]
            print p
            L.append([sum([Id[:,j]*Id[(n-1),:]+Id[:,(n-1)]*Id[j,:] for j in range(Integer(n/2))])+sum([Id[:,j]*Id[j-p[j-(n/2)],:]+Id[:,j-p[j-(n/2)]]*Id[j,:] for j in range(Integer(n/2),n-1)]),[(n-1)-j for j in range(1,Integer(n/2))]+[p[j] for j in range(Integer(n/2)-1)]+[n-1]])
            Graph(sum([Id[:,j]*Id[(n-1),:]+Id[:,(n-1)]*Id[j,:] for j in range(Integer(n/2))])+sum([Id[:,j]*Id[j-p[j-(n/2)],:]+Id[:,j-p[j-(n/2)]]*Id[j,:] for j in range(Integer(n/2),n-1)])).plot().save(str([(n-1)-j for j in range(1,Integer(n/2))]+[p[j] for j in range(Integer(n/2)-1)]+[n-1]).replace(', ','_').replace('[','').replace(']','')+'.png')
        return L
   
def HalfPermutationGeneratorUndirectedAdjacencyMatrixDrawingsII(n):
    """
    Returns adjacency matrices associated with the
    half permutation construction. The input must be 
    an odd positive integer.

    EXAMPLES:
    ::
        sage: HalfPermutationGeneratorUndirectedAdjacencyMatrixDrawingsII(3)
        []


    AUTHORS:
    - Edinah K. Gnang
    """
    if n%2 == 1:
        # Initialization of the list of Permutations of elements from 1 to ((n-1)-floor(n/2)-1)
        P=Permutations((n-1)-floor(n/2)-1)
        # Initialization of the identity matrix
        Id=identity_matrix(n)
        # Initializing the list of adjacency matrices.
        L=[]
        # Loop constructing the adjacency matrices.
        for q in P:
            p=[q[i] for i in range((n-1)-floor(n/2)-1)]
            TmpA=sum([Id[:,j]*Id[(n-1),:]+Id[:,(n-1)]*Id[j,:] for j in range(ceil(n/2))])+sum([Id[:,j]*Id[j-p[j-ceil(n/2)],:]+Id[:,j-p[j-ceil(n/2)]]*Id[j,:] for j in range(ceil(n/2),n-1)])
            if (TmpA*ones_matrix(TmpA.ncols(),1))[0,0]!=1 and (TmpA*ones_matrix(TmpA.ncols(),1))[TmpA.nrows()-1,0]!=1:
                L.append([TmpA, [(n-1)-j for j in range(1,ceil(n/2))]+[p[j] for j in range((n-1)-floor(n/2)-1)]+[n-1]])
                Graph(sum([Id[:,j]*Id[(n-1),:]+Id[:,(n-1)]*Id[j,:] for j in range(ceil(n/2))])+sum([Id[:,j]*Id[j-p[j-ceil(n/2)],:]+Id[:,j-p[j-ceil(n/2)]]*Id[j,:] for j in range(ceil(n/2),n-1)])).plot().save(str([(n-1)-j for j in range(1,ceil(n/2))]+[p[j] for j in range((n-1)-floor(n/2)-1)]+[n-1]).replace(', ','_').replace('[','').replace(']','')+'.png')
        return L
    else:
        # Initialization of the list of Permutations of elements from 1 to ((n-1)-floor(n/2)-1)
        P=Permutations(Integer(n/2)-1)
        # Initialization of the identity matrix
        Id=identity_matrix(n)
        # Initializing the list of adjacency matrices.
        L=[]
        # Loop constructing the adjacency matrices.
        for q in P:
            p=[q[i] for i in range(Integer(n/2)-1)]
            TmpA=sum([Id[:,j]*Id[(n-1),:]+Id[:,(n-1)]*Id[j,:] for j in range(Integer(n/2))])+sum([Id[:,j]*Id[j-p[j-(n/2)],:]+Id[:,j-p[j-(n/2)]]*Id[j,:] for j in range(Integer(n/2),n-1)])
            if (TmpA*ones_matrix(TmpA.ncols(),1))[0,0]!=1 and (TmpA*ones_matrix(TmpA.ncols(),1))[TmpA.nrows()-1,0]!=1:
                L.append([TmpA,[(n-1)-j for j in range(1,Integer(n/2))]+[p[j] for j in range(Integer(n/2)-1)]+[n-1]])
                Graph(sum([Id[:,j]*Id[(n-1),:]+Id[:,(n-1)]*Id[j,:] for j in range(Integer(n/2))])+sum([Id[:,j]*Id[j-p[j-(n/2)],:]+Id[:,j-p[j-(n/2)]]*Id[j,:] for j in range(Integer(n/2),n-1)])).plot().save(str([(n-1)-j for j in range(1,Integer(n/2))]+[p[j] for j in range(Integer(n/2)-1)]+[n-1]).replace(', ','_').replace('[','').replace(']','')+'.png')
        return L

def HalfPermutationUndirectedCharPolyAdjacencyMatrix(n, x0=0, x1=1):
    """
    Returns adjacency matrices associated with the
    half permutation construction. The input must be 
    an odd positive integer.

    EXAMPLES:
    ::
        sage: HalfPermutationUndirectedCharPolyAdjacencyMatrix(3)
        [[
        [0 0 1]                      
        [0 0 1]                      
        [1 1 0], x^3 - 2*x, [1, 2, 3]
        ]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Permutations of elements from 1 to (n-1)
    P=Permutations(floor(n/2))
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initializing the list of adjacency matrices.
    L=[]
    # Loop constructing the adjacency matrices.
    for q in P:
        p=[q[i] for i in range(floor(n/2))]
        # Initialization of the adjacency matrix
        TmpA=sum([Id[:,j]*Id[j+p[j-1],:]+Id[:,j+p[j-1]]*Id[j,:] for j in range(1,ceil(n/2))])+sum([Id[:,j]*Id[0,:]+Id[:,0]*Id[j,:] for j in range(ceil(n/2),n)])
        A=TmpA*x1+(ones_matrix(n,n)-TmpA)*x0
        rpt=False
        for M in L:
            if M[1]==A.charpoly():
                rpt=True
                break
        if rpt==False:
            L.append([TmpA, A.charpoly(), [p[j-1] for j in range(1,ceil(n/2))]+range(ceil(n/2),n+1)])
    return L

def HalfPermutationUndirectedCharPolyAdjacencyMatrixII(n, x0=0, x1=1):
    """
    Returns adjacency matrices associated with the
    half permutation construction. The input must be 
    an odd positive integer.

    EXAMPLES:
    ::
        sage: HalfPermutationUndirectedCharPolyAdjacencyMatrixII(3)
        [[
        [0 0 1]                   
        [0 0 1]                   
        [1 1 0], x^3 - 2*x, [1, 2]
        ]]


    AUTHORS:
    - Edinah K. Gnang
    """
    if n%2 == 1:
        # Initialization of the list of Permutations of elements from 1 to ((n-1)-floor(n/2)-1)
        P=Permutations((n-1)-floor(n/2)-1)
        # Initialization of the identity matrix
        Id=identity_matrix(n)
        # Initializing the list of adjacency matrices.
        L=[]
        # Loop constructing the adjacency matrices.
        for q in P:
            p=[q[i] for i in range((n-1)-floor(n/2)-1)]
            TmpA=sum([Id[:,j]*Id[(n-1),:]+Id[:,(n-1)]*Id[j,:] for j in range(ceil(n/2))])+sum([Id[:,j]*Id[j-p[j-ceil(n/2)],:]+Id[:,j-p[j-ceil(n/2)]]*Id[j,:] for j in range(ceil(n/2),n-1)])
            A=TmpA*x1+(ones_matrix(n,n)-TmpA)*x0
            rpt=False
            for M in L:
                if M[1]==A.charpoly():
                    rpt=True
                    break
            if rpt==False:
                L.append([TmpA, A.charpoly(), [(n-1)-j for j in range(1,ceil(n/2))]+[p[j] for j in range((n-1)-floor(n/2)-1)]+[n-1]])
        return L
    else:
        # Initialization of the list of Permutations of elements from 1 to ((n-1)-floor(n/2)-1)
        P=Permutations(Integer(n/2)-1)
        # Initialization of the identity matrix
        Id=identity_matrix(n)
        # Initializing the list of adjacency matrices.
        L=[]
        # Loop constructing the adjacency matrices.
        for q in P:
            p=[q[i] for i in range(Integer(n/2)-1)]
            TmpA=sum([Id[:,j]*Id[(n-1),:]+Id[:,(n-1)]*Id[j,:] for j in range(Integer(n/2))])+sum([Id[:,j]*Id[j-p[j-(n/2)],:]+Id[:,j-p[j-(n/2)]]*Id[j,:] for j in range(Integer(n/2),n-1)])
            A=TmpA*x1+(ones_matrix(n,n)-TmpA)*x0
            rpt=False
            for M in L:
                if M[1]==A.charpoly():
                    rpt=True
                    break
            if rpt==False:
                L.append([TmpA, A.charpoly(), [(n-1)-j for j in range(1,Integer(n/2))]+[p[j] for j in range(Integer(n/2)-1)]+[n-1]])
        return L

def HalfPermutationUndirectedCharPolyRepresentativeDrawings(n, x0=0, x1=1):
    """
    Returns adjacency matrices associated with the
    half permutation construction. The input must be 
    an odd positive integer.

    EXAMPLES:
    ::
        sage: HalfPermutationUndirectedCharPolyRepresentativeDrawings(3)
        [[x^3 - 2*x, [1, 2, 3]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Permutations of elements from 1 to (n-1)
    P=Permutations(floor(n/2))
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initializing the list of adjacency matrices.
    L=[]
    # Loop constructing the adjacency matrices.
    for q in P:
        p=[q[i] for i in range(floor(n/2))]
        # Initialization of the adjacency matrix
        TmpA=sum([Id[:,j]*Id[j+p[j-1],:]+Id[:,j+p[j-1]]*Id[j,:] for j in range(1,ceil(n/2))])+sum([Id[:,j]*Id[0,:]+Id[:,0]*Id[j,:] for j in range(ceil(n/2),n)])
        A=TmpA*x1+(ones_matrix(n,n)-TmpA)*x0
        rpt=False
        for M in L:
            if M[0]==A.charpoly():
                rpt=True
                break
        if rpt==False:
            L.append([A.charpoly(), [p[j-1] for j in range(1,ceil(n/2))]+range(ceil(n/2),n+1)])
            Graph(TmpA).plot().save(str(p).replace(', ','_').replace('[','').replace(']','')+'.png')
    return L

def HalfPermutationUndirectedCharPolyRepresentativeDrawingsII(n, x0=0, x1=1):
    """
    Returns adjacency matrices associated with the
    half permutation construction. The input must be 
    an odd positive integer.

    EXAMPLES:
    ::
        sage: HalfPermutationUndirectedCharPolyRepresentativeDrawingsII(3,2,3)
        [[x^3 - 6*x^2 - 10*x, [1, 2]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    if n%2 == 1:
        # Initialization of the list of Permutations of elements from 1 to ((n-1)-floor(n/2)-1)
        P=Permutations((n-1)-floor(n/2)-1)
        # Initialization of the identity matrix
        Id=identity_matrix(n)
        # Initializing the list of adjacency matrices.
        L=[]
        # Loop constructing the adjacency matrices.
        for q in P:
            p=[q[i] for i in range((n-1)-floor(n/2)-1)]
            TmpA=sum([Id[:,j]*Id[(n-1),:]+Id[:,(n-1)]*Id[j,:] for j in range(ceil(n/2))])+sum([Id[:,j]*Id[j-p[j-ceil(n/2)],:]+Id[:,j-p[j-ceil(n/2)]]*Id[j,:] for j in range(ceil(n/2),n-1)])
            A=TmpA*x1+(ones_matrix(n,n)-TmpA)*x0
            rpt=False
            for M in L:
                if M[0]==A.charpoly():
                    rpt=True
                    break
            if rpt==False:
                L.append([A.charpoly(), [(n-1)-j for j in range(1,ceil(n/2))]+[p[j] for j in range((n-1)-floor(n/2)-1)]+[n-1]])
                Graph(TmpA).plot().save(str([(n-1)-j for j in range(1,ceil(n/2))]+[p[j] for j in range((n-1)-floor(n/2)-1)]+[n-1]).replace(', ','_').replace('[','').replace(']','')+'.png')
        return L
    else:
        # Initialization of the list of Permutations of elements from 1 to ((n-1)-floor(n/2)-1)
        P=Permutations(Integer(n/2)-1)
        # Initialization of the identity matrix
        Id=identity_matrix(n)
        # Initializing the list of adjacency matrices.
        L=[]
        # Loop constructing the adjacency matrices.
        for q in P:
            p=[q[i] for i in range(Integer(n/2)-1)]
            TmpA=sum([Id[:,j]*Id[(n-1),:]+Id[:,(n-1)]*Id[j,:] for j in range(Integer(n/2))])+sum([Id[:,j]*Id[j-p[j-(n/2)],:]+Id[:,j-p[j-(n/2)]]*Id[j,:] for j in range(Integer(n/2),n-1)])
            A=TmpA*x1+(ones_matrix(n,n)-TmpA)*x0
            rpt=False
            for M in L:
                if M[0]==A.charpoly():
                    rpt=True
                    break
            if rpt==False:
                L.append([A.charpoly(), [(n-1)-j for j in range(1,Integer(n/2))]+[p[j] for j in range(Integer(n/2)-1)]+[n-1]])
                Graph(TmpA).plot().save(str([(n-1)-j for j in range(1,Integer(n/2))]+[p[j] for j in range(Integer(n/2)-1)]+[n-1]).replace(', ','_').replace('[','').replace(']','')+'.png')
        return L

def GracefulUndirectedNonTreeAdjacencyMatrix(n):
    """
    Goes through all the list of functions 
    derived from graceful permutations and
    construct directed adjacency matrices 
    of binary trees.

    EXAMPLES:
    ::
        sage: GracefulUndirectedNonTreeAdjacencyMatrix(3)
        [[[], [1, 2]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of Functions.
    L=GracefulFunctions(n)
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the result
    Rslt=[]
    for l in L:
        Rslt.append([[]])
        for f in l[0]:
            # Initializing the adjacency matrix
            A = sum([Id[:,j]*Id[f[j-1],:]+Id[:,f[j-1]]*Id[j,:]  for j in range(1,n)])
            # Testing binary treeness
            if not is_Tree(A):
                Rslt[len(Rslt)-1][0].append(A)
        Rslt[len(Rslt)-1].append(l[len(l)-1])
    return Rslt

def GracefulNonTreeDrawings(n):
    """
    Goes through all the permutation of n>1 elements constructs
    the functions for constructing graceful graphs and saves 
    in the working directory the png files associated with
    drawings of non trees gracefully labeled graphs. This 
    implementation is much better creating images because it 
    never store a list it just saves the image as we go


    EXAMPLES:

    ::

        sage: GracefulNonTreeDrawings(7)
        ['2_4_1_3_5_6_0__1.png',
         '3_4_2_1_5_6_0__0.png',
         '5_1_2_4_3_6_0__3.png',
         '5_3_1_4_2_6_0__0.png']
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of names
    List_of_Names=[]
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the list of Graceful Permutations.
    P=Permutations(n-1)
    # Initialization of the sets.
    S1=Set(range(n)); S2=Set(range(2,n))
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Loop going through the GracefulPermutations.
    for q in P:
        # Appending 0 at end of the permutation
        p=[q[i] for i in range(n-1)]+[0]
        # Initialization of our boolean variable
        bl=True
        for j in range(n-1):
            # Verifying that the graceful criteria is met
            if ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                bl=False
                break
        if bl:
            L.append(p)
            # Initialization of the list of functions.
            c=[]
            for j in range(n-1):
                # testing that only the first criteria is met
                if ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[abs((j+1)-p[j])]
                # testing that only the second criteria is met
                elif ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)+p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)+p[j]]
                # testing that all two criterias are met
                elif ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                        c.append([(j+1)+p[j]])
                    else:
                        d=copy(c)
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)-p[j]]
                        for i in range(len(d)):
                            d[i]=d[i]+[(j+1)+p[j]]
                        c=c+d
            # Put here the stuff that is needed.
            for t in range(len(c)):
                # Initialization of the adjacency matrix.
                A=sum([Id[:,k]*Id[c[t][k-1],:]+Id[:,c[t][k-1]]*Id[k,:]  for k in range(1,n)])
                if not is_Tree(A):
                    Graph(A).plot().save(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
                    List_of_Names.append(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
    return List_of_Names
 
def GracefulTreeDrawingsII(n):
    """
    Goes through all the permutation of n>1 elements constructs
    the functions for constructing graceful graphs and saves 
    in the working directory the png files associated with
    drawings of gracefully labeled trees. This implemntation 
    never stores a list it just saves the images on the go

    EXAMPLES:

    ::
        sage: GracefulTreeDrawingsII(3)
        ['1_2_0__0.png', '1_2_0__1.png']


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of names
    List_of_Names=[]
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the list of Graceful Permutations.
    P=Permutations(n-1)
    # Initialization of the sets.
    S1=Set(range(n)); S2=Set(range(2,n))
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Loop going through the GracefulPermutations.
    for q in P:
        # Appending 0 at end of the permutation
        p=[q[i] for i in range(n-1)]+[0]
        # Initialization of our boolean variable
        bl=True
        for j in range(n-1):
            # Verifying that the graceful criteria is met
            if ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                bl=False
                break
        if bl:
            L.append(p)
            # Initialization of the list of functions.
            c=[]
            for j in range(n-1):
                # testing that only the first criteria is met
                if ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[abs((j+1)-p[j])]
                # testing that only the second criteria is met
                elif ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)+p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)+p[j]]
                # testing that all two criterias are met
                elif ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                        c.append([(j+1)+p[j]])
                    else:
                        d=copy(c)
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)-p[j]]
                        for i in range(len(d)):
                            d[i]=d[i]+[(j+1)+p[j]]
                        c=c+d
            # Put here the stuff that is needed.
            for t in range(len(c)):
                # Initialization of the adjacency matrix.
                A=sum([Id[:,k]*Id[c[t][k-1],:]+Id[:,c[t][k-1]]*Id[k,:]  for k in range(1,n)])
                if is_Tree(A):
                    Graph(A).plot().save(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
                    List_of_Names.append(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
    return List_of_Names

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

def GracefulGeneratorGraphList(n):
    """
    Goes through all the permutation of n>1 elements constructs
    the functions for constructing graceful graphs and saves 
    in the working directory the png files associated with
    drawings of gracefully labeled trees. This implemntation 
    never stores a list it just saves the images on the go

    EXAMPLES:
    ::
        sage: L=GracefulGeneratorGraphList(3)
        sage: L
        []
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the list of Graceful Permutations.
    P=Permutations(n-1)
    # Initialization of the sets.
    S1=Set(range(n)); S2=Set(range(2,n))
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Initializaing the list storing the matrices
    Lm=[]
    # Loop going through the GracefulPermutations.
    for q in P:
        # Appending 0 at end of the permutation
        p=[q[i] for i in range(n-1)]
        # Initialization of our boolean variable
        bl=True
        for j in range(n-1):
            # Verifying that the graceful criteria is met
            if ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                bl=False
                break
        if bl:
            L.append(p)
            # Initialization of the list of functions.
            c=[]
            for j in range(n-1):
                # testing that only the first criteria is met
                if ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[abs((j+1)-p[j])]
                # testing that only the second criteria is met
                elif ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)+p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)+p[j]]
                # testing that all two criterias are met
                elif ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                        c.append([(j+1)+p[j]])
                    else:
                        d=copy(c)
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)-p[j]]
                        for i in range(len(d)):
                            d[i]=d[i]+[(j+1)+p[j]]
                        c=c+d
            # Put here the stuff that is needed.
            for t in range(len(c)):
                # Initialization of the adjacency matrix.
                TmpA=sum([Id[:,k]*Id[c[t][k-1],:]+Id[:,c[t][k-1]]*Id[k,:]  for k in range(1,n)])
                if (TmpA*ones_matrix(TmpA.ncols(),1))[0,0]!=1 and (TmpA*ones_matrix(TmpA.ncols(),1))[TmpA.nrows()-1,0]!=1:
                    Lm.append([TmpA, p])
    return Lm

def GracefulGeneratorTreeDrawings(n):
    """
    Goes through all the permutation of n>1 elements constructs
    the functions for constructing graceful graphs and saves 
    in the working directory the png files associated with
    drawings of gracefully labeled trees. This implemntation 
    never stores a list it just saves the images on the go

    EXAMPLES:
    ::
        sage: GracefulGeneratorTreeDrawings(3)
        []
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the list of Graceful Permutations.
    P=Permutations(n-1)
    # Initialization of the sets.
    S1=Set(range(n)); S2=Set(range(2,n))
    # Intialization of the list collecting the graceful permutations.
    L=[]
    # Initializaing the list storing the matrices
    Lm=[]
    # Loop going through the GracefulPermutations.
    for q in P:
        # Appending 0 at end of the permutation
        p=[q[i] for i in range(n-1)]
        # Initialization of our boolean variable
        bl=True
        for j in range(n-1):
            # Verifying that the graceful criteria is met
            if ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                bl=False
                break
        if bl:
            L.append(p)
            # Initialization of the list of functions.
            c=[]
            for j in range(n-1):
                # testing that only the first criteria is met
                if ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[abs((j+1)-p[j])]
                # testing that only the second criteria is met
                elif ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)+p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)+p[j]]
                # testing that all two criterias are met
                elif ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                        c.append([(j+1)+p[j]])
                    else:
                        d=copy(c)
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)-p[j]]
                        for i in range(len(d)):
                            d[i]=d[i]+[(j+1)+p[j]]
                        c=c+d
            # Put here the stuff that is needed.
            for t in range(len(c)):
                # Initialization of the adjacency matrix.
                TmpA=sum([Id[:,k]*Id[c[t][k-1],:]+Id[:,c[t][k-1]]*Id[k,:]  for k in range(1,n)])
                if is_Tree(TmpA) and (TmpA*ones_matrix(TmpA.ncols(),1))[0,0]!=1 and (TmpA*ones_matrix(TmpA.ncols(),1))[TmpA.nrows()-1,0]!=1:
                    Lm.append([TmpA, p])
                    Graph(TmpA).plot().save(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'__'+str(len(Lm))+'.png')
    return Lm

def GracefulGeneratorGraphDrawings(n):
    """
    Goes through all the permutation of n>1 elements constructs
    the functions for constructing graceful graphs and saves 
    in the working directory the png files associated with
    drawings of gracefully labeled graph. This implemntation 
    never stores a list it just saves the images on the go

    EXAMPLES:
    ::
        sage: GracefulGeneratorGraphDrawings(3)
        []
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the list of Graceful Permutations.
    P=Permutations(n-1)
    # Initialization of the sets.
    S1=Set(range(n)); S2=Set(range(2,n))
    # Initializaing the list storing the matrices
    Lm=[]
    # Loop going through the GracefulPermutations.
    for q in P:
        # Appending 0 at end of the permutation
        p=[q[i] for i in range(n-1)]
        # Initialization of our boolean variable
        bl=True
        for j in range(n-1):
            # Verifying that the graceful criteria is met
            if ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                bl=False
                break
        if bl:
            # Initialization of the list of functions.
            c=[]
            for j in range(n-1):
                # testing that only the first criteria is met
                if ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[abs((j+1)-p[j])]
                # testing that only the second criteria is met
                elif ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)+p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)+p[j]]
                # testing that all two criterias are met
                elif ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                        c.append([(j+1)+p[j]])
                    else:
                        d=copy(c)
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)-p[j]]
                        for i in range(len(d)):
                            d[i]=d[i]+[(j+1)+p[j]]
                        c=c+d
            # Put here the stuff that is needed.
            for t in range(len(c)):
                # Initialization of the adjacency matrix.
                TmpA=sum([Id[:,k]*Id[c[t][k-1],:]+Id[:,c[t][k-1]]*Id[k,:]  for k in range(1,n)])
                if (TmpA*ones_matrix(TmpA.ncols(),1))[0,0]!=1 and (TmpA*ones_matrix(TmpA.ncols(),1))[TmpA.nrows()-1,0]!=1:
                    Lm.append([TmpA, p])
                    Graph(TmpA).plot().save(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'__'+str(len(Lm))+'.png')
    return Lm

def FastGracefulGeneratorTreeDrawings(n):
    """
    Goes through all the permutation of n>1 elements constructs
    the functions for constructing graceful graphs and saves 
    in the working directory the png files associated with
    drawings of gracefully labeled trees. This implemntation 
    never stores a list it just saves the images on the go

    EXAMPLES:
    ::
        sage: FastGracefulGeneratorTreeDrawings(5)
        ['1_2_3_4__2.png', '1_2_3_4__3.png', '3_2_1_4__0.png', '3_2_1_4__2.png']
 

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of names
    List_of_Names=[]
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the list of Graceful Permutations.
    P=Permutations(n-1)
    # Initialization of the sets.
    S1=Set(range(n)); S2=Set(range(2,n))
    # Intialization of the list collecting the graceful permutations.
    # Loop going through the GracefulPermutations.
    for q in P:
        # Appending 0 at end of the permutation
        p=[q[i] for i in range(n-1)]
        # Initialization of our boolean variable
        bl=True
        for j in range(n-1):
            # Verifying that the graceful criteria is met
            if ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                bl=False
                break
        if bl:
            # Initialization of the list of functions.
            c=[]
            for j in range(n-1):
                # testing that only the first criteria is met
                if ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[abs((j+1)-p[j])]
                # testing that only the second criteria is met
                elif ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)+p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)+p[j]]
                # testing that all two criterias are met
                elif ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                        c.append([(j+1)+p[j]])
                    else:
                        d=copy(c)
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)-p[j]]
                        for i in range(len(d)):
                            d[i]=d[i]+[(j+1)+p[j]]
                        c=c+d
            # Put here the stuff that is needed.
            for t in range(len(c)):
                # Initialization of the adjacency matrix.
                TmpA=sum([Id[:,k]*Id[c[t][k-1],:]+Id[:,c[t][k-1]]*Id[k,:]  for k in range(1,n)])
                if is_Tree(TmpA) and (TmpA*ones_matrix(TmpA.ncols(),1))[0,0]!=1 and (TmpA*ones_matrix(TmpA.ncols(),1))[TmpA.nrows()-1,0]!=1:
                    Graph(TmpA).plot().save(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
                    List_of_Names.append(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
    return List_of_Names

def FastGracefulGeneratorGraphDrawings(n):
    """
    Goes through all the permutation of n>1 elements constructs
    the functions for constructing graceful graphs and saves 
    in the working directory the png files associated with
    drawings of gracefully labeled graphs. This implementation 
    never stores a list it just saves the images on the go

    EXAMPLES:

    ::

        sage: FastGracefulGeneratorGraphDrawings(5)
        ['1_2_3_4__2.png', '1_2_3_4__3.png', '3_2_1_4__0.png', '3_2_1_4__2.png']


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the list of names
    List_of_Names=[]
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Initialization of the list of Graceful Permutations.
    P=Permutations(n-1)
    # Initialization of the sets.
    S1=Set(range(n)); S2=Set(range(2,n))
    # Loop going through the GracefulPermutations.
    for q in P:
        # Appending 0 at end of the permutation
        p=[q[i] for i in range(n-1)]
        # Initialization of our boolean variable
        bl=True
        for j in range(n-1):
            # Verifying that the graceful criteria is met
            if ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                bl=False
                break
        if bl:
            # Initialization of the list of functions.
            c=[]
            for j in range(n-1):
                # testing that only the first criteria is met
                if ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) not in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[abs((j+1)-p[j])]
                # testing that only the second criteria is met
                elif ( ((j+1)-p[j]) not in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)+p[j]])
                    else:
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)+p[j]]
                # testing that all two criterias are met
                elif ( ((j+1)-p[j]) in S1 ) and ( ((j+1)+p[j]) in S2 ):
                    if len(c)==0:
                        c.append([(j+1)-p[j]])
                        c.append([(j+1)+p[j]])
                    else:
                        d=copy(c)
                        for i in range(len(c)):
                            c[i]=c[i]+[(j+1)-p[j]]
                        for i in range(len(d)):
                            d[i]=d[i]+[(j+1)+p[j]]
                        c=c+d
            # Put here the stuff that is needed.
            for t in range(len(c)):
                # Initialization of the adjacency matrix.
                TmpA=sum([Id[:,k]*Id[c[t][k-1],:]+Id[:,c[t][k-1]]*Id[k,:]  for k in range(1,n)])
                if (TmpA*ones_matrix(TmpA.ncols(),1))[0,0]!=1 and (TmpA*ones_matrix(TmpA.ncols(),1))[TmpA.nrows()-1,0]!=1:
                    Graph(TmpA).plot().save(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
                    List_of_Names.append(str(p).replace(', ','_').replace('[','').replace(']','')+'__'+str(t)+'.png')
    return List_of_Names

def TreeAdjacencyMatrixList(n):
    """
    Computes the list of adjacency matrices of trees on n
    vertices via the Matrix Tree theorem.

    EXAMPLES:
    ::
        sage: TreeAdjacencyMatrixList(3)
        [[
        [0 1 1]            
        [1 0 0]            
        [1 0 0], [a02, a10]
        ],
         [
        [0 1 0]            
        [1 0 1]            
        [0 1 0], [a01, a12]
        ],
         [
        [0 0 1]            
        [0 0 1]            
        [1 1 0], [a02, a12]
        ]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the symbolic matrix.
    TmpA = Matrix(SR, n, n, [var('a'+str(i)+str(j)) for i in range(n) for j in range(n)])
    # Initialization of the submatrix.
    A = ((diagonal_matrix((TmpA*ones_matrix(n,1)).list())-TmpA)[0:n-1,0:n-1])
    # Initializing the list of permutations.
    P = Permutations(range(A.nrows()))
    # Initialization of the edge list
    L = sum([Permutation([p[i]+1 for i in range(len(p))]).signature()*prod([A[k][p[k]] for k in range(A.nrows())]) for p in P]).expand().operands()
    # Initialization of the final list.
    Lt = []
    for f in L:
        # Initializing the list of edges.
        l=f.operands()
        # Initialization of the corresponding adjacency matrix.
        B=Matrix(SR, zero_matrix(n,n))
        for i in range(n):
            for j in range(n):
                if not Set([TmpA[i,j]]).intersection(Set(l)).is_empty():
                    B[i,j]=1; B[j,i]=1
        Lt.append([B, l]) 
    return Lt
