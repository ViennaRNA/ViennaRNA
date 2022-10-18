"""
Minimal RNA clustering module.

Provides basic utility functions to analyze the output of
the DIANA hierarchical clustering algorithm, and chose the
best number of clusters for of RNA secondary structures
obtained from sampling methods.
"""

import sys, math, csv


def parseFile(fpath):
    """
    Parses a FASTA-ish file and extract a list of secondary structures.

    Args:
        fpath (string): Path to file in FASTA-ish notation, featuring a list
        of Vienna-formatted secondary structures.

    Returns:
        list: List of secondary structures, each represented as a pair `(v,bps)`,
        where `v` is the Vienna notation and `bps` is a set of base-pairs.
    """
    res = []
    for l in open(fpath):
        if not l.startswith("#") and not l.startswith(";") and not l.startswith(">"):
            data = l.split()
            if len(data)==1:
                ss = data[0]
                res.append((ss.strip(),parseSS(ss.strip())))
    return res

def parseSS(ss):
    """
    Parses a secondary structure denoted by a well parenthesized expression.

    Args:
        ss (string): Vienna-formatted secondary structure.

    Returns:
        list: Set of base-pairs in secondary structures.
    """
    p = []
    res = []
    for i,c in enumerate(ss):
        if c == "(":
            p.append(i)
        elif  c == ")":
            j = p.pop()
            res.append((j,i))
    return set(res)

def computeDistance(ss1,ss2):
    """
    Computes the base-pair distance between to secondary structures.

    Args:
        ss1 (list): Secondary structure represented as a list of parenthesis.
        ss2 (list): Secondary structure represented as a list of parenthesis.

    Returns:
        int: Distance between `ss1` and `ss2`.
    """
    return len(ss1^ss2)

def CHIndex(clusters):
    """Computes the CH-Index associated with a bunch of clusters.

    Args:
        clusters (list): List of secondary structure clusters.

    Returns:
        real: CH-Index for the proposed clustering.

    """
    n = sum([len(clusters[i]) for i in clusters])
    k = len(clusters)
    return (B(clusters)/(k-1.))/(W(clusters)/(n-k))

def centroid(cluster):
    """
    Computes the centroid structure (~the median) for a given structure.

    Args:
        cluster (list): List of secondary structures.

    Returns:
        set: Set of base-pairs defining the centroid structure.
    """    
    allBPs = set()
    for ss in cluster:
        allBPs |= ss
    counts = {b:0 for b in allBPs}
    for ss in cluster:
        for b in ss:
            counts[b] += 1
    return set([b for b in counts if counts[b]>len(cluster)/2.])

def Wm(c):
    """Computes the 'within-ness' measure for a cluster, computed
    with respect to the centroid to avoid pairwise comparisons.

    Args:
        c (list): List of secondary structures.

    Returns:
        set: The within-ness measure.
    """
    cent = centroid(c)
    return sum([math.pow(computeDistance(ss,cent),2) for ss in c])
    
def W(clusters):
    """
    Computes the total 'within-ness' for a clustering .

    Args:
        clusters (list): List of clusters, ie lists of secondary structures.

    Returns:
        set: The within-ness measure for the proposed clustering.
    """
    res = 0.
    for m in clusters:
        c = clusters[m]
        res += Wm(c)
    return res

def B(clusters):
    """
    Computes the 'between-ness' measure for a cluster, computed
    with respect to the centroids to avoid pairwise comparisons.

    Args:
        clusters (list): List of secondary structures.

    Returns:
        set: The between-ness measure.
    """
    res = 0.
    for m1 in clusters:
        c1 = clusters[m1]
        cent1 = centroid(c1)                
        for m2 in clusters:
            c2 = clusters[m2]
            cent2 = centroid(c2)                
            if m1<m2:
                res += sum([math.pow(computeDistance(ss1,cent2),2) for ss1 in c1])
                res += sum([math.pow(computeDistance(ss2,cent1),2) for ss2 in c2])
    return res

def loadClusters(fclusters,structs):
    """
    Loads clusters computed externally using R's implementation of the
    DIANA algorithm in the cluster package for hierarchical clustering,
    and explicitly builds the associated clusters.

    Args:
        fclusters (string): Path to output of R script.
        structs (list): List of secondary structure represented as pairs `(v,bps)`,
            where `v` is the Vienna notation and `bps` is a set of base-pairs.

    Returns:
        list: A list of clustering propositions (ordered increasingly on number of classes).
    """    
    with open(fclusters, 'rb') as csvfile:
        rr = csv.reader(csvfile, delimiter=',', quotechar='"')
        tab = [[x for x in row][1:] for row in rr][1:]
        maxNumClust = len(tab[0])
        clusters = [{i:[] for i in range(1,k+1)} for k in range(1,maxNumClust+1)]
        for i,r in enumerate(tab):
            for j in range(maxNumClust):
                v = r[j]
                clusters[j][int(v)].append(structs[i][1])
        return clusters

def refinementQuality(cl1,cl2):
    """
    Decides whether or not to accept the subdivision proposed by the clustering.
    To that purpose, compares the scores of the former and current clustering, and
    compares the (necessary) observed improvement to that expected in a random model.

    Args:
        cl1 (list): Former clustering, list of lists of secondary structures.
        cl2 (list): Current clustering, list of lists of secondary structures.

    Returns:
        Boolean: True, if the improvement in quality is greater than expected at
        random, False otherwise.
    """    
    size1 = {}
    size2 = {}
    for i in cl1:
        c = cl1[i]
        for ss in c:
            size1[tuple(ss)] = (len(c),i)
    for i in cl2:
        c = cl2[i]
        for ss in c:
            size2[tuple(ss)] = (len(c),i)
    oldClusts=set()
    newClusts=set()
    
    for ss in size1:
        l1,i1 = size1[ss]
        l2,i2 = size2[ss]
        if l1!=l2:
            oldClusts.add(i1)
            newClusts.add(i2)
    old = list(oldClusts)[0]
    a,b = tuple(newClusts)
    return (Wm(cl2[a])+Wm(cl2[b]))/(Wm(cl1[old]))<(pow(len(cl2[a]),2.)+pow(len(cl2[b]),2.))/(pow(len(cl1[old]),2.))


def convertBack(ss,n):
    """
    Renders a secondary structure, denoted by a list of base-pairs and total length,
    as a Vienna-style well parenthesized expression.

    Args:
        ss (list): List of base-pairs.
        n (int): Length of secondary structure.

    Returns:
        String: Dot-Bracket (Vienna) notation for secondary structure.
    """
    res = ['.' for c in range(n)]
    for a,b in ss:
        res[a]='('
        res[b]=')'
    return "".join(res)

def displayBestCluster(cl,n):
    """
    Displays a clustering as a list of clusters (lists of structures), each
    preceded by its centroid. Each structure/centroid is rendered using the DBN
    notation.

    Args:
        cl (list): List of clusters.
        n (int): Length of secondary structures.
    """    
    for num in cl:
        print "ClusterID:",num
        print convertBack(centroid(cl[num]),n)
        for ss in cl[num]:
            print convertBack(ss,n)
        
def doClustering(structs,fclusters):
    """
    Computes and displays the best clustering of a set of RNA secondary
    structures, determining the best possible number of clusters.
    A hierarchical clustering must have been computed externally using R's
    implementation of the DIANA algorithm in the `cluster` package.

    Args:
        fpath (string): Path to list of secondary structures.
        fclusters (string): Path to output of R script.
    """        
    n = sum([len(s) for s,ss in structs])/len(structs)
    clusters = loadClusters(fclusters,structs)
    maxCH = -sys.maxint
    cl1 = None
    for i,cl in enumerate(clusters):
        if i>0 and i<(len(clusters)-1):
            locCH = CHIndex(cl)
            print i+1, locCH
            if cl1 is not None:
                cl2 = cl
                keepGoing = refinementQuality(cl1,cl2)
                if not keepGoing:
                    displayBestCluster(clusters[i-1],n)
                    return clusters[i-1]
            cl1 = cl
    return clusters[-1]
    


if __name__ == "__main__":
    fpath = sys.argv[1]
    fclusters = sys.argv[2]
    structs = parseFile(fpath)
    doClustering(structs,fclusters)
    
