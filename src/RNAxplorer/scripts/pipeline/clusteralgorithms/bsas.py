"""
BSAS Clustering module for clustering RNA structures according to their basepairdistance.
The cluster distance function computes the distance between the centroids of the clusters.
"""

import sys, RNA


def parseFile(fpath):
    """
    Parses a FASTA-ish file and extract a list of secondary structures.

    @param fpath (string): Path to file in FASTA-ish notation, featuring a list
        of Vienna-formatted secondary structures.

    @return list: List of secondary structures, each represented as a pair `(v,bps)`,
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

    @param ss (string): Vienna-formatted secondary structure.
    @return list: Set of base-pairs in secondary structures.
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

    @param ss1 (list): Secondary structure represented as a list of parenthesis.
    @param ss2 (list): Secondary structure represented as a list of parenthesis.

    @return int: Distance between `ss1` and `ss2`.
    """
    return len(ss1^ss2)

def centroid(cluster):
    """
    Computes the centroid structure (~the median) for a given structure.

    @param cluster (list): List of secondary structures.
    @return set: Set of base-pairs defining the centroid structure.
    """    
    allBPs = set()
    for ss in cluster:
        allBPs |= ss
    counts = {b:0 for b in allBPs}
    for ss in cluster:
        for b in ss:
            counts[b] += 1
    return set([b for b in counts if counts[b]>len(cluster)/2.])

def convertBack(ss,n):
    """
    Renders a secondary structure, denoted by a list of base-pairs and total length,
    as a Vienna-style well parenthesized expression.

    @param ss (list): List of base-pairs.
    @param n (int): Length of secondary structure.
    """
    res = ['.' for c in range(n)]
    for a,b in ss:
        res[a]='('
        res[b]=')'
    return "".join(res)


class Cluster:
    def __init__(self):
        self.structures = []
        self.representative = set()
    


def selectClosestClusterIndex(clusters, structure):
    """
    @param clusters = a dictionary with integers as indices and structureIndices as values
    @param structure = structure as basepair list.
    """
    keys = clusters.keys()
    minKey = ""
    minDist = sys.maxint
    for k in keys:
        maxDistInCluster=0
        for s in clusters[k].structures:
            dist = RNA.bp_distance(structure, s)
            if maxDistInCluster < dist:
                maxDistInCluster = dist
        dist2 = maxDistInCluster
        if dist2 < minDist:
            minDist = dist2
            minKey = k
    return (minKey,minDist)

class BSAS:
    @staticmethod             
    def doClustering(structs, threshold, maxClusters):
        """
        Computes the basic sequential algorithmic scheme clustering
        
        @param structs - list - secondary structures in dot-bracket notation
        @param threshold - int - minimal basepairdistance for belonging to a cluster.
        @param maxClusters - int - maximum number of clusters, that the algorithm will create.
        @return list of lists - clusters
        """
        #1. m = 1; Cm = {x1}; // Init first cluster = first sample     
        #2. for every sample x from 2 to N      
        #   a. find cluster Ck such that min d(x, Ck)     
        #   b. if d(x, Ck) > threshold AND (m < maxNumberOfClusters)
        #     i. m = m + 1; Cm = {x} // Create a new cluster      
        #   c. else      
        #     i. Ck = Ck + {x} // Add sample to the nearest cluster      
        #     ii. Update representative if needed      
        #3. end algorithm
    
        structs = list(structs)
        #start BSAS
        allClusters = {}
        #1. initialization
        m = 1
        CM = Cluster()
        firstStructure = structs[0]
        CM.structures.append(firstStructure)
        #CM.representative = firstStructure
        allClusters[m] = CM
        #2.
        for i in range(1,len(structs)):
            #a.
            structure = structs[i]
            clusterIndex, minDist = selectClosestClusterIndex(allClusters,structure)
            #b.
            if minDist > threshold and (m < maxClusters):
                #i.
                m = m + 1
                allClusters[m] = Cluster()
                allClusters[m].structures.append(structure)
            #c.
            else:
                #i.
                allClusters[clusterIndex].structures.append(structure)
        #end BSAS
            
        
        #convert to list of lists
        resultClusters = []
        keys = allClusters.keys()
        for k in keys:
            c = allClusters[k]
            resultClusters.append(c.structures)
            
        return resultClusters
    

if __name__ == "__main__":
    structureFileName = sys.argv[1]
    structs = parseFile(structureFileName)
    threshold = 4 #int(math.floor(len(structs[0][1])*0.3))
    maxClusters = sys.maxint
    print "# Threshold:",threshold,"Max Clusters:",maxClusters
    Bsas.doClustering(structs, threshold, maxClusters)












