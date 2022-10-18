"""
Diana Clustering module for clustering RNA structures according to their basepairdistance.
The cluster distance function computes the distance between the centroids of the clusters.

! This algorithm may present different results than the algorithm in R. 
The reason is that the maximal diameter (criterion for selecting cluster to split) is not unique !
A second reason is that the object with maximal average distance is not unique. 

"""

import sys, math, csv, RNA, numpy,re


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
            match = re.search('([\.\(\)]+)', l)
            if match:
                sequence = match.group(1)
                res.append(sequence.strip())
    return res


class Cluster:
    def __init__(self):
        self.structures = []
        self.representative = ""
        self.childNodes = []

def computeBasePairDistanceMatrix(structs):
    matrix = numpy.zeros((len(structs),len(structs)),dtype=numpy.int)
    for i in range(len(structs)):
        for j in range(i+1,len(structs)):
            dist = RNA.bp_distance(structs[i], structs[j])
            matrix[i][j] = dist
            matrix[j][i] = dist
    return matrix

def averageDissimilarity(index,cluster,dissMatrix):
    """
    computes the average dissimilarity
    index is a structure index, which is conform the the index in the dissMatrix.
    cluster is a list of structure indices.
    """
    avgDiss = 0.0
    subtract = 0.0
    if index in cluster:
        subtract = 1.0
    for i in cluster:
        avgDiss += float(dissMatrix[index][i])
    NumberOfElementsWithoutIndexElement = float(len(cluster)-subtract)
    if NumberOfElementsWithoutIndexElement > 0.0:
        avgDiss = avgDiss / NumberOfElementsWithoutIndexElement
    else:
        avgDiss = 0.0
    return avgDiss
    
def objectIndexWithHighestDissimilarity(cluster,dissMatrix):
    maxAvgDiss = 0.0
    maxIndex = -1
    for i1 in cluster:
        avgDiss = averageDissimilarity(i1, cluster, dissMatrix)
        if (avgDiss > maxAvgDiss):
            maxAvgDiss = avgDiss
            maxIndex = i1
    return maxIndex     
    
def diameter(cluster, dissMatrix):
    maxDist = 0
    for i in range(len(cluster)):
        for j in range(i+1,len(cluster)):
            dist = dissMatrix[cluster[i]][cluster[j]]
            if dist > maxDist:
                maxDist = dist
    return maxDist

def WmAll(c,dissMatrix):
    """Computes the 'within-ness' measure for a cluster, with all pairwise comparisons.

    Args:
        c (list): List of secondary structures.

    Returns:
        set: The within-ness measure.
    """
    sum = 0.0
    for i in range(0,len(c)):
        for j in range(i+1,len(c)):
            sum = sum + math.pow(dissMatrix[c[i]][c[j]],2)
    return sum

def varianceRatio(parentCluster,clusterA, clusterB, dissMatrix):
    """
    compute the variance of the old clusters and the new cluster.
    Args:
        cl1 (list): Former clustering, list of lists of secondary structures.
        cl2 (list): Current clustering, list of lists of secondary structures.
    """
    
    wmA = WmAll(clusterA.structures,dissMatrix)
    wmB = WmAll(clusterB.structures,dissMatrix)
    wmParent = WmAll(parentCluster.structures,dissMatrix)
    
    withinRatio = (wmA+wmB)/wmParent
    randomRatio = (pow(len(clusterA.structures),2.)+pow(len(clusterB.structures),2.))/(pow(len(parentCluster.structures),2.))
    
    return withinRatio, randomRatio

def maxDiameterInLeafNodes(clusterTree,dissMatrix):
        maxDiameter = -1
        maxCluster = None
        if len(clusterTree.childNodes) > 0:
            for cn in clusterTree.childNodes:
                dia, cluster = maxDiameterInLeafNodes(cn,dissMatrix)
                if dia > maxDiameter:
                    maxDiameter = dia
                    maxCluster = cluster
        else:
            maxDiameter = diameter(clusterTree.structures,dissMatrix)
            maxCluster = clusterTree
        
        return maxDiameter, maxCluster


def createClusterTree(c_root, dissMatrix):
    """
    The core of the diana algorithm (recursive function).
    c_root = the rootnode of the clusterTree. It contains the main cluster as childnode.
    """
    #1. select cluster with the largest diameter from all leafnodes.
    maxDiameter, c_m = maxDiameterInLeafNodes(c_root, dissMatrix)
    
    #if maxDiameter <=8:
    #    return      
            
    if c_m == None:
        return
    if len(c_m.structures) > 1:
        #2. object with highest dissimilarity to all others defines the new cluster (c_newA).
        o = objectIndexWithHighestDissimilarity(c_m.structures,dissMatrix)
        c_newA = Cluster()
        c_newA.structures.append(o)
        # B contains all other elements
        c_newB = Cluster()
        c_newB.structures.extend([x for x in c_m.structures if x != o])
        
        #3. select best elements for each cluster (move closest elements from B to A).
        #for each object outside the splinter group:
        positiveDiExists = True
        while(positiveDiExists):
            positiveDiExists = False
            largestD_i = 0.0
            bestElement = None
            for j in c_newB.structures:
                d_i = averageDissimilarity(j, c_newB.structures, dissMatrix) - averageDissimilarity(j, c_newA.structures, dissMatrix)
                if (d_i > largestD_i):
                    largestD_i = d_i
                    bestElement = j
            if (bestElement != None) & (len(c_newB.structures)>1):
                positiveDiExists = True
                c_newA.structures.append(bestElement)
                c_newB.structures.remove(bestElement)
        
        ratios = varianceRatio(c_m, c_newA, c_newB, dissMatrix)
        criterion = ratios[0] < ratios[1]
        print ratios[0], ratios[1]
        if(criterion):
            c_m.childNodes.append(c_newB)
            c_m.childNodes.append(c_newA)
            
        #build the subtrees for each child
        createClusterTree(c_root, dissMatrix)

def generateID():
    k = 0
    while True:
       k += 1
       yield k
ID = generateID().next

def printLeafNodes(clusterTree,structs):
        if len(clusterTree.childNodes) > 0:
            for cn in clusterTree.childNodes:
                printLeafNodes(cn,structs)
        else:
            index = ID()
            print "ClusterID:",index
            for i in clusterTree.structures:
                print structs[i]
 
def doClustering(structs, threshold = None):
    """
    Computes the DIANA clustering
    
    threshold = minimal basepairdistance for belonging to a cluster.
    """

    
    #start DIANA
    #1. initialization
    dissMatrix = computeBasePairDistanceMatrix(structs)
    c_root = Cluster()
    c_root.structures.extend([x for x in range(0,len(structs))])
    
    #start the real DIANA algorithm.
    createClusterTree(c_root, dissMatrix)
    #end DIANA
          
    #print the clusters (all leafs)
    printLeafNodes(c_root,structs)
    

if __name__ == "__main__":
    structureFileName = sys.argv[1]
    structs = parseFile(structureFileName)
    doClustering(structs)












