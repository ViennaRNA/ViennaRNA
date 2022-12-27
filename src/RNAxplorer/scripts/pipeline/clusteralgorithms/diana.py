#!/usr/bin/python3
"""
Diana Clustering module for clustering RNA structures according to their basepairdistance.
The cluster distance function computes the distance between the centroids of the clusters.

! This algorithm may present different results than the algorithm in R. 
The reason is that the maximal diameter (criterion for selecting cluster to split) is not unique !
A second reason is that the object with maximal average distance is not unique. 

"""

import sys, math, RNA, numpy, argparse, re, resource
from multiprocessing import Pool

sys.setrecursionlimit(int(pow(2,31)-1))
resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))


class Cluster:
    """
    diana internal tree-like cluster structure.
    """
    def __init__(self):
        self.structures = []
        self.representative = ""
        self.childNodes = []

def computeBasePairDistanceMatrix(structs):
    len_structs = len(structs)
    matrix = numpy.zeros((len_structs, len_structs), dtype=numpy.int)
    pool = Pool(4)
    res_list = []
    for i in range(len_structs):
        for j in range(i + 1, len_structs):
            func = pool.apply_async(RNA.bp_distance, args=(structs[i], structs[j]))
            res_list.append(func)
    res_index = 0
    for i in range(len_structs):
        for j in range(i + 1, len_structs):
            dist = res_list[res_index].get()
            #dist = RNA.bp_distance(structs[i], structs[j])
            matrix[i][j] = dist
            matrix[j][i] = dist
            res_index += 1
    return matrix

def averageDissimilarity(index, cluster, dissMatrix):
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
    NumberOfElementsWithoutIndexElement = float(len(cluster) - subtract)
    if NumberOfElementsWithoutIndexElement > 0.0:
        avgDiss = avgDiss / NumberOfElementsWithoutIndexElement
    else:
        avgDiss = 0.0
    return avgDiss
    
def objectIndexWithHighestDissimilarity(cluster, dissMatrix):
    maxAvgDiss = 0.0
    maxIndex = -1
    for i1 in cluster:
        avgDiss = averageDissimilarity(i1, cluster, dissMatrix)
        if (avgDiss > maxAvgDiss):
            maxAvgDiss = avgDiss
            maxIndex = i1
    return maxIndex     
    
def diameter(cluster, dissMatrix):
    """
    computes the largest basepairdistance between all pairs of structures in the cluster.
    """
    maxDist = 0
    for i in range(len(cluster)):
        for j in range(i + 1, len(cluster)):
            dist = dissMatrix[cluster[i]][cluster[j]]
            if dist > maxDist:
                maxDist = dist
    return maxDist

def db_structure_to_bp_set(ss):
    p = []
    res = set()
    for i,c in enumerate(ss):
        if c == "(":
            p.append(i)
        elif  c == ")":
            j = p.pop()
            res.add((j,i))
    return res

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

def centroid(cluster, structures_bp_sets):
    """
    Computes the centroid structure (~the median) for a given structure.

    Args:
        cluster (list): List of secondary structures.

    Returns:
        set: Set of base-pairs defining the centroid structure.
    """    
    allBPs = set()
    for si in cluster:
        allBPs |= structures_bp_sets[si]
    counts = {b:0 for b in allBPs}
    for si in cluster:
        ss = structures_bp_sets[si]
        for b in ss:
            counts[b] += 1
    return set([b for b in counts if counts[b]>len(cluster)/2.])

def Wm(c, cent, structures_bp_sets):
    """Computes the 'within-ness' measure for a cluster, computed
    with respect to the centroid to avoid pairwise comparisons.

    Args:
        c (list): List of secondary structures.

    Returns:
        set: The within-ness measure.
    """
    #cent = centroid(c, structures_bp_sets)
    sum_dist = 0
    for si in c:
        sum_dist += computeDistance(structures_bp_sets[si],cent)
    return sum_dist

def pooled_within_group_sum_of_squares(clusters, structures_bp_sets, len_clusters, cluster_centroids):
    pwgss = 0
    if clusters == None:
        return pwgss
    for i in range(len(clusters)):
        c = clusters[i]
        len_c = len_clusters[i]
        cent = cluster_centroids[i]
        pwgss += Wm(c, cent, structures_bp_sets) / len_c
        """
        faster alternatice is could be sum_i_N(pow(bpdist(centroid_structure(c), structure(c_i)),2))
        """
    return pwgss

def pooled_between_group_sum_of_squares(clusters, structures_bp_sets, dataset_centroid, len_clusters, cluster_centroids):
    pbgss = 0
    if clusters == None:
        return pbgss
    
    for i in range(len(clusters)):
        c = clusters[i]
        len_c = len_clusters[i]
        cluster_centoid = cluster_centroids[i]
        pbgss += computeDistance(cluster_centoid, dataset_centroid) * len_c
        """
        faster alternatice is could be sum_i_N(pow(bpdist(centroid_structure(dataset), centroid_structure(c_i)),2)) /float(len(c))
        """
    return pbgss

def calinski_harabasz_index(clusters, structures_bp_sets, dataset_centroid):
    ch = 0
    len_clusters = [ float(len(x)) for x in clusters ]
    cluster_centroids = [ centroid(c, structures_bp_sets) for c in clusters ]
    pbgss = pooled_between_group_sum_of_squares(clusters, structures_bp_sets, dataset_centroid, len_clusters, cluster_centroids)
    pwgss = pooled_within_group_sum_of_squares(clusters, structures_bp_sets, len_clusters, cluster_centroids)
    number_of_clusters = len(clusters)
    number_of_structures = len(structures_bp_sets)
    if (number_of_clusters <= 1):
        ch = 0
    elif (pwgss == 0):
        ch = float('inf')
    else:
        ch =  ((pbgss) / (pwgss)) * ((number_of_structures - number_of_clusters) / float(number_of_clusters -1))
    return ch
    

def maxDiameterInLeafNodes(clusterTree, dissMatrix):
    maxDiameter = -1
    maxCluster = None
    if len(clusterTree.childNodes) > 0:
        for cn in clusterTree.childNodes:
            dia, cluster = maxDiameterInLeafNodes(cn, dissMatrix)
            if dia > maxDiameter:
                maxDiameter = dia
                maxCluster = cluster
    else:
        maxDiameter = diameter(clusterTree.structures, dissMatrix)
        maxCluster = clusterTree
    
    return maxDiameter, maxCluster

def averageDiameterInLeafNodes(clusterTree, dissMatrix):
    """
    also known as Divisive Coefficient (DC).
    """
    sumDiameters = 0
    numberOfClusters = 0
    stackChildNodes = []
    stackChildNodes.append(clusterTree)
    while(len(stackChildNodes) > 0):
        cn = stackChildNodes.pop()
        if len(cn.childNodes) == 0:
            #is leaf node/cluster --> sum
            sumDiameters += diameter(cn.structures, dissMatrix)
            numberOfClusters += 1
        else:    
            for c in cn.childNodes:
                stackChildNodes.append(c)
                
    averageDiameter = sumDiameters / float(numberOfClusters)
    return averageDiameter

def createClusterTree(c_root, dissMatrix, maxDiameterThreshold, maxAverageDiameterThreshold, do_ch_first_local_min = False, prev_ch = None, structures_bp_sets = None, dataset_centroid = None,  hierarchy=1):
    """
    The core of the diana algorithm (recursive function).
    c_root = the rootnode of the clusterTree. It contains the main cluster as childnode.
    """
    # 1. select cluster with the largest diameter from all leafnodes.
    dc = averageDiameterInLeafNodes(c_root, dissMatrix)
    
    if dc < maxAverageDiameterThreshold:
        return
    if maxDiameterThreshold >= 0:
        maxDiameter, c_m = maxDiameterInLeafNodes(c_root, dissMatrix)
        if maxDiameter <= maxDiameterThreshold:
            return      
            
    if c_m == None:
        return
    
    #print(hierarchy, dc)
    if len(c_m.structures) > 1:
        hierarchy+=1
        # 2. object with highest dissimilarity to all others defines the new cluster (c_newA).
        o = objectIndexWithHighestDissimilarity(c_m.structures, dissMatrix)
        c_newA = Cluster()
        c_newA.structures.append(o)
        # B contains all other elements
        c_newB = Cluster()
        c_newB.structures.extend([x for x in c_m.structures if x != o])
        
        # 3. select best elements for each cluster (move closest elements from B to A).
        # for each object outside the splinter group:
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
            if (bestElement != None) & (len(c_newB.structures) > 1):
                positiveDiExists = True
                c_newA.structures.append(bestElement)
                c_newB.structures.remove(bestElement)
        
        c_m.childNodes.append(c_newB)
        c_m.childNodes.append(c_newA)

        if do_ch_first_local_min:
            cis = convertTreeToListOfClusters_indices(c_root, [])
            #print(cis, len(cis))
            if len(cis) > 1:
                ch = calinski_harabasz_index(cis, structures_bp_sets, dataset_centroid)
                if prev_ch != None and ch < prev_ch:
                    if len(c_m.childNodes) > 0:
                        c_m.childNodes = []
                        print(hierarchy, "ch", "{:10.3f}".format(ch))
                    return # ch index > last ch index --> first local min --> break
                print(hierarchy, "ch", "{:10.3f}".format(ch))
                prev_ch = ch
        
        #remove c_m.structures because it is an inner node
        c_m.structures.clear()
        
        # build the subtrees for each child
        createClusterTree(c_root, dissMatrix, maxDiameterThreshold, maxAverageDiameterThreshold, do_ch_first_local_min, prev_ch, structures_bp_sets, dataset_centroid, hierarchy)


def convertTreeToListOfClusters_indices(clusterTree, clusterList=[]):
    if len(clusterTree.childNodes) > 0:
        for cn in clusterTree.childNodes:
            convertTreeToListOfClusters_indices(cn, clusterList)
    else:
        cluster = [ c for c in clusterTree.structures]
        clusterList.append(cluster)
    return clusterList 
         
class DIANA:
    @staticmethod
    def doClustering(structures, maxDiameter, maxAverageDiameter, do_ch_first_local_min = False):
        """
        Computes the DIANA clustering
        
        Args:
            maxDiameter = threshold for clustering
            maxAverageDiameter = "
        """
    
        structures = list(structures)
        # start DIANA
        # 1. initialization
        dissMatrix = computeBasePairDistanceMatrix(structures)
        c_root = Cluster()
        c_root.structures.extend([x for x in range(0, len(structures))])
        
        # start the real DIANA algorithm.
        structures_bp_sets = None
        #print(c_root.structures,structures_bp_sets)
        dataset_centroid = None
        if do_ch_first_local_min:
            structures_bp_sets = [ db_structure_to_bp_set(ss) for ss in structures ]
            dataset_centroid = centroid(c_root.structures, structures_bp_sets)
        #print("dc",dataset_centroid)
        createClusterTree(c_root, dissMatrix, maxDiameter, maxAverageDiameter, do_ch_first_local_min, None, structures_bp_sets, dataset_centroid, 1)
        # end DIANA
              
        clusters = convertTreeToListOfClusters_indices(c_root, [])
        
        #cis = convertTreeToListOfClusters_indices(c_root, structs)
        #if len(cis) > 1:
        #    ch = calinski_harabasz_index(cis, dissMatrix, len(structs))
        #    print("ch", ch)
        #else:
        #    print("ch", 0)
        return clusters
    
    @staticmethod
    def printClusters(clusters, structures_db):
        """
        print a list of lists.
        """
        if clusters == None:
            return
        cid = 0
        for c in clusters:
            cid +=1
            print("ClusterID:",cid)
            for i in c:
                s = structures_db[i]
                print(s)
    
    @staticmethod
    def printClustersMFEs(sequence, clusters, structures_db):
        if clusters == None:
            return
        cid = 0
        for c in clusters:
            cid +=1
            print("ClusterID:",cid)
            
            mfe = sys.float_info.max
            mfe_str = None
            for i in c:
                s = structures_db[i]
                energy = RNA.energy_of_struct(sequence, s)
                if energy < mfe:
                    mfe = energy
                    mfe_str = s
            print(s, "{:10.2f}".format(mfe))
    
    @staticmethod
    def printClustersCentroids(clusters, structures_db):
        structures_bp_sets = [ db_structure_to_bp_set(ss) for ss in structures_db ]
        if clusters == None:
            return
        cid = 0
        str_len = len(structures_db[0])
        for c in clusters:
            cid +=1
            print("ClusterID:",cid)
            cent = centroid(c, structures_bp_sets)
            cent_list = ["."] * str_len
            for bp in list(cent):
                x, y = sorted(list(bp))
                cent_list[x] = "("
                cent_list[y] = ")"
            print("".join(cent_list))
        

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Implememntation of the DIANA divisive clustering for RNA secondary structures.')
    parser.add_argument("-f", "--file", type=str, required=True, help="Fasta file with secondary structures")
    parser.add_argument("-d", "--diameter-threshold", type=int, default=0, required=False, help="Cluster diameter threshold (max bp distance within a cluster)")
    parser.add_argument("-m", "--average-diameter-threshold", type=float, default=0, required=False, help="Average diameter threshold")
    parser.add_argument("-c", "--calinski-harabasz-threshold", action='store_true', default=False, required=False, help="Abort if the Calinski-Harabasz-Index reaches the first local minimum")
    parser.add_argument("--mfe-output", action='store_true', default=False, required=False, help="Reduced output with mfe representatives. Requres an RNA sequence as input")
    parser.add_argument("-s", "--sequence", type=str, required=False, help="RNA sequence (ACGU)")
    parser.add_argument("--centroid-output", action='store_true', default=False, required=False, help="Reduced output with centroid representatives")
    args = parser.parse_args()
    structureFileName = args.file #sys.argv[1]
    structs = parseFile(structureFileName)
    maxDiameter = args.diameter_threshold
    maxAverageDiameter = args.average_diameter_threshold
    d = DIANA()
    clusters = d.doClustering(structs, maxDiameter, maxAverageDiameter, args.calinski_harabasz_threshold)
    if args.mfe_output:
        d.printClustersMFEs(args.sequence, clusters, structs)
    elif args.centroid_output:
        d.printClustersCentroids(clusters, structs)
    else:
        d.printClusters(clusters, structs)
    








