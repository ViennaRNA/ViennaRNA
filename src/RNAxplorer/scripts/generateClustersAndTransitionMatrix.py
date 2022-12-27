"""
1. Generate clusters with the DIANA algorithm
2. Compute a distance matrix: use centroids or mfe-structure and use RNAxplorer to compute the paths between all cluster-representatives

"""

import sys, re, RNA, subprocess, os, operator
import numpy as np
from generateSamples import *
from tinyRNAClustering import *

def exportDisMatrix(m, fp):
    f = open(fp, "w")
    for i in range(len(m)):
        for j in range(len(m[i])):
            if j != 0:
                f.write(",")
            f.write("%s" % m[i][j])
        
        f.write("\n")

def computeDisMatrix(structs):
    Nstructs = len(structs)
    matrix = np.zeros((Nstructs, Nstructs))
    
    for i in range(Nstructs):
        for j in range (i+1, Nstructs):
            dist = RNA.bp_distance(structs[i], structs[j])
            matrix[i][j] = dist
            matrix[j][i] = dist
    
    return matrix

def runDIANA(matPath):
    os.system("Rscript ../R/pw_distance.R %s > dummy"%matPath)
    return matPath+".clusters.csv"


def generateClusters(structures):
    disMatPath = "distMat.dat"
    mat = computeDisMatrix(structures)
    exportDisMatrix(mat,disMatPath)
    clustersPath = runDIANA(disMatPath)
    fclusters = "distMat.dat.clusters.csv"
    structs = []
    for ss in structures:
        structs.append((ss.strip(),parseSS(ss.strip())))
    #returns the best hierarchy of clusters.
    cl = doClustering(structs,fclusters)
    # convert back.
    n = sum([len(s) for s,ss in structs])/len(structs)
    clusters = []
    for num in cl:
        cluster = []
        for ss in cl[num]:
            s = convertBack(ss,n)
            cluster.append(s)
        clusters.append(cluster)
    return clusters

def extractClusterRepresentatives(seq, clusters):
    """
    returns the minima and their energies.
    """
    representatives = []
    for c in clusters:
        minE = sys.maxint
        minS = ""
        for s in c:
            e = RNA.energy_of_struct(seq, s)
            if e < minE:
                minE = e
                minS = s
        representatives.append((minS, minE))
    return representatives 

class Rates:
    #gas constant in kcal/(mol*K)
    K = 0.00198717
    #temperature (37 Celsius in Kelvin)
    T = 273.15 + 37
    
    def arrheniusRates(self, sequence, structures):
        Nstructs = len(structures)
        matrix = np.zeros((Nstructs, Nstructs))
        for i in range(Nstructs):
            for j in range (i+1, Nstructs):
                s1, e1 = structures[i]
                s2, e2 = structures[j]
                #Method 2 (using RNAxplorer):
                #call RNAxplorer to compute the path and extract the barrier
                to_call = "echo -e \""+sequence+"\n"+s1+"\n"+s2+"\" | RNAxplorer -M GW | grep -F 'barrier' | grep -o -P '\-?\d*\.\d*'"
                barrier = subprocess.check_output(to_call, shell=True)
                barrier = float(barrier)
                matrix[i][j] = math.exp(-(barrier-e1)/(self.K*self.T))
                matrix[j][i] = math.exp(-(barrier-e2)/(self.K*self.T))
        return matrix

def writeMatrixToFile(m,fp):
    f = open(fp,"w")
    for i in range(len(m)):
        for j in range(len(m[i])):
            f.write("%10.4g"%m[i][j]+" ")            
        f.write("\n")
    f.close()

def writeRepresentativesToFile(seq, representatives, pathToFile):
    """
    representatives are the structures and their free energies.
    """
    f = open(pathToFile,"w")
    f.write(" "+seq+"\n")
    for i in range(len(representatives)):
        f.write(" "+str(i+1)+" "+representatives[i][0]+" "+str(representatives[i][1])+"\n")
    f.close()

if __name__ == "__main__":
    #inputfile = sys.argv[1]
    #records = readFasta(inputfile)
    
    #records = readRNAxplorerFile(inputfile)
    seq = "GGGAAUUAUUGUUCCCUGAGAGCGGUAGUUCUC"
    (mfe_struct, mfe) = RNA.fold(seq)
    print mfe
    print mfe_struct
    sss = mainloop(seq, mfe_struct, mfe, 1)
    sss = uniq_list(sss)
    
    #cluster secondary-structure-set
    clusters = generateClusters(sss)   
    #extract cluster representatives (could be centroid or mfe; here we chose mfe)
    representatives = extractClusterRepresentatives(seq, clusters)
    
    #openchain = "." * len(seq)
    #energy = RNA.energy_of_struct(seq,openchain)
    #representatives.insert(0, (openchain,energy))
    
    #sort according to the free energy and if it is equal, lexographically.
    representatives.sort(key=operator.itemgetter(1,0))
    #approximate the rates. An alternative is to compute paths between representatives, extract the barrier and compute rates
    r = Rates()
    rateMatrix = r.arrheniusRates(seq, representatives)
    
    structuresFile = "structureRepresentatives.txt"
    writeRepresentativesToFile(seq, representatives, structuresFile)
    ratesFile = "rateMatrix.txt"
    writeMatrixToFile(rateMatrix,ratesFile)
    
    #call treekin
    kineticFile = "treekin-kinetic.txt"
    startStructureID = 1 #TODO: find the "right" start structure.
    startTime = "0.001"
    endTime = "1e20"
    to_call = "cat " + structuresFile + " | treekin --p0 " + str(startStructureID) + "=1 --ratesfile " + ratesFile + " --t0 " + startTime + " --t8 " + endTime + " -m I > "+kineticFile
    os.system(to_call)
    
    
    #visualize the kinetic
    #--> maybe in d3 js.
    
    
    
    
    
    
    
    
    
    
    
    
    