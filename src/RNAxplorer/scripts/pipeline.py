from tinyRNAClustering import *
import sys,os

RNAFOLD_SAMPLING = "rnasubopt -p 10000 -s <"
R_EXEC = "RScript"

def unique(structs):
    res = {}
    for s,ss in structs:
        res[s] = ss
    return [(x,res[x]) for x in res]

def generateStructs(seq):
    outpath = "structs.faa"
    outfile = open("tmp.faa","w")
    outfile.write(seq+"\n")
    outfile.close()
    os.system(RNAFOLD_SAMPLING+ " tmp.faa > "+outpath)
    structs = sorted(unique(parseFile(outpath)))
    print "\n".join([a for a,b in structs])
    return structs

def computeDisMatrix(structs):
    return [[computeDistance(ss1,ss2) for j,(s2,ss2) in enumerate(structs)] for i,(s1,ss1) in enumerate(structs)]
    
def exportDisMatrix(m,fp):
    f = open(fp,"w")
    for i in range(len(m)):
        for j in range(len(m[i])):
            if j!=0:
                f.write(",") 
            f.write("%s"%m[i][j])            
        f.write("\n")
    f.close()

def runDIANA(matPath):
    os.system(R_EXEC+ " ..\R\pw_distance.R %s > dummy"%matPath)
    return matPath+".clusters.csv"

def createClusters(seq):
    disMatPath = "distMat.dat"
    structs = generateStructs(seq)
    #print structs
    mat = computeDisMatrix(structs)
    exportDisMatrix(mat,disMatPath)
    clustersPath = runDIANA(disMatPath)
    doClustering(structs,clustersPath)
    



class Matrix2D:
    def __init__(self,structs,ss1,ss2):
        self.ss1 = ss1
        self.ss2 = ss2
        self.structs = structs
        self.doProjection()

    def doProjection(self):
        # [Stub] Something magical happens, and then...
        pass

    def estimateEnergyBarrier(self, s1, s2):
        # Estimates the energy between s1 and s2 based on 2D projection
        return 0.
    
    def createClusters(self):
        res = [] # Sequence of StructureSet objects
        # [Stub] Something magical happens, and then...
        return res
        

class KineticLandscape:

    def __init__(self, seq):
        self.rna = seq
        self.allStructs = StructureSet()
        self.energyBarriers = {}
        self.states = {}

    def sampleDistortedLandscape(self, ss1, ss2):
        """
            Attempt to generate samples along the path from ss1 to ss2
        """
        res = StructureSet()
        # [Stub] Something magical happens, and then...
        return res
        
    def buildLandscape(self):
        # Start from MFE and unfolded state
        mfe = getMFEFolding(self.seq)
        unfolded = getUnfoldedStructure(self.seq)
        
        # Initialize Clusters and interaction stacks
        clusters = {0:StructureSet([unfolded]), 1:StructureSet([mfe])}
        barriers2D = {}
        interactionStack = [(0,1)]
        
        # Construct graph
        while len(interactionsStack)>0:
            (i,j) = interactionsStack.pop()
            c1,c2 = clusters[i],clusters[j]
            ss1,ss2 = c1.getCentroid(),c2.getCentroid()
            
            sample = self.sampleDistortedLandscape(ss1,ss2)
            self.allStructs = self.allStructs + sample
            
            m2D = Matrix2D(self.allStructs,ss1,ss2)
            clusters2D = m2D.doClustering2D()

            # First we cluster based on isocurves in 2D projection
            for c2D in clusters2D:
                # Then we cluster the structures (which may correspond to different "depth")
                for c in c2D.doClustering():                
                    # Discard "garbage" clusters
                    if c.isSufficientlyTight():
                        # Try to merge current cluster with existing ones
                        target = self.mergingCandidate(clusters)
                        if target is not None:
                            # There is such a cluster, add current cluster to it
                            target.selectiveAdd(c)
                        else:
                            # Otherwise add cluster to list of cluster, and schedule
                            # the exploration of other interactions.
                            x = len(clusters)
                            clusters[i] = c
                            interactionStack += [(x,y) for y in clusters if y != x]

        # Now we've got a complete graph, build approximate energy barriers
        for i in range(len(clusters)):
            for j in range(len(clusters)):
                if i != j:
                    c1,c2 = clusters[i],clusters[j]
                    ss1,ss2 = c1.getCentroid(),c2.getCentroid()
                    m2D = self.allStructs.project(ss1,ss2)
                    barriers[(i,j)] = m2D.estimateEnergyBarrier(ss1,ss2)
        
        
        
class StructureSet:
    
    def __init__(self, structs = []):
        self.structs = set(structs)
        
    def project(self,ss1,ss2):
        res = Matrix2D()
        # [Stub] Something magical happens, and then...
        return res
    
    def doClustering(self):
        res = [] # Sequence of StructureSet objects
        # [Stub] Something magical happens, and then...
        return res

    def getCentroid(self):
        return centroid(self.structs)

    def hasStructure(self,ss):
        return ss in self.structs

    def isSufficientlyTight(self):
        
        return True


if __name__=="__main__":
    createClusters(sys.argv[1])
