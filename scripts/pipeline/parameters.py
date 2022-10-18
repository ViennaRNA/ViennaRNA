
class TwoDSamplingParameters:
    def __init__(self):
        self.Sequence = ""
        self.Reference_one = ""
        self.Reference_two = ""
        
        self.MaxXplorerSamples = 100 
        #iterations for calls of RNAxplorer
        self.SamplingIterations = 0

class Parameters:
    def __init__(self):
        self.TwoDSamplingParam = TwoDSamplingParameters()
        
        # input treekin
        self.StartStructureID = 1
        self.StartTime = "0.001"
        self.EndTime = "1e20"    
        
        # output
        self.RatesFilePath = "./ratesFile.txt"
        self.StructuresFilePath = "./structuresFile.txt"
        self.KineticFile = "./treekin-kinetic.txt"
