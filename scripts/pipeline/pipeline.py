#!/usr/bin/env python

import RNA
import samplegenerator
from clusteralgorithms.diana import DIANA
from rates import Rates
import filewriter
from filereader import FastaFileReader
import os
import operator
import sys
import math
import argparse
from parameters import Parameters, TwoDSamplingParameters

def mfeStructure(seq, cluster):
    mfe = sys.float_info.max
    mfe_s = ""
    for s in cluster:
        e = RNA.energy_of_struct(seq, s)
        if e < mfe:
            mfe = e
            mfe_s = s
    return mfe_s

def selectRepresentatives(seq, clusters):
    """
    extract mfe structures from clusters
    @param seq - string - the RNA sequence.
    @return list - structures in dot-bracket notation
    """
    representatives = []
    for c in clusters:
        # select mfe or TODO: centroid structure.
        s = mfeStructure(seq, c)
        representatives.append(s)
    return representatives
    

if __name__ == "__main__":
   
    parser = argparse.ArgumentParser(description='Script to execute the iterative 2D-sampling and clustering of the RNA folding landscape.', \
                                     usage='At first Treekin and RNAxplorer have to be installed. After that you can run "python %(prog)s -I <input_file> [options]"')
    parser.add_argument("-s", "--sequence", action="store", type=str, default="", required=False, help="The RNA sequence (overwrites -I).")
    parser.add_argument("-r1", "--ref_one", action="store", type=str, default="", required=False, help="The first reference structure (overwrites -I). Default: mfe")
    parser.add_argument("-r2", "--ref_two", action="store", type=str, default="", required=False, help="The second reference structure (overwrites -I).Default: open chain")
    parser.add_argument("-I", "--infile", action="store", type=str, required=False, help="Input file with a sequence in FASTA format and to initial reference structures in dot-bracket format.")
    
    parser.add_argument("--samplingIterations", action="store", type=int, default=5, required=False, help="2D-sampling iterations.")
    parser.add_argument("--maxRNAxplorerSamples", action="store", type=int, default=500, required=False, help="The maximal number of samples produced by RNAxplorer")
    
    parser.add_argument("--maxDiameterThreshold", action="store", type=int, default=0, required=False, help="Maximal allowed diameter for hierarchical clustering.")
    parser.add_argument("--maxAverageDiameterThreshold", action="store", type=int, default=6, required=False, help="Maximal allowed average diameter for hierarchical clustering.")
    
    parser.add_argument("--treekinStart", action="store", type=float, default=0.001, required=False, help="Treekin start time.")
    parser.add_argument("--treekinEnd", action="store", type=float, default=1e20, required=False, help="Treekin end time.")
    
    parser.add_argument("--ratesFile", action="store", type=str, default="./ratesFile.txt", required=False, help="Path to the file where the rate matrix will be stored.")
    parser.add_argument("--structuresFile", action="store", type=str, default="./structuresFile.txt", required=False, help="Path to the file where the structures will be stored.")
    parser.add_argument("--kineticFile", action="store", type=str, default="./treekin-kinetic.txt", required=False, help="Path to the file where the treekin output will be stored.")
    # TODO:   read treekin start structure (mfe or oc) or initial distribution.
    
    
    args = parser.parse_args()
    #set minimal required parameters
    seq = ""
    ref_one = ""
    ref_two = ""
    if args.infile:
        seq, ref_one , ref_two = FastaFileReader.readFile(args.infile)
    if args.sequence:
        seq = args.sequence
    if not seq:
        raise Exception("Error: no sequence given as input!")
    if args.ref_one:
        ref_one = args.ref_one
    if args.ref_two:
        ref_two = args.ref_two
    if not ref_one:
        (mfe_struct, mfe) = RNA.fold(seq)
        ref_one = mfe_struct
    if not ref_two:
        openchain = "." * len(seq)
        ref_two = openchain
    
    param = Parameters()
    twoDSamplingParameters = TwoDSamplingParameters()
    twoDSamplingParameters.Sequence = seq
    #2D sampling
    twoDSamplingParameters.Reference_one = ref_one
    twoDSamplingParameters.Reference_two = ref_two
    twoDSamplingParameters.MaxXplorerSamples = args.maxRNAxplorerSamples
    twoDSamplingParameters.SamplingIterations = args.samplingIterations
    param.TwoDSamplingParam = twoDSamplingParameters
    # input treekin
    param.StartTime = args.treekinStart
    param.EndTime = args.treekinEnd
    # output
    param.RatesFilePath = args.ratesFile
    param.StructuresFilePath = args.structuresFile
    param.KineticFile = args.kineticFile
    #Diana Threshold (max basepair difference in cluster)
    maxDiameterThreshold = args.maxDiameterThreshold
    maxAverageDiameterThreshold = args.maxAverageDiameterThreshold
    
    # call samplegenerator which calls RNAxplorer (in matrix2d objects) and executes a clustering (2D or DIANA, etc.) only to select 
    # new structure representatives for samplegeneration.
    structures = samplegenerator.mainloop(twoDSamplingParameters)

    # call final clusteralg.
    clusters = DIANA.doClustering(structures, maxDiameterThreshold, maxAverageDiameterThreshold)
    DIANA.printClusters(clusters)
    representatives = selectRepresentatives(seq, clusters)
    
    # prepare representatives and sort them, to make the comparison of the output easier.
    repsAndEnergies = [ (x, RNA.energy_of_struct(seq, x)) for x in representatives]
    repsAndEnergies.sort(key=operator.itemgetter(1, 0))
    
    # rate computation.
    rateMatrix = Rates.arrheniusRates(seq, repsAndEnergies)
    
    # write files for treekin (matrix and structures)
    fw = filewriter.FileWriter()
    fw.writeMatrixToFile(rateMatrix, param.RatesFilePath)
    fw.writeRepresentativesToFile(seq, repsAndEnergies, param.StructuresFilePath)
    
    # run treekin
    to_call = "cat " + param.StructuresFilePath + " | treekin --p0 " + str(param.StartStructureID) + "=1 --ratesfile " + param.RatesFilePath + " --t0 " + \
            str(param.StartTime) + " --t8 " + str(param.EndTime) + " -m I > " + param.KineticFile
    os.system(to_call)
    
    
    
    
