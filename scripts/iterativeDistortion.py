#!/usr/bin/python

"""
Generate structures with iterative distortions.
the formula for changing the distortion parameter x (for the first reference) is:
data->x += relaxFactor * orig_x / maxSteps;
"""

import RNA
import RNAxplorer
import re
import argparse

parser = argparse.ArgumentParser(description='Iterative 2-D sampling with different distortion parameters.\n'\
                                 'The formula data->x += relaxFactor * orig_x / maxSteps is computed from 0 to maxSteps.', \
                                     usage='run "python %(prog)s -I inputfile [options]"')
#(vc, grid, relaxFactor, relax, shift, shift_to_first, verbose, maxIterations, maxSteps)
parser.add_argument("-f", "--inputFile", action="store", type=str, default="", required=True, help="File with structure\n reference one\n reference two")
parser.add_argument("-B", "--shiftBoth", action="store_true", required=False, help="Shift iterative to first and then to second reference.")
parser.add_argument("-r", "--relax", action="store_true", required=False, help="Changes the sign (increment or decrement).")
parser.add_argument("-rf", "--relaxfactor", action="store", type=float, default=1.0, required=False, help="Multiply with the given relaxfactor.")
parser.add_argument("-s", "--shift", action="store_true", required=False, help="Shift.")
parser.add_argument("-sf", "--shiftToFirst", action="store_true", required=False, help="Shift to first if True, else to second.")
parser.add_argument("-v", "--verbose", action="store_true", required=False, help="More output.")
parser.add_argument("-i", "--maxIterations", action="store", type=int, default=1, required=False, help="Max number of samples.")
parser.add_argument("-ms", "--maxSteps", action="store", type=int, default=10, required=False, help="Max steps for iterative distortion.")
parser.add_argument("-b", "--betaScale", action="store", type=float, default=1, required=False, help="Factor for scaling the temperature.")
args = parser.parse_args()

def mapGridToStructureList(gridlist):
    """
    write k,l for each structure.
    """
    structures = {}
    for entry in gridlist:
        for s in entry[3]:
            # k, l, energy, structure
            structures[s] = (entry[0], entry[1], entry[2], s)
    return structures

def createFoldCompound(seq):
    """
    RNA.cvar.uniq_ML =1 #use local var if available
    md = RNA.md("global")
    md = RNA.md()
    md.circ     = 0
    #md.uniq_ML  = 1 #TODO: use this if available
    # in case we need M1 arrays
    md.compute_bpp = 0
    md.betaScale = args.betaScale
    """
    md = RNAxplorer.createModelDetails(0, 1, 0, args.betaScale)
    vc = RNA.fold_compound(seq, md, RNA.OPTION_MFE | RNA.OPTION_PF)
    return vc

if __name__ == "__main__":
    
    seq = ""
    s1=''
    s2=''
    
    file = open(args.inputFile,"r")
    for l in file:
        l.rstrip('\n')
        # start actual parsing
        match = re.search('([ACGUTacgutnN]+)', l)
        if match:
            seq = match.group(1)
        match = re.search('([\.\)\(]+)', l)
        if match:
            if not s1:
                s1 = match.group(1)
            elif not s2:
                s2 = match.group(1)
    file.close()
    
    print seq
    print s1, RNA.energy_of_struct(seq, s1)
    print s2, RNA.energy_of_struct(seq, s2)
    
    relaxFactor = args.relaxfactor
    relax = int(args.relax)
    both = int(args.shiftBoth)
    shift = int(args.shift)
    shift_to_first = int(args.shiftToFirst)
    verbose = int(args.verbose)
    maxSteps = args.maxSteps
    maxIterations = args.maxIterations
    
    
    vc = createFoldCompound(seq)
    grid = RNAxplorer.initLandscape(seq, s1, s2)
    distortion_x, distortion_y = RNAxplorer.computeInitialDistortion(vc, s1, s2)
    print "distortion", distortion_x, distortion_y
    
    #rescale
    mfe_struct, mmfe = RNA.fold(seq) #TODO: use extended fold with modelDetails in release
    bp_dist_mfe_ref1 = RNA.bp_distance(mfe_struct, s1)
    bp_dist_mfe_ref2 = RNA.bp_distance(mfe_struct, s2)
    rescale = mmfe + (bp_dist_mfe_ref1 * distortion_x) + (bp_dist_mfe_ref2 * distortion_y)
    RNAxplorer.rescaleEnergy(vc, rescale) #TODO: use RNAlib rescale if available
    
    RNAxplorer.addSoftconstraints(vc, s1, s2, distortion_x, distortion_y)

    print relaxFactor,relax,shift,shift_to_first,verbose,maxIterations,maxSteps
    if both:
        RNAxplorer.fillGridStepwiseBothRef(vc, grid, relaxFactor, relax, shift, shift_to_first, verbose, maxIterations, maxSteps)
        RNAxplorer.printLandscape(grid, vc)
        gridlist = RNAxplorer.convertGrid_toList(grid)
        structs = mapGridToStructureList(gridlist)
        print "both",len(structs.values())
        #print structs.keys()
    elif shift_to_first:
        grid = RNAxplorer.initLandscape(seq, s1, s2)
        RNAxplorer.fillGridStepwiseFirstRef(vc, grid, relaxFactor, relax, verbose, maxIterations, maxSteps)
        RNAxplorer.printLandscape(grid, vc)
        gridlist = RNAxplorer.convertGrid_toList(grid)
        structs = mapGridToStructureList(gridlist)
        print "first ref",s1,len(structs.values())
        #print structs.keys()
    elif shift:
        grid = RNAxplorer.initLandscape(seq, s1, s2)
        RNAxplorer.fillGridStepwiseSecondRef(vc, grid, relaxFactor, relax, verbose, maxIterations, maxSteps)
        RNAxplorer.printLandscape(grid, vc)
        gridlist = RNAxplorer.convertGrid_toList(grid)
        structs = mapGridToStructureList(gridlist)
        print "second ref",s2,len(structs.values())
    else:
        # shift to first and to second
        grid = RNAxplorer.initLandscape(seq, s1, s2)
        RNAxplorer.fillGridStepwiseFirstRef(vc, grid, relaxFactor, relax, verbose, maxIterations, maxSteps)
        RNAxplorer.fillGridStepwiseSecondRef(vc, grid, relaxFactor, relax, verbose, maxIterations, maxSteps)
        RNAxplorer.printLandscape(grid, vc)
        gridlist = RNAxplorer.convertGrid_toList(grid)
        structs = mapGridToStructureList(gridlist)
        print "first and second ref",s2,len(structs.values())
        #print structs.keys()


    




