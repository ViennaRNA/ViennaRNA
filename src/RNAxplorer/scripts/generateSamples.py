"""
General pipeline for diverse structure sample set generation
using RNAxplorer
"""

import sys, re, RNA, subprocess
from flooding import *


def readFasta(filename):
    """
    Read a FASTA file and extract records as pairs of
    id, and sequence
    """
    res = []
    counter = 0
    is_fasta  = 0
    record_id = "seq_" + `counter`
    sequence = ""
    for l in open(filename):
        l.rstrip('\n')

        # start actual parsing
        if l.startswith(">"):
            match = re.search('>\s*(\S+)', l)
            if match:
                is_fasta  = 1
                record_id = match.group(1)
                sequence  = ""
        elif not l.startswith("#")  and not l.startswith(";"):
            match = re.search('([ACGUTacgutnN]+)', l)
            if match:
                sequence += match.group(1)

        if not is_fasta:
            res.append((record_id, sequence))
            sequence = ""
            counter += 1
            record_id = "seq_" + `counter`

    return res

class XplorerData(object):
    Sequence = ""
    Ref_struc1 = ""
    Ref_struc2 = ""
    Structures = []

    def __init__(self, seq, str1, str2, structures):
        self.Sequence = seq
        self.Ref_struc1 = str1
        self.Ref_struc2 = str2
        self.Structures = structures

def readRNAxplorerFile(filename):
    """
    Read a RNAxplorer output file and extract the references and structures.
    """
    structures = [];
    ref_struct1 = "";
    ref_struct2 = "";
    
    for l in open(filename):
        l.rstrip('\n')

        # start actual parsing
        match = re.search('(^[ACGUTacgutnN]+$)', l)
        if match:
            sequence = match.group(1)
        
        match = re.search('^([\(\)\.]+)\s+(\-?\d+\.\d+)$', l)
        if match:
            structure = match.group(1)
            if(ref_struct1 == ""):
                ref_struct1 = structure
            else:
                if (ref_struct2 == ""):
                    ref_struct2 = structure
            
        match = re.search('(^(\d+)\s+(\d+)\s+(\-?\d+\.\d+)\s+\(\d+\)\s+([\(\)\.]+)$)', l)
        if match:
            structures.append(match.group(5))
    
    result = XplorerData(sequence, ref_struct1, ref_struct2, structures)
    return result


def callRNAxplorer(seq, ref_struct1, ref_struct2, n=100):
    """
    Call RNAxplorer using ref_struct1, ref_struct2 as
    structures that are made equally probable within
    the ensemble. Tell RNAxplorer to generate n samples

    At this point, we might have redundancy in the sample
    set which needs to be removed later
    """
    structures = []

    RNAxplorer="./RNAxplorer -M SM -e \"MBSF\" --betaScale=1.2 -i "+str(n)
    #Unique=" | sort -k5 | uniq"
    sequenceAndStructures=seq+"\n"+ref_struct1+"\n"+ref_struct2
    # Run RNAxplorer
    result=subprocess.check_output("echo -e \""+sequenceAndStructures+"\" | "+RNAxplorer, shell=True)
    result = result.splitlines()
    for line in result:   
        match = re.search('(^(\d+)\s+(\d+)\s+(\-?\d+\.\d+)\s+\(\d+\)\s+([\(\)\.]+)$)', line)
        if match:
            structures.append(match.group(5))
    
    print "Generating %d samples with RNAxplorer using\n%s\n%s\n%s" % (n, seq, ref_struct1, ref_struct2)
    
    return structures


def addToClusters(clusters, structures):
    """
    Group the provided structures into clusters,
    and for each cluster compute a representative,
    e.g. the centroid, or MFE, and possibly its
    ensemble diversity and partition function,
    i.e. whatever property we might be interested in

    We can either utilize the hierarchical clustering
    method that we employed before, or, as Maria
    suggested, create gradient basins at this stage.
    """
    new_clusters = clusters

    # this loop is just a placeholder
    # we actually need to assign each structure to
    # a corresponding cluster instead
    for s in structures:
        new_clusters.append(s)

    return new_clusters


def getInterestingClusters(clusters):
    """
    Extract a list of interesting structures that
    represent the provided clusters. These structures
    are then used in the next sampling iteration to
    emphasize a particular region in the landscape.

    This could be a cluster with low sample number,
    or one that has the highest distance to all the
    others. Open for suggestions here...
    """


def get2DBasins(samples, sequence, ref1, ref2):
    """
    Given two reference structures ref1, and ref2, generate
    the 2D projection of the provided sample set of structures.
    Using the MFE representatives of each distance class,
    determine the basins of attraction in this projection,
    and associate each structure to such basin.
    Basin detection might be implemented using a flooding
    technique that stops as soon as a saddle point in the
    projection is encountered
    
    This function returns a list of lists of (s,idx) pairs, where
    s is a secondary structure in dot-bracket notation, and idx
    is the corresponding index in the samples list
    """

    twoDlandscapeData = []
    #generate k,l neighborhoods.
    for structure in samples:
        k = RNA.bp_distance(structure,ref1)
        l = RNA.bp_distance(structure,ref2)
        energy = RNA.energy_of_struct(sequence, structure)        
        twoDlandscapeData.append( [k, l, energy, structure] )

    return watershed_flooding(twoDlandscapeData)

def uniq_list(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def generateSamples(seq, mfe_struct, reference_stack):
    """
    Generate a diverse structure sample set by calling
    RNAxplorer for each pair of reference structures
    that is supplied

    After calling RNAxplorer we extract 2D basins from
    the projection and retrieve a set of new interesting
    reference structures from them.
    """

    # The global storage of all structures sampled in this run.
    sample_set = []

    # stack to hold interesting new reference structures
    new_reference_stack = []


    # now process all the structures in the reference_stack
    # to obtain more structure samples. 
    while len(reference_stack) > 0:
        # at this point, pop two structures at a time to get
        # new references and only default to using mfe_struct
        # as first reference if we have just a single new_ref left
        new_ref1 = reference_stack.pop()
        if len(reference_stack) > 0:
            new_ref2 = reference_stack.pop()
        else:
            new_ref2 = mfe_struct

        # call RNAxplorer with our reference structures
        new_samples = callRNAxplorer(seq, new_ref1, new_ref2)
        new_samples = uniq_list(new_samples)

        # determine which structures belong to the same region
        # of interest, according to current 2D projection.
        basins_2D = get2DBasins(new_samples, seq, new_ref1, new_ref2)

        # cluster the structures within each basin to
        # obtain potentially interesting new reference
        # structures
        for i, b in enumerate(basins_2D):
            print i
            for (s,idx) in b:
                print s

            #basin_clusters = getClusters(b)
            #s = getInterestingClusters(basin_clusters)
            #for c in s:
            #    new_reference_stack.append(c)

        # finally, we add the new_samples to our sample_set
        for s in new_samples:
            sample_set.append(s)

    return sample_set, new_reference_stack


def mainloop(seq, mfe_struct, mfe, iterations=0):
    """
    This is the mainloop that iteratively calls
    generateSamples() to grow the global_structures list
    """
    # the global storage of all structures sampled so far.
    # this variable is called clusters, however, we actually
    # don't require the structures to be clustered here.
    global_structures = []

    # Initial step:
    # start with mfe_struct and open chain as references
    # to retrieve first sample set
    ref_structs = []
    openchain = "." * len(seq)
    ref_structs.append(mfe_struct)
    ref_structs.append(openchain)

    while iterations >= 0 and len(ref_structs) > 0:
        # get new samples and potentially new reference structures
        samples, new_ref_structs = generateSamples(seq, mfe_struct, ref_structs)
        # add new samples to set of global structues
        global_structures = addToClusters(global_structures, samples)
        # assign new reference structures for next iteration
        ref_structs = new_ref_structs
        iterations = iterations - 1

    return global_structures


if __name__ == "__main__":
    inputfile = sys.argv[1]
    #records = readFasta(inputfile)
    
    #records = readRNAxplorerFile(inputfile)
    seq = "GGGAAUUAUUGUUCCCUGAGAGCGGUAGUUCUC"
    (mfe_struct, mfe) = RNA.fold(seq)
    print mfe
    print mfe_struct
    sss = mainloop(seq, mfe_struct, mfe, 1)

    print sss

    #ref1 = "................................."
    #ref2 = "((((((((((((((.....))))))))))))))"
    #samples = callRNAxplorer(seq, ref1, ref2, 10)

    #result = XplorerData(seq, ref1, ref2, samples)
    #xplorerData = readRNAxplorerFile(inputfile)
    #seq = xplorerData.Sequence
    #ref1 = xplorerData.Ref_struc1
    #ref2 = xplorerData.Ref_struc2
    #samples = xplorerData.Structures
    #basins = get2DBasins(samples, seq, ref1, ref2)
    #print basins
    
    
"""
    for fasta_id, seq in records:
        print "...processing new Record:\n\nID: %s\n%s" % (fasta_id, seq)

        # compute MFE and corresponding structure
        mfe_struct, mfe = RNA.fold(seq)
        print "%s (%6.2f)" % (mfe_struct, mfe)

        # generate diverse sample set
        clusters = generateSamples(seq,mfe_struct,mfe)


        print "\n...done\n"
"""
