import re
from samplegenerator import Matrix2D

class FastaFileReader:
    @staticmethod
    def readFile(filename):
        """
        Read a FASTA file and extract the sequence and two reference structures
        """
        sequence = ""
        ref_one = ""
        ref_two = ""
        for l in open(filename):
            l.rstrip('\n')
            # start actual parsing
            match = re.search('([ACGUTacgutnN]+)', l)
            if match:
                sequence += match.group(1)
            match = re.search('([\.\)\(]+)', l)
            if match:
                if not ref_one:
                    ref_one = match.group(1)
                elif not ref_two:
                    ref_two = match.group(1)
        return sequence, ref_one, ref_two

class xplorerFileReader:
    @staticmethod
    def readFile(filename):
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
        
        result = Matrix2D(sequence, structures, ref_struct1, ref_struct2)
        return result
