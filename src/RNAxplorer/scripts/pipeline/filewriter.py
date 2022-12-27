
class FileWriter:
    @staticmethod
    def writeMatrixToFile(m, fp):
        f = open(fp, "w")
        for i in range(len(m)):
            for j in range(len(m[i])):
                f.write("%10.4g" % m[i][j] + " ")            
            f.write("\n")
        f.close()
    @staticmethod    
    def writeRepresentativesToFile(seq, representatives, pathToFile):
        """
        representatives are the structures and their free energies.
        """
        f = open(pathToFile, "w")
        f.write(" " + seq + "\n")
        for i in range(len(representatives)):
            f.write(" " + str(i + 1) + " " + representatives[i][0] + " " + str(representatives[i][1]) + "\n")
        f.close()
