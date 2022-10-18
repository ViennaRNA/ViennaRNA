import numpy as np
import subprocess
import math
import select

import RNAxplorer
import RNA

class Rates:
    # gas constant in kcal/(mol*K)
    K = 0.00198717
    # temperature (37 Celsius in Kelvin)
    T = 273.15 + 37
    
    @staticmethod
    def computeBarriers(sequence, structsAndEneries):
        """
        compute all barriers in one process to save time
        """
        md = RNA.md()
        md.circ     = 0
        md.uniq_ML  = 1 # in case we need M1 arrays
        md.compute_bpp = 0
        md.betaScale = 1
        
        n = len(structsAndEneries)
        barriersMatrix = np.zeros((n, n))
        for i in range(n):
            s1 = structsAndEneries[i][0]
            for j in range(i + 1, n):
                s2 = structsAndEneries[j][0]
                #barrierEnergy = RNAxplorer.barrier_estimate_2D(sequence, md, s1, s2, 10, 10)
                maxKeep = 10
                foldingPath = RNA.get_path(sequence, s1, s2, maxKeep);
                saddlePointPathType = RNAxplorer.getSaddlePoint(foldingPath);
                spList = RNAxplorer.convertVRNA_PATH_toList(saddlePointPathType)
                barrierEnergy = spList[0][1] 
                #print barrierEnergy
                barriersMatrix[i][j] = barrierEnergy
                barriersMatrix[j][i] = barrierEnergy
        return barriersMatrix

    
    @staticmethod
    def computeBarriersOld(sequence, structsAndEneries):
        """
        compute all barriers in one process to save time
        """
        cmd = subprocess.Popen(['bash'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        poll = select.poll()
        poll.register(cmd.stdout.fileno(), select.POLLIN)

        n = len(structsAndEneries)
        barriersMatrix = np.zeros((n, n))
        for i in range(n):
            s1 = structsAndEneries[i][0]
            for j in range(i + 1, n):
                s2 = structsAndEneries[j][0]
                to_call = "echo -e \"" + sequence + "\n" + s1 + "\n" + s2 + "\" | RNAxplorer -M GW | grep -F 'barrier' | grep -o -P '\-?\d*\.\d*'"
                cmd.stdin.write(to_call)
                cmd.stdin.write('\n')
                cmd.stdin.flush()  # Must include this to ensure data is passed to child process
                ready = poll.poll()
                if ready:
                    result = cmd.stdout.readline()
                    result = result.strip()
                    result = float(result)
                    barriersMatrix[i][j] = result
                    barriersMatrix[j][i] = result
        cmd.kill()
        return barriersMatrix
    
    @staticmethod
    def arrheniusRates(sequence, structsAndEneries):
        barriers = Rates.computeBarriers(sequence, structsAndEneries)
        
        Nstructs = len(structsAndEneries)
        matrix = np.zeros((Nstructs, Nstructs))
        for i in range(Nstructs):
            s1, e1 = structsAndEneries[i]
            for j in range (i + 1, Nstructs):            
                s2, e2 = structsAndEneries[j]
                barrier = barriers[i][j]
                matrix[i][j] = math.exp(-(barrier - e1) / (Rates.K * Rates.T))
                matrix[j][i] = math.exp(-(barrier - e2) / (Rates.K * Rates.T))
        return matrix

                       
                       
