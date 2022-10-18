"""
Script for 2D flooding algorithms and helper functions.
"""

import operator  # for sort with 2 criteria.
import sys
import math

"""
Objects of the RNADfoldColumn class represent a column in a
2D projection with respect to some reference structures

In essence, this just collects a number of secondary structures
in dot-bracket notation, and stores its mfe representative. Thus,
it can be used for any projection, no matter the number of dimensions
"""
class RNA2DfoldColumn(object):

    def __init__(self, structure, e, idx):
        self.structures = []
        self.mfe = 10000.0
        self.mfe_struct = ""

        self.addStructure(structure, e, idx)

    def addStructure(self, s, e, idx):
        self.structures.append((s, idx))
        if self.mfe > e:
            self.mfe = e
            self.mfe_struct = s

class WatershedFlooder:
    @staticmethod
    def doflooding(landscapeData):
        """"
        @param landscapeData =[ [2,5,-3,"...)"], [2,3,-3,".).)"],...] ref1 dist, ref2 dist, free energy, structure
        @return list of lists - clustered structures
        """
        
        if len(landscapeData) <= 0 or len(landscapeData[0]) < 4:
            print "Error: landscapeData has the wrong structure."
            return []
    
        # put data into a 2D grid (actually a dict with (k,l) keys)
        grid2D = {}
        for i in range(0, len(landscapeData)):
            k = landscapeData[i][0]
            l = landscapeData[i][1]
            e = landscapeData[i][2]
            s = landscapeData[i][3]
            if (k, l) in grid2D:
                grid2D[(k, l)].addStructure(s, e, i)
            else:
                grid2D[(k, l)] = RNA2DfoldColumn(s, e, i)
    
        # create a to-do list of columns to process
        todo = []
        for k in grid2D.keys():
            e = grid2D[k].mfe
            todo.append((k, e))
    
        # sort the to-do list by free energy
        todo.sort(key=operator.itemgetter(1))
    
        # list of already processed columns
        done = []
    
        basins = []
    
        for i in range(0, len(todo)):
            basin = []
            todo_stack = []
    
            ((k, l), e) = todo[i]
            
            todo_stack.append((k, l))
    
            while len(todo_stack) > 0:
                # cell is a (k,l) tuple
                cell = todo_stack.pop()
    
                if not cell in done:
                    basin.append(cell)
    
                    cell_e = grid2D[cell].mfe
    
                    # determine neighbor grid cells of current cell
                    neighbors = []
                    nn = [(k - 1, l - 1), (k + 1, l - 1), (k - 1, l + 1), (k + 1, l + 1)]
                    for n in nn:
                        if n in grid2D:
                            neighbors.append(n)
    
                    for n in neighbors:
                        n_e = grid2D[n].mfe
                        if n_e >= cell_e and n not in basin:
                            todo_stack.append(n)
    
            if len(basin) > 0:
                basins.append(basin)
    
            for b in basin:
                done.append(b)
    
        # now collect all the structures again
        bbb = []
    
        def combineBasins(b):
            c = []
            for n in b:
                map(lambda v: c.append(v), grid2D[n].structures)
            bbb.append(c)
    
        map(combineBasins, basins)
    
    
        return bbb

    @staticmethod
    def doFloodingAndExtractInterestingClusters(landscapeData):
        """"
        2D flooding that returns both, 2D clusters and interesting clusters which are close to the diagonal.
        
        @param landscapeData =[ [2,5,-3,"...)"], [2,3,-3,".).)"],...] ref1 dist, ref2 dist, free energy, structure
        @return array with clusters
        """
        if len(landscapeData) <= 0 or len(landscapeData[0]) < 4:
            print "Error: landscapeData has the wrong structure."
            return []
    
        # put data into a 2D grid (actually a dict with (k,l) keys)
        grid2D = {}
        for i in range(0, len(landscapeData)):
            k = landscapeData[i][0]
            l = landscapeData[i][1]
            e = landscapeData[i][2]
            s = landscapeData[i][3]
            if (k, l) in grid2D:
                grid2D[(k, l)].addStructure(s, e, i)
            else:
                grid2D[(k, l)] = RNA2DfoldColumn(s, e, i)
    
        # create a to-do list of columns to process
        todo = []
        for k in grid2D.keys():
            e = grid2D[k].mfe
            todo.append((k, e))
    
        # sort the to-do list by free energy
        todo.sort(key=operator.itemgetter(1))
    
        # list of already processed columns
        done = []
    
        basins = []
    
        for i in range(0, len(todo)):
            basin = []
            todo_stack = []
    
            ((k, l), e) = todo[i]
            
            todo_stack.append((k, l))
    
            while len(todo_stack) > 0:
                # cell is a (k,l) tuple
                cell = todo_stack.pop()
    
                if not cell in done:
                    basin.append(cell)
    
                    cell_e = grid2D[cell].mfe
    
                    # determine neighbor grid cells of current cell
                    neighbors = []
                    nn = [(k - 1, l - 1), (k + 1, l - 1), (k - 1, l + 1), (k + 1, l + 1)]
                    for n in nn:
                        if n in grid2D:
                            neighbors.append(n)
    
                    for n in neighbors:
                        n_e = grid2D[n].mfe
                        if n_e >= cell_e and n not in basin:
                            todo_stack.append(n)
    
            if len(basin) > 0:
                basins.append(basin)
    
            for b in basin:
                done.append(b)
    
        """
        select interesting basins that are close to the diagonal
        and also closest to the origin and farthest away from the origin
        """
        
        
        def minCellDistanceFromDiagonal(basin):
            minDist = sys.maxint
            for cell in basin:
                dist = math.fabs(cell[0] - cell[1]) 
                if dist < minDist:
                    minDist = dist
            return  minDist
        
        def closestToDiagonalComparison(basin1, basin2): 
            dist1 = minCellDistanceFromDiagonal(basin1)
            dist2 = minCellDistanceFromDiagonal(basin2)
            if dist1 < dist2:
                return -1
            elif dist1 > dist2:
                return 1
            else:
                return 0
        
        sortedBasinsCTD = sorted(basins, cmp=closestToDiagonalComparison)
        
        def cellClosestToOrigin(basin):
            minDist = sys.float_info.max
            for cell in basin:
                dist = math.sqrt(cell[0]*cell[0]+cell[1]*cell[1])
                if dist < minDist:
                    minDist = dist
            return minDist
        
        def closestToOriginComparison(basin1, basin2):
            dist1 = cellClosestToOrigin(basin1)
            dist2 = cellClosestToOrigin(basin2)
            if dist1 < dist2:
                return -1
            elif dist1 > dist2:
                return 1
            else:
                return 0
        
        sortedBasinsCTO = sorted(sortedBasinsCTD, cmp=closestToOriginComparison)
        firstBasin = sortedBasinsCTO[0]
        
        def cellFarthestFromOrigin(basin):
            maxDist = 0.0
            for cell in basin:
                dist = math.sqrt(cell[0]*cell[0]+cell[1]*cell[1])
                if dist > maxDist:
                    maxDist = dist
            return maxDist
        
        def farthestFromOrigin(basin1, basin2):
            dist1 = cellFarthestFromOrigin(basin1)
            dist2 = cellFarthestFromOrigin(basin2)
            if dist1 > dist2:
                return -1
            elif dist1 < dist2:
                return 1
            else:
                return 0
        
        sortedBasinsFFO = sorted(sortedBasinsCTD, cmp=farthestFromOrigin)
        secondBasin = sortedBasinsFFO[0]
        #sortedBasinsFFO.reverse()
        #sortedBasinsCTD.reverse()
        #sortedBasinsFFO = sortedBasinsCTD 
        #sortedBasinsCTO.sort(cmp=closestToDiagonalComparison)
        #sortedBasinsFFO.sort(cmp=closestToDiagonalComparison)
        
        mostInterestingBasins = [firstBasin, secondBasin]
        
        """
        """
        
        
    
        # now collect all the structures again
        bbb = []
    
        def combineBasins(b):
            c = []
            for n in b:
                map(lambda v: c.append(v), grid2D[n].structures)
            bbb.append(c)
    
        #map(combineBasins, basins)
        map(combineBasins, mostInterestingBasins)
    
        return bbb
    
    
    
    
    
    