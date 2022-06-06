#!/usr/bin/python3

import math as math
from viewshedMath import viewshedMath as vsm

# init
maxRow = maxCol = 5
srcRow = srcCol = maxRow // 2

# "direction" : [rowInc, colInc]
dirIncEdge = {"N" : [-1,0], "NE" : [-1,1],
          "E" : [0,1],  "SE" : [1,1],
          "S" : [1,0],  "SW" : [1,-1],
          "W" : [0,-1], "NW" : [-1,-1]}


# "direction" : [xi, xj, yi, yj]
# i is row, col is j
dirIncSlice = {"NW_N" : [1, 1, 1, 0], 
               "N_NE": [1, 0, 1, -1],
               "NE_E" : [1, -1, 0, -1], 
               "E_SE": [0, -1, -1, -1],
               "SE_S" : [-1, -1, -1, 0],  
               "S_SW": [-1, 0, -1, 1],
               "SW_W" : [-1, 1, 0, 1], 
               "W_NW": [0, 1, 1, 1]}



# if maxRow = maxCol = 3
# [1, 2, 3]
# [1, 2, 3]
# [1, 2, 3]
#DEM = [[i+1 for i in range(maxCol)] for j in range(maxRow)]
DEM = [[1, 1, 1, 1, 1],
       [1, 9, 9, 9, 1],
       [1, 9, 1, 9, 1],
       [1, 9, 9, 9, 1],
       [1, 1, 1, 1, 1]]

# Line below copis ref to DEM into AuxGrid
# therefore modifying AuxGrid with DEM does nothing?
# 
# Either way, below way does not allow proper modifying of AuxGrid
#AuxGrid = DEM

# The below way makes a deep copy, works as intended
AuxGrid = [[DEM[i][j] for i in range(len(DEM[j]))] for j in range(len(DEM))]
TrackerGrid = [[0 for i in range(maxCol)] for j in range(maxRow)] 
vizScore = [[0 for i in range(maxCol)] for j in range(maxRow)]
vizViews = [[0 for i in range(maxCol)] for j in range(maxRow)]


# 3D point = (row, col, height)
# viewPoint = (2, 2, DEM[2][2])


# decides if some point is in grid
def inGrid(currRow, currCol):
    rowIn = (currRow >= 0) and (currRow < maxRow)
    colIn = (currCol >= 0) and (currCol < maxCol)
    return rowIn and colIn



# gets all potential indices that need to be tested from a slice
def getSliceIndices(direction, initRow, initCol):
    indices = []
    indices.append((initRow, initCol))
    i = j = 0
    tracker = 0

    if direction == "NW_N":
        while True:
            j = 0
            i -= 1
            tracker = 0
            while j > (i - 1):
                if inGrid(initRow+i, initCol+j) == True:
                    newIndex = (initRow+i, initCol+j)
                    indices.append(newIndex)
                    tracker += 1
                j -= 1
            if tracker == 0: break # no more indices
        
    elif direction == "N_NE":
        while True: # to keep incrementing j and searching
            j = 0
            i -= 1
            tracker = 0
            while j < (-1 * (i - 1)):
                if inGrid(initRow+i, initCol+j) == True:
                    newIndex = (initRow+i, initCol+j)
                    indices.append(newIndex)
                    tracker += 1
                j += 1
            if tracker == 0: break # no more indices
    
    elif direction == "NE_E":
        while True: # to keep incrementing j and searching
            j += 1
            i = 0
            tracker = 0
            while i > (-1 * (j + 1)):
                if inGrid(initRow+i, initCol+j) == True:
                    newIndex = (initRow+i, initCol+j)
                    indices.append(newIndex)
                    tracker += 1
                i -= 1
            if tracker == 0: break # no more indices

    elif direction == "E_SE":
        while True:
            j += 1
            i = 0
            tracker = 0
            while i < (j + 1):

                #print(f"row: {initRow+i}")
                #print(f"col: {initCol+j}")

                if inGrid(initRow+i, initCol+j) == True:
                    newIndex = (initRow+i, initCol+j)
                    indices.append(newIndex)
                    tracker += 1
                i += 1
            if tracker == 0: break # no more indices
    
    elif direction == "SE_S":
        while True:
            j = 0
            i += 1
            tracker = 0
            while j < (i + 1):
                if inGrid(initRow+i, initCol+j) == True:
                    newIndex = (initRow+i, initCol+j)
                    indices.append(newIndex)
                    tracker += 1
                j += 1
            if tracker == 0: break # no more indices
    
    elif direction == "S_SW":
        while True:
            j = 0
            i += 1
            tracker = 0
            while j > (-1 * (i + 1)):
                if inGrid(initRow+i, initCol+j) == True:
                    newIndex = (initRow+i, initCol+j)
                    indices.append(newIndex)
                    tracker += 1
                j -= 1
            if tracker == 0: break # no more indices
    
    elif direction == "SW_W":
        while True:
            j -= 1
            i = 0
            tracker = 0
            while i < (-1 * (j - 1)):
                if inGrid(initRow+i, initCol+j) == True:
                    newIndex = (initRow+i, initCol+j)
                    indices.append(newIndex)
                    tracker += 1
                i += 1
            if tracker == 0: break # no more indices
    
    elif direction == "W_NW":
        while True:
            j -= 1
            i = 0
            tracker = 0
            while i > (j - 1):
                if inGrid(initRow+i, initCol+j) == True:
                    newIndex = (initRow+i, initCol+j)
                    indices.append(newIndex)
                    tracker += 1
                i -= 1
            if tracker == 0: break # no more indices

    return indices


  
# processing sector boundaries
# 8 different lines
def processEdge(calc, vpRow, vpCol, direction):
    rowInc = dirIncEdge[direction][0]
    colInc = dirIncEdge[direction][1]

    #print(f"vpRow: {vpRow}")
    #print(f"vpCol: {vpCol}")

    # getting viewpoint
    vp = (vpRow, vpCol, AuxGrid[vpRow][vpCol])
    
    # getting row and col of destination
    destRow = vpRow + rowInc
    destCol = vpCol + colInc

    # Setting up previous node used for vert angle
    prev = (0,0,0)
    visible = True

    # while points accessed are in the grid, keep going
    while inGrid(destRow, destCol) == True:
        if TrackerGrid[destRow][destCol] == 1:
            destRow += rowInc
            destCol += colInc 
            continue
        else:
            TrackerGrid[destRow][destCol] = 1

        #print(f"destRow: {destRow}, destCol: {destCol}")
        #print(f"destRow: {destRow}, destCol: {destCol}")
        destPoint = (destRow, destCol, AuxGrid[destRow][destCol])
        
        # <vp dest> = dest - vp
        viewVec = calc.vecSub(destPoint, vp)      # sample vector
        viewAng = calc.vecAngle((0,0,1), viewVec) # zvs angle
        
        if prev[2] != 0:
            # <vp, prev> = prev - vp
            prevVec = calc.vecSub(prev, vp)           # prev vector 
            prevAng = calc.vecAngle((0,0,1), prevVec) # zvp angle

            # condition for visibility (zvp > zvs)
            if prevAng >= viewAng:
                visible = True
            # if not visible, store minimum height of inter point
            else:
                visible = False
                # line of destPoint expand into inifinite Z-axis
                aboveDestPoint = (destPoint[0], destPoint[1], destPoint[2]+1)

                # Min height for visibility
                projectedHeight = calc.getLargerHeight(vp, prev, destPoint, aboveDestPoint)
                AuxGrid[destPoint[0]][destPoint[1]] = projectedHeight

        # setting prev to current desitination
        prev = (destRow, destCol, AuxGrid[destRow][destCol])
        
        # updating visibilty structs
        if visible:
            vizScore[vpRow][vpCol] += 1
            vizViews[destRow][destCol] += 1

        # update row and col to get new dest
        destRow += rowInc
        destCol += colInc


# Processing each slice
def processSlice(calc, vpRow, vpCol, direction):
    # get viewpoint from auxgrid
    vp = (vpRow, vpCol, AuxGrid[vpRow][vpCol])
    #sampleRow = dirIncSlice[direction][0]
    #sampleCol = dirIncSlice[direction][1] 
    prev1RowInc = dirIncSlice[direction][0] 
    prev1ColInc = dirIncSlice[direction][1] 
    prev2RowInc = dirIncSlice[direction][2] 
    prev2ColInc = dirIncSlice[direction][3]

    #if inGrid(sampleRow, sampleCol) == True:
    indices = getSliceIndices(direction, vpRow, vpCol)
    
    # for each possible index, get the previous pts (prev1, prev2)
    for i in indices:
        if TrackerGrid[i[0]][i[1]] == 1:
            continue
        else:
            TrackerGrid[i[0]][i[1]] = 1
        
        # use incs to get prev1
        prev1Row = i[0] + prev1RowInc
        prev1Col = i[1] + prev1ColInc
        prev1Point = (prev1Row, prev1Col, AuxGrid[prev1Row][prev1Col])
        
        # use incs to get prev2
        prev2Row = i[0] + prev2RowInc
        prev2Col = i[1] + prev2ColInc
        prev2Point = (prev2Row, prev2Col, AuxGrid[prev2Row][prev2Col])

        # make plane from (vp, prev1, prev2))

        # defining current point
        currPoint = (i[0], i[1], AuxGrid[i[0]][i[1]])

        # a point directly above currPoint
        aboveCurrPoint = (currPoint[0], currPoint[1], currPoint[2]+1)
        
        # intersection of plane and currPoint going infinite into Z
        interPoint = calc.linePlaneIntersect(currPoint, aboveCurrPoint, 
                                             vp, prev1Point, prev2Point)
        
        projectedHeightOnPlane = interPoint[2]
        actualHeight = currPoint[2]

        # Not visibile
        if projectedHeightOnPlane > actualHeight:
            AuxGrid[i[0]][i[1]] = projectedHeightOnPlane
        # visible!
        else: # projectedHeight < actualHeight
            #AuxGrid[i[0]][i[1]] = actualHeight
            vizScore[vpRow][vpCol] += 1
            vizViews[i[0]][i[1]] += 1

# Process all boundaries
def processAllEdges(calc, vpRow, vpCol):
    for key in dirIncEdge:
        print(f"{key}")
        processEdge(calc, vpRow, vpCol, key)


# Process all boundaries
def processAllSlices(calc, vpRow, vpCol):
    for key in dirIncSlice:
        print(f"{key}")
        processSlice(calc, vpRow, vpCol, key)


# zero out everything
def refreshZero(ll):
    for i in range(0,len(ll)):
        for j in range(0,len(ll[i])):
            ll[i][j] = 0

# match DEM
# Currently unused
def refreshAux(ll):
    for i in range(0,len(ll)):
        for j in range(0,len(ll[i])):
            ll[i][j] = DEM[i][j]


# match DEM
def refreshAuxAlt():
    for i in range(0,len(AuxGrid)):
        for j in range(0,len(AuxGrid[i])):
            AuxGrid[i][j] = DEM[i][j]
        


def run(case): 
    myvsm = vsm()
    
    refreshAuxAlt()
    refreshZero(TrackerGrid)
    refreshZero(vizViews)
    refreshZero(vizScore)

    ################
    if case == "oneEdge":
        print("oneEdge")
        processEdge(myvsm, srcRow, srcCol, "N")
    ################
    elif case == "NW_N_Slice":
        print("NW_N_Slice")
        processSlice(myvsm, srcRow, srcCol, "NW_N")
    ################
    elif case == "N_NE_Slice":
        print("N_NE_Slice")
        processSlice(myvsm, srcRow, srcCol, "N_NE")
    ################
    elif case == "allSlices":
        print("allSlices")
        processAllSlices(myvsm, srcRow, srcCol) 
    ################
    elif case == "allEdges":
        print("allEdges")
        processAllEdges(myvsm, srcRow, srcCol)
    
    print("AuxGrid: ")
    for ll in AuxGrid:
        print(ll)

    print("\nvizViews: ")
    for ll in vizViews:
        print(ll)

    print("\nvizScore: ")
    for ll in vizScore:
        print(ll)
    
    print("\nTrackerGrid: ")
    for ll in TrackerGrid:
        print(ll)




#############
# /V\ain
# * global variables can be modified in main()
# using functions works on global vars as well
# Pass the global vars into functions to modify them reliably
#############
def main():
    print("hello")
    t1 = (1, 2, 3)
    
    t2 = (4, 5, 6)
    myvsm = vsm()
    
    myvsm.vecAdd(t1, t2)
    print(myvsm.rlist)
    
    myvsm.rlist = [4, 5, 6]
    print(myvsm.rlist)
    
    #print(f"{vecAdd(t1, t2)}")
    #print(f"{vecSub(t1, t2)}")
    #print(f"{scalarMult(-1, t1)}")
    #print(f"{vecAngle(t1, t2)}")
    
    run("allSlices") 
    run("allEdges")
    

main()



