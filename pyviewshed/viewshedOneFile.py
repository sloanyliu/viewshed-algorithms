#!/usr/bin/python3

import math as math

# init
maxRow = maxCol = 5
srcRow = srcCol = maxRow // 2


# "direction" : [rowInc, colInc]
dirIncEdge = {"N" : [-1,0], "NE" : [-1,1],
          "E" : [0,1],  "SE" : [1,1],
          "S" : [1,0],  "SW" : [1,-1],
          "W" : [0,-1], "NW" : [-1,-1]}


# "direction" : [testi, testj, xi, xj, yi, yj]
# i is row, col is j

#dirIncSlice = {"NW_N" : [srcRow-2, srcCol-1, 1, 1, 1, 0], 
#               "N_NE": [srcRow-2, srcCol+1, 1, 0, 1, -1],
#               "NE_E" : [srcRow-1, srcCol+2, 1, -1, 0, -1], 
#               "E_SE": [srcRow+1, srcCol+2, 0, -1, -1, -1],
#               "SE_S" : [srcRow+2, srcCol+1, -1, -1, -1, 0],  
#               "S_SW": [srcRow+2, srcCol-1, -1, 0, -1, 1],
#               "SW_W" : [srcRow+1, srcCol-2, -1, 1, 0, 1], 
#               "W_NW": [srcRow-1, srcCol-2, 0, 1, 1, 1]}

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
#viewPoint = (2, 2, DEM[2][2])


# decides if some point is in grid
def inGrid(currRow, currCol):
    rowIn = (currRow >= 0) and (currRow < maxRow)
    colIn = (currCol >= 0) and (currCol < maxCol)
    return rowIn and colIn


# (x1, y1, z1) - (x2, y2, z2)
def vecSub(v1, v2):
    #print(f"v1: {v1}")
    #print(f"v2: {v2}")
    return (v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2])


# (x1, y1, z1) - (x2, y2, z2)
def vecAdd(v1, v2):
    return (v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2])

# v1 dot v2
def dotProd(v1, v2):
    return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

# c * (x, y, z)
def scalarMult(c, v):
    return (c*v[0], c*v[1], c*v[2])


# finds angle between v1 and v2
def vecAngle(v1, v2):
    x1s = v1[0] ** 2
    y1s = v1[1] ** 2
    z1s = v1[2] ** 2

    x2s = v2[0] ** 2
    y2s = v2[1] ** 2
    z2s = v2[2] ** 2

    dotPr = dotProd(v1, v2)
    v1Mag = math.sqrt(x1s + y1s + z1s)
    v2Mag = math.sqrt(x2s + y2s + z2s)

    angle = math.acos(dotPr / (v1Mag * v2Mag))
    return angle


# used in intersection finding
def getD(a, b, c, d):
    partX = (a[0] - b[0]) * (c[0] - d[0])
    partY = (a[1] - b[1]) * (c[1] - d[1])
    partZ = (a[2] - b[2]) * (c[2] - d[2])
    return partX + partY + partZ


# finds shortest distance
# one -> vp ViewPoint
# two -> prev prevPoint
# three -> dp DestinationPoint
# four -> idz (0, 0, 1)
#
#      |/
#      * <-- min height neede for dest to be visible
#     /|
#    / |
#   /  |
# vp  dest
#
# Technically this algo finds the shortestDistanceBetween two 3D Lines
# - we take the largest height, of the two points that make up the line
# - this is the minimum height needed for viz
#
# If the two lines do touch:
# - distance of lines would be zero
# - and point on either line would be the same
def shortestDistBetweenLines(one, two, three, four):
    d_1343 = float(getD(one, three, four, three))
    d_4321 = float(getD(four, three, two, one))
    d_1321 = float(getD(one, three, two, one))
    d_4343 = float(getD(four, three, four, three))
    d_2121 = float(getD(two, one, two, one))
    
    mu_a_numer = float((d_1343 * d_4321) - (d_1321 * d_4343))
    mu_a_denom = float((d_2121 * d_4343) - (d_4321 * d_4321))
    #print(f"{mu_a_numer} / {mu_a_denom}")
    if mu_a_numer == 0.0:
        mu_a = 0.0
    else:
        mu_a = float(mu_a_numer) / float(mu_a_denom)

    mu_b_numer = float(d_1343 + (mu_a * d_4321))
    mu_b_denom = d_4343
    #print(f"{mu_b_numer} / {mu_b_denom}")
    if mu_b_numer == 0.0:
        mu_b = 0.0
    else:
        mu_b = float(mu_b_numer) / float(mu_b_denom)
    
    pa = vecAdd(scalarMult(mu_a, vecSub(two, one)), one)
    pb = vecAdd(scalarMult(mu_b, vecSub(four, three)), three)

    return (pa, pb)


# Gets the minimum height of visibility
# = largest heighyt between two points of shortest line between 2 3D lines
def getMinHeight(one, two, three, four):
    # line pa -- pb is the shortest line between two 3D lines
    # We return the largest height
    (pa, pb) = shortestDistBetweenLines(one, two, three, four)
    mh = pa[2] if (pa[2] > pb[2]) else pb[2]
    return mh
    

# This technically find the intersection of two 2D lines
#
# But this works for us since we do not care if they intersect in 3D
# We just want to find the heigh at which the two vectors cross paths:
#
# To do this, we ditch the axis that is same amongst them
# then treat what remains as a 2 2D lines
def lineLineIntersect(one, two, three, four):
    print(f"one: {one}")
    print(f"two: {two}")
    print(f"three: {three}")
    print(f"four: {four}")



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

# finds coefficients of a plane given three points
# Ax + By + Cz + D = 0 (A, B, C, D)
def findPlaneCoeffs(p1, p2, p3):
    x1 = p1[0]
    y1 = p1[1]
    z1 = p1[2]
    
    x2 = p2[0]
    y2 = p2[1]
    z2 = p2[2]
    
    x3 = p3[0]
    y3 = p3[1]
    z3 = p3[2]

    a = (y1*(z2 - z3)) + (y2*(z3 - z1)) + (y3*(z1 - z2))
    b = (z1*(x2 - x3)) + (z2*(x3 - x1)) + (z3*(x1 - x2))
    c = (x1*(y2 - y3)) + (x2*(y3 - y1)) + (x3*(y1 - y2))
    
    d1 = x1*((y2*z3) - (y3*z2))
    d2 = x2*((y3*z1) - (y1*z3))
    d3 = x3*((y1*z2) - (y2*z1))
    d = -1*(d1 + d2 + d3)

    return (a, b, c, d)


# uses findPlaneCoeffs to find Normal Vector of Plane
def findPlaneNormal(p1, p2, p3):
    cfs = findPlaneCoeffs(p1, p2, p3)
    return (cfs[0], cfs[1], cfs[2])


# finds the intersection point of line and plane
def linePlaneIntersect(l1, l2, p1, p2, p3):
    norm = findPlaneNormal(p1, p2, p3)
    
    u_numer = dotProd(norm, vecSub(p3, l1))
    u_denom = dotProd(norm, vecSub(l2, l1))
   
    #print(f"numer: {u_numer}")
    #print(f"denom: {u_denom}")

    # if denom is 0, that means line is on the plane
    # -> which means we are visibile (inifinite solns)
    # returning -1 to ensure projectedHeight < actualHeight == visible
    if u_denom == 0:
        return (-1, -1, -1)

    u = float(u_numer / u_denom)
    
    # intersect pt = p1 + u (p2 - p1)
    # given p1 and p2 lie on the line
    p_a = scalarMult(u, vecSub(l2, l1))
    p = vecAdd(l1, p_a)
    return p

  
# processing sector boundaries
# 8 different lines
def processEdge(vpRow, vpCol, direction):
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
        viewVec = vecSub(destPoint, vp)      # sample vector
        viewAng = vecAngle((0,0,1), viewVec) # zvs angle
        
        if prev[2] != 0:
            # <vp, prev> = prev - vp
            prevVec = vecSub(prev, vp)           # prev vector 
            prevAng = vecAngle((0,0,1), prevVec) # zvp angle

            # condition for visibility (zvp > zvs)
            if prevAng >= viewAng:
                visible = True
            # if not visible, store minimum height of inter point
            else:
                visible = False
                # line of destPoint expand into inifinite Z-axis
                aboveDestPoint = (destPoint[0], destPoint[1], destPoint[2]+1)

                # Min height for visibility
                projectedHeight = getMinHeight(vp, prev, destPoint, aboveDestPoint)
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
def processSlice(vpRow, vpCol, direction):
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
        interPoint = linePlaneIntersect(currPoint, aboveCurrPoint, 
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
def processAllEdges(vpRow, vpCol):
    for key in dirIncEdge:
        print(f"{key}")
        processEdge(vpRow, vpCol, key)


# Process all boundaries
def processAllSlices(vpRow, vpCol):
    for key in dirIncSlice:
        print(f"{key}")
        processSlice(vpRow, vpCol, key)


# zero our everything
def refreshZero(ll):
    for i in range(0,len(ll)):
        for j in range(0,len(ll[i])):
            ll[i][j] = 0

# match DEM
def refreshAux(ll):
    for i in range(0,len(ll)):
        for j in range(0,len(ll[i])):
            ll[i][j] = DEM[i][j]
        


def run(case): 
    ################
    if case == "oneEdge":
        processEdge(srcRow, srcCol, "N")
    ################
    elif case == "NW_N_Slice":
        processSlice(srcRow, srcCol, "NW_N")
    ################
    elif case == "N_NE_Slice":
        processSlice(srcRow, srcCol, "N_NE")
    ################
    elif case == "N_NE_Slice":
        processSlice(srcRow, srcCol, "N_NE")
    ################
    elif case == "allSlices":
        processAllSlices(srcRow,srcCol) 
    ################
    elif case == "allEdges":
        processAllEdges(srcRow,srcCol)
    
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
#############
def main():
    print("hello")
    t1 = (1, 2, 3)
    t2 = (4, 5, 6)
    #run("oneSlice")
    #run("oneEdge")
    
    refreshAux(AuxGrid)
    refreshZero(TrackerGrid)
    refreshZero(vizViews)
    refreshZero(vizScore)
    run("allSlices")
    
    refreshAux(AuxGrid)
    #AuxGrid = DEM
    refreshZero(TrackerGrid)
    refreshZero(vizViews)
    refreshZero(vizScore)
    run("allEdges")
    
    #print(f"{vecAdd(t1, t2)}")
    #print(f"{vecSub(t1, t2)}")
    #print(f"{scalarMult(-1, t1)}")
    #print(f"{vecAngle(t1, t2)}")


main()

#refreshAux(AuxGrid)
#refreshZero(TrackerGrid)
#refreshZero(vizViews)
#refreshZero(vizScore)
#run("NW_N_Slice")
#run("N_NE_Slice")
#run("allEdges")
#run("allSlices")

#refreshAux(AuxGrid)
#refreshZero(vizViews)
#refreshZero(vizScore)
#run("allEdges")



