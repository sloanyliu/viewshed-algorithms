#!/usr/bin/python3

import math as math


class viewshedMath:
    # dummy variable for sake of oop testing
    rlist = [1, 2, 4] 

    def __init__(self):
        pass

    # (x1, y1, z1) - (x2, y2, z2)
    def vecSub(self, v1, v2):
        #print(f"v1: {v1}")
        #print(f"v2: {v2}")
        return (v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2])

    def vecMag(self, v1):
        return math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)

    # (x1, y1, z1) - (x2, y2, z2)
    # This function does have self, 
    # so we must call it from instance of class like so:
    # >> import viewshedMath as vsmath
    #    v1 = (1, 2, 3)
    #    v2 = (4, 5, 6)
    #    vsInst = vsmath.viewshedMath()
    #    res = vsInst.vecAdd(v1, v2)  
    #
    def vecAdd(self, v1, v2):
        return (v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2])
    

    # Since this function does not have a self as an arg,
    # we can't call it with an instance of this class
    # we need to call it directly like so: 
    # >> import viewshedMath as vsmath
    #    v1 = (1, 2, 3)
    #    v2 = (4, 5, 6)
    #    res = vsmath.viewshedMath.vecAdd(v1, v2)
    #
    #def vecAdd(v1, v2):
    #    return (v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2])

    # v1 dot v2
    def dotProd(self, v1, v2):
        return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

    # c * (x, y, z)
    def scalarMult(self, c, v):
        return (c*v[0], c*v[1], c*v[2])


    # finds angle between v1 and v2
    def vecAngle(self, v1, v2):
        x1s = v1[0] ** 2
        y1s = v1[1] ** 2
        z1s = v1[2] ** 2

        x2s = v2[0] ** 2
        y2s = v2[1] ** 2
        z2s = v2[2] ** 2

        dotPr = self.dotProd(v1, v2)
        #v1Mag = math.sqrt(x1s + y1s + z1s)
        #v2Mag = math.sqrt(x2s + y2s + z2s)
        v1Mag = self.vecMag(v1)
        v2Mag = self.vecMag(v2)

        angle = math.acos(dotPr / (v1Mag * v2Mag))
        return angle


    # used in intersection finding
    def getD(self, a, b, c, d):
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
    def shortestDistBetweenLines(self, one, two, three, four):
        d_1343 = float(self.getD(one, three, four, three))
        d_4321 = float(self.getD(four, three, two, one))
        d_1321 = float(self.getD(one, three, two, one))
        d_4343 = float(self.getD(four, three, four, three))
        d_2121 = float(self.getD(two, one, two, one))
        
        mu_a = 0.0
        mu_b = 0.0
        mu_a_numer = float((d_1343 * d_4321) - (d_1321 * d_4343))
        mu_a_denom = float((d_2121 * d_4343) - (d_4321 * d_4321))
        #print(f"{mu_a_numer} / {mu_a_denom}")
        if mu_a_numer != 0.0:
            mu_a = float(mu_a_numer) / float(mu_a_denom)

        mu_b_numer = float(d_1343 + (mu_a * d_4321))
        mu_b_denom = d_4343
        #print(f"{mu_b_numer} / {mu_b_denom}")
        if mu_b_numer != 0.0:
            mu_b = float(mu_b_numer) / float(mu_b_denom)
        
        pa1 = self.scalarMult(mu_a, self.vecSub(two, one))
        pb1 = self.scalarMult(mu_b, self.vecSub(four, three))
        pa = self.vecAdd(pa1, one)
        pb = self.vecAdd(pb1, three)

        return (pa, pb)


    # Gets largest height between 2 points of shortest line between 2 3D lines
    def getLargerHeight(self, one, two, three, four):
        # line pa -- pb is the shortest line between two 3D lines
        # We return the largest height
        (pa, pb) = self.shortestDistBetweenLines(one, two, three, four)
        mh = pa[2] if (pa[2] > pb[2]) else pb[2]
        return mh
    

    # finds coefficients of a plane given three points
    # Ax + By + Cz + D = 0 (A, B, C, D)
    def findPlaneCoeffs(self, p1, p2, p3):
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
    def findPlaneNormal(self, p1, p2, p3):
        cfs = self.findPlaneCoeffs(p1, p2, p3)
        return (cfs[0], cfs[1], cfs[2])


    # finds the intersection point of line and plane
    def linePlaneIntersect(self, l1, l2, p1, p2, p3):
        norm = self.findPlaneNormal(p1, p2, p3)
        
        u_numer = self.dotProd(norm, self.vecSub(p3, l1))
        u_denom = self.dotProd(norm, self.vecSub(l2, l1))
       
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
        p_a = self.scalarMult(u, self.vecSub(l2, l1))
        p = self.vecAdd(l1, p_a)
        return p



