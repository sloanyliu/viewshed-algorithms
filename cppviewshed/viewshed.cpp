//#include <iostream>
#include <cmath>
//#include <set>
//#include <map>
//#include <vector>

#include "viewshed.hpp"

// ------------------
// Ctor(dummy)
ViewShed::ViewShed(float vpRow, float vpCol)
{
  VPRow = vpRow;
  VPCol = vpCol;
}


// ------------------
// Dtor
ViewShed::~ViewShed()
{
}


// ------------------
// Dummy hello function
void ViewShed::hello(void)
{
  std::cout << "hello!\n";
}

/*
// Below implementation works but its verbose, does not need all the ViewShed::'s
// Only in the return type of function
ViewShed::testS ViewShed::testFunc(ViewShed::testS tempS, ViewShed::testE tempE)
{
  ViewShed::testS ret = {tempS.a, tempS.b};
  std::cout << "tempE: " << static_cast<int>(tempE) << std::endl;
  std::cout << "tempS: (" << tempS.a << "," << tempS.b << ")" << std::endl;
  return ret;
}

// This one is less verbose
ViewShed::testS ViewShed::testFunc(testS tempS, testE tempE)
{
  testS ret = {tempS.a, tempS.b};
  std::cout << "tempE: " << static_cast<int>(tempE) << std::endl;
  std::cout << "tempS: (" << tempS.a << "," << tempS.b << ")" << std::endl;
  return ret;
}


// Below implementation works but its verbose, does not need all the ViewShed::'s
// Only in the return type of function
void ViewShed::privStructs(void)
{
  ViewShed::testS tq = {1, 2};
  ViewShed::testE tw = ViewShed::testE::b;
  ViewShed::testS tt = testFunc(tq, tw);
}

// This one is less verbose
void ViewShed::privStructs(void)
{
  testS tq = {1, 2};
  testE tw = testE::b;
  testS tt = testFunc(tq, tw);
}
*/


void ViewShed::printXYZ(XYZ p)
{
  std::cout << "<" << p.X << ", " 
            << p.Y << ", "
            << p.Z << ">" << std::endl;
}

// ------------------
// Vector Add function
ViewShed::XYZ ViewShed::vecAdd(XYZ p1, XYZ p2)
{
  XYZ res;
  res.X = p1.X + p2.X;
  res.Y = p1.Y + p2.Y;
  res.Z = p1.Z + p2.Z;
  return res;
}

// ------------------
// Vector Sub function
ViewShed::XYZ ViewShed::vecSub(XYZ p1, XYZ p2)
{
  XYZ res;
  res.X = p1.X - p2.X;
  res.Y = p1.Y - p2.Y;
  res.Z = p1.Z - p2.Z;
  return res;
} // p1 - p2

// ------------------
// Vector Magnitude function
float ViewShed::vecMag(XYZ p1)
{
  float xs = p1.X * p1.X;
  float ys = p1.Y * p1.Y;
  float zs = p1.Z * p1.Z;
  float res = sqrt(xs + ys + zs);  
  return res;
} 

// ------------------
// Dot Product function
float ViewShed::dotProduct(XYZ p1, XYZ p2)
{
  float tx = p1.X * p2.X;
  float ty = p1.Y * p2.Y;
  float tz = p1.Z * p2.Z;
  return tx + ty + tz;
} // p1 dot p2

// ------------------
// Scalar Multiplication
ViewShed::XYZ ViewShed::scalarMult(float c, XYZ p1)
{
  XYZ res;
  res.X = p1.X * c;
  res.Y = p1.Y * c;
  res.Z = p1.Z * c;
  return res;
} // c * p1

// ------------------
// Vector Angle function
float ViewShed::vecAngle(XYZ p1, XYZ p2)
{
  float p1_xs = p1.X * p1.X;
  float p1_ys = p1.Y * p1.Y;
  float p1_zs = p1.Z * p1.Z;

  float p2_xs = p2.X * p2.X;
  float p2_ys = p2.Y * p2.Y;
  float p2_zs = p2.Z * p2.Z;
  
  float dprod = dotProduct(p1, p2);
  float p1_mag = vecMag(p1);
  float p2_mag = vecMag(p2);

  float mag_prod = p1_mag * p2_mag;
  float angle = acos(dprod / mag_prod);
  return angle;
} // angle between p1 and p2

// ------------------
// Used for dist between 2 3D lines
float ViewShed::getD(XYZ a, XYZ b, XYZ c, XYZ d)
{
  float px = (a.X - b.X) * (c.X - d.X);
  float py = (a.Y - b.Y) * (c.Y - d.Y);
  float pz = (a.Z - b.Z) * (c.Z - d.Z);
  return px + py + pz;
} 

// ------------------
// Finds normal vector of a plane
ViewShed::XYZ ViewShed::planeNormal(XYZ p1, XYZ p2, XYZ p3)
{
  float x1 = p1.X;
  float y1 = p1.Y;
  float z1 = p1.Z;
  float x2 = p2.X;
  float y2 = p2.Y;
  float z2 = p2.Z;
  float x3 = p3.X;
  float y3 = p3.Y;
  float z3 = p3.Z;
  
  XYZ norm;
  norm.X = (y1*(z2-z3)) + (y2*(z3-z1)) + (y3*(z1-z2));
  norm.Y = (z1*(x2-x3)) + (z2*(x3-x1)) + (z3*(x1-x2));
  norm.Z = (x1*(y2-y3)) + (x2*(y3-y1)) + (x3*(y1-y2));

  return norm;
} // Gets normal vector of a plane

// ------------------
// LinePlane intersection point
ViewShed::XYZ ViewShed::linePlaneIntersect(XYZ l1, XYZ l2, XYZ p1, XYZ p2, XYZ p3)
{
  XYZ norm = planeNormal(p1, p2, p3);

  XYZ sub1 = vecSub(l2, l1);
  XYZ sub2 = vecSub(p3, l1);
  float u_denom = dotProduct(norm, sub1);
  float u_numer = dotProduct(norm, sub2);

  if (u_denom == 0) {
    XYZ ret = {-1, -1, -1};
    return ret;
  }

  float u = u_numer / u_denom;

  XYZ inter_a = scalarMult(u, sub1);
  XYZ inter = vecAdd(l1, inter_a);

  return inter;
} 

// ------------------
// Minimum Height needed for visibility from viewpoint
// Finds the z value of the dest point projectd onto the prev vector
// Finds the interstion of (destpoint-destpoint to Z)
float ViewShed::minVisHeight(XYZ p1, XYZ p2, XYZ p3, XYZ p4)
{
  float d_1343 = getD(p1, p3, p4, p3);
  float d_4321 = getD(p4, p3, p2, p1);
  float d_1321 = getD(p1, p3, p2, p1);
  float d_4343 = getD(p4, p3, p4, p3);
  float d_2121 = getD(p2, p1, p2, p1);

  float mu_a = 0;
  float mu_a_numer = (d_1343 * d_4321) - (d_1321 * d_4343);
  float mu_a_denom = (d_2121 * d_4343) - (d_4321 * d_4321);
  if (mu_a_numer != 0) {
    mu_a = mu_a_numer / mu_a_denom;
  }

  float mu_b = 0;
  float mu_b_numer = d_1343 + (mu_a * d_4321);
  float mu_b_denom = d_4343;
  if (mu_b_numer != 0) {
    mu_b = mu_b_numer / mu_b_denom;
  }

  XYZ suba = vecSub(p2, p1);
  XYZ ia_1 = scalarMult(mu_a, suba);
  XYZ ia = vecAdd(ia_1, p1);

  XYZ subb = vecSub(p4, p3);
  XYZ ib_1 = scalarMult(mu_b, subb);
  XYZ ib = vecAdd(ib_1, p3);
  
  float res = (ib.Z > ia.Z) ? ib.Z : ia.Z;
  return res;
} 


// ------------------
// Load in DEM
// AuxGrid = DEM
// Sotre all zeros in all other grids
void ViewShed::bootStrap(std::vector<std::vector<float>> dem)
{
  MaxRow = dem.size();
  MaxCol = dem[0].size();

  for (int i = 0; i < MaxRow; i++) {
    AuxGrid.push_back(dem[i]);
  }

  for (int i = 0; i < MaxRow; i++) {
    std::vector<int> temp;
    for (int j = 0; j < MaxCol; j++) {
      temp.push_back(0);
    }
    vizViews.push_back(temp);
  }

  for (int i = 0; i < MaxRow; i++) {
    std::vector<int> temp;
    for (int j = 0; j < MaxCol; j++) {
      temp.push_back(0);
    }
    vizScore.push_back(temp);
  }
  
  for (int i = 0; i < MaxRow; i++) {
    std::vector<int> temp;
    for (int j = 0; j < MaxCol; j++) {
      temp.push_back(0);
    }
    TrackerGrid.push_back(temp);
  }
}


// ------------------
// Sets AuxGrid to DEM
// Sets rest to all 0's
void ViewShed::refreshAllGrids(std::vector<std::vector<float>> dem)
{
  refreshAG(dem);
  refreshVS();
  refreshVV();
  refreshTG();
}



void ViewShed::refreshVV(void)
{
  for (int i = 0; i < vizViews.size(); i++) {
    for (int j = 0; vizViews[i].size(); j++) {
      (vizViews[i])[j] = 0;
    }
  }

}

void ViewShed::refreshVS(void)
{
  for (int i = 0; i < vizScore.size(); i++) {
    for (int j = 0; vizScore[i].size(); j++) {
      (vizScore[i])[j] = 0;
    }
  }

}

void ViewShed::refreshTG(void)
{
  for (int i = 0; i < TrackerGrid.size(); i++) {
    for (int j = 0; TrackerGrid[i].size(); j++) {
      (TrackerGrid[i])[j] = 0;
    }
  }

}

void ViewShed::refreshAG(std::vector<std::vector<float>> dem)
{
  for (int i = 0; i < AuxGrid.size(); i++) {
    for (int j = 0; AuxGrid[i].size(); j++) {
      (AuxGrid[i])[j] = (dem[i])[j];
    }
  }
}




// ------------------
// Prints AuxGrid
void ViewShed::printAuxGrid(void)
{
  for (int i = 0; i < MaxRow; i++) {
    for (int j = 0; j < AuxGrid[i].size(); j++) {
      std::cout << (AuxGrid[i])[j] << " ";
    }
    std::cout << std::endl;
  }
}


// ------------------
// Prints vizViews
void ViewShed::printVizViews(void)
{
  for (int i = 0; i < MaxRow; i++) {
    for (int j = 0; j < vizViews[i].size(); j++) {
      std::cout << (vizViews[i])[j] << " ";
    }
    std::cout << std::endl;
  }
}


// ------------------
// Prints vizScore
void ViewShed::printVizScore(void)
{
  for (int i = 0; i < MaxRow; i++) {
    for (int j = 0; j < vizScore[i].size(); j++) {
      std::cout << (vizScore[i])[j] << " ";
    }
    std::cout << std::endl;
  }
}


// ------------------
// Prints TrackerGrid
void ViewShed::printTrackerGrid(void)
{
  for (int i = 0; i < MaxRow; i++) {
    for (int j = 0; j < TrackerGrid[i].size(); j++) {
      std::cout << (TrackerGrid[i])[j] << " ";
    }
    std::cout << std::endl;
  }
}


// ------------------
// Decides whether a point is within the grid coordinates
bool ViewShed::inGrid(float cRow, float cCol)
{
  bool rowIn = (cRow >= 0) && (cRow < MaxRow);
  bool colIn = (cCol >= 0) && (cCol < MaxCol);
  return rowIn && colIn;
}


// ------------------
// All the indices in a slice plus an edge
// Will use TrackerGrid to skip those already analyzed
std::vector<ViewShed::XYZ> ViewShed::getSliceIndices(GSlice slc, float iRow, float iCol)
{
  std::vector<XYZ> indices;
  XYZ ipt = {iRow, iCol, 0};
  indices.push_back(ipt);

  int i = 0, j = 0, tracker = 0;
  
  // Slice NW_N
  if (slc == GSlice::NW_N) { 
    while (true) {
      j=0;
      i--;
      tracker = 0;
      while (j > (i - 1)) {
        if (inGrid(iRow + i, iCol + j) == true) {
          XYZ newIndex = {iRow + i, iCol + j, 0};
          indices.push_back(newIndex);
          tracker++;
        }
        j--;
      }
      if (tracker == 0) break;
    }
  // Slice N_NE
  } else if (slc == GSlice::N_NE) {
    while (true) {
      j=0;
      i--;
      tracker = 0;
      while (j < (-1 * (i - 1))) {
        if (inGrid(iRow + i, iCol + j) == true) {
          XYZ newIndex = {iRow + i, iCol + j, 0};
          indices.push_back(newIndex);
          tracker++;
        }
        j++;
      }
      if (tracker == 0) break;
    }
  // Slice NE_E
  } else if (slc == GSlice::NE_E) {
    while (true) {
      j++;
      i=0;
      tracker = 0;
      while (i > (-1 * (j + 1))) {
        if (inGrid(iRow + i, iCol + j) == true) {
          XYZ newIndex = {iRow + i, iCol + j, 0};
          indices.push_back(newIndex);
          tracker++;
        }
        i--;
      }
      if (tracker == 0) break;
    }
  // Slice E_SE
  } else if (slc == GSlice::E_SE) {
    while (true) {
      j++;
      i=0;
      tracker = 0;
      while (i < (j + 1)) {
        if (inGrid(iRow + i, iCol + j) == true) {
          XYZ newIndex = {iRow + i, iCol + j, 0};
          indices.push_back(newIndex);
          tracker++;
        }
        i++;
      }
      if (tracker == 0) break;
    }
  // Slice SE_S
  } else if (slc == GSlice::SE_S) {
    while (true) {
      j=0;
      i++;
      tracker = 0;
      while (j < (i + 1)) {
        if (inGrid(iRow + i, iCol + j) == true) {
          XYZ newIndex = {iRow + i, iCol + j, 0};
          indices.push_back(newIndex);
          tracker++;
        }
        j++;
      }
      if (tracker == 0) break;
    }
  // Slice S_SW
  } else if (slc == GSlice::S_SW) {
    while (true) {
      j=0;
      i++;
      tracker = 0;
      while (j > (-1 * (i + 1))) {
        if (inGrid(iRow + i, iCol + j) == true) {
          XYZ newIndex = {iRow + i, iCol + j, 0};
          indices.push_back(newIndex);
          tracker++;
        }
        j--;
      }
      if (tracker == 0) break;
    }
  // Slice SW_W
  } else if (slc == GSlice::SW_W) {
    while (true) {
      j--;
      i=0;
      tracker = 0;
      while (i < (-1 * (j - 1))) {
        if (inGrid(iRow + i, iCol + j) == true) {
          XYZ newIndex = {iRow + i, iCol + j, 0};
          indices.push_back(newIndex);
          tracker++;
        }
        i++;
      }
      if (tracker == 0) break;
    }
  // Slice W_NW
  } else if (slc == GSlice::W_NW) {
    while (true) {
      j--;
      i=0;
      tracker = 0;
      while (i > (j - 1)) {
        if (inGrid(iRow + i, iCol + j) == true) {
          XYZ newIndex = {iRow + i, iCol + j, 0};
          indices.push_back(newIndex);
          tracker++;
        }
        i--;
      }
      if (tracker == 0) break;
    }
  }

  return indices;
}


// ------------------
// Processes a single edge
void ViewShed::processEdge(GEdge se, float vpRow, float vpCol)
{
  int rowInc = (mapIncEdge.at(se)).rowInc1;
  int colInc = (mapIncEdge.at(se)).colInc1;

  int vRow = static_cast<int>(vpRow);
  int vCol = static_cast<int>(vpCol);

  XYZ vp = {vpRow, vpCol, (AuxGrid[vRow])[vCol]};

  int dRow = vpRow + rowInc;
  int dCol = vpCol + colInc;

  XYZ prev = {0, 0, 0};
  bool visible = true;

  while (inGrid(dRow, dCol) == true) {
    if ((TrackerGrid[dRow])[dCol] == 1) {
      dRow += rowInc;
      dCol += colInc;
    } else {
      (TrackerGrid[dRow])[dCol] = 1;
    }

    XYZ destPoint = {static_cast<float>(dRow), 
                     static_cast<float>(dCol),
                     (AuxGrid[dRow])[dCol]};

    XYZ viewVec = vecSub(destPoint, vp);
    XYZ zVec = {0, 0, 1};
    float viewAng = vecAngle(zVec, viewVec);

    if (prev.Z != 0) {
      XYZ prevVec = vecSub(prev, vp);
      float prevAng = vecAngle(zVec, prevVec);

      if (prevAng >= viewAng) {
        visible = true;
      } else {
        visible = false;
        XYZ aboveDest = {destPoint.X, destPoint.Y, destPoint.Z + 1};
        float projHeight = minVisHeight(vp, prev, destPoint, aboveDest);
        (AuxGrid[destPoint.X])[destPoint.Y] = projHeight;
      }
    }
    prev = {destPoint.X, destPoint.Y, destPoint.Z};

    if (visible) {
      (vizScore[vRow])[vCol]++;
      (vizViews[dRow])[dCol]++;
    }

    dRow += rowInc;
    dCol += colInc;
  }
}


// ------------------
// Processes a single slice
void ViewShed::processSlice(GSlice slc, float vpRow, float vpCol)
{
  int prev1RowInc = (mapIncSlice.at(slc)).rowInc1;
  int prev1ColInc = (mapIncSlice.at(slc)).colInc1;
  int prev2RowInc = (mapIncSlice.at(slc)).rowInc2;
  int prev2ColInc = (mapIncSlice.at(slc)).colInc2;

  int vRow = static_cast<int>(vpRow);
  int vCol = static_cast<int>(vpCol);

  XYZ vp = {vpRow, vpCol, (AuxGrid[vRow])[vCol]};
  
  std::vector<XYZ> indices = getSliceIndices(slc, vpRow, vpCol);

  for (int i = 0; i < indices.size(); i++) {
    XYZ currIndex = indices[i];
    int cRow = static_cast<int>(currIndex.X);
    int cCol = static_cast<int>(currIndex.Y);

    if ((TrackerGrid[cRow])[cCol] == 1) {
      continue;
    } else {
      (TrackerGrid[cRow])[cCol] = 1;
    }

    int prev1Row = cRow + prev1RowInc;
    int prev1Col = cCol + prev1ColInc;
    XYZ prev1Point = {static_cast<float>(prev1Row), 
                      static_cast<float>(prev1Col), 
                      (AuxGrid[prev1Row])[prev1Col]};
    
    int prev2Row = cRow + prev2RowInc;
    int prev2Col = cCol + prev2ColInc;
    XYZ prev2Point = {static_cast<float>(prev2Row), 
                      static_cast<float>(prev2Col), 
                      (AuxGrid[prev2Row])[prev2Col]};
    
    XYZ currPoint = {currIndex.X, currIndex.Y, (AuxGrid[currIndex.X])[currIndex.Y]};

    XYZ aboveCurr = {currPoint.X, currPoint.Y, currPoint.Z + 1};

    XYZ interPoint = linePlaneIntersect(currPoint, aboveCurr, vp, prev1Point, prev2Point);

    float projHeight = interPoint.Z;
    float realHeight = currPoint.Z;

    if (projHeight > realHeight) {
      (AuxGrid[cRow])[cCol] = projHeight;
    } else {
      (vizScore[vRow])[vCol]++;
      (vizViews[cRow])[cCol]++;
    }
  }
}


// ------------------
// Processes all edges based on a view point
void ViewShed::processAllEdges(float vpRow, float vpCol)
{
  for (int i = 0; i < allEdges.size(); i++) {
    GEdge currEdge = allEdges.at(i);
    processEdge(currEdge, vpRow, vpCol);
  }
}


// -----------------
// Processes all slices based on a view point-
void ViewShed::processAllSlices(float vpRow, float vpCol)
{
  for (int i = 0; i < allSlices.size(); i++) {
    GSlice currSlice = allSlices.at(i);
    processSlice(currSlice, vpRow, vpCol);
  }
}






