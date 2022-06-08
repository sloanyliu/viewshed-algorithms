#ifndef VIEWSHED_H
#define VIEWSHED_H

#include <iostream>
#include <set>
#include <map>
#include <vector>


class ViewShed 
{
  public:
    
    //ViewShed(int vpRow, int vpCol);
    ViewShed(int oRow, int oCol);
    ~ViewShed();

    void hello(void);

    
    void bootStrap(std::vector<std::vector<float>> dem);
    void refreshAllGrids(void);
    
    void refreshAG(void);
    void refreshVS(void);
    void refreshVV(void);
    void refreshTG(void);

    void printAuxGrid(void);
    void printVizViews(void);
    void printVizScore(void);
    void printTrackerGrid(void);
    
    void processAllEdges(int vpRow, int vpCol);
    void processAllSlices(int vpRow, int vpCol);  
    void processAllPoints(void);

    bool queryVisibility(int srcRow, int srcCol, int destRow, int destCol);
    //void privStructs(void);
  
  private:

    typedef std::vector<int> vec1d;
    typedef std::vector<vec1d> vec2d;

    int ORow;
    int OCol;
    int rowAdj;
    int colAdj;

    struct XY
    {
      int X;
      int Y;
    };
    
    struct XYZ 
    {
      float X;
      float Y;
      float Z;
    };

    struct IncInfo2
    {
      int rowInc1;
      int colInc1;
    };

    struct IncInfo4 
    {
      int rowInc1;
      int colInc1;
      int rowInc2;
      int colInc2;
    };

    enum class GEdge 
    {
      N, NE, E, SE, S, SW, W, NW
    };

    enum class GSlice 
    {
      NW_N, N_NE, NE_E, E_SE, SE_S, S_SW, SW_W, W_NW
    };
    
    /*
    struct testS
    {
      int a;
      int b;
    };

    enum class testE
    {
      a, b, c
    };

    testS testFunc(testS tempS, testE tempE);
    */

    //const int some = 20;
    //int other;
    int MaxRow;
    int MaxCol;
    //first [r][c] is for look up
    //at each [r][c], there is 2d vector that is vizViews of [r][c]
    std::vector<std::vector<vec2d>> lut;
 
    std::vector<std::vector<float>> DEM;
    std::vector<std::vector<float>> AuxGrid;
    std::vector<std::vector<int>> vizScore;
    std::vector<std::vector<int>> vizViews;
    std::vector<std::vector<int>> TrackerGrid;

    void printXYZ(XYZ p);
    XYZ vecAdd(XYZ p1, XYZ p2);                // p1 + p2
    XYZ vecSub(XYZ p1, XYZ p2);                // p2 - p2
    float vecMag(XYZ p1);                      // |p1|
    float dotProduct(XYZ p1, XYZ p2);          // p1 dot p2
    XYZ scalarMult(float c, XYZ p1);           // c * p1
    float vecAngle(XYZ p1, XYZ p2);            // angle between p1 and p2
    float getD(XYZ a, XYZ b, XYZ c, XYZ d);    // Used for dist between 2 3D lines
    XYZ planeNormal(XYZ p1, XYZ p2, XYZ p3);   // Gets normal vector of a plane
    
    // intersection pt of line and plane
    XYZ linePlaneIntersect(XYZ l1, XYZ l2, XYZ p1, XYZ p2, XYZ p3); 
    
    // Gets min visible height when processing edges
    float minVisHeight(XYZ p1, XYZ p2, XYZ p3, XYZ p4); 
    bool inGrid(int cRow, int cCol);
    std::vector<XY> getSliceIndices(GSlice slc, int iRow, int iCol);
    void processEdge(GEdge se, int vpRow, int vpCol);
    void processSlice(GSlice slc, int vpRow, int vpCol);

    const std::vector<GEdge> allEdges = {GEdge::N, GEdge::NE, GEdge::E, GEdge::SE, 
                                         GEdge::S, GEdge::SW, GEdge::W, GEdge::NW};

    const std::vector<GSlice> allSlices = {GSlice::NW_N, GSlice::N_NE, 
                                           GSlice::NE_E, GSlice::E_SE, 
                                           GSlice::SE_S, GSlice::S_SW, 
                                           GSlice::SW_W, GSlice::W_NW};

    const std::map<GEdge, IncInfo2> mapIncEdge = {{GEdge::N, {-1, 0}},
                                                 {GEdge::NE, {-1, 1}},
                                                 {GEdge::E, {0, 1}},
                                                 {GEdge::SE, {1, 1}},
                                                 {GEdge::S, {1, 0}},
                                                 {GEdge::SW, {1, -1}},
                                                 {GEdge::W, {0, -1}},
                                                 {GEdge::NW, {-1, -1}}};
    
    const std::map<GSlice, IncInfo4> mapIncSlice = {{GSlice::NW_N, {1, 1, 1, 0}},
                                                   {GSlice::N_NE, {1, 0, 1, -1}},
                                                   {GSlice::NE_E, {1, -1, 0, -1}},
                                                   {GSlice::E_SE, {0, -1, -1, -1}},
                                                   {GSlice::SE_S, {-1, -1, -1, 0}},
                                                   {GSlice::S_SW, {-1, 0, -1, 1}},
                                                   {GSlice::SW_W, {-1, 1, 0, 1}},
                                                   {GSlice::W_NW, {0, 1, 1, 1}}};

};


#endif
