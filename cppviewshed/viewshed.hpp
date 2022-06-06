#ifndef VIEWSHED_H
#define VIEWSHED_H

#include <iostream>
#include <set>
#include <map>
#include <vector>


class ViewShed 
{
  public:
    float VPRow;
    float VPCol;
    
    ViewShed(float vpRow, float vpCol);
    ~ViewShed();

    void hello(void);

    
    void bootStrap(std::vector<std::vector<float>> dem);
    void refreshAllGrids(std::vector<std::vector<float>> dem);
    
    void refreshAG(std::vector<std::vector<float>> dem);
    void refreshVS(void);
    void refreshVV(void);
    void refreshTG(void);

    void printAuxGrid(void);
    void printVizViews(void);
    void printVizScore(void);
    void printTrackerGrid(void);
    
    void processAllEdges(float vpRow, float vpCol);
    void processAllSlices(float vpRow, float vpCol);  

    //void privStructs(void);
  
  private:
    
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
    float MaxRow;
    float MaxCol;
 
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
    bool inGrid(float cRow, float cCol);
    std::vector<XYZ> getSliceIndices(GSlice slc, float iRow, float iCol);
    void processEdge(GEdge se, float vpRow, float vpCol);
    void processSlice(GSlice slc, float vpRow, float vpCol);

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
