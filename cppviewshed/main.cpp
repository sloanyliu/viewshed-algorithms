#include <iostream>
#include "viewshed.hpp"


int main(int argc, char* argv[])
{
 /*
  XYZ testDot;
  testDot.X = testDot.Y = testDot.Z = 4;
  std::cout << testDot.X << std::endl;
  
  ViewShed dummy = ViewShed(1, 1); 
  dummy.hello();
  
  XYZ testDot1 = {1, 2, 3};
  XYZ testDot2 = {3, 2, 1};
  XYZ testDot3 = dummy.vecAdd(testDot2, testDot1); 
  dummy.printXYZ(testDot3);

  //std::cout << "<" << testDot3.X << ", " 
  //         << testDot3.Y << ", "
  //          << testDot3.Z << ">" << std::endl;
  
  XYZ testDot4 = dummy.vecSub(testDot3, testDot1); 
  dummy.printXYZ(testDot4);

  float testAng1 = dummy.vecAngle(testDot1, testDot2);
  std::cout << "testAng1: " << testAng1 << std::endl;

  float testAng2 = dummy.vecAngle(testDot3, testDot1);
  std::cout << "testAng2: " << testAng2 << std::endl;

  XYZ ta1 = {1, 1, 1};
  XYZ ta2 = {3, 2, 7};
  float testAng3 = dummy.vecAngle(ta1, ta2);
  std::cout << "testAng3: " << testAng3 << std::endl;
  
  ta1 = {1, 0, 0};
  ta2 = {0, 1, 0};
  float testAng4 = dummy.vecAngle(ta1, ta2);
  std::cout << "testAng4: " << testAng4 << std::endl;
*/
  //=====================
  std::vector<std::vector<float>> tdem1 = {{1, 2, 3, 4}, 
                                           {1, 2, 3, 4}, 
                                           {1, 2, 3, 4}, 
                                           {1, 2, 3, 4}};
  
  std::vector<std::vector<float>> tdem2 = {{1, 1, 1, 1, 1}, 
                                           {1, 9, 9, 9, 1}, 
                                           {1, 9, 1, 9, 1}, 
                                           {1, 9, 9, 9, 1},
                                           {1, 1, 1, 1, 1}};
  
  std::vector<std::vector<float>> tdem3 = {{1, 2, 3, 4, 5}, 
                                           {1, 2, 3, 4, 5}, 
                                           {1, 2, 3, 4, 5}, 
                                           {1, 2, 3, 4, 5},
                                           {1, 2, 3, 4, 5}};

  std::vector<std::vector<float>> testerDEM = tdem2;

  //ViewShed tester = ViewShed(testerDEM.size() / 2, testerDEM.size() / 2);
  ViewShed tester = ViewShed(0, tdem3.size() - 1);
  tester.bootStrap(testerDEM); 
  tester.processAllEdges(1, 0);
  tester.processAllSlices(1, 0);
  tester.printAuxGrid();
  tester.printVizViews();
  tester.printVizScore();
  tester.printTrackerGrid();
  //tester.processAllEdges(tester.VPRow, tester.VPCol);
  //tester.processAllSlices(tester.VPRow, tester.VPCol);
  
  /*
  tester.processAllEdges(2, 2);
  tester.processAllSlices(2, 2);
  
  */
  // FIXME: function below is broken
  //tester.processAllPoints();

  return 0;
}


