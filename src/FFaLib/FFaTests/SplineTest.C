// $Id$

#include "FFaLib/FFaAlgebra/FFaSpline.H"
#include <fstream>


int spline_test (const std::string& fname, size_t nSample,
                 std::vector<FaVec3>* samples)
{
  std::ifstream inp(fname,std::ios::in);
  std::vector<FaVec3> points;

  FaVec3 XYZ;
  inp >> XYZ;
  while (inp.good())
  {
    points.push_back(XYZ);
    inp >> XYZ;
  }
  inp.close();
  std::cout <<"#Read "<< points.size() <<" points from "<< fname << std::endl;
  if (points.size() < 2) return -1;

  FFaSpline spline(points);
  double t0 = spline.start();
  double t1 = spline.stop();
  double dt = (t1-t0)/(nSample-1);
  for (double t = t0; t <= t1; t += dt)
    if (samples)
      samples->push_back(spline(t));
    else
      std::cout << spline(t) << std::endl;

  return samples ? samples->size() : 0;
}


int SplineTest (const char* fname, size_t nSample)
{
  return spline_test(fname,nSample,NULL);
}
