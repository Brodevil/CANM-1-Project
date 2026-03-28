#pragma once
#include <functional>

class Rootfinder
{
public:
  Rootfinder(double tolerance=1e-6, bool verbose=false);
  double Root(std::function<double(double)> f);
  double Bisection_method(std::function<double(double)> f, double higher, double lower);
  double Newton_Raphson(std::function<double(double)> f,double x0=0.0);
  double Differential_evolution(std::function<double(double)> f);
private:
  double FDMderivative(std::function<double(double)> f,double x, double h=1e-6);
  double tol;
  bool verbose;
};