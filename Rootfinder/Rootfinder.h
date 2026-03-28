#pragma once
#include <functional>

class Rootfinder
{
public:
  Rootfinder(double tolerance=1e-6, bool verbose=false);
  double Root(std::function<double(double)> f);
  double Bisection_method(std::function<double(double)> f, double higher, double lower,int num_iter=120);
  double Newton_Raphson(std::function<double(double)> f,double x0=0.0,int num_iter=40);
  double Secant_Method(std::function<double(double)> f,double x0=0.0,int num_iter=100);
  double IQI(std::function<double(double)> f,double x0=0.0,int num_iter=100);
private:
  double FDMderivative(std::function<double(double)> f,double x, double h=1e-6);
  double tol;
  bool verbose;
};