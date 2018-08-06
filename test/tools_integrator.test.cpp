#include <iostream>
#include <nuSQuIDS/tools.h>

int main (int argc, char const *argv[])
{
  using namespace nusquids;

  double integral = integrate([](double x){ return x*x; }, 0.0, 2.0);
  if( fabs(integral - (2.0*2.0*2.0)/3.) > 1.e-5 ) std::cout << "Integral tolerance exceed" << std::endl;

  return 0;
}
