#include "field.h"

int main(int,char**)
{
  using namespace boost::numeric;
  field<double> m (4,8);

  double x = 0.1;
  for ( auto i =0 ; i<m.size1() ; ++i )
    for ( auto j =0 ; j<m.size2() ; ++j )
    {
      m(i,j) = i+1+0.1*(double)j;
    }
  std::cout << m << std::endl;
  std::cout << m.density() << std::endl;

  std::cout << m.step.dv << std::endl;
  return 0;
}
