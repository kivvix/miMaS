#include "field.h"
#include <boost/numeric/ublas/matrix_proxy.hpp>

int main(int,char**)
{
  using namespace boost::numeric;
  field<double> m (16,8);

  double x = 0.1;
  for ( auto i =0 ; i<m.size1() ; ++i )
    for ( auto j =0 ; j<m.size2() ; ++j )
    {
      m(i,j) = i+0.01*(double)j;
    }
  std::cout << m << std::endl;
  std::cout << m.density() << std::endl;

  std::cout << m.step.dv << std::endl;
  std::cout << m.flux()[0] << std::endl;
  std::cout << m.flux()[1] << std::endl;
  std::cout << m.flux()[2] << std::endl;
  std::cout <<  m.stencil<direction::v>(0,0) << std::endl;
  auto tmp = m.stencil<direction::v>(0,0);
  for (ublas::matrix_vector_slice_periodic<ublas::matrix<double>,direction::v>::const_iterator it=tmp.cbegin();it<tmp.cend();++it)
    std::cout << *it << std::endl;
  std::cout << "-- " << tmp.size() <<std::endl;

  for (auto k=0;k<tmp.size();++k)
  {std::cout << k << " " ; std::cout<< tmp[k] << " " << tmp(k) << std::endl;}
  return 0;
}

