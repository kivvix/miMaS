#include "field.h"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "weno.h"

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
  std::cout <<  m.stencil<direction::x>(0,0) << std::endl;
  auto tmp = m.stencil<direction::x>(0,0);
  for (ublas::matrix_vector_slice_periodic<ublas::matrix<double>,direction::x>::const_iterator it=tmp.cbegin();it<tmp.cend();++it)
    std::cout << *it << std::endl;
  std::cout << "-- " << tmp.size() <<std::endl;

  for (auto k=0;k<tmp.size();++k)
  {std::cout << k << " " ; std::cout<< tmp[k] << " " << tmp(k) << std::endl;}


  field<double> c(100,100);
  c.range.x_min=0; c.range.x_max=2.*boost::math::constants::pi<double>();
  c.range.v_min=0; c.range.v_max=2.*boost::math::constants::pi<double>();
  c.step.dx = (c.range.x_max-c.range.x_min)/c.size1();
  c.step.dv = (c.range.v_max-c.range.v_min)/c.size2();

  for(auto i=0;i<c.size1();++i)
    for(auto k=0;k<c.size2();++k)
      c(i,k) = cos(i*c.step.dx);

  ublas::vector<double> v(c.size2());
  ublas::vector<double> E(c.size1());

  for (auto k=0;k<v.size();++k) { v(k) = 1.; }
  for (auto i=0;i<E.size();++i) { E(i) = 0.; }

  weno::trp2D(c,v,E);
  return 0;
}

