#include <iostream>
#include <algorithm>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;
template <typename T>
struct field : public ublas::matrix<T>
{
  struct step
  {
    T dx,dv;
    step() = default;
    step(step&) = default;
    step(step&&) = default;
    step(T dx_,T dv_) : dx(dx_) , dv(dv_) {}
  } step;

  field ()
    : ublas::matrix<T>()
  {}
  
  field ( typename ublas::matrix<T>::size_type size1 , typename ublas::matrix<T>::size_type size2 )
    : ublas::matrix<T>(size1,size2) , step(1.,1.)
  {}

  field ( const field &m )
    : ublas::matrix<T>(m) , step(m.step)
  {}

  template<class AE>
  field ( const ublas::matrix_expression<AE> &ae )
    : ublas::matrix<T>(ae) , step(1.,1.)
  {}

  field ( typename ublas::matrix<T>::size_type size1 , typename ublas::matrix<T>::size_type size2 ,
      T dx , T dv )
    : ublas::matrix<T>(size1,size2) , step(dx,dv)
  {}

  auto
  density () const
  {
    static const ublas::scalar_vector<T> v(this->size2(),step.dv);
    return prod(*this,v);
  }
};

