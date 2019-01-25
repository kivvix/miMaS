#ifndef _FIELD_H_
#define _FIELD_H_

#include <array>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multi_array.hpp>

#include "direction.h"

using namespace boost::numeric;
namespace math = boost::math::constants;

template <typename _T,std::size_t NumDims>
class field
{
  public:
    typedef typename boost::multi_array<_T,NumDims>::value_type             value_type;
    typedef typename boost::multi_array<_T,NumDims>::reference              reference;
    typedef typename boost::multi_array<_T,NumDims>::const_reference        const_reference;
    typedef typename boost::multi_array<_T,NumDims>::iterator               iterator;
    typedef typename boost::multi_array<_T,NumDims>::const_iterator         const_iterator;
    typedef typename boost::multi_array<_T,NumDims>::reverse_iterator       reverse_iterator;
    typedef typename boost::multi_array<_T,NumDims>::element                element;
    typedef typename boost::multi_array<_T,NumDims>::size_type              size_type; 
    typedef typename boost::multi_array<_T,NumDims>::difference_type        difference_type; 
    typedef typename boost::multi_array<_T,NumDims>::index                  index; 
    typedef typename boost::multi_array<_T,NumDims>::extent_range           extent_range;
    typedef typename boost::multi_array<_T,NumDims>::const_reverse_iterator const_reverse_iterator;
    
    template <std::size_t NDims>
    struct const_array_view {
      typedef boost::detail::multi_array::const_multi_array_view<_T,NDims> type;
    };
    template <std::size_t NDims>
    struct array_view {
      typedef boost::detail::multi_array::multi_array_view<_T,NDims> type;
    };

  protected:
    boost::multi_array<_T,NumDims> m_data;

  public:
    struct step_t {
      _T dx=0.1 , dv=0.1;
    } step;
    struct range_t {
      _T x_min=0.  , x_max=1.;
      _T v_min=-1. , v_max=1.;
    } range;


    field ()
      : m_data() , step() , range()
    {}

    field ( field const& rhs )
      : m_data(rhs.m_data) , step(rhs.step) , range(rhs.range)
    {}

    explicit
    field ( boost::detail::multi_array::extent_gen<NumDims> const& ranges )
      : m_data(ranges) , step() , range()
    {}

    template < typename ExtentList >
    field ( ExtentList const& sizes )
      : m_data(sizes) , step() , range()
    {}
    template < typename ExtentList >
    field ( ExtentList const& sizes , step_t s , range_t r )
      : m_data(sizes) , step(s) , range(r)
    {}

    ~field ()
    {}

    field &
    operator = ( field const& rhs )
    {
      if ( this != &rhs ) {
        m_data = rhs.m_data;
        step  = rhs.step;
        range = rhs.range;
      }
      return *this;
    }

    inline size_type
    size ( size_type n ) const
    { return m_data.shape()[n]; }

    const auto
    operator[] ( size_type i ) const
    { return m_data[i]; }
    auto
    operator[] ( size_type i )
    { return m_data[i]; }

    auto
    shape () const
    { return m_data.shape(); }

    template < direction::direction D >
    typename std::enable_if< D == direction::x ,
             typename const_array_view<1>::type >::type
    stencil ( size_type i , size_type k ) const
    { return m_data[ boost::indices[typename boost::multi_array<_T,NumDims>::index_range(i-2,i+3)][k] ]; }
    template < direction::direction D >
    typename std::enable_if< D == direction::x ,
             typename array_view<1>::type >::type
    stencil ( size_type i , size_type k )
    { return m_data[ boost::indices[typename boost::multi_array<_T,NumDims>::index_range(i-2,i+3)][k] ]; }

    template < direction::direction D >
    typename std::enable_if< D == direction::v ,
             typename const_array_view<1>::type >::type
    stencil ( size_type i , size_type k ) const
    { return m_data[ boost::indices[i][typename boost::multi_array<_T,NumDims>::index_range(k-2,k+3)] ]; }
    template < direction::direction D >
    typename std::enable_if< D == direction::v ,
             typename array_view<1>::type >::type
    stencil ( size_type i , size_type k )
    { return m_data[ boost::indices[i][typename boost::multi_array<_T,NumDims>::index_range(k-2,k+3)] ]; }

    // stencil sur le bord
    template < direction::direction D >
    typename std::enable_if< D == direction::x ,
             std::array<_T,5> >::type
    stencil_border ( size_type i , size_type k ) const
    {
      size_type im2 = (i-2+m_data.shape()[0])%(m_data.shape()[0]);
      size_type im1 = (i-1+m_data.shape()[0])%(m_data.shape()[0]);
      size_type ip1 = (i+1)%(m_data.shape()[0]);
      size_type ip2 = (i+2)%(m_data.shape()[0]);
      i = i%(m_data.shape()[0]); // si appel pour stencil en $i+1$ avec `i=size1-1`
      
      return { m_data[im2][k] , m_data[im1][k] , m_data[i][k] , m_data[ip1][k] , m_data[ip2][k] };
    }
    template < direction::direction D >
    typename std::enable_if< D == direction::v ,
             std::array<_T,5> >::type
    stencil_border ( size_type i , size_type k ) const
    {
      size_type km2 = (k-2+m_data.shape()[1])%(m_data.shape()[1]);
      size_type km1 = (k-1+m_data.shape()[1])%(m_data.shape()[1]);
      size_type kp1 = (k+1)%(m_data.shape()[1]);
      size_type kp2 = (k+2)%(m_data.shape()[1]);
      k = k%(m_data.shape()[1]); // si appel pour stencil en $i+1$ avec `i=size1-1`
      
      return { m_data[i][km2] , m_data[i][km1] , m_data[i][k] , m_data[i][kp1] , m_data[i][kp2] };
    }


    auto
    density () const
    {
      ublas::vector<_T> rho(size(0),0);
      size_type i=0;
      for ( auto it=m_data.begin() ; it!=m_data.end() ; ++it , ++i ) {
        rho(i) = std::accumulate( it->begin() , it->end() , _T(0.) );
      }
      return rho;
    }

    auto
    moments () const
    {
      _T vk;
      std::array<ublas::vector<_T>,3> U = { ublas::vector<_T>(size(0),0.) ,
                                            ublas::vector<_T>(size(0),0.) ,
                                            ublas::vector<_T>(size(0),0.) };
      for ( size_type i=0 ; i<size(0) ; ++i ) {
        for (size_type k=0 ; k<size(1) ; ++k ) {
          vk = k*step.dv+range.v_min;
          U[0][i] += m_data[i][k];
          U[1][i] += vk*m_data[i][k];
          U[2][i] += 0.5*vk*vk*m_data[i][k];
        }
      }
      return U;
    }

    auto
    flux () const
    {
      _T vk;
      std::array<ublas::vector<_T>,3> U = { ublas::vector<_T>(size(0),0.) ,
                                            ublas::vector<_T>(size(0),0.) ,
                                            ublas::vector<_T>(size(0),0.) };
      for ( size_type i=0 ; i<size(0) ; ++i ) {
        for (size_type k=0 ; k<size(1) ; ++k ) {
          vk = k*step.dv+range.v_min;
          U[0][i] += vk*m_data[i][k];
          U[1][i] += vk*vk*m_data[i][k];
          U[2][i] += 0.5*vk*vk*vk*m_data[i][k];
        }
      }
      return U;
    }
    
};

#endif

