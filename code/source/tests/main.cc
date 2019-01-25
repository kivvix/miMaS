#include <algorithm>
#include <iostream>
#include <iterator>

#include <fstream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multi_array.hpp>
/*
namespace direction { enum direction { x,v }; }

using namespace boost::numeric;
namespace math = boost::math::constants;

template <typename _T,std::size_t NumDims>
class field
  : public boost::multi_array<_T,NumDims>
{
  typedef boost::multi_array<_T,NumDims> mom_type;
public:
  typedef typename mom_type::value_type             value_type;
  typedef typename mom_type::reference              reference;
  typedef typename mom_type::const_reference        const_reference;
  typedef typename mom_type::iterator               iterator;
  typedef typename mom_type::const_iterator         const_iterator;
  typedef typename mom_type::reverse_iterator       reverse_iterator;
  typedef typename mom_type::element                element;
  typedef typename mom_type::size_type              size_type; 
  typedef typename mom_type::difference_type        difference_type; 
  typedef typename mom_type::index                  index; 
  typedef typename mom_type::extent_range           extent_range;
  typedef typename mom_type::const_reverse_iterator const_reverse_iterator;
  
  template <std::size_t NDims>
  struct const_array_view {
    typedef boost::detail::multi_array::const_multi_array_view<_T,NDims> type;
  };
  template <std::size_t NDims>
  struct array_view {
    typedef boost::detail::multi_array::multi_array_view<_T,NDims> type;
  };

  struct step {
    _T dx,dv;
  } step;
  struct range {
    _T x_min,x_max;
    _T v_min,v_max;
  } range;

  explicit field ()
    : mom_type()
  {}

  explicit field ( const boost::detail::multi_array::extent_gen<NumDims> & ranges )
    : mom_type(ranges)
  {}

  field ( const field & rhs )
    : mom_type(rhs)
  {}

  template < direction::direction D >
  typename std::enable_if< D == direction::x ,
           typename const_array_view<1>::type >::type
  stencil ( size_type i , size_type k ) const
  { return this->operator[]( boost::indices[mom_type::index_range(i-2,i+2)][k] ); }
  template < direction::direction D >
  typename std::enable_if< D == direction::x ,
           typename array_view<1>::type >::type
  stencil ( size_type i , size_type k )
  { return this->operator[]( boost::indices[mom_type::index_range(i-2,i+2)][k] ); }

  template < direction::direction D >
  typename std::enable_if< D == direction::v ,
           typename const_array_view<1>::type >::type
  stencil ( size_type i , size_type k ) const
  { return this->operator[]( boost::indices[i][mom_type::index_range(k-2,k+2)] ); }
  template < direction::direction D >
  typename std::enable_if< D == direction::v ,
           typename array_view<1>::type >::type
  stencil ( size_type i , size_type k )
  { return this->operator[]( boost::indices[i][mom_type::index_range(k-2,k+2)] ); }

  auto
  density () const
  {
    ublas::vector<_T> rho(mom_type::shape()[0],0);
    size_type i=0;
    for ( auto it=mom_type::begin() ; it!=mom_type::end() ; ++it , ++i ) {
      rho(i) = std::accumulate( it->begin() , it->end() , _T(0.0) );
    }
    return rho;
  }

  auto
  moments () const
  {
    _T vk;
    std::array<ublas::vector<_T>,3> U = { ublas::vector<_T>(mom_type::shape()[0],0) ,
                                          ublas::vector<_T>(mom_type::shape()[0],0) ,
                                          ublas::vector<_T>(mom_type::shape()[0],0) };
    
    for ( size_type i=0 ; i<mom_type::shape()[0] ; ++i ) {
      for ( size_type k=0 ; k<mom_type::shape()[1] ; ++k ) {
        vk = k*step.dv+range.v_min;
        U[0][i] += this->operator[](i)[k];
        U[1][i] += vk*this->operator[](i)[k];
        U[2][i] += 0.5*vk*vk*this->operator[](i)[k];
      }
    }

    return U;
  }

  auto
  flux () const
  {
    _T vk;
    std::array<std::vector<_T>,3> U = { ublas::vector<_T>(mom_type::shape()[0],0) ,
                                        ublas::vector<_T>(mom_type::shape()[0],0) ,
                                        ublas::vector<_T>(mom_type::shape()[0],0) };
    
    for ( size_type i=0 ; i<mom_type::shape()[0] ; ++i ) {
      for ( size_type k=0 ; k<mom_type::shape()[1] ; ++k ) {
        vk = k*step.dv+range.v_min;
        U[0][i] += vk*this->operator[](i)[k];
        U[1][i] += vk*vk*this->operator[](i)[k];
        U[2][i] += 0.5*vk*vk*vk*this->operator[](i)[k];
      }
    }

    return U;
  }

};

template <typename _T>
struct variabili
{
  ublas::vector<_T> rho;
  ublas::vector<_T> u;
  ublas::vector<_T> T;

  variabili ( std::array<ublas::vector<_T>,3> && U )
    : rho(std::move(U[0])) , u(std::move(U[1])) , T(U[2].size())
  {
    for ( std::size_t i=0 ; i<rho.size() ; ++i ) {
      u(i) /= rho(i);
      T(i)  = (U[2](i) - U[1](i)*U[1](i)/(2.*rho(i)))/(0.5*rho(i));
    }
  }
};

template <typename _T,std::size_t NumDims>
struct maxwellian
{
  maxwellian ( typename field<_T,NumDims>::step s , typename field<_T,NumDims>::range r )
    : step(s) , range(r) , v((r.v_max-r.v_min)/s.dv ) ,
      data(boost::extents[(r.x_max-r.x_min)/s.dx][(r.v_max-r.v_min)/s.dv])
  {
    for ( auto i=0 ; i<v.size() ; ++i )
      { v(i) = i*step.dv + range.v_min; }
  }

  void
  update ( variabili<_T> const& var )
  { update(var.rho,var.u,var.T); }

  void
  update ( ublas::vector<_T> const& rho , ublas::vector<_T> const& u , ublas::vector<_T> const& T )
  {
    for ( std::size_t i=0 ; i<data.shape()[0] ; ++i ) {
      for ( std::size_t k=0 ; k<data.shape()[1] ; ++k ) {
        data[i][k] = rho(i)/std::sqrt(2.*math::pi<_T>*T(i))*std::exp(-0.5*(v(k)-u(i))*(v(k)-u(i))/T(i));
      }
    }
  }
  
  _T const&
  operator () ( std::size_t i , std::size_t k ) const
  { return data[i][k]; }
  _T &
  operator () ( std::size_t i , std::size_t k )
  { return data[i][k]; }

  typename field<_T,NumDims>::step step;
  typename field<_T,NumDims>::range range;
  ublas::vector<_T> v;
  boost::multi_array<_T,NumDims> data;
};

namespace weno {

using namespace boost::numeric;
#define SQ(X) ((X)*(X))

_T
local_flux_p ( const _T * f_loc )
{
  static const _T epso = 1e-6;
  _T w0 = (13./12.)*SQ( f_loc[0] - 2.*f_loc[1] + f_loc[2] ) + 0.25*SQ(    f_loc[0] - 4.*f_loc[1] + 3.*f_loc[2] );
  _T w1 = (13./12.)*SQ( f_loc[1] - 2.*f_loc[2] + f_loc[3] ) + 0.25*SQ(    f_loc[1]               -    f_loc[3] );
  _T w2 = (13./12.)*SQ( f_loc[2] - 2.*f_loc[3] + f_loc[4] ) + 0.25*SQ( 3.*f_loc[2] - 4.*f_loc[3] +    f_loc[4] );

  w0 = 0.1/SQ(w0+epsi);
  w1 = 0.6/SQ(w1+epsi);
  w2 = 0.3/SQ(w2+epsi);

  sum_w = w0+w1+w2;
  w0 /= sum_w;
  w1 /= sum_w;
  w2 /= sum_w;

  return w0*( (2./6.)*f_loc[0] - (7./6.)*f_loc[1] + (11./6.)*f_loc[2] )
       + w1*(-(1./6.)*f_loc[1] + (5./6.)*f_loc[2] +  (2./6.)*f_loc[3] )
       + w2*( (2./6.)*f_loc[2] + (5./6.)*f_loc[3] +  (1./6.)*f_loc[4] );
}

_T
local_flux_m ( const _T * f_loc )
{
  static const _T epso = 1e-6;
  _T w0 = (13./12.)*SQ( f_loc[2] - 2.*f_loc[3] + f_loc[4] ) + 0.25*SQ( 3.*f_loc[2] - 4.*f_loc[3] +    f_loc[4] );
  _T w1 = (13./12.)*SQ( f_loc[1] - 2.*f_loc[2] + f_loc[3] ) + 0.25*SQ(    f_loc[1]               -    f_loc[3] );
  _T w2 = (13./12.)*SQ( f_loc[0] - 2.*f_loc[1] + f_loc[2] ) + 0.25*SQ(    f_loc[0] - 4.*f_loc[1] + 3.*f_loc[2] );

  w0 = 0.1/SQ(w0+epsi);
  w1 = 0.6/SQ(w1+epsi);
  w2 = 0.3/SQ(w2+epsi);

  sum_w = w0+w1+w2;
  w0 /= sum_w;
  w1 /= sum_w;
  w2 /= sum_w;

  return w2*(-(1./6.)*f_loc[0] + (5./6.)*f_loc[1] + (2./6.)*f_loc[2] )
       + w1*( (2./6.)*f_loc[1] + (5./6.)*f_loc[2] - (1./6.)*f_loc[3] )
       + w2*((11./6.)*f_loc[2] - (7./6.)*f_loc[3] + (2./6.)*f_loc[4] );
}

#undef SQ

template < direction::direction D> inline int& index_velocity ( int &i , int &k );
template <> inline int& index_velocity<direction::x> (int &i , int &k ) { return k; }
template <> inline int& index_velocity<direction::v> (int &i , int &k ) { return i; }

template < direction::direction D > inline int index_i_m ( int i );
template <> inline int index_i_m<direction::x> ( int i ) { return i+1; }
template <> inline int index_i_m<direction::v> ( int i ) { return i;   }
template < direction::direction D > inline int index_k_m ( int k );
template <> inline int index_k_m<direction::x> ( int k ) { return k;   }
template <> inline int index_k_m<direction::v> ( int k ) { return k+1; }

template < typename _T >
struct array_view
{
  typedef _T                                    value_type;
  typedef _T *                                  iterator;
  typedef const _T *                            const_iterator;
  typedef std::reverse_iterator<iterator>       reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
  typedef _T &                                  reference;
  typedef const _T &                            const_reference;
  typedef _T *                                  pointer;
  typedef std::ptrdiff_t                        difference_type;
  typedef std::size_t                           size_type;

  const_iterator begin () const { return data;      }
  iterator       begin ()       { return data;      }
  const_iterator end   () const { return data+size; }
  iterator       end   ()       { return data+size; }

  const_reverse_iterator rbegin () const { return std::make_reverse_iterator(end());   }
  reverse_iterator       rbegin ()       { return std::make_reverse_iterator(end());   }
  const_reverse_iterator rend   () const { return std::make_reverse_iterator(begin()); }
  reverse_iterator       rend   ()       { return std::make_reverse_iterator(begin()); }

  _T * data;
  std::size_t size;
};

template < typename _T , std::size_t NumDims , direction::direction D >
auto
flux_p ( field<_T,NumDims> const& u , ublas::vector<_T> const& v )
{
  boost::multi_array<_T,NumDims> fip12(array_view<_T>(u.shape(),NumDims)); // TODO: regarder les constructeur de multi_array (ou plus exactement de boost::extents) pour pouvoir copier la shape d'un autre multi_array

  for ( std::size_t i=0 ; i<u.shape()[0] ; ++i ) {
    for ( std::size_t k=0 ; k<u.shape()[1] ; ++j ) {
      fip12[i][k] = v(index_velocity<D>(i,k))*local_flux_p( u.field<_T,NumDims>::template stencil<D>(i,k) );
    }
  }

  return fip12;
}

template < typename _T , std::size_t NumDims , direction::direction D >
auto
flux_m ( field<_T,NumDims> const& u , ublas::vector<_T> const& v )
{
  boost::multi_array<_T,NumDims> fim12(u.shape()); // TODO: regarder les constructeur de multi_array (ou plus exactement de boost::extents) pour pouvoir copier la shape d'un autre multi_array

  for ( std::size_t i=0 ; i<u.shape()[0] ; ++i ) {
    for ( std::size_t k=0 ; k<u.shape()[1] ; ++j ) {
      fim12[i][k] = v(index_velocity<D>(i,k))*local_flux_m( u.field<_T,NumDims>::template stencil<D>(index_i_m<D>(i),index_k_m<D>(k)) );
    }
  }

  return fim12;
}

template < typename _T , std::size_t NumDims >
auto
trp2D ( field<_T,NumDims> const& u , ublas::vector<_T> const& v , ublas::vector<_T> const& E )
{
  
}

} // namespace weno
*/

#include "physics.h"
#include "direction.h"
#include "field.h"
#include "weno.h"

/*
int main(int,char**)
{
  field<double,2> f(boost::extents[10][10]);
  for (field<double,2>::size_type i=0 ; i<f.size<0>() ; ++i ) {
    for (field<double,2>::size_type k=0 ; k<f.size<1>() ; ++k ) {
      f[i][k] = i*k;
    }
  }

  auto rho = f.density();
  std::cout << rho << std::endl;
  auto U = f.moments();
  std::cout << U[0] << std::endl;
  std::cout << U[1] << std::endl;
  std::cout << U[2] << std::endl;
  U = f.flux();
  std::cout << U[0] << std::endl;
  std::cout << U[1] << std::endl;
  std::cout << U[2] << std::endl;
  
  return 0;
}
*/

int main(int,char**)
{
	field<double,2> f(boost::extents[51][51]);
	f.range.x_min = 0.;  f.range.x_max = 1.;
	f.range.v_min = -1.; f.range.v_max = 1.;
	f.step.dx = (f.range.x_max - f.range.x_min)/50. ; f.step.dv = (f.range.v_max - f.range.v_min)/50.;
	for (field<double,2>::size_type i=0 ; i<f.size(0) ; ++i) {
    for (field<double,2>::size_type k=0 ; k<f.size(1) ; ++k ) {
    	//std::cout << i << " " << k << "\r";
      f[i][k] = cos(2.*boost::math::constants::pi<double>()*f.step.dv*k);
    }
	}
	std::ofstream of_init; of_init.open("init.dat");
	for (field<double,2>::size_type i=0 ; i<f.size(0) ; ++i) {
    for (field<double,2>::size_type k=0 ; k<f.size(1) ; ++k ) {
    	of_init << f.step.dx*(double)i+f.range.x_min << " " << f.step.dv*(double)k+f.range.v_min << " " << f[i][k] << "\n";
    }
    of_init << std::endl;
	}
	of_init.close();

	ublas::vector<double> v(51),E(51);
	std::generate(v.begin(),v.end(),[](){return 0.;});
	std::generate(E.begin(),E.end(),[](){return -1.;});

	field<double,2> f1(boost::extents[51][51]);
	field<double,2> f2(boost::extents[51][51]);
	f1.step = f.step; f2.step = f.step;
	double dt = 0.5*f.step.dv;

	//std::cout << f.step.dx << " " << f.step.dv << " " << f.size(1) << " " << f.size(2) << std::endl;

	for ( auto t=0 ; t<1 ; ++t ) {
		//std::cout << t << std::flush << "\r";
		auto L = weno::trp2D(f,v,E);
		for (field<double,2>::size_type i=0 ; i<f.size(0) ; ++i) {
	    for (field<double,2>::size_type k=0 ; k<f.size(1) ; ++k ) {
	    	f1[i][k] = f[i][k] - dt*L[i][k];
	    }
	  }
	  L = weno::trp2D(f1,v,E);
		for (field<double,2>::size_type i=0 ; i<f.size(0) ; ++i) {
	    for (field<double,2>::size_type k=0 ; k<f.size(1) ; ++k ) {
	    	f2[i][k] = 0.75*f[i][k] + 0.25*f1[i][k] - 0.25*dt*L[i][k];
	    }
	  }
	  L = weno::trp2D(f2,v,E);
		for (field<double,2>::size_type i=0 ; i<f.size(0) ; ++i) {
	    for (field<double,2>::size_type k=0 ; k<f.size(1) ; ++k ) {
	    	f[i][k] = (1./3.)*f[i][k] + (2./3.)*f2[i][k] - (2./3.)*dt*L[i][k];
	    }
	  }
	}
	//std::cout << std::endl;

	std::ofstream of_vp; of_vp.open("vp.dat");
	for (field<double,2>::size_type i=0 ; i<f.size(0) ; ++i) {
    for (field<double,2>::size_type k=0 ; k<f.size(1) ; ++k ) {
    	of_vp << f.step.dx*(double)i+f.range.x_min << " " << f.step.dv*(double)k+f.range.v_min << " " << f[i][k] << "\n";
    }
    of_vp << std::endl;
	}
	of_vp.close();

	return 0;
}
