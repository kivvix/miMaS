#ifndef _FIELD_H_
#define _FIELD_H_

#include <array>
#include <fstream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array.hpp>

#include <boost/iterator/zip_iterator.hpp>

using namespace boost::numeric;

/*
template < typename _T >
struct step
{
  _T dx = 0.1;
  _T dv = 0.1; // idéalement il faudrait remplacer ceci par un tableau de taille NumDimsV pour avoir un dv différents dans chaque direction
};
template < typename _T >
struct range
{
  _T x_min=0.  , x_max=1.;
  _T v_min=-1. , v_max=1.; // même remarque que pour step
};
*/

template < typename _T , std::size_t NumDimsV >
class field
{
  public:
    typedef typename boost::multi_array< _T , NumDimsV+1 >::value_type             value_type;
    typedef typename boost::multi_array< _T , NumDimsV+1 >::reference              reference;
    typedef typename boost::multi_array< _T , NumDimsV+1 >::const_reference        const_reference;
    typedef typename boost::multi_array< _T , NumDimsV+1 >::iterator               iterator;
    typedef typename boost::multi_array< _T , NumDimsV+1 >::const_iterator         const_iterator;
    typedef typename boost::multi_array< _T , NumDimsV+1 >::reverse_iterator       reverse_iterator;
    typedef typename boost::multi_array< _T , NumDimsV+1 >::element                element;
    typedef typename boost::multi_array< _T , NumDimsV+1 >::size_type              size_type; 
    typedef typename boost::multi_array< _T , NumDimsV+1 >::difference_type        difference_type; 
    typedef typename boost::multi_array< _T , NumDimsV+1 >::index                  index; 
    typedef typename boost::multi_array< _T , NumDimsV+1 >::extent_range           extent_range;
    typedef typename boost::multi_array< _T , NumDimsV+1 >::const_reverse_iterator const_reverse_iterator;

    typedef typename boost::multi_array< _T , NumDimsV+1 >::index_gen   index_gen;
    typedef typename boost::multi_array< _T , NumDimsV+1 >::index_range index_range;
    
    template <std::size_t NDims>
    struct const_array_view {
      typedef boost::detail::multi_array::const_multi_array_view<_T,NDims> type;
    };
    template <std::size_t NDims>
    struct array_view {
      typedef boost::detail::multi_array::multi_array_view<_T,NDims> type;
    };

  protected:
    // stored data
    boost::multi_array<_T,NumDimsV+1> m_data;

  public:
    struct step {
      _T dx=0.1 , dv=0.1; // default value
    } step;
    struct range {
      _T x_min=0.  , x_max=1.; // default value
      _T v_min=-1. , v_max=1.; // default value

      _T len_x () const { return x_max - x_min; }
      _T len_v () const { return v_max - v_min; }
    } range;

    //// CONSTRUCTORS ///////////////////////
    field ()
      : m_data()
    {}

    field ( field const& rhs )
      : m_data(rhs.m_data) , step(rhs.step) , range(rhs.range)
    {}

    explicit
    field ( boost::detail::multi_array::extent_gen<NumDimsV+1> const& ranges )
      : m_data(ranges)
    {}

    template < typename ExtentList >
    field ( ExtentList const& sizes )
      : m_data(sizes)
    {}

    template < typename ExtentList >
    field ( ExtentList const& sizes , struct step s , struct range r )
      : m_data(sizes) , step(s) , range(r)
    {}

    //// DESTRUCTORS ////////////////////////
    ~field()
    {}


    //// OPERATORS //////////////////////////
    field &
    operator = ( field const& rhs )
    {
      if ( this != &rhs ) {
        m_data = rhs.m_data;
        step = rhs.step;
        range = rhs.range;
      }
      return *this;
    }

    //// CAPACITY ///////////////////////////
    inline size_type
    size ( size_type n ) const
    { return m_data.shape()[n]; }

    inline size_type
    size_x () const
    { return m_data.shape()[NumDimsV]; }

    auto
    shape () const
    { return m_data.shape(); }

    //// ELEMENT ACCESS /////////////////////
    const auto
    operator [] ( size_type k ) const
    { return m_data[k]; }
    auto
    operator [] ( size_type k )
    { return m_data[k]; }

    template <int NDims>
    const auto
    operator [] ( boost::detail::multi_array::index_gen<NumDimsV+1,NDims> const& indices ) const
    { return m_data[indices]; }

    template <int NDims>
    auto
    operator [] ( boost::detail::multi_array::index_gen<NumDimsV+1,NDims> const& indices )
    { return m_data[indices]; }

    //// ITERATORS //////////////////////////
    const auto
    origin () const
    { return m_data.origin(); }
    auto
    origin ()
    { return m_data.origin(); }

    const auto
    num_elements () const
    { return m_data.num_elements(); }

    const auto
    begin () const
    { return m_data.begin(); }
    const auto
    cbegin () const
    { return m_data.begin(); }
    auto
    begin ()
    { return m_data.begin(); }

    const auto
    end () const
    { return m_data.end(); }
    const auto
    cend () const
    { return m_data.end(); }
    auto
    end ()
    { return m_data.end(); }

    const auto
    begin ( size_type k ) const
    { return m_data[k].begin(); }
    const auto
    cbegin ( size_type k ) const
    { return m_data[k].begin(); }
    auto
    begin ( size_type k )
    { return m_data[k].begin(); }

    const auto
    end ( size_type k ) const
    { return m_data[k].end(); }
    const auto
    cend ( size_type k ) const
    { return m_data[k].end(); }
    auto
    end ( size_type k )
    { return m_data[k].end(); }

    // FONCTIONS DÉPENDANT DE NumDimsV
    // ici NumDimsV = 1
    const auto
    begin_stencil ( size_type k ) const
    { return boost::make_zip_iterator(boost::make_tuple( m_data[k-2].begin() , m_data[k-1].begin() , m_data[k].begin() , m_data[k+1].begin() , m_data[k+2].begin() , m_data[k+3].begin() )); }
    const auto
    cbegin_stencil ( size_type k ) const
    { return boost::make_zip_iterator(boost::make_tuple( m_data[k-2].begin() , m_data[k-1].begin() , m_data[k].begin() , m_data[k+1].begin() , m_data[k+2].begin() , m_data[k+3].begin() )); }
    auto
    begin_stencil ( size_type k )
    { return boost::make_zip_iterator(boost::make_tuple( m_data[k-2].begin() , m_data[k-1].begin() , m_data[k].begin() , m_data[k+1].begin() , m_data[k+2].begin() , m_data[k+3].begin() )); }

    const auto
    end_stencil ( size_type k ) const
    { return boost::make_zip_iterator(boost::make_tuple( m_data[k-2].end() , m_data[k-1].end() , m_data[k].end() , m_data[k+1].end() , m_data[k+2].end() , m_data[k+3].end() )); }
    const auto
    cend_stencil ( size_type k ) const
    { return boost::make_zip_iterator(boost::make_tuple( m_data[k-2].end() , m_data[k-1].end() , m_data[k].end() , m_data[k+1].end() , m_data[k+2].end() , m_data[k+3].end() )); }
    auto
    end_stencil ( size_type k )
    { return boost::make_zip_iterator(boost::make_tuple( m_data[k-2].end() , m_data[k-1].end() , m_data[k].end() , m_data[k+1].end() , m_data[k+2].end() , m_data[k+3].end() )); }

    const auto
    begin_border_stencil ( size_type k ) const
    {
      std::size_t km2=(k-2+size(0))%size(0) , km1=(k-1+size(0))%size(0) , kp1=(k+1)%size(0) , kp2=(k+2)%size(0) , kp3=(k+3)%size(0);
      return boost::make_zip_iterator(boost::make_tuple( m_data[km2].begin() , m_data[km1].begin() , m_data[k].begin() , m_data[kp1].begin() , m_data[kp2].begin() , m_data[kp3].begin() ));
    }
    const auto
    cbegin_border_stencil ( size_type k ) const
    {
      std::size_t km2=(k-2+size(0))%size(0) , km1=(k-1+size(0))%size(0) , kp1=(k+1)%size(0) , kp2=(k+2)%size(0) , kp3=(k+3)%size(0);
      return boost::make_zip_iterator(boost::make_tuple( m_data[km2].begin() , m_data[km1].begin() , m_data[k].begin() , m_data[kp1].begin() , m_data[kp2].begin() , m_data[kp3].begin() ));
    }
    auto
    begin_border_stencil ( size_type k )
    {
      std::size_t km2=(k-2+size(0))%size(0) , km1=(k-1+size(0))%size(0) , kp1=(k+1)%size(0) , kp2=(k+2)%size(0) , kp3=(k+3)%size(0);
      return boost::make_zip_iterator(boost::make_tuple( m_data[km2].begin() , m_data[km1].begin() , m_data[k].begin() , m_data[kp1].begin() , m_data[kp2].begin() , m_data[kp3].begin() ));
    }

    const auto
    end_border_stencil ( size_type k ) const
    {
      std::size_t km2=(k-2+size(0))%size(0) , km1=(k-1+size(0))%size(0) , kp1=(k+1)%size(0) , kp2=(k+2)%size(0) , kp3=(k+3)%size(0);
      return boost::make_zip_iterator(boost::make_tuple( m_data[km2].end() , m_data[km1].end() , m_data[k].end() , m_data[kp1].end() , m_data[kp2].end() , m_data[kp3].end() ));
    }
    const auto
    cend_border_stencil ( size_type k ) const
    {
      std::size_t km2=(k-2+size(0))%size(0) , km1=(k-1+size(0))%size(0) , kp1=(k+1)%size(0) , kp2=(k+2)%size(0) , kp3=(k+3)%size(0);
      return boost::make_zip_iterator(boost::make_tuple( m_data[km2].end() , m_data[km1].end() , m_data[k].end() , m_data[kp1].end() , m_data[kp2].end() , m_data[kp3].end() ));
    }
    auto
    end_border_stencil ( size_type k )
    {
      std::size_t km2=(k-2+size(0))%size(0) , km1=(k-1+size(0))%size(0) , kp1=(k+1)%size(0) , kp2=(k+2)%size(0) , kp3=(k+3)%size(0);
      return boost::make_zip_iterator(boost::make_tuple( m_data[km2].end() , m_data[km1].end() , m_data[k].end() , m_data[kp1].end() , m_data[kp2].end() , m_data[kp3].end() ));
    }

    //// OTHER METHODS //////////////////////
    void
    compute_steps ()
    {
      step.dv = (range.v_max - range.v_min)/size(0);
      step.dx = (range.x_max - range.x_min)/size_x();
    }

    auto
    density () const
    {
      ublas::vector<_T> rho(size_x(),0.); // x direction is le last one
      for ( auto k=0 ; k<size(0) ; ++k ) {
        for ( auto i=0 ; i<size(NumDimsV) ; ++i ) {
          rho(i) += m_data[k][i]*step.dv;
        }
      }
      return rho;
    }

    auto
    courant () const
    {
      _T vk;
      ublas::vector<_T> J(size_x(),0.); // x direction is le last one
      for ( auto k=0 ; k<size(0) ; ++k ) {
        vk = k*step.dv+range.v_min;
        for ( auto i=0 ; i<size(NumDimsV) ; ++i ) {
          J(i) += vk*m_data[k][i]*step.dv;
        }
      }
      return J;
    }

    auto
    moments () const
    {
      _T vk;
      std::array<ublas::vector<_T>,3> U = { ublas::vector<_T>(size_x(),0.) ,
                                            ublas::vector<_T>(size_x(),0.) ,
                                            ublas::vector<_T>(size_x(),0.) };
      for ( size_type k=0 ; k<size(0) ; ++k ) {
        vk = k*step.dv+range.v_min;
        for ( size_type i=0 ; i<size_x() ; ++i ) {
          U[0][i] += m_data[k][i]*step.dv;
          U[1][i] += vk*m_data[k][i]*step.dv;
          U[1][i] += 0.5*vk*vk*m_data[k][i]*step.dv;
        }
      }
      return U;
    }

    auto
    flux () const
    {
      _T vk;
      std::array<ublas::vector<_T>,3> U = { ublas::vector<_T>(size_x(),0.) ,
                                            ublas::vector<_T>(size_x(),0.) ,
                                            ublas::vector<_T>(size_x(),0.) };
      for ( size_type k=0 ; k<size(0) ; ++k ) {
        vk = k*step.dv+range.v_min;
        for ( size_type i=0 ; i<size_x() ; ++i ) {
          U[0][i] += vk*m_data[k][i]*step.dv;
          U[1][i] += vk*vk*m_data[k][i]*step.dv;
          U[1][i] += 0.5*vk*vk*vk*m_data[k][i]*step.dv;
        }
      }
      return U;
    }
    
    void
    write ( std::string const& filename ) const
    {
      std::ofstream f(filename);
      for ( size_type k=0 ; k<size(0) ; ++k ) {
        for ( size_type i=0 ; i<size_x() ; ++i ) {
          f << i*step.dx+range.x_min << " " <<  k*step.dv+range.v_min << " " << m_data[k][i] << "\n";
        }
        f << std::endl;
      }
      f.close();
    }
};

template < typename _T , std::size_t NumDimsV >
std::ostream &
operator << ( std::ostream & os , field<_T,NumDimsV> const& f )
{
  for ( typename field<_T,NumDimsV>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for ( typename field<_T,NumDimsV>::size_type i=0 ; i<f.size_x() ; ++i ) {
      f << i*f.step.dx+f.range.x_min << " " <<  k*f.step.dv+f.range.v_min << " " << f[k][i] << "\n";
    }
    f << std::endl;
  }
  return os;
}

#define SQ(X) ((X)*(X))

template < typename _T , std::size_t NumDimsV >
_T
energy ( field<_T,NumDimsV> const& f , ublas::vector<_T> const& E )
{
  _T H = double{0.};

  for ( typename field<_T,NumDimsV>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for ( typename field<_T,NumDimsV>::size_type i=0 ; i<f.size_x() ; ++i ) {
      H += SQ( (_T(k)*f.step.dv+f.range.v_min) ) * f[k][i] * f.step.dv*f.step.dx;
    }
  }

  for ( auto it=E.begin() ; it!=E.end() ; ++it ) {
    H += SQ(*it)*f.step.dx;
  }

  return H;
}

template < typename _T , std::size_t NumDimsV >
_T
kinetic_energy ( field<_T,NumDimsV> const& f )
{
  _T Ec = double{0.};

  for ( typename field<_T,NumDimsV>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for ( typename field<_T,NumDimsV>::size_type i=0 ; i<f.size_x() ; ++i ) {
      Ec += SQ( (_T(k)*f.step.dv+f.range.v_min) ) * f[k][i] * f.step.dv*f.step.dx;
    }
  }

  return Ec;
}

#undef SQ

#endif
