#ifndef _MATRIX_PROXY_H_
#define _MATRIX_PROXY_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace direction { enum direction { x,v }; }

namespace boost { namespace numeric { namespace ublas {

  template <class M,direction::direction D>
  class matrix_vector_slice_periodic:
    public matrix_vector_slice<M>
  {
    public:
      BOOST_UBLAS_INLINE
      matrix_vector_slice_periodic ( typename ublas::matrix_vector_slice<M>::matrix_data &data, const typename ublas::matrix_vector_slice<M>::slice_type &s1, const typename ublas::matrix_vector_slice<M>::slice_type &s2)
       : matrix_vector_slice<M>(data,s1,s2)
      {}

#ifndef BOOST_UBLAS_PROXY_CONST_MEMBER
      BOOST_UBLAS_INLINE
      typename matrix_vector_slice<M>::const_reference
      operator () (typename matrix_vector_slice<M>::size_type i) const;
      BOOST_UBLAS_INLINE
      typename matrix_vector_slice<M>::reference
      operator () (typename matrix_vector_slice<M>::size_type i);
#else
      BOOST_UBLAS_INLINE
      typename matrix_vector_slice<M>::reference
      operator () (typename matrix_vector_slice<M>::size_type i) const;
#endif
  };

  template <class M>
  class matrix_vector_slice_periodic<M,direction::x>:
    public matrix_vector_slice<M>
  {
    public:
      BOOST_UBLAS_INLINE
      matrix_vector_slice_periodic ( M &data, const typename ublas::matrix_vector_slice<M>::slice_type &s1, const typename ublas::matrix_vector_slice<M>::slice_type &s2)
       : matrix_vector_slice<M>(data,s1,s2)
      {}

      //// CECI FONCTIONNE ENFIN UNIQUEMENT POUR L'OPÉRATEUR (), par contre je viens de me rendre compte que la majorité des autres fonctions liées au matrix_vector_slice utilisent les itérateurs qu'il peut y avoir dessus (itérateurs n'utilisant pas cet opérateur () mais leur propre accès aux données, ce qui est logique quand on y pense car on ne viendrait pas penser que quelqu'un souhaiterait accéder aux données autrement...)
#ifndef BOOST_UBLAS_PROXY_CONST_MEMBER
      BOOST_UBLAS_INLINE
      typename matrix_vector_slice<M>::const_reference
      operator () (typename matrix_vector_slice<M>::size_type i) const
      {
        return ublas::matrix_vector_slice<M>::data()( (this->start1()+i*this->stride1())%this->data().size1() , (this->start2()+i*this->stride2()) );
      }
      BOOST_UBLAS_INLINE
      typename matrix_vector_slice<M>::reference
      operator () (typename matrix_vector_slice<M>::size_type i)
      {
        return ublas::matrix_vector_slice<M>::data()( (this->start1()+i*this->stride1())%this->data().size1() , (this->start2()+i*this->stride2()) );
      }
#else
      BOOST_UBLAS_INLINE
      typename matrix_vector_slice<M>::reference
      operator () (typename matrix_vector_slice<M>::size_type i) const
      {
        return ublas::matrix_vector_slice<M>::data()( (this->start1()+i*this->stride1())%this->data().size1() , (this->start2()+i*this->stride2()) );
      }
#endif
  };

  template <class M>
  class matrix_vector_slice_periodic<M,direction::v>:
    public matrix_vector_slice<M>
  {
    public:
      BOOST_UBLAS_INLINE
      matrix_vector_slice_periodic ( M &data, const typename ublas::matrix_vector_slice<M>::slice_type &s1, const typename ublas::matrix_vector_slice<M>::slice_type &s2)
       : matrix_vector_slice<M>(data,s1,s2)
      {}

#ifndef BOOST_UBLAS_PROXY_CONST_MEMBER
      BOOST_UBLAS_INLINE
      typename matrix_vector_slice<M>::const_reference
      operator () (typename matrix_vector_slice<M>::size_type i) const
      {
        return ublas::matrix_vector_slice<M>::data()( (this->start1()+i*this->stride1()) , (this->start2()+i*this->stride2())%this->data().size2()  );
        //return this->operator()((i+this->data().size2())%this->data().size2());
      }
      BOOST_UBLAS_INLINE
      typename matrix_vector_slice<M>::reference
      operator () (typename matrix_vector_slice<M>::size_type i)
      {
        return ublas::matrix_vector_slice<M>::data()( (this->start1()+i*this->stride1()) , (this->start2()+i*this->stride2())%this->data().size2()  );
        //return this->operator()((i+this->data().size2())%this->data().size2());
      }
#else
      BOOST_UBLAS_INLINE
      typename matrix_vector_slice<M>::reference
      operator () (typename matrix_vector_slice<M>::size_type i) const
      {
        return ublas::matrix_vector_slice<M>::data()( (this->start1()+i*this->stride1()) , (this->start2()+i*this->stride2())%this->data().size2()  );
        //return this->operator()((i+this->data().size2())%this->data().size2());
      }
#endif
  };

}}}

#endif

