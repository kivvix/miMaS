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

    BOOST_UBLAS_INLINE
    typename matrix_vector_slice<M>::const_reference
    operator [] (typename matrix_vector_slice<M>::size_type i) const
    { return (*this) (i); }
    BOOST_UBLAS_INLINE
    typename matrix_vector_slice<M>::reference
    operator [] (typename matrix_vector_slice<M>::size_type i)
    { return (*this) (i); }
#else
    BOOST_UBLAS_INLINE
    typename matrix_vector_slice<M>::reference
    operator () (typename matrix_vector_slice<M>::size_type i) const
    {
      return ublas::matrix_vector_slice<M>::data()( (this->start1()+i*this->stride1())%this->data().size1() , (this->start2()+i*this->stride2()) );
    }

    BOOST_UBLAS_INLINE
    typename matrix_vector_slice<M>::reference
    operator [] (typename matrix_vector_slice<M>::size_type i) const
    { return (*this) (i); }
#endif

    class const_iterator :
      virtual public matrix_vector_slice<M>::const_iterator,
      public container_const_reference<matrix_vector_slice_periodic>
    {
      public:
        BOOST_UBLAS_INLINE
        const_iterator () {}
        BOOST_UBLAS_INLINE
        const_iterator ( const matrix_vector_slice_periodic &mvs , const typename ublas::matrix_vector_slice<M>::slice_type::const_iterator &it1 , const typename ublas::matrix_vector_slice<M>::slice_type::const_iterator &it2 )
          : container_const_reference<matrix_vector_slice_periodic>(mvs) , matrix_vector_slice<M>::const_iterator(mvs,it1,it2)
        {}
        BOOST_UBLAS_INLINE
        const_iterator ( const typename matrix_vector_slice_periodic::iterator &it )
          : container_const_reference<matrix_vector_slice_periodic>(it()) , matrix_vector_slice<M>::const_iterator(it)
        {}
        BOOST_UBLAS_INLINE
        const_iterator ( const matrix_vector_slice_periodic &mvs , const typename ublas::matrix_vector_slice<M>::const_iterator &it )
          : container_const_reference<matrix_vector_slice_periodic>(mvs) , matrix_vector_slice<M>::const_iterator(it)
        {}

        typename matrix_vector_slice<M>::const_iterator::reference
        operator * () const
        { return this->container_const_reference<matrix_vector_slice_periodic>::operator()()( this->index()%( this->container_const_reference<matrix_vector_slice_periodic>::operator()() .data().size1()) ); }
    };

    BOOST_UBLAS_INLINE
    const_iterator
    begin () const
    { return const_iterator(*this,matrix_vector_slice<M>::begin()); }
    BOOST_UBLAS_INLINE
    const_iterator
    cbegin () const
    { return begin(); }

    BOOST_UBLAS_INLINE
    const_iterator
    end () const
    { return const_iterator(*this,matrix_vector_slice<M>::end()); }
    BOOST_UBLAS_INLINE
    const_iterator
    cend () const
    { return end(); }

    class iterator :
      virtual public matrix_vector_slice<M>::iterator,
      public container_reference<matrix_vector_slice_periodic>
    {
      public:
        BOOST_UBLAS_INLINE
        iterator () {}
        BOOST_UBLAS_INLINE
        iterator ( const matrix_vector_slice_periodic &mvs , const typename ublas::matrix_vector_slice<M>::slice_type::const_iterator &it1 , const typename ublas::matrix_vector_slice<M>::slice_type::const_iterator &it2 )
          : container_reference<matrix_vector_slice_periodic>(mvs) , matrix_vector_slice<M>::iterator(mvs,it1,it2)
        {}
        BOOST_UBLAS_INLINE
        iterator ( const typename matrix_vector_slice_periodic::iterator &it )
          : container_reference<matrix_vector_slice_periodic>(it()) , matrix_vector_slice<M>::iterator(it)
        {}
        BOOST_UBLAS_INLINE
        iterator ( const matrix_vector_slice_periodic &mvs , const typename ublas::matrix_vector_slice<M>::iterator &it )
          : container_reference<matrix_vector_slice_periodic>(mvs) , matrix_vector_slice<M>::iterator(it)
        {}

        typename matrix_vector_slice<M>::iterator::reference
        operator * () const
        { return this->container_reference<matrix_vector_slice_periodic>::operator()()( this->index()%( this->container_reference<matrix_vector_slice_periodic>::operator()() .data().size1()) ); }
    };

    BOOST_UBLAS_INLINE
    iterator
    begin ()
    { return iterator(*this,matrix_vector_slice<M>::begin()); }
    
    BOOST_UBLAS_INLINE
    iterator
    end ()
    { return iterator(*this,matrix_vector_slice<M>::end()); }
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
    }
    BOOST_UBLAS_INLINE
    typename matrix_vector_slice<M>::reference
    operator () (typename matrix_vector_slice<M>::size_type i)
    {
      return ublas::matrix_vector_slice<M>::data()( (this->start1()+i*this->stride1()) , (this->start2()+i*this->stride2())%this->data().size2()  );
    }

    BOOST_UBLAS_INLINE
    typename matrix_vector_slice<M>::const_reference
    operator [] (typename matrix_vector_slice<M>::size_type i) const
    { return (*this) (i); }
    BOOST_UBLAS_INLINE
    typename matrix_vector_slice<M>::reference
    operator [] (typename matrix_vector_slice<M>::size_type i)
    { return (*this) (i); }
#else
    BOOST_UBLAS_INLINE
    typename matrix_vector_slice<M>::reference
    operator () (typename matrix_vector_slice<M>::size_type i) const
    {
      return ublas::matrix_vector_slice<M>::data()( (this->start1()+i*this->stride1()) , (this->start2()+i*this->stride2())%this->data().size2()  );
    }

    BOOST_UBLAS_INLINE
    typename matrix_vector_slice<M>::reference
    operator [] (typename matrix_vector_slice<M>::size_type i) const
    { return (*this) (i); }
#endif


    class const_iterator :
      virtual public matrix_vector_slice<M>::const_iterator,
      public container_const_reference<matrix_vector_slice_periodic>
    {
      public:
        BOOST_UBLAS_INLINE
        const_iterator () {}
        BOOST_UBLAS_INLINE
        const_iterator ( const matrix_vector_slice_periodic &mvs , const typename ublas::matrix_vector_slice<M>::slice_type::const_iterator &it1 , const typename ublas::matrix_vector_slice<M>::slice_type::const_iterator &it2 )
          : container_const_reference<matrix_vector_slice_periodic>(mvs) , matrix_vector_slice<M>::const_iterator(mvs,it1,it2)
        {}
        BOOST_UBLAS_INLINE
        const_iterator ( const typename matrix_vector_slice_periodic::iterator &it )
          : container_const_reference<matrix_vector_slice_periodic>(it()) , matrix_vector_slice<M>::const_iterator(it)
        {}
        BOOST_UBLAS_INLINE
        const_iterator ( const matrix_vector_slice_periodic &mvs , const typename ublas::matrix_vector_slice<M>::const_iterator &it )
          : container_const_reference<matrix_vector_slice_periodic>(mvs) , matrix_vector_slice<M>::const_iterator(it)
        {}

        typename matrix_vector_slice<M>::const_iterator::reference
        operator * () const
        { return this->container_const_reference<matrix_vector_slice_periodic>::operator()()( this->index()%( this->container_const_reference<matrix_vector_slice_periodic>::operator()() .data().size2()) ); }
    };

    BOOST_UBLAS_INLINE
    const_iterator
    begin () const
    { return const_iterator(*this,matrix_vector_slice<M>::begin()); }
    BOOST_UBLAS_INLINE
    const_iterator
    cbegin () const
    { return begin(); }

    BOOST_UBLAS_INLINE
    const_iterator
    end () const
    { return const_iterator(*this,matrix_vector_slice<M>::end()); }
    BOOST_UBLAS_INLINE
    const_iterator
    cend () const
    { return end(); }

    class iterator :
      virtual public matrix_vector_slice<M>::iterator,
      public container_reference<matrix_vector_slice_periodic>
    {
      public:
        BOOST_UBLAS_INLINE
        iterator () {}
        BOOST_UBLAS_INLINE
        iterator ( const matrix_vector_slice_periodic &mvs , const typename ublas::matrix_vector_slice<M>::slice_type::const_iterator &it1 , const typename ublas::matrix_vector_slice<M>::slice_type::const_iterator &it2 )
          : container_reference<matrix_vector_slice_periodic>(mvs) , matrix_vector_slice<M>::iterator(mvs,it1,it2)
        {}
        BOOST_UBLAS_INLINE
        iterator ( const typename matrix_vector_slice_periodic::iterator &it )
          : container_reference<matrix_vector_slice_periodic>(it()) , matrix_vector_slice<M>::iterator(it)
        {}
        BOOST_UBLAS_INLINE
        iterator ( const matrix_vector_slice_periodic &mvs , const typename ublas::matrix_vector_slice<M>::iterator &it )
          : container_reference<matrix_vector_slice_periodic>(mvs) , matrix_vector_slice<M>::iterator(it)
        {}

        typename matrix_vector_slice<M>::iterator::reference
        operator * () const
        { return this->container_reference<matrix_vector_slice_periodic>::operator()()( this->index()%( this->container_reference<matrix_vector_slice_periodic>::operator()() .data().size2()) ); }
    };

    BOOST_UBLAS_INLINE
    iterator
    begin ()
    { return iterator(*this,matrix_vector_slice<M>::begin()); }
    
    BOOST_UBLAS_INLINE
    iterator
    end ()
    { return iterator(*this,matrix_vector_slice<M>::end()); }

};

}}}

#endif

