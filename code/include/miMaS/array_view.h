#ifndef cqder_ARRAY_VIEW_H
#define cqder_ARRAY_VIEW_H

#include <utility>
#include <iterator>
#include <stdexcept>
#include <sstream>


/**
struct:
  name: array_view<T>
  brief: this class present a pointer like an array with facilities to iterrate over it and safety access to memory
  details: >
    This class doesn't allocate or free memory, this is just a wrapper to use a
    pointer like a C++ container with iterators and safety access to memory, it
    could be usefull to examge data with a C or Fortran library.

**/

/**
~~~[TEST]
  TEST(constructors)
  {
    std::size_t n = 42;
    double * x = new double[n];
    tools::array_view<double> vx {x,n};

    expect_eq( n , vx.size() );

    std::size_t m=0;
    for ( const auto & elm : vx )
    { ++m; }
    
    expect_eq( n , m );
  }
~~~
**/
namespace tools
{

template < typename T >
struct array_view
{
  // this class doesn't allocate or free memory, this is just a wrapper to use a
  // pointer like a C++ container with iterators and [] operator

  // defined type like other containers
  typedef T                                     value_type;
  typedef std::size_t                           size_type;
  typedef std::ptrdiff_t                        difference_type;
  typedef value_type &                          reference;
  typedef value_type const&                     const_reference;
  typedef value_type *                          pointer;
  typedef value_type const*                     const_pointer;
  typedef pointer                               iterator;
  typedef const_pointer                         const_iterator;
  typedef std::reverse_iterator<iterator>       reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

  private:
  // attributs -----------------------------------------------------------------
  // just a pointer on data
  // size of array (can't be determine only with a pointer)
    pointer   m_data;
    size_type m_size;

  protected:
    // safety check used only from at()
    void
    range_check ( size_t pos ) const
    {
      if ( pos >= m_size )
      {
        std::stringstream ss;
        ss << "array_view::range_check: pos (which is "<< pos << ") >= this->size() (which is " << m_size << ")";
        throw std::out_of_range(ss.str());
      }
    }

  public:
  // constructors --------------------------------------------------------------
    // don't allocate memory here! this class is just a wrapper not a memory manager
    // default constructor (useless in main cases)
    array_view ()
      : m_data(nullptr) , m_size(0)
    {}

    // copy constructor (absurd but why not)
    array_view ( array_view const& rhs )
      : m_data(rhs.m_data) , m_size(rhs.m_size)
    {}

    // constructor from a pair of attributs types
    array_view ( std::pair<pointer,size_type> t_pair )
      : m_data(t_pair.first) , m_size(t_pair.second)
    {}

    // constructor with two parameter, one for each attribut
    array_view ( pointer t_data , size_type t_size )
      : m_data(t_data) , m_size(t_size)
    {}

    array_view ( std::initializer_list<T> il )
      : m_data(il.begin()) , m_size(il.size())
    {}

  // destructor
    // don't free memory here! this class is just a wrapper not a memory manager
    virtual
    ~array_view ()
    {}

  // standard operator
    inline array_view &
    operator = ( array_view const& rhs )
    {
      m_data = rhs.m_data;
      m_size = rhs.m_size;
      return *this;
    }

  // iterator api
    // begin and end member functions for iterator and const_iterator
    inline const_iterator
    begin () const
    { return const_iterator(m_data); }
    inline const_iterator
    cbegin () const
    { return const_iterator(m_data); }

    inline iterator
    begin ()
    { return iterator(m_data); }

    inline const_iterator
    end () const
    { return const_iterator(m_data + m_size); }
    inline const_iterator
    cend () const
    { return const_iterator(m_data + m_size); }

    inline iterator
    end ()
    { return iterator(m_data + m_size); }

    // rbegin and rend member functions for reverse_iterator and const_reverse_iterator
    inline const_reverse_iterator
    rbegin () const
    { return std::reverse_iterator<const_iterator>(end()); }
    inline const_reverse_iterator
    crbegin () const
    { return std::reverse_iterator<const_iterator>(end()); }
    inline reverse_iterator
    rbegin ()
    { return std::reverse_iterator<iterator>(end()); }

    inline const_reverse_iterator
    rend () const
    { return std::reverse_iterator<const_iterator>(begin()); }
    inline const_reverse_iterator
    crend () const
    { return std::reverse_iterator<const_iterator>(begin()); }
    inline reverse_iterator
    rend ()
    { return std::reverse_iterator<iterator>(begin()); }

  // array api
    // operator [] to use it like an array (doesn't check limits)
    inline const_reference
    operator [] ( size_type n ) const
    { return *(m_data + n); }
    inline reference
    operator [] ( size_type n )
    { return *(m_data + n); }

    // safety aceess to memory
    inline const_reference
    at ( size_type pos ) const
    {
      range_check(pos);
      return (*this)[pos];
    }

    inline reference
    at ( size_type pos )
    {
      range_check(pos);
      return (*this)[pos];
    }

    // return the size of array
    inline size_type
    size () const
    { return m_size; }

    inline size_type
    max_size () const
    { return m_size; }

    inline pointer
    data ()
    { return m_data; }
    inline const_pointer
    data () const
    { return m_data; }

    inline reference
    front ()
    { return *m_data; }
    inline const_reference
    front () const
    { return *m_data; }

    inline reference
    back ()
    { return *(m_data+m_size-1); }
    inline const_reference
    back () const
    { return *(m_data+m_size-1); }

    inline array_view
    queue () const
    { return array_view( (m_size>1)?m_data+1:nullptr , m_size-1 ); }

    inline bool
    empty () const
    { return m_size == 0; }

    // becareful of double free
    inline void
    clear ()
    { delete[] m_data; m_data = nullptr; m_size = 0; }

    inline void
    swap ( array_view & rhs )
    {
      std::swap(m_size,rhs.m_size);
      std::swap(m_data,rhs.m_data);
    }

  // array_view api

    // return a subarray of `array_view`
    inline array_view
    subarray ( size_type pos ) const
    {
      range_check(pos);
      return array_view( m_data+pos , m_size-pos );
    }

    // return a subarray of lenght `count` of `array_view`
    inline array_view
    subarray ( size_type pos , std::size_t count ) const
    {
      range_check(pos);
      return array_view( m_data+pos , std::min(count,m_size-pos) );
    }

};

} // namespace tools

#endif
