#ifndef _COMPLEX_FIELD_H_
#define _COMPLEX_FIELD_H_

#include <complex>

#include <boost/multi_array.hpp>

using namespace boost::numeric;

template < typename _T , std::size_t NumDimsV >
using complex_field = boost::multi_array< std::complex<_T> , NumDimsV+1 >;

#endif
