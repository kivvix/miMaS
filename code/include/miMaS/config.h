#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <utility>

#include <filesystem>

struct convertor {
  std::map<std::string,std::string> map_config;

  convertor ()
  { ; }

  template < typename _T >
  _T operator () ( std::string && key , _T && default_value )
  { return std::move(default_value); }

};

#define converter( type ) template <> \
type \
convertor::operator () ( std::string && , type && )

converter(int);
converter(long);
converter(unsigned long);
converter(unsigned long long);
converter(float);
converter(double);
converter(long double);

#undef converter

struct config {
  std::size_t Nx;
  std::size_t Nv;
  double Tc;
  double Tf;
  double alpha;
  double tol;
  std::filesystem::path output_dir;

  config ( std::filesystem::path );
};

std::ostream &
operator << ( std::ostream & , const config & );

#endif
