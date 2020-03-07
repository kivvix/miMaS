#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <utility>
#include <algorithm>
#include <iterator>

#include <filesystem>

struct convertor {
  std::map<std::string,std::string> map_config;

  convertor ()
  { ; }
  template < typename OPENABLE >
  convertor ( OPENABLE input )
  {
    std::ifstream ifs(input);
    std::string key,value;
    while ( ifs >> key >> value ) {
      map_config[key] = value;
    }
    ifs.close();
  }

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
converter(std::string);

#undef converter

struct config {
  std::string name;
  std::size_t Nx;
  std::size_t Nv;
  double Tc;
  double Tf;
  double ui;
  double alpha;
  double tol;
  std::filesystem::path output_dir;

  template <typename OPENABLE>
  config ( OPENABLE path_config )
    : name("")
  {
    convertor convert(path_config);

    Nx = convert("Nx",135);
    Nv = convert("Nv",128);
    Tc = convert("Tc",0.01);
    Tf = convert("Tf",10.0);
    ui = convert("ui",3.4);
    alpha = convert("alpha",0.2);
    tol = convert("tol",1e-5);
    output_dir = convert("output_dir",std::move(std::string{"."}));
  }

  bool
  create_output_directory () const;
};

std::ostream &
operator << ( std::ostream & , const config & );

namespace monitoring
{

template < typename Container , typename Writer >
struct data
{
  std::filesystem::path file;
  const Container * dat;
  Writer writer;

  data ( std::filesystem::path && _file , const Container & _dat ,  Writer _writer )
    : file(_file) , dat(&_dat) , writer(_writer)
  { ; }
};

} // namespace monitoring

template < typename Container , typename Writer >
void
operator << ( const config & c , const monitoring::data<Container,Writer> & dat )
{
  std::ofstream of( c.output_dir / dat.file );
  std::transform( std::begin(*dat.dat) , std::end(*dat.dat) , std::ostream_iterator<std::string>(of,"\n") , dat.writer );
  of.close();
}

#endif
