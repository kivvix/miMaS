#include "miMaS/config.h"

#define converter( type , stotype ) template <> \
type \
convertor::operator () ( std::string && key , type && default_value ) \
{ \
  auto it = map_config.find(key); \
  if ( it != map_config.end() ) { return stotype(it->second); }\
  return std::move(default_value); \
}

  converter(int,std::stoi)
  converter(long,std::stol)
  converter(unsigned long,std::stoul)
  converter(unsigned long long,std::stoull)
  converter(float,std::stof)
  converter(double,std::stod)
  converter(long double, std::stold)
#undef converter

template <>
std::string
convertor::operator () ( std::string && key , std::string && default_value )
{
  auto it = map_config.find(key);
  if ( it != map_config.end() ) { return it->second; }
  return std::move(default_value);
}

bool
config::create_output_directory () const
{ return std::filesystem::create_directories(output_dir); }

std::ostream &
operator << ( std::ostream & os , const config & c )
{
  os << "Nx " << c.Nx << "\n"
     << "Nv " << c.Nv << "\n"
     << "Tc " << c.Tc << "\n"
     << "Tf " << c.Tf << "\n"
     << "ui " << c.ui << "\n" 
     << "alpha " << c.alpha << "\n"
     << "tol " << c.tol << "\n"
     << "output_dir " << c.output_dir.string();
  return os;
}
