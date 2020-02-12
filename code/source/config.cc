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

config::config ( std::filesystem::path p_config )
{
  convertor convert;

  std::ifstream if_config(p_config);
  std::string key,value;
  while( if_config >> key >> value ) {
    convert.map_config[key] = value;
  }
  if_config.close();

  Nx = convert("Nx",135);
  Nv = convert("Nv",128);
  Tc = convert("Tc",0.01);
  Tf = convert("Tf",10.0);
  alpha = convert("alpha",0.2);
  tol = convert("tol",1e-7);
  output_dir = convert("output_dir",std::string{"."});
}

std::ostream &
operator << ( std::ostream & os , const config & c )
{
  os << "Nx " << c.Nx << "\n"
     << "Nv " << c.Nv << "\n"
     << "Tc " << c.Tc << "\n"
     << "Tf " << c.Tf << "\n"
     << "alpha " << c.alpha << "\n"
     << "tol " << c.tol << "\n"
     << "output_dir " << c.output_dir.string();
  return os;
}
