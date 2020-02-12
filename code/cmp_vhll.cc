#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <fstream>
#include <valarray>
#include <sstream>
#include <filesystem>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include "miMaS/field.h"
#include "miMaS/complex_field.h"
#include "miMaS/weno.h"
#include "miMaS/fft.h"
#include "miMaS/array_view.h"
#include "miMaS/poisson.h"
#include "miMaS/rk.h"
#include "miMaS/config.h"

namespace o2 {
  template < typename _T , std::size_t NumDimsV >
  auto
  trp_v ( field<_T,NumDimsV> const & u , ublas::vector<_T> const& E )
  {
    field<_T,NumDimsV> trp(tools::array_view<const std::size_t>(u.shape(),NumDimsV+1));

    { auto k=0, km1=trp.size(0)-1;
      for ( auto i=0 ; i<trp.size(1) ; ++i ) {
        trp[k][i] = ( E(i)*(u[k+1][i]-u[km1][i])/(2.*u.step.dv) );
      }
    }
    for ( auto k=1 ; k<trp.size(0)-1 ; ++k ) {
      for ( auto i=0 ; i<trp.size(1) ; ++i ) {
        trp[k][i] = ( E(i)*(u[k+1][i]-u[k-1][i])/(2.*u.step.dv) );
      }
    }
    { auto k=trp.size(0)-1, kp1=0;
      for ( auto i=0 ; i<trp.size(1) ; ++i ) {
        trp[k][i] = ( E(i)*(u[kp1][i]-u[k-1][i])/(2.*u.step.dv) );
      }
    }

    return trp;
  }
}

struct iter_s {
  int iter;
  double dt;
  double current_time;
  double L;
};


namespace math = boost::math::constants;
const std::complex<double> & I = std::complex<double>(0.,1.);

#define SQ(X) ((X)*(X))
#define Xi(i) (i*fh.step.dx+fh.range.x_min)
#define Vk(k) (k*fh.step.dv+fh.range.v_min)


auto
maxwellian ( double rho , double u , double T ) {
  return [=](double x,double v){ return rho/(std::sqrt(2.*math::pi<double>()*T))*std::exp( -0.5*SQ(v-u)/T ); };
}

int
main(int argc, char const *argv[])
{
  std::filesystem::path p("config.init");
  if ( argc > 1 ) {
    p = argv[1];
  }
  auto c = config(p);

  std::filesystem::create_directories(c.output_dir);
  std::ofstream oconfig( c.output_dir / "config.init" );
  oconfig << c << std::endl;
  oconfig.close();

/* --------------------------------------------------------------- */
  field<double,1> fh(boost::extents[c.Nv][c.Nx]);
  complex_field<double,1> hfh(boost::extents[c.Nv][c.Nx]);

  fh.range.v_min = -8.; fh.range.v_max = 8.;
  fh.step.dv = (fh.range.v_max-fh.range.v_min)/c.Nv;

  double Kx = 0.5;
  fh.range.x_min = 0.; fh.range.x_max = 2./Kx*math::pi<double>();
  fh.step.dx = (fh.range.x_max-fh.range.x_min)/c.Nx;

  ublas::vector<double> v (c.Nv,0.);
  for ( std::size_t k=0 ; k<c.Nv ; ++k ) { v[k] = Vk(k); }

  ublas::vector<double> kx(c.Nx);
  {
    double l = fh.range.len_x();
    for ( auto i=0 ; i<c.Nx/2 ; ++i ) { kx[i]      = 2.*math::pi<double>()*i/l; }
    for ( int i=-c.Nx/2 ; i<0 ; ++i ) { kx[c.Nx+i] = 2.*math::pi<double>()*i/l; }
  }

  double alpha = c.alpha , ui = 2.;
  auto tb_M1    = maxwellian( 0.5*alpha ,  ui , 1.   ) , tb_M2  = maxwellian( 0.5*alpha , -ui , 1.  );

  for (field<double,2>::size_type k=0 ; k<fh.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<fh.size(1) ; ++i ) {
      // tb
      fh[k][i] = ( tb_M1(Xi(i),Vk(k)) + tb_M2(Xi(i),Vk(k)) )*(1. + 0.01*std::cos(Kx*Xi(i)));
      //fh[k][i] = ( tb_M1(Xi(i),Vk(k)) + tb_M2(Xi(i),Vk(k)) )*(1.);
    }
    fft::fft(&(fh[k][0]),&(fh[k][c.Nx-1])+1,&(hfh[k][0]));
  }
  fh.write( c.output_dir / "init_vhll.dat" );

  unsigned int i_t = 0;
  double current_time = 0.;
  double dt = 1.*fh.step.dv;


  ublas::vector<double> uc(c.Nx,0.);
  ublas::vector<double> E (c.Nx,0.);

  std::vector<double> ee;   ee.reserve(int(std::ceil(c.Tf/dt))+1);
  std::vector<double> Emax; Emax.reserve(int(std::ceil(c.Tf/dt))+1);
  std::vector<double> H;    H.reserve(int(std::ceil(c.Tf/dt))+1);

  std::vector<double> times; times.reserve(int(std::ceil(c.Tf/dt))+1);
  std::vector<double> dts;   dts.reserve(int(std::ceil(c.Tf/dt))+1);
  std::vector<iter_s> iterations; iterations.reserve(int(std::ceil(c.Tf/dt))+1);


  const double rho_c = 1.-alpha;
  const double sqrt_rho_c = std::sqrt(rho_c);
  {
    poisson<double> poisson_solver(c.Nx,fh.range.len_x());
    ublas::vector<double> rho(c.Nx,0.);
    rho = fh.density(); // compute density from init data
    for ( auto i=0 ; i<c.Nx ; ++i ) { rho[i] += (1.-alpha); } // add (1-alpha) for cold particules
    E = poisson_solver(rho);

    Emax.push_back( std::abs(*std::max_element( E.begin() , E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} )) );
    double electric_energy = 0.;
    for ( const auto & ei : E ) { electric_energy += ei*ei*fh.step.dx; }
    ee.push_back( std::sqrt(electric_energy) );

    double total_energy = energy(fh,E);
    total_energy += 0.; // sum(rho_c*u_c*u_c) = 0 because u_c = 0 at time 0
    H.push_back( total_energy );
    times.push_back( 0. );
  }


  // initialize memory for all temporary variables
  ublas::vector<double> J(c.Nx);
  fft::spectrum_ d(c.Nx);
  field<double,1> Edvf(tools::array_view<const std::size_t>(fh.shape(),2));
  ublas::vector<double> uc1(c.Nx) , uc2(c.Nx) , uc3(c.Nx) , uc4(c.Nx) , uc5(c.Nx) , uc6(c.Nx) , uc7(c.Nx),
                        E1(c.Nx)  , E2(c.Nx)  , E3(c.Nx)  , E4(c.Nx)  , E5(c.Nx)  , E6(c.Nx)  , E7(c.Nx);
  complex_field<double,1> hfh1(boost::extents[c.Nv][c.Nx]),hfh2(boost::extents[c.Nv][c.Nx]),hfh3(boost::extents[c.Nv][c.Nx]),hfh4(boost::extents[c.Nv][c.Nx]),hfh5(boost::extents[c.Nv][c.Nx]),hfh6(boost::extents[c.Nv][c.Nx]),hfh7(boost::extents[c.Nv][c.Nx]);

  while (  current_time < c.Tf ) {
    std::cout<<"\r ["<<std::setw(6)<<i_t<<"] "<< std::setw(8) << current_time << " (" << std::setw(9) << dt << ")"<<std::flush;
    

/////////////////////////////////////////////////////////////////////
/**
    // RK33

    // STAGE 1
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh[k][0]),&(hfh[k][0])+c.Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E);

      double c1 = std::cos(dt*sqrt_rho_c), s1 = std::sin(dt*sqrt_rho_c);
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        uc1[i] =  uc[i]*c1 + E[i]*s1/sqrt_rho_c - dt*J[i]*s1/sqrt_rho_c;
        E1[i]  = -uc[i]*s1*sqrt_rho_c + E[i]*c1 - dt*J[i]*c1;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh1[k][i] = hfh[k][i]*std::exp(-I*kx[i]*Vk(k)*dt) - dt*d[i]*std::exp(-I*kx[i]*Vk(k)*dt);
        }
      }
    } // end stage 1

    // STAGE 2
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh1[k][0]),&(hfh1[k][0])+c.Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E1);

      double c2 = std::cos(0.5*dt*sqrt_rho_c), s2 = std::sin(0.5*dt*sqrt_rho_c);
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        uc2[i] =  0.75*uc[i]*c2 + 0.75*E[i]*s2/sqrt_rho_c + 0.25*uc1[i]*c2 - 0.25*E1[i]*s2/sqrt_rho_c + 0.25*dt*J[i]*s2/sqrt_rho_c;
        E2[i]  = -0.75*uc[i]*s2*sqrt_rho_c + 0.75*E[i]*c2 + 0.25*uc1[i]*s2*sqrt_rho_c + 0.25*E1[i]*c2 - 0.25*dt*J[i]*c2;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh2[k][i] = 0.75*hfh[k][i]*std::exp(-I*kx[i]*Vk(k)*0.5*dt) + 0.25*hfh1[k][i]*std::exp(I*kx[i]*Vk(k)*0.5*dt) - 0.25*dt*d[i]*std::exp(I*kx[i]*Vk(k)*0.5*dt);
        }
      }
    } // end stage 2

    // STAGE 3
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh2[k][0]),&(hfh2[k][0])+c.Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E2);

      double c1 = std::cos(dt*sqrt_rho_c), s1 = std::sin(dt*sqrt_rho_c);
      double c2 = std::cos(0.5*dt*sqrt_rho_c), s2 = std::sin(0.5*dt*sqrt_rho_c);
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        double tmp_uc =  (1./3.)*uc[i]*c1 + (1./3.)*E[i]*s1/sqrt_rho_c + (2./3.)*uc2[i]*c2 + (2./3.)*E2[i]*s2/sqrt_rho_c - (2./3.)*dt*J[i]*s2/sqrt_rho_c;
        E[i]          = -(1./3.)*uc[i]*s1*sqrt_rho_c + (1./3.)*E[i]*c1 - (2./3.)*uc2[i]*s2*sqrt_rho_c + (2./3.)*E2[i]*c2 - (2./3.)*dt*J[i]*c2;
        uc[i] = tmp_uc;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh[k][i] = (1./3.)*hfh[k][i]*std::exp(-I*kx[i]*Vk(k)*dt) + (2./3.)*hfh2[k][i]*std::exp(-I*kx[i]*Vk(k)*0.5*dt) - (2./3.)*dt*d[i]*std::exp(-I*kx[i]*Vk(k)*0.5*dt);
        }
      }
    } // end stage 3

/////////////////////////////////////////////////////////////////////
/**/
/////////////////////////////////////////////////////////////////////

    //DP5
    // STAGE 1
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh[k][0]),&(hfh[k][0])+c.Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E);

      double c02 = std::cos(0.2*dt*sqrt_rho_c) , s02 = std::sin(0.2*dt*sqrt_rho_c);
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        uc1[i] =  uc[i]*c02 + E[i]*s02/sqrt_rho_c - 0.2*dt*J[i]*s02/sqrt_rho_c;
        E1[i]  = -uc[i]*s02*sqrt_rho_c + E[i]*c02 - 0.2*dt*J[i]*c02;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh1[k][i] = hfh[k][i]*std::exp(-0.2*I*kx[i]*Vk(k)*dt) - 0.2*dt*d[i]*std::exp(-0.2*I*kx[i]*Vk(k)*dt);
        }
      }
    } // end stage 1

    // STAGE 2
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh1[k][0]),&(hfh1[k][0])+c.Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E1);

      double c01 = std::cos(0.1*sqrt_rho_c*dt) , s01 = std::sin(0.1*sqrt_rho_c*dt) , c03 = std::cos(0.3*sqrt_rho_c*dt) , s03 = std::sin(0.3*sqrt_rho_c*dt);
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        uc2[i] = 0.625*(  uc[i]*c03 + E[i]*s03/sqrt_rho_c ) + 0.375*( uc1[i]*c01 + E1[i]*s01/sqrt_rho_c ) - 0.225*dt*J[i]*s01/sqrt_rho_c;
        E2[i]  = 0.625*( -uc[i]*s03*sqrt_rho_c + E[i]*c03 ) + 0.375*( -uc1[i]*sqrt_rho_c*s01 + E1[i]*c01 ) - 0.225*dt*J[i]*c01;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh2[k][i] = 0.625*hfh[k][i]*std::exp(-0.3*I*kx[i]*Vk(k)*dt) + 0.375*hfh1[k][i]*std::exp(-0.1*I*kx[i]*Vk(k)*dt) - 0.225*dt*d[i]*std::exp(-0.1*I*kx[i]*Vk(k));
        }
      }
    } // end stage 2

    // STAGE 3
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh2[k][0]),&(hfh2[k][0])+c.Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E2);

      double c08 = std::cos(0.8*dt*sqrt_rho_c) , s08 = std::sin(0.8*dt*sqrt_rho_c) ,
             c06 = std::cos(0.6*dt*sqrt_rho_c) , s06 = std::sin(0.6*dt*sqrt_rho_c) ,
             c05 = std::cos(0.5*dt*sqrt_rho_c) , s05 = std::sin(0.5*dt*sqrt_rho_c) ;
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        uc3[i] = (175./27.)*(  uc[i]*c08 + E[i]*s08/sqrt_rho_c ) + (100./9.)*(  uc1[i]*c06 + E1[i]*s06/sqrt_rho_c ) - (448./27.)*(  uc2[i]*c05 + E2[i]*s05/sqrt_rho_c ) - (32./9.)*dt*J[i]*s05/sqrt_rho_c;
        E3[i]  = (175./27.)*( -uc[i]*s08*sqrt_rho_c + E[i]*c08 ) + (100./9.)*( -uc1[i]*s06*sqrt_rho_c + E1[i]*c06 ) - (448./27.)*( -uc2[i]*s05*sqrt_rho_c + E2[i]*c05 ) - (32./9.)*dt*J[i]*c05;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh3[k][i] = (175./27.)*hfh[k][i]*std::exp(-0.8*I*kx[i]*Vk(k)*dt) + (100./9.)*hfh1[k][i]*std::exp(-0.6*I*kx[i]*Vk(k)*dt) - (448./27.)*hfh2[k][i]*std::exp(-0.5*I*kx[i]*Vk(k)*dt) - (32./9.)*dt*d[i]*std::exp(-0.5*I*kx[i]*Vk(k));
        }
      }
    } // end stage 3

    // STAGE 4
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh3[k][0]),&(hfh3[k][0])+c.Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E3);

      double c89    = std::cos((8./9.)*sqrt_rho_c*dt)    , s89    = std::sin((8./9.)*sqrt_rho_c*dt)    ,
             c3145  = std::cos((31./45.)*sqrt_rho_c*dt)  , s3145  = std::sin((31./45.)*sqrt_rho_c*dt)  ,
             c5390  = std::cos((53./90.)*sqrt_rho_c*dt)  , s5390  = std::sin((53./90.)*sqrt_rho_c*dt)  ,
             c44827 = std::cos((448./27.)*sqrt_rho_c*dt) , s44827 = std::sin((448./27.)*sqrt_rho_c*dt) ,
             c445   = std::cos((4./45.)*sqrt_rho_c*dt)   , s445   = std::sin((4./45.)*sqrt_rho_c*dt)   ;
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        uc4[i] = (3551./6561.)*(  uc[i]*c89 + E[i]*s89/sqrt_rho_c ) + (7420./2187.)*(  uc1[i]*c3145 + E1[i]*s3145/sqrt_rho_c ) - (37376./6561.)*(  uc2[i]*c5390 + E2[i]*s5390/sqrt_rho_c ) + (2014./729.)*(  uc3[i]*c445 + E3[i]*s445/sqrt_rho_c ) + (212./729.)*dt*J[i]*s445/sqrt_rho_c;
        E4[i]  = (3551./6561.)*( -uc[i]*s89*sqrt_rho_c + E[i]*c89 ) + (7420./2187.)*( -uc1[i]*s3145*sqrt_rho_c + E1[i]*c3145 ) - (37376./6561.)*( -uc2[i]*s5390*sqrt_rho_c + E2[i]*c5390 ) + (2014./729.)*( -uc3[i]*s445*sqrt_rho_c + E3[i]*c445 ) + (212./729.)*dt*J[i]*c445;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh4[k][i] = (3551./6561.)*hfh[k][i]*std::exp(-(8./9.)*I*kx[i]*Vk(k)*dt) + (7420./2187.)*hfh1[k][i]*std::exp(-(31./45.)*I*kx[i]*Vk(k)*dt) - (37376./6561.)*hfh2[k][i]*std::exp(-(53./90.)*I*kx[i]*Vk(k)*dt) + (2014./729.)*hfh3[k][i]*std::exp(-(4./45.)*I*kx[i]*Vk(k)*dt) + (212./729.)*dt*d[i]*std::exp(-(4./45.)*kx[i]*Vk(k)*dt);
        }
      }
    } // end stage 4

    // STAGE 5
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh4[k][0]),&(hfh4[k][0])+c.Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E4);

      double c1  = std::cos(sqrt_rho_c*dt)         , s1  = std::sin(sqrt_rho_c*dt)         ,
             c08 = std::cos(0.8*sqrt_rho_c*dt)     , s08 = std::sin(0.8*sqrt_rho_c*dt)     ,
             c07 = std::cos(0.7*sqrt_rho_c*dt)     , s07 = std::sin(0.7*sqrt_rho_c*dt)     ,
             c02 = std::cos(0.2*sqrt_rho_c*dt)     , s02 = std::sin(0.2*sqrt_rho_c*dt)     ,
             c19 = std::cos((1./9.)*sqrt_rho_c*dt) , s19 = std::sin((1./9.)*sqrt_rho_c*dt) ;
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        uc5[i] = (313397./335808.)*(  uc[i]*c1 + E[i]*s1/sqrt_rho_c ) + (424025./55968.)*(  uc1[i]*c08 + E1[i]*s08/sqrt_rho_c ) - (61400./5247.)*(  uc2[i]*c07 + E2[i]*s07/sqrt_rho_c ) + (96075./18656.)*(  uc3[i]*c02 + E3[i]*s02/sqrt_rho_c ) - (35721./37312.)*(  uc4[i]*c19 + E4[i]*s19/sqrt_rho_c ) + (5103./18656.)*dt*J[i]*s19/sqrt_rho_c;
        E5[i]  = (313397./335808.)*( -uc[i]*s1*sqrt_rho_c + E[i]*c1 ) + (424025./55968.)*( -uc1[i]*s08*sqrt_rho_c + E1[i]*c08 ) - (61400./5247.)*( -uc2[i]*s07*sqrt_rho_c + E2[i]*c07 ) + (96075./18656.)*( -uc3[i]*s02*sqrt_rho_c + E3[i]*c02 ) - (35721./37312.)*( -uc4[i]*s19*sqrt_rho_c + E4[i]*c19 ) + (5103./18656.)*dt*J[i]*c19;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh5[k][i] = (313397./335808.)*hfh[k][i]*std::exp(-I*kx[i]*Vk(k)*dt) + (424025./55968.)*hfh1[k][i]*std::exp(-0.8*I*kx[i]*Vk(k)*dt) - (61400./5247.)*hfh2[k][i]*std::exp(-0.7*I*kx[i]*Vk(k)*dt) + (96075./18656.)*hfh3[k][i]*std::exp(-0.2*I*kx[i]*Vk(k)*dt) - (35721./37312.)*hfh4[k][i]*std::exp(-(1./9.)*I*kx[i]*Vk(k)*dt) + (5103./18656.)*dt*d[i]*std::exp(-(1./9.)*I*kx[i]*Vk(k)*dt);
        }
      }
    } // end stage 5

    // STAGE 6
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh5[k][0]),&(hfh5[k][0])+c.Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E5);

      double c1  = std::cos(sqrt_rho_c*dt)         , s1  = std::sin(sqrt_rho_c*dt)         ,
             c08 = std::cos(0.8*sqrt_rho_c*dt)     , s08 = std::sin(0.8*sqrt_rho_c*dt)     ,
             c07 = std::cos(0.7*sqrt_rho_c*dt)     , s07 = std::sin(0.7*sqrt_rho_c*dt)     ,
             c02 = std::cos(0.2*sqrt_rho_c*dt)     , s02 = std::sin(0.2*sqrt_rho_c*dt)     ,
             c19 = std::cos((1./9.)*sqrt_rho_c*dt) , s19 = std::sin((1./9.)*sqrt_rho_c*dt) ;
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        uc6[i] = -(563./3456.)*(  uc[i]*c1 + E[i]*s1/sqrt_rho_c ) - (575./252.)*(  uc1[i]*c08 + E1[i]*s08/sqrt_rho_c ) + (31400./10017.)*(  uc2[i]*c07 + E2[i]*s07/sqrt_rho_c ) + (325./1344.)*(  uc3[i]*c02 + E3[i]*s02/sqrt_rho_c ) - (7533./6784.)*(  uc4[i]*c19 + E4[i]*s19/sqrt_rho_c ) + (33./28.)*uc5[i];
        E6[i]  = -(563./3456.)*( -uc[i]*s1*sqrt_rho_c + E[i]*c1 ) - (575./252.)*( -uc1[i]*s08*sqrt_rho_c + E1[i]*c08 ) + (31400./10017.)*( -uc2[i]*s07*sqrt_rho_c + E2[i]*c07 ) + (325./1344.)*( -uc3[i]*s02*sqrt_rho_c + E3[i]*c02 ) - (7533./6784.)*( -uc4[i]*s19*sqrt_rho_c + E4[i]*c19 ) + (33./28.)*E5[i] - (11./84.)*dt*J[i];
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh6[k][i] = -(563./3456.)*hfh[k][i]*std::exp(-I*kx[i]*Vk(k)*dt) - (575./252.)*hfh1[k][i]*std::exp(-0.8*I*kx[i]*Vk(k)*dt) + (31400./10017.)*hfh2[k][i]*std::exp(-0.7*I*kx[i]*Vk(k)*dt) + (325./1344.)*hfh3[k][i]*std::exp(-0.2*I*kx[i]*Vk(k)*dt) - (7533./6784.)*hfh4[k][i]*std::exp(-(1./9.)*I*kx[i]*Vk(k)*dt) + (33./28.)*hfh5[k][i] - (11./84.)*dt*d[i];
        }
      }
    } // end stage 6

    // STAGE 7
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh6[k][0]),&(hfh6[k][0])+c.Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E6);

      double c1  = std::cos(sqrt_rho_c*dt)         , s1  = std::sin(sqrt_rho_c*dt)         ,
             c08 = std::cos(0.8*sqrt_rho_c*dt)     , s08 = std::sin(0.8*sqrt_rho_c*dt)     ,
             c07 = std::cos(0.7*sqrt_rho_c*dt)     , s07 = std::sin(0.7*sqrt_rho_c*dt)     ,
             c02 = std::cos(0.2*sqrt_rho_c*dt)     , s02 = std::sin(0.2*sqrt_rho_c*dt)     ,
             c19 = std::cos((1./9.)*sqrt_rho_c*dt) , s19 = std::sin((1./9.)*sqrt_rho_c*dt) ;
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        uc7[i] = (8813./172800.)*(  uc[i]*c1 + E[i]*s1/sqrt_rho_c ) + (41./180.)*(  uc1[i]*c08 + E1[i]*s08/sqrt_rho_c ) - (4294./10017.)*(  uc2[i]*c07 + E2[i]*s07/sqrt_rho_c ) + (263./384.)*(  uc3[i]*c02 + E3[i]*s02/sqrt_rho_c ) - (137781./339200.)*(  uc4[i]*c19 + E4[i]*s19/sqrt_rho_c ) + (803./4200)*uc5[i] + (17./25.)*uc6[i];
        E7[i]  = (8813./172800.)*( -uc[i]*s1*sqrt_rho_c + E[i]*c1 ) + (41./180.)*( -uc1[i]*s08*sqrt_rho_c + E1[i]*c08 ) - (4294./10017.)*( -uc2[i]*s07*sqrt_rho_c + E2[i]*c07 ) + (263./384.)*( -uc3[i]*s02*sqrt_rho_c + E3[i]*c02 ) - (137781./339200.)*( -uc4[i]*s19*sqrt_rho_c + E4[i]*c19 ) + (803./4200)*E5[i] + (17./25.)*E6[i] - (1./40.)*dt*J[i];
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh7[k][i] = (8813./172800.)*hfh[k][i]*std::exp(-I*kx[i]*Vk(k)*dt) + (41./180.)*hfh1[k][i]*std::exp(-0.8*I*kx[i]*Vk(k)*dt) - (4294./10017.)*hfh2[k][i]*std::exp(-0.7*I*kx[i]*Vk(k)*dt) + (263./384.)*hfh3[k][i]*std::exp(-0.2*I*kx[i]*Vk(k)*dt) - (137781./339200.)*hfh4[k][i]*std::exp(-(1./9.)*I*kx[i]*Vk(k)*dt) + (803./4200)*hfh5[k][i] + (17./25.)*hfh6[k][i] - (1./40.)*dt*d[i];
        }
      }
    } // end stage 7
/**/
    // end of time loop

    // CHECK local error for compute next time step
    //double L_E = std::sqrt( std::inner_product( E7.begin() , E7.end() , E6.begin() , 0. ,
    //  [] ( const auto & a  , const auto & b   ) { return a+b; } ,
    //  [&]( const auto & en , const auto & en1 ) { return SQ( en1 - en )*fh.step.dx; } ) );
    double L_E = 0.;
    for ( auto i=0 ; i<c.Nx ; ++i )
      { L_E += SQ( E7[i] - E6[i] )*fh.step.dx; }
    L_E = std::sqrt(L_E);

    double L_hfh = 0.;
    for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
      for ( auto i=0 ; i<c.Nx ; ++ i ) {
        L_hfh += SQ( std::abs(hfh6[k][i] - hfh7[k][i]) )*fh.step.dv*fh.step.dx; // ??? dx n'a pas de sens ici...
      }
    }
    L_hfh = std::sqrt(L_hfh);

    if ( std::isnan(L_E) ) { L_E = 100.*c.tol; }

    std::cout << " -- " << std::setw(10) << L_hfh << "        " << std::flush;

    if ( std::abs(L_hfh - c.tol) < c.tol )
    {
      // SAVE TIME STEP
      std::copy( uc7.begin()  , uc7.end()  , uc.begin()  );
      std::copy( E7.begin()   , E7.end()   , E.begin()   );
      std::copy( hfh7.begin() , hfh7.end() , hfh.begin() );

      // MONITORING

      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh[k][0]),&(hfh[k][0])+c.Nx,&(fh[k][0])); }
      Emax.push_back( std::abs(*std::max_element( E.begin() , E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} )) );
      double electric_energy = 0.;
      for ( const auto & ei : E ) { electric_energy += ei*ei*fh.step.dx; }
      ee.push_back( std::sqrt(electric_energy) );

      double total_energy = energy(fh,E);
      {
        auto rhoh = fh.density();
        fft::spectrum_ hrhoh(c.Nx); hrhoh.fft(&rhoh[0]);
        fft::spectrum_ hE(c.Nx); hE.fft(&E[0]);
        fft::spectrum_ hrhoc(c.Nx);
        hrhoc[0] = I*kx[0]*hE[0] - hrhoh[0] + 1.;
        for ( auto i=1 ; i<c.Nx ; ++i ) {
          hrhoc[i] = I*kx[i]*hE[i] - hrhoh[i];
        }
        ublas::vector<double> rhoc (c.Nx,0.); hrhoc.ifft(&rhoc[0]);

        for ( auto i=0 ; i<c.Nx ; ++i ) {
          total_energy += rhoc[i]*uc[i]*uc[i];
        }
      }
      H.push_back( total_energy );

      // increment time
      current_time += dt;
      times.push_back( current_time );
    }
    iterations.push_back( { i_t , dt , current_time , L_hfh } );

    ++i_t;
    //std::cout << i_t << " " << current_time << " (" << dt << ") " << c.tol/L_E << "\n";
    dt = std::pow( c.tol/L_hfh , 0.2 )*dt; // compute optimum next time step

  } // while (  i_t*dt < Tf )
  std::cout<<" ["<<std::setw(5)<<i_t<<"] "<<i_t*dt <<std::endl;

  for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh[k][0]),&(hfh[k][0])+c.Nx,&(fh[k][0])); }
  fh.write( c.output_dir / "vp_vhll.dat" );


  std::ofstream of;

  of.open( c.output_dir / "iterations.dat" );
  std::transform( iterations.begin() , iterations.end() , std::ostream_iterator<std::string>(of,"\n") ,
    [](auto const& it) mutable { std::stringstream ss; ss << it.iter << " " << it.dt << " " << it.current_time << " " << it.L; return ss.str(); } );
  of.close();

  auto dt_y = [&,count=0](auto const& y) mutable { std::stringstream ss; ss<<times[count++]<<" "<<y; return ss.str(); };

  of.open( c.output_dir / "ee_vhll.dat" );
  std::transform( ee.begin() , ee.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();
  
  of.open( c.output_dir / "Emax_vhll.dat" );
  std::transform( Emax.begin() , Emax.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();
  
  of.open( c.output_dir / "H_vhll.dat" );
  std::transform( H.begin() , H.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();

  auto dx_y = [&,count=0](auto const& y) mutable { std::stringstream ss; ss<< fh.step.dx*(count++) <<" "<<y; return ss.str(); };

  of.open( c.output_dir / "E_vhll.dat" );
  std::transform( E.begin() , E.end() , std::ostream_iterator<std::string>(of,"\n") , dx_y );
  of.close();

  ublas::vector<double> rho  (c.Nx,0.);
  ublas::vector<double> rhoc (c.Nx,0.);
  {
    auto rhoh = fh.density();
    fft::spectrum_ hrhoh(c.Nx); hrhoh.fft(&rhoh[0]);
    fft::spectrum_ hE(c.Nx); hE.fft(&E[0]);
    fft::spectrum_ hrhoc(c.Nx);
    hrhoc[0] = I*kx[0]*hE[0] - hrhoh[0] + 1.;
    for ( auto i=1 ; i<c.Nx ; ++i ) {
      hrhoc[i] = I*kx[i]*hE[i] - hrhoh[i];
    }
    hrhoc.ifft(&rhoc[0]);
    rho = rhoc + rhoh;
  }
  of.open( c.output_dir / "rho_vhll.dat" );
  std::transform( rho.begin() , rho.end() , std::ostream_iterator<std::string>(of,"\n") , dx_y );
  of.close();

  of.open( c.output_dir / "uc_vhll.dat" );
  std::transform( uc.begin() , uc.end() , std::ostream_iterator<std::string>(of,"\n") , dx_y );
  of.close();

  auto Jh = fh.courant();
  ublas::vector<double> Jc (c.Nx,0.), Jtmp (c.Nx,0.);
  for ( auto i=0 ; i<c.Nx ; ++i ) {
    Jc[i] = rhoc[i]*uc[i];
  }
  Jtmp = Jc + Jh;
  of.open( c.output_dir / "J_vhll.dat" );
  std::transform( Jtmp.begin() , Jtmp.end() , std::ostream_iterator<std::string>(of,"\n") , dx_y );
  of.close();

  return 0;
}
