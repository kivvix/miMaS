#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include "miMaS/field.h"
#include "miMaS/weno.h"
#include "miMaS/fft.h"
#include "miMaS/array_view.h"
#include "miMaS/poisson.h"

namespace math = boost::math::constants;
const std::complex<double> & I = std::complex<double>(0.,1.);

#define SQ(X) ((X)*(X))
#define Xi(i) (i*f.step.dx+f.range.x_min)
#define Vk(k) (k*f.step.dv+f.range.v_min)

void
order_Hdt( double dt )
{
  std::size_t Nx = 42, Nv = 1024;
  field<double,1> f(boost::extents[Nv][Nx]);

  f.range.v_min = -12.; f.range.v_max = 12;
  f.step.dv = (f.range.v_max-f.range.v_min)/Nv;
  f.range.x_min = 0.; f.range.x_max = 20.*math::pi<double>();
  f.step.dx = (f.range.x_max-f.range.x_min)/Nx;
  
  ublas::vector<double> v (Nv,0.);
  ublas::vector<double> E (Nx,0.),rho(Nx);
  for ( std::size_t k=0 ; k<Nv ; ++k ) { v[k] = Vk(k); }

  const double l = f.range.x_max-f.range.x_min;
  ublas::vector<double> kx(Nx);
  for ( auto i=0 ; i<Nx/2 ; ++i ) { kx[i]    = 2.*math::pi<double>()*i/l; }
  for ( int i=-Nx/2 ; i<0 ; ++i ) { kx[Nx+i] = 2.*math::pi<double>()*i/l; }
  
  double np = 0.9 , nb = 0.2 , ui = 4.5;
  for (field<double,2>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<f.size(1) ; ++i ) {
      f[k][i] = ( std::exp(-0.5*SQ(Vk(k)))*np/std::sqrt(2.*math::pi<double>()) + nb/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-ui)/0.25) )*(1.+0.04*std::cos(0.3*Xi(i)));
      //f[k][i] = ( std::exp(-0.5*SQ(Vk(k)))*np/std::sqrt(2.*math::pi<double>()) + nb/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-ui)/0.25) )*(1.+0.01*std::cos(0.3*Xi(i)));
    }
  }


  std::vector<double> H; //(int(std::ceil(Tf/dt)+1),0.);
  std::vector<double> Emax; //(int(std::ceil(Tf/dt)+1),0.);
  std::vector<double> ee;

  poisson<double> poisson_solver(Nx,l);
  rho = f.density();
  E = poisson_solver(rho);
  H.push_back( energy(f,E) );

  //const double Tf = 3.8;
  const double Tf = 10.;
  int i_t=0;

  ee.push_back(0.);
  for ( auto i=0 ; i<Nx ; ++i ) {
    ee[i_t] += SQ(E(i))*f.step.dx;
  }
/*
  std::cout << "dt " << dt << std::endl;
  std::cout << "dx " << f.step.dx << std::endl;
  std::cout << "dv " << f.step.dv << std::endl;
  std::cout << "Tf " << Tf << std::endl;
*/

  field<double,1> f1(tools::array_view<const std::size_t>(f.shape(),2)),
                  f2(tools::array_view<const std::size_t>(f.shape(),2)),
                  f3(tools::array_view<const std::size_t>(f.shape(),2)),
                  f4(tools::array_view<const std::size_t>(f.shape(),2)),
                  f5(tools::array_view<const std::size_t>(f.shape(),2)),
                  f6(tools::array_view<const std::size_t>(f.shape(),2)),
                  f7(tools::array_view<const std::size_t>(f.shape(),2));
  fft::spectrum_ hf(Nx),hf1(Nx),hf2(Nx),hf3(Nx),hf4(Nx),hf5(Nx),hf6(Nx),hf7(Nx),
                 hEdvf(Nx);

#define L (-v(k)*I*kx[i])
  while (  i_t*dt < Tf ) {
    std::cout<<" \r"<<i_t<<" / "<<int(Tf/dt)<<std::flush;


    /**
    // RK(3,3)
    // RK(3,3) eq19
    E = poisson_solver(f.density());
    field<double,1> Edvf = weno::trp_v(f,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        //hf1[i] = std::exp(L*dt)*( hf[i]-dt*hEdvf[i] );
        hf1[i] = 0.5*std::exp((2./3.)*L*dt)*hf[i] + 0.5*std::exp((2./3.)*dt*L)*( hf[i] - (4./3.)*dt*hEdvf[i] );
      }
      hf1.ifft(&(f1[k][0]));
    }

    E = poisson_solver(f1.density());
    Edvf = weno::trp_v(f1,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        //hf2[i] = 0.75*std::exp(0.5*L*dt)*hf[i] + 0.25*std::exp(-0.5*L*dt)*( hf1[i]-dt*hEdvf[i] );
        hf2[i] = (2./3.)*std::exp((2./3.)*dt*L)*hf[i] + (1./3.)*( hf1[i] - (4./3.)*dt*hEdvf[i] );
      }
      hf2.ifft(&(f2[k][0]));
    }

    E = poisson_solver(f2.density());
    Edvf = weno::trp_v(f2,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        //hf[i] = (1./3.)*std::exp(L*dt)*hf[i] + (2./3.)*std::exp(0.5*L*dt)*( hf2[i]-dt*hEdvf[i] );
        hf[i] = (59./128.)*std::exp(L*dt)*hf[i] + (15./128.)*std::exp(L*dt)*( 2.*hf1[i]*std::exp(-(2./3.)*L*dt) - hf[i] ) + (27./64.)*std::exp((1./3.)*dt*L)*( hf2[i] - (4./3.)*dt*hEdvf[i] );
      }
      hf.ifft(&(f[k][0]));
    }
    **/
    /**
    // RK(4,4)
    // RK(4,4) 3/8 rule
    E = poisson_solver( f.density() );
    field<double,1> Edvf = weno::trp_v(f,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        //hf1[i] = std::exp(0.5*L*dt)*( hf[i] - 0.5*dt*hEdvf[i] );
        hf1[i] = std::exp((1./3.)*L*dt)*( hf[i] - (1./3.)*dt*hEdvf[i] );
      }

      hf1.ifft(&(f1[k][0]));
    }

    E = poisson_solver( f1.density() );
    Edvf = weno::trp_v(f1,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        //hf2[i] = std::exp(0.5*L*dt)*hf[i] - 0.5*dt*hEdvf[i] ;
        hf2[i] = 2.*std::exp((2./3.)*L*dt)*hf[i] - std::exp((1./3.)*L*dt)*hf1[i] - dt*std::exp((1./3.)*L*dt)*hEdvf[i];
      }
      
      hf2.ifft(&(f2[k][0]));
    }

    E = poisson_solver( f2.density() );
    Edvf = weno::trp_v(f2,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        //hf3[i] = std::exp(L*dt)*hf[i] - dt*std::exp(0.5*L*dt)*hEdvf[i];
        hf3[i] = 2.*std::exp((2./3.)*L*dt)*hf1[i] - std::exp((1./3.)*L*dt)*hf2[i] - dt*std::exp((1./3.)*L*dt)*hEdvf[i];
      }

      hf3.ifft(&(f3[k][0]));
    }

    E = poisson_solver( f3.density() );
    Edvf  = weno::trp_v(f3,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hf3.fft(&(f3[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        //hf[i] = -(1./3.)*std::exp(L*dt)*hf[i] + (1./3.)*std::exp(0.5*L*dt)*hf1[i] + (2./3.)*std::exp(0.5*L*dt)*hf2[i] + (1./3.)*hf3[i] - (1./6.)*dt*hEdvf[i];
        hf[i] = -(1./8.)*std::exp(L*dt)*hf[i] + 0.75*std::exp((1./3.)*L*dt)*hf2[i] + (3./8.)*hf3[i] - (1./8.)*dt*hEdvf[i];
      }

      hf.ifft(&(f[k][0]));
    }

    **/
    /**
    // RK (5,3)
    E = poisson_solver(f.density());
    field<double,1> Edvf = weno::trp_v(f,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf1[i] = std::exp((1./7.)*L*dt)*hf[i] - (1./7.)*dt*std::exp((1./7.)*L*dt)*hEdvf[i];
      }
      hf1.ifft(&(f1[k][0]));
    }

    E = poisson_solver(f1.density());
    Edvf = weno::trp_v(f1,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf2[i] = std::exp((3./16.)*L*dt)*hf[i] - (3./16.)*dt*std::exp((5./112.)*L*dt)*hEdvf[i];
      }
      hf2.ifft(&(f2[k][0]));
    }

    E = poisson_solver(f2.density());
    Edvf = weno::trp_v(f2,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf3[i] = std::exp((1./3.)*L*dt)*hf[i] - (1./3.)*dt*std::exp((7./48.)*L*dt)*hEdvf[i];
      }
      hf3.ifft(&(f3[k][0]));
    }

    E = poisson_solver(f3.density());
    Edvf = weno::trp_v(f3,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf4[i] = std::exp((2./3.)*L*dt)*hf[i] - (2./3.)*dt*std::exp((1./3.)*L*dt)*hEdvf[i];
      }
      hf4.ifft(&(f4[k][0]));
    }

    E = poisson_solver(f4.density());
    Edvf = weno::trp_v(f4,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf[i] = -0.75*std::exp(L*dt)*hf[i] + 1.75*std::exp((6./7.)*L*dt)*hf1[i] - 0.75*dt*std::exp((1./3.)*L*dt)*hEdvf[i];
      }
      hf.ifft(&(f[k][0]));
    }
    **/
    /**
    // DP5
    E = poisson_solver(f.density());
    field<double,1> Edvf = weno::trp_v(f,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf1[i] = std::exp((1./5.)*L*dt)*hf[i] - (1./5.)*dt*std::exp((1./5.)*L*dt)*hEdvf[i];
      }
      hf1.ifft(&(f1[k][0]));
    }

    E = poisson_solver(f1.density());
    Edvf = weno::trp_v(f1,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf2[i] = (5./8.)*std::exp((3./10.)*L*dt)*hf[i] + (3./8.)*std::exp((1./10.)*L*dt)*hf1[i] - (9./40.)*dt*std::exp((1./10.)*L*dt)*hEdvf[i];
      }
      hf2.ifft(&(f2[k][0]));
    }

    E = poisson_solver(f2.density());
    Edvf = weno::trp_v(f2,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf3[i] = (175./27.)*std::exp((4./5.)*L*dt)*hf[i] + (100./9.)*std::exp((3./5.)*L*dt)*hf1[i] - (448./27.)*std::exp(0.5*L*dt)*hf2[i] - (32./9.)*dt*std::exp(0.5*L*dt)*hEdvf[i];
      }
      hf3.ifft(&(f3[k][0]));
    }

    E = poisson_solver(f3.density());
    Edvf = weno::trp_v(f3,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hf3.fft(&(f3[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf4[i] = (3551./6561.)*std::exp((8./9.)*L*dt)*hf[i] + (7420./2187.)*std::exp((31./45.)*L*dt)*hf1[i] - (37376./6561.)*std::exp((53./90.)*L*dt)*hf2[i] + (2014./729.)*std::exp((4./45.)*L*dt)*hf3[i] + (212./729.)*dt*std::exp((4./45.)*L*dt)*hEdvf[i];
      }
      hf4.ifft(&(f4[k][0]));
    }

    E = poisson_solver(f4.density());
    Edvf = weno::trp_v(f4,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hf3.fft(&(f3[k][0]));
      hf4.fft(&(f4[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf5[i] = (313397./335808.)*std::exp(L*dt)*hf[i] + (424025./55968.)*std::exp((4./5.)*L*dt)*hf1[i] - (61400./5247.)*std::exp((7./10.)*L*dt)*hf2[i] + (96075./18656.)*std::exp((1./5.)*L*dt)*hf3[i] - (35721./37312.)*std::exp((1./9.)*L*dt)*hf4[i] + (5103./18656.)*dt*std::exp((1./9.)*L*dt)*hEdvf[i];
      }
      hf5.ifft(&(f5[k][0]));
    }

    E = poisson_solver(f5.density());
    Edvf = weno::trp_v(f5,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hf3.fft(&(f3[k][0]));
      hf4.fft(&(f4[k][0]));
      hf5.fft(&(f5[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf[i] = -(563./3456.)*std::exp(L*dt)*hf[i] - (575./252.)*std::exp((4./5.)*L*dt)*hf1[i] + (31400./10017.)*std::exp((7./10.)*L*dt)*hf2[i] + (325./1344.)*std::exp((1./5.)*L*dt)*hf3[i] - (7533./6784.)*std::exp((1./9.)*L*dt)*hf4[i] + (33./28.)*hf5[i] - (11./84.)*dt*hEdvf[i];
      }
      hf.ifft(&(f[k][0]));
    }

    **/
    /**/
    // RK(8,6)
    E = poisson_solver(f.density());
    field<double,1> Edvf = weno::trp_v(f,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf1[i] = std::exp((1./9.)*L*dt)*hf[i] - (1./9.)*dt*std::exp((1./9.)*dt*L)*hEdvf[i];
      }
      hf1.ifft(&(f1[k][0]));
    }

    E = poisson_solver(f1.density());
    Edvf = weno::trp_v(f1,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf2[i] = (5./8.)*std::exp((1./6.)*L*dt)*hf[i] + (3./8.)*std::exp((1./18.)*L*dt)*hf1[i] - (1./8.)*dt*std::exp((1./18.)*L*dt)*hEdvf[i];
      }
      hf2.ifft(&(f2[k][0]));
    }

    E = poisson_solver(f2.density());
    Edvf = weno::trp_v(f2,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf3[i] = 2.*std::exp((1./3.)*L*dt)*hf[i] + 3.*std::exp((2./9.)*L*dt)*hf1[i] -4.*std::exp((1./6.)*L*dt)*hf2[i] - (2./3.)*dt*std::exp((1./6.)*L*dt)*hEdvf[i];
      }
      hf3.ifft(&(f3[k][0]));
    }

    E = poisson_solver(f3.density());
    Edvf = weno::trp_v(f3,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hf3.fft(&(f3[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf4[i] = (305./1268.)*std::exp(0.5*L*dt)*hf[i] + (2817./1268.)*std::exp((7./18.)*L*dt)*hf1[i] - (927./317.)*std::exp((1./3.)*L*dt)*hf2[i] + (927./634.)*std::exp((1./6.)*L*dt)*hf3[i] - (321./1268.)*dt*std::exp((1./6.)*L*dt)*hEdvf[i];
      }
      hf4.ifft(&(f4[k][0]));
    }

    E = poisson_solver(f4.density());
    Edvf = weno::trp_v(f4,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hf3.fft(&(f3[k][0]));
      hf4.fft(&(f4[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf5[i] = (3191./321.)*std::exp((2./3.)*L*dt)*hf[i] - (2436./107.)*std::exp((5./9.)*L*dt)*hf1[i] - (2404./107.)*std::exp(0.5*L*dt)*hf2[i] + (12330./107.)*std::exp((1./3.)*L*dt)*hf3[i] - (25340./321.)*std::exp((1./6.)*L*dt)*hf4[i] - 8.*dt*std::exp((1./6.)*L*dt)*hEdvf[i];
      }
      hf5.ifft(&(f5[k][0]));
    }

    E = poisson_solver(f5.density());
    Edvf = weno::trp_v(f5,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hf3.fft(&(f3[k][0]));
      hf4.fft(&(f4[k][0]));
      hf5.fft(&(f5[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf6[i] = -(15130159./6286464.)*std::exp((5./6.)*L*dt)*hf[i] + (2014319./349248.)*std::exp((13./18.)*L*dt)*hf1[i] + (1194095./523872.)*std::exp((2./3.)*L*dt)*hf2[i] - (1471057./116416.)*std::exp(0.5*L*dt)*hf3[i] + (12601453./1571616.)*std::exp((1./3.)*L*dt)*hf4[i] - (433./19584.)*std::exp((1./6.)*L*dt)*hf5[i] - (33./1088.)*dt*std::exp((1./6.)*L*dt)*hEdvf[i];
      }
      hf6.ifft(&(f6[k][0]));
    }

    E = poisson_solver(f6.density());
    Edvf = weno::trp_v(f6,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hf3.fft(&(f3[k][0]));
      hf4.fft(&(f4[k][0]));
      hf5.fft(&(f5[k][0]));
      hf6.fft(&(f6[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf7[i] = (805187./78966.)*std::exp(L*dt)*hf[i] - (2263766./48257.)*std::exp((8./9.)*L*dt)*hf1[i] + (2745422./144771.)*std::exp((5./6.)*L*dt)*hf2[i] + (2271108./48257.)*std::exp((2./3.)*L*dt)*hf3[i] - (13115270./434313.)*std::exp(0.5*L*dt)*hf4[i] - (227./2706.)*std::exp((1./3.)*L*dt)*hf5[i] + (888./451.)*std::exp((1./6.)*L*dt)*hf6[i] - (36./41.)*dt*std::exp((1./6.)*L*dt)*hEdvf[i];
      }
      hf7.ifft(&(f7[k][0]));
    }

    E = poisson_solver(f7.density());
    Edvf = weno::trp_v(f7,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hf3.fft(&(f3[k][0]));
      hf4.fft(&(f4[k][0]));
      hf5.fft(&(f5[k][0]));
      hf6.fft(&(f6[k][0]));
      hf7.fft(&(f7[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf[i] = -(193999./179760.)*std::exp(L*dt)*hf[i] + (2487363./329560.)*std::exp((8./9.)*L*dt)*hf1[i] - (847909./164780.)*std::exp((5./6.)*L*dt)*hf2[i] - (1600251./329560.)*std::exp((2./3.)*L*dt)*hf3[i] + (362713./98868.)*std::exp(0.5*L*dt)*hf4[i] + (109./1232.)*std::exp((1./3.)*L*dt)*hf5[i] + (186./385.)*std::exp((1./6.)*L*dt)*hf6[i] + (41./140.)*hf7[i] - (41./840.)*dt*hEdvf[i];
      }
      hf.ifft(&(f[k][0]));
    }
    /**/

    E = poisson_solver(f.density());
    ee.push_back(0.);
    for ( auto i=0 ; i<Nx ; ++i ) {
      ee[i_t] += SQ(E(i))*f.step.dx;
    }
    H.push_back( energy(f,E) );
    //Emax.push_back( std::abs(*std::max_element( E.begin() , E.end() , [](double a,double b){return std::abs(a)<std::abs(b);} )) );
    ++i_t;
  }

  //E = poisson_solver(f.density());
  //H.push_back( energy(f,E) );
  //Emax(i_t) = std::abs(*std::max_element( E.begin() , E.end() , [](double a,double b){return std::abs(a)<std::abs(b);} ));

#undef L
  double h = std::abs(*std::max_element( H.begin() , H.end() , [&](double a,double b){return ( std::abs((a-H[0])/std::abs(H[0])) < std::abs((b-H[0])/std::abs(H[0])) );} ));
  //std::cout << dt << " " << std::abs((h-H[0])/std::abs(H[0])) << std::endl;

  std::ofstream of;
  std::stringstream ss; ss << "H.dat";
  of.open(ss.str());
  std::size_t count=0;
  std::transform(H.begin(),H.end(),std::ostream_iterator<std::string>(of,"\n"),[&,count=0](auto h)mutable{std::stringstream ss; ss<<dt*count++<<" "<<h;return ss.str();});
  of.close();
  ss.str(""); ss << "ee.dat";
  of.open(ss.str());
  std::transform(ee.begin(),ee.end(),std::ostream_iterator<std::string>(of,"\n"),[&,count=0](auto e)mutable{std::stringstream ss; ss<<dt*count++<<" "<<e;return ss.str();});
  of.close();
}

int
main (int,char**)
{
  // dt = 1.606*24/1024/Emax
  // 1.606: min CFL
  // 24 lenght of v domain
  // 1024 number of points
  // Emax = (0.04/0.3)=0.133 ou plus petite perturbation (0.01/0.3)=0.035
  // dt = 0.28 ou 1.12
  //double dt = 0.9;
  //double dt = 0.7;
  //double dt = 0.95;
  double dt = 0.05;
  //double dt = 0.23525;
  for (int i=1;i<=1;++i) {
    order_Hdt(dt/i);
  }

  return 0;
}

