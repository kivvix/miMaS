#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <fstream>
#include <valarray>
#include <sstream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include "miMaS/field.h"
#include "miMaS/weno.h"
#include "miMaS/fft.h"
#include "miMaS/array_view.h"
#include "miMaS/poisson.h"

#ifndef SIGMA
#define SIGMA (2.*std::sqrt(2.))
#endif
#define E_MAX (0.6)

/*
template <unsigned int i>
struct phi
{
  static std::complex<double>
  operator () ( std::complex<double> const & z ) {
    static std::valarray<std::complex<double>> coeff(i);
    coeff[0] = 1.;

    for ( unsigned int k=1 ; k<coeff.size() ; ++k ) {
      coeff[k] = coeff[k-1] * z / (double(k));
    }

    return (std::exp(z) - std::accumulate( std::begin(coeff) , std::end(coeff) , std::complex<double>(0.,0.) ))/(std::pow(z,i));
  }
};
*/
template <unsigned int i>
std::complex<double>
phi ( std::complex<double> const & _z )
{
  std::valarray<std::complex<double>> coeff(i);
  coeff[0] = 1.;

  std::complex<double> z = _z;
  if ( _z == 0. ) { z = std::complex<double>(1.,0.); }

  for ( unsigned int k=1 ; k<coeff.size() ; ++k ) {
    coeff[k] = coeff[k-1] * z / (double(k));
  }
  //std::copy(std::begin(coeff),std::end(coeff),std::ostream_iterator<std::complex<double>>(std::cout," . "));
  //std::cout << std::endl;

  if ( z != 0. ) {
    return (std::exp(z) - std::accumulate( std::begin(coeff) , std::end(coeff) , std::complex<double>(0.,0.) ))/(std::pow(z,i));  
  }
  return coeff[i-1];
}

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

namespace math = boost::math::constants;
const std::complex<double> & I = std::complex<double>(0.,1.);

#define SQ(X) ((X)*(X))
#define Xi(i) (i*f.step.dx+f.range.x_min)
#define Vk(k) (k*f.step.dv+f.range.v_min)

int main(int,char**)
{
	std::size_t Nx = 135, Nv = 256 , Nb_iter=10;
	field<double,1> f(boost::extents[Nv][Nx]);

	f.range.v_min = -8.; f.range.v_max = 8.;
	f.step.dv = (f.range.v_max-f.range.v_min)/Nv;

  const double Kx = 0.3;
  f.range.x_min = 0.; f.range.x_max = 2./Kx*math::pi<double>();
	//f.range.x_min = 0.; f.range.x_max = 20.*math::pi<double>();
  //f.range.x_min = 0.; f.range.x_max = 4.0*math::pi<double>();
  //f.range.x_min = -8.; f.range.x_max = 8.;
	f.step.dx = (f.range.x_max-f.range.x_min)/Nx;

  // SIGMA is the CFL number 0.45 is E_max in our test case
	const double dt = 0.05;//SIGMA*f.step.dv/E_MAX; //1.606*f.step.dv/0.6; //0.005; //1.606*f.step.dv/0.6; //0.5*6.*math::pi<double>()/(Nv*f.range.v_max);

  //field<double,1> f_sol = f;
  //field<double,1> f_ini = f;
	
	ublas::vector<double> v (Nv,0.);
  ublas::vector<double> E (Nx,0.),rho(Nx);
  for ( std::size_t k=0 ; k<Nv ; ++k ) { v[k] = Vk(k); }
  //for ( std::size_t i=0 ; i<Nx ; ++i ) { E[i] = -Xi(i); }

	const double l = f.range.x_max-f.range.x_min;
	ublas::vector<double> kx(Nx);
	//for ( auto i=0 ; i<Nx/2+1 ; ++i )   { kx[i] = 2.*math::pi<double>()*i/l; }
	//for ( auto i=0 ; i<((Nx/2)) ; ++i ) { kx[i+Nx/2+1] = -kx[Nx/2-i]; }
  for ( auto i=0 ; i<Nx/2 ; ++i ) { kx[i]    = 2.*math::pi<double>()*i/l; }
  for ( int i=-Nx/2 ; i<0 ; ++i ) { kx[Nx+i] = 2.*math::pi<double>()*i/l; }
	

  double ui = 3.4;
  double alpha = 0.1;
  double Tc = 0.0001;

  //double np = 0.9 , nb = 0.2 , ui = 4.5;
  for (field<double,2>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<f.size(1) ; ++i ) {
      //f[k][i] = ( std::exp(-0.5*SQ(Vk(k)))*np/std::sqrt(2.*math::pi<double>()) + nb/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-ui)/0.25) )*(1.+0.04*std::cos(0.5*Xi(i)));
      
      // Bump on Tail
      f[k][i] = ( std::exp(-0.5*SQ(Vk(k)))*np/std::sqrt(2.*math::pi<double>()) + nb/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-ui)/0.25) )*(1.+0.04*std::cos(0.3*Xi(i)));
      // Landau dumpping test
      //f[k][i] = (1./std::sqrt(2*math::pi<double>()))*std::exp(-0.5*SQ(Vk(k)))*(1.+0.001*std::cos(0.5*Xi(i)));
      
      //f[k][i] = std::exp( -SQ(Xi(i)-6) );
      //f[k][i] = std::exp(-SQ(Xi(i)-3)/0.5 - SQ(Vk(k))/2.);
      //f_sol[k][i] = std::exp(-SQ(Xi(i)-3)/0.5 - SQ(Vk(k))/2.);
      //f[k][i] = SQ(v(k))*std::exp(-0.5*SQ(Vk(k)))/std::sqrt(2.*math::pi<double>())*(1.+0.01*std::cos(0.5*Xi(i)));
      //f_ini[k][i] = f[k][i];
      //f_sol[k][i] =  std::exp(-SQ(Xi(i)-3.-v[k]*Nb_iter*dt)/0.5 - SQ(Vk(k)-E[i]*Nb_iter*dt)/2.);

      // triple bump
      //f[k][i] = ( 0.5*alpha/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-ui)) + 0.5*alpha/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)+ui)) )*(1.+0.04*std::cos(Kx*Xi(i))) + ((1-alpha)/(std::sqrt(2.*math::pi<double>()*Tc))*std::exp(-0.5*SQ(Vk(k))/Tc));
    }
  }
  f.write("vphl/kin/init.dat");

  poisson<double> poisson_solver(Nx,l);
  rho = f.density();
  E = poisson_solver(rho);
  

  //double Tf = 60.;//2*math::pi<double>();
  double Tf = 40.;
  int i_t=0;

  std::cout << "Nx: " << Nx << "\n";
  std::cout << "Nv: " << Nv << "\n";
  std::cout << "v_min: " << f.range.v_min << "\n";
  std::cout << "v_max: " << f.range.v_max << "\n";
  std::cout << "x_min: " << f.range.x_min << "\n";
  std::cout << "x_max: " << f.range.x_max << "\n";
  std::cout << "dt: " << dt << "\n";
  std::cout << "dx: " << f.step.dx << "\n";
  std::cout << "dv: " << f.step.dv << "\n";
  std::cout << "Tf: " << Tf << "\n";
  std::cout << "f_0: " << "\"bot\"" << "\n";
  std::cout << std::endl;

  std::ofstream info("info.yaml");

  info << "Nx: " << Nx << "\n";
  info << "Nv: " << Nv << "\n";
  info << "v_min: " << f.range.v_min << "\n";
  info << "v_max: " << f.range.v_max << "\n";
  info << "x_min: " << f.range.x_min << "\n";
  info << "x_max: " << f.range.x_max << "\n";
  info << "dt: " << dt << "\n";
  info << "dx: " << f.step.dx << "\n";
  info << "dv: " << f.step.dv << "\n";
  info << "Tf: " << Tf << "\n";
  info << "f_0: " << "\"bot\"" << "\n";
  info << std::endl;
  info.close();

  ublas::vector<double> ee(int(std::ceil(Tf/dt))+1,0.);
  ublas::vector<double> Emax(int(std::ceil(Tf/dt))+1,0.);
  ublas::vector<double> H(int(std::ceil(Tf/dt))+1,0.);
  //ublas::vector<double> ee(Nb_iter);
  //ublas::vector<double> Emax(Nb_iter);
  //ublas::vector<double> H(Nb_iter);

  field<double,1> f1(tools::array_view<const std::size_t>(f.shape(),2)),f2(tools::array_view<const std::size_t>(f.shape(),2)),f3(tools::array_view<const std::size_t>(f.shape(),2)),f4(tools::array_view<const std::size_t>(f.shape(),2)),f5(tools::array_view<const std::size_t>(f.shape(),2));
  fft::spectrum_ hf(Nx),hf1(Nx),hf2(Nx),hf3(Nx),hf4(Nx),hf5(Nx),hEdvf(Nx),hEdvf1(Nx),hEdvf2(Nx),hEdvf3(Nx),hEdvf4(Nx);


  rho = f.density();
  E = poisson_solver(rho);
  //for ( int i=0 ; i<E.size() ; ++i ) { E[i] = 1.; }
  //for ( int i=0 ; i<E.size() ; ++i ) { std::cout << E[i] << " "; }
  //  std::cout << std::endl;

  //ee(i_t) = 0.;
  //for ( auto i=0 ; i<Nx ; ++i ) {
  //  ee(i_t) += SQ(E(i))*f.step.dx;
  //}
  //ee(i_t) = std::sqrt(ee(i_t));
  H(i_t) = energy(f,E);


#define L (-v(k)*I*kx[i])
//#define L (0.)
  while (  i_t*dt < Tf ) {
  //while (  i_t < Nb_iter ) {
  	//if (t%32==0) { std::cout<<"\r"<<t<<" "<<std::flush ; }
    std::cout<<" ["<<std::setw(5)<<i_t<<"] "<<i_t*dt<<"\r"<<std::flush;
    Emax(i_t) = std::abs(*std::max_element( E.begin() , E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} ));
    ee(i_t) = 0.;
    for ( auto i=0 ; i<Nx ; ++i ) { ee(i_t) += SQ(E(i))*f.step.dx; }
    ee(i_t) = std::sqrt(ee(i_t));

    /**
    // exprk(2,2) =============================================================
    #define SCHEME "expRK22"
    // SIGMA = 0.551 (10^-2)
    E = poisson_solver(f.density());
    field<double,1> Edvf = o2::trp_v(f,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf1[i] = std::exp(-L*dt)*hf[i] + dt*phi<1>(dt*L)*hEdvf[i];
      }
      hf1.ifft(&(f1[k][0]));
    }

    E = poisson_solver(f1.density());
    field<double,1> Edvf1 = o2::trp_v(f1,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&f[k][0]);
      hf1.fft(&f1[k][0]);
      hEdvf.fft(&(Edvf[k][0]));
      hEdvf1.fft(&(Edvf1[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf[i] = std::exp(-L*dt)*hf[i] + dt*( (phi<1>(dt*L)-phi<2>(dt*L))*hEdvf[i] + phi<2>(dt*L)*hEdvf1[i] );
      }

      hf.ifft(&(f[k][0]));
    }

    **/
    /**
    // Cox-Matthews ===========================================================
    #define SCHEME "CM"
    // SIGMA = 0.450 (10^-2)
    E = poisson_solver( f.density() );
    field<double,1> Edvf = o2::trp_v(f,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf1[i] = std::exp(-0.5*L*dt)*hf[i] + 0.5*dt*phi<1>(-0.5*L*dt)*hEdvf[i];
      }

      hf1.ifft(&(f1[k][0]));
    }

    E = poisson_solver( f1.density() );
    field<double,1> Edvf1 = o2::trp_v(f1,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hEdvf.fft(&(Edvf[k][0]));
      hEdvf1.fft(&(Edvf1[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf2[i] = std::exp(-0.5*L*dt)*hf[i] + 0.5*dt*phi<1>(-0.5*L*dt)*hEdvf1[i];
      }
      
      hf2.ifft(&(f2[k][0]));
    }

    E = poisson_solver( f2.density() );
    field<double,1> Edvf2 = o2::trp_v(f2,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hEdvf.fft(&(Edvf[k][0]));
      hEdvf1.fft(&(Edvf1[k][0]));
      hEdvf2.fft(&(Edvf2[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf3[i] = std::exp(-L*dt)*hf[i] + 0.5*dt*phi<1>(-0.5*L*dt)*(std::exp(-0.5*L*dt)-1.)*hEdvf[i] + dt*phi<1>(-0.5*L*dt)*hEdvf2[i];
      }

      hf3.ifft(&(f3[k][0]));
    }

    E = poisson_solver( f3.density() );
    field<double,1> Edvf3  = o2::trp_v(f3,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hf3.fft(&(f3[k][0]));
      hEdvf.fft(&(Edvf[k][0]));
      hEdvf1.fft(&(Edvf1[k][0]));
      hEdvf2.fft(&(Edvf2[k][0]));
      hEdvf3.fft(&(Edvf3[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf4[i] = std::exp(-L*dt)*hf[i] + dt*( (phi<1>(-L*dt)-3.*phi<2>(-L*dt)+4.*phi<3>(-L*dt))*hEdvf[i]
                                             + (2.*phi<2>(-L*dt)-4.*phi<3>(-L*dt))*(hEdvf1[i]+hEdvf2[i])
                                             + (-phi<2>(-L*dt)+4.*phi<3>(-L*dt))*hEdvf3[i] );
      }

      hf4.ifft(&(f[k][0]));
    }

    **/
    /**
    // Krogstad ===============================================================
    #define SCHEME "K"
    // SIGMA = 0.200 (10^-2)
    E = poisson_solver( f.density() );
    field<double,1> Edvf = o2::trp_v(f,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf1[i] = std::exp(-0.5*L*dt)*hf[i] + 0.5*dt*phi<1>(-0.5*L*dt)*hEdvf[i];
      }

      hf1.ifft(&(f1[k][0]));
    }

    E = poisson_solver( f1.density() );
    field<double,1> Edvf1 = o2::trp_v(f1,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hEdvf.fft(&(Edvf[k][0]));
      hEdvf1.fft(&(Edvf1[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf2[i] = std::exp(-0.5*L*dt)*hf[i] + dt*(0.5*phi<1>(-0.5*L*dt)-phi<2>(-0.5*L*dt))*hEdvf[i] + dt*phi<2>(-0.5*L*dt)*hEdvf1[i];
      }
      
      hf2.ifft(&(f2[k][0]));
    }

    E = poisson_solver( f2.density() );
    field<double,1> Edvf2 = o2::trp_v(f2,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hEdvf.fft(&(Edvf[k][0]));
      hEdvf1.fft(&(Edvf1[k][0]));
      hEdvf2.fft(&(Edvf2[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf3[i] = std::exp(-L*dt)*hf[i] + dt*(phi<1>(-L*dt)-2.*phi<2>(-0.5*L*dt))*hEdvf[i] + 2.*dt*phi<2>(-L*dt)*hEdvf2[i];
      }

      hf3.ifft(&(f3[k][0]));
    }

    E = poisson_solver( f3.density() );
    field<double,1> Edvf3  = o2::trp_v(f3,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hf3.fft(&(f3[k][0]));
      hEdvf.fft(&(Edvf[k][0]));
      hEdvf1.fft(&(Edvf1[k][0]));
      hEdvf2.fft(&(Edvf2[k][0]));
      hEdvf3.fft(&(Edvf3[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf4[i] = std::exp(-L*dt)*hf[i] + dt*( (phi<1>(-L*dt)-3.*phi<2>(-L*dt)+4.*phi<3>(-L*dt))*hEdvf[i]
                                             + (2.*phi<2>(-L*dt)-4.*phi<3>(-L*dt))*(hEdvf1[i]+hEdvf2[i])
                                             + (-phi<2>(-L*dt)+4.*phi<3>(-L*dt))*hEdvf3[i] );
      }

      hf4.ifft(&(f[k][0]));
    }

    **/
    /**
    // Hochbruck-Ostermann ====================================================
    #define SCHEME "HO"
    // SIGMA = 0.501 (10^-2)
    E = poisson_solver( f.density() );
    field<double,1> Edvf = o2::trp_v(f,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf1[i] = std::exp(-0.5*L*dt)*hf[i] + 0.5*dt*phi<1>(-0.5*L*dt)*hEdvf[i];
      }

      hf1.ifft(&(f1[k][0]));
    }

    E = poisson_solver( f1.density() );
    field<double,1> Edvf1 = o2::trp_v(f1,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hEdvf.fft(&(Edvf[k][0]));
      hEdvf1.fft(&(Edvf1[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf2[i] = std::exp(-0.5*L*dt)*hf[i] + dt*(0.5*phi<1>(-0.5*L*dt)-phi<2>(-0.5*L*dt))*hEdvf[i] + dt*phi<2>(-0.5*L*dt)*hEdvf1[i];
      }
      
      hf2.ifft(&(f2[k][0]));
    }

    E = poisson_solver( f2.density() );
    field<double,1> Edvf2 = o2::trp_v(f2,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hEdvf.fft(&(Edvf[k][0]));
      hEdvf1.fft(&(Edvf1[k][0]));
      hEdvf2.fft(&(Edvf2[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf3[i] = std::exp(-L*dt)*hf[i] + dt*(phi<1>(-L*dt)-2.*phi<2>(-L*dt))*hEdvf[i] + dt*phi<2>(-L*dt)*hEdvf1[i] + dt*phi<2>(-L*dt)*hEdvf2[i];
      }

      hf3.ifft(&(f3[k][0]));
    }

    E = poisson_solver( f3.density() );
    field<double,1> Edvf3  = o2::trp_v(f3,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hf3.fft(&(f3[k][0]));
      hEdvf.fft(&(Edvf[k][0]));
      hEdvf1.fft(&(Edvf1[k][0]));
      hEdvf2.fft(&(Edvf2[k][0]));
      hEdvf3.fft(&(Edvf3[k][0]));

#define a52 (0.5*phi<2>(-0.5*L*dt)-phi<3>(-L*dt)+0.25*phi<2>(-L*dt)-0.5*phi<3>(-0.5*L*dt))
#define a54 (0.25*phi<2>(-0.5*L*dt)-a52)
      for ( auto i=0 ; i<Nx ; ++i ) {
        hf4[i] = std::exp(-0.5*L*dt)*hf[i] + dt*(0.5*phi<1>(-0.5*L*dt)-2.*a52-a54)*hEdvf[i] + dt*a52*(hEdvf1[i]+hEdvf2[i]) + dt*(0.25*phi<2>(-0.5*L*dt)-a52)*hEdvf3[i];
      }
#undef a54
#undef a52

      hf4.ifft(&(f4[k][0]));
    }

    E = poisson_solver( f4.density() );
    field<double,1> Edvf4  = o2::trp_v(f4,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hf3.fft(&(f3[k][0]));
      hf4.fft(&(f4[k][0]));
      hEdvf.fft(&(Edvf[k][0]));
      hEdvf1.fft(&(Edvf1[k][0]));
      hEdvf2.fft(&(Edvf2[k][0]));
      hEdvf3.fft(&(Edvf3[k][0]));
      hEdvf4.fft(&(Edvf4[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf5[i] = std::exp(-L*dt)*hf[i] + dt*( (phi<1>(-L*dt)-3.*phi<2>(-L*dt)+4.*phi<3>(-L*dt))*hEdvf[i]
                                             + (-phi<2>(-L*dt)+4.*phi<3>(-L*dt))*hEdvf3[i]
                                             + (4.*phi<2>(-L*dt)-8.*phi<3>(-L*dt))*hEdvf4[i] );
      }

      hf5.ifft(&(f[k][0]));
    }
    **/
    /**
    // RK(3,2) best ===========================================================
    #define SCHEME "RK32"
    // SIGMA = 2. (y_max)
    // SIGMA = 1.344 (WENO)
    E = poisson_solver(f.density());
    field<double,1> Edvf = o2::trp_v(f,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf1[i] = std::exp(0.5*L*dt)*( hf[i]-0.5*dt*hEdvf[i] );
      }
      hf1.ifft(&(f1[k][0]));
    }

    E = poisson_solver(f1.density());
    Edvf = o2::trp_v(f1,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      //hf1.fft(&(f1[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf2[i] = std::exp(0.5*L*dt)*hf[i] - 0.5*dt*hEdvf[i];
      }
      hf2.ifft(&(f2[k][0]));
    }

    E = poisson_solver(f2.density());
    Edvf = o2::trp_v(f2,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      //hf1.fft(&(f1[k][0]));
      //hf2.fft(&(f2[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf[i] = std::exp(L*dt)*hf[i] - dt*std::exp(0.5*L*dt)*hEdvf[i];
      }
      hf.ifft(&(f[k][0]));
    }
    **/
    /**/
    // RK(3,3) ================================================================
    // RK(3,3) eq19 ===========================================================
    #define SCHEME "RK33"
    //#define SCHEME "RK33_eq19"
    // SIGMA = std::sqrt(3) (y_max)
    // SIGMA = 1.433 (WENO)
    E = poisson_solver(f.density());
    
    //std::cout << "\n" << L << " " << *std::max_element(E.begin(),E.end()) << " " << *std::min_element(E.begin(),E.end()) << std::endl;

    field<double,1> Edvf = weno::trp_v(f,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf1[i] = std::exp(L*dt)*( hf[i]-dt*hEdvf[i] );
        //hf1[i] = 0.5*std::exp((2./3.)*L*dt)*hf[i] + 0.5*std::exp((2./3.)*dt*L)*( hf[i] - (4./3.)*dt*hEdvf[i] );
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
        hf2[i] = 0.75*std::exp(0.5*L*dt)*hf[i] + 0.25*std::exp(-0.5*L*dt)*( hf1[i]-dt*hEdvf[i] );
        //hf2[i] = (2./3.)*std::exp((2./3.)*dt*L)*hf[i] + (1./3.)*( hf1[i] - (4./3.)*dt*hEdvf[i] );
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
        hf[i] = (1./3.)*std::exp(L*dt)*hf[i] + (2./3.)*std::exp(0.5*L*dt)*( hf2[i]-dt*hEdvf[i] );
        //hf[i] = (59./128.)*std::exp(L*dt)*hf[i] + (15./128.)*std::exp(L*dt)*( 2.*hf1[i]*std::exp(-(2./3.)*L*dt) - hf[i] ) + (27./64.)*std::exp((1./3.)*dt*L)*( hf2[i] - (4./3.)*dt*hEdvf[i] );
      }
      hf.ifft(&(f[k][0]));
    }
    /**/
    /**
    // RK(4,4) ================================================================
    // RK(4,4) 3/8 rule =======================================================
    #define SCHEME "RK44"
    //#define SCHEME "RK44_38"
    // SIGMA = 2.*std::sqrt(2.) (y_max)
    // SIGMA = 1.731 (WENO)
    E = poisson_solver( f.density() );
    field<double,1> Edvf = weno::trp_v(f,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf1[i] = std::exp(0.5*L*dt)*( hf[i] - 0.5*dt*hEdvf[i] );
        //hf1[i] = std::exp((1./3.)*L*dt)*( hf[i] - (1./3.)*dt*hEdvf[i] );
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
        hf2[i] = std::exp(0.5*L*dt)*hf[i] - 0.5*dt*hEdvf[i] ;
        //hf2[i] = 2.*std::exp((2./3.)*L*dt)*hf[i] - std::exp((1./3.)*L*dt)*hf1[i] - dt*std::exp((1./3.)*L*dt)*hEdvf[i];
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
        hf3[i] = std::exp(L*dt)*hf[i] - dt*std::exp(0.5*L*dt)*hEdvf[i];
        //hf3[i] = 2.*std::exp((2./3.)*L*dt)*hf1[i] - std::exp((1./3.)*L*dt)*hf2[i] - dt*std::exp((1./3.)*L*dt)*hEdvf[i];
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
        hf[i] = -(1./3.)*std::exp(L*dt)*hf[i] + (1./3.)*std::exp(0.5*L*dt)*hf1[i] + (2./3.)*std::exp(0.5*L*dt)*hf2[i] + (1./3.)*hf3[i] - (1./6.)*dt*hEdvf[i];
        //hf[i] = -(1./8.)*std::exp(L*dt)*hf[i] + 0.75*std::exp((1./3.)*L*dt)*hf2[i] + (3./8.)*hf3[i] - (1./8.)*dt*hEdvf[i];
      }

      hf.ifft(&(f[k][0]));
    }
    **/
    /**
    // RK(5,3) ================================================================
    #define SCHEME "RK53"
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
    // DP5  ===================================================================
    #define SCHEME "DP5"
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
    /**
    // RK(8,6) ================================================================
    #define SCHEME"RK86"
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
    **/

    // end of time loop
    ++i_t;
    rho = f.density();
    E = poisson_solver(rho);
    //ee(i_t) = 0.;
    //for ( auto i=0 ; i<Nx ; ++i ) {
    //  ee(i_t) += SQ(E(i))*f.step.dx;
    //}
    //ee(i_t) = std::sqrt(ee(i_t));
    H(i_t) = energy(f,E);


//#define FOLDER "lukas/vp/"
#define FOLDER "vphl/kin/"
#define SPACE_SCHEME "weno"
/*
    if ( i_t == int(15./dt) ) {
      std::stringstream ss; ss << FOLDER << "vp_" << SCHEME << "_" << SPACE_SCHEME << "_15.dat";
      f.write(ss.str());
      std::cout << std::endl;
    }
    if ( i_t == int(20./dt) ) {
      std::stringstream ss; ss << FOLDER << "vp_" << SCHEME << "_" << SPACE_SCHEME << "_20.dat";
      f.write(ss.str());
      std::cout << std::endl;
    }
    if ( i_t == int(25./dt) ) {
      std::stringstream ss; ss << FOLDER << "vp_" << SCHEME << "_" << SPACE_SCHEME << "_25.dat";
      f.write(ss.str());
      std::cout << std::endl;
    }
    if ( i_t == int(30./dt) ) {
      std::stringstream ss; ss << FOLDER << "vp_" << SCHEME << "_" << SPACE_SCHEME << "_30.dat";
      f.write(ss.str());
      std::cout << std::endl;
    }
    if ( i_t == int(35./dt) ) {
      std::stringstream ss; ss << FOLDER << "vp_" << SCHEME << "_" << SPACE_SCHEME << "_35.dat";
      f.write(ss.str());
      std::cout << std::endl;
    }
*/
	} // while (  i_t*dt < Tf )
#undef L
  std::cout<<" ["<<std::setw(5)<<i_t<<"] "<<i_t*dt<<"   \r"<<std::endl;

  std::stringstream ss; ss << FOLDER << "vp_" << SCHEME << "_" << SPACE_SCHEME << "_dt" << dt << ".dat";
  f.write(ss.str());

  rho = f.density();
  auto dx_y = [&,count=0](auto const& y) mutable { std::stringstream ss; ss<<Xi(count++)<<" "<<y; return ss.str(); };
  ss.str(std::string());
  ss << FOLDER << "rho_" << SCHEME << "_" << SPACE_SCHEME << "_" << Tf << ".dat";
  std::ofstream of;
  of.open(ss.str()); ss.str(std::string());
  std::transform( rho.begin() , rho.end() , std::ostream_iterator<std::string>(of,"\n") , dx_y );
  of.close();

  ss.str(std::string());
  //std::ofstream of;
  std::size_t count = 0;
  auto dt_y = [&,count=0](auto const& y) mutable { std::stringstream ss; ss<<(count++)*dt<<" "<<y; return ss.str(); };
  ss << FOLDER << "ee_" << SCHEME << "_" << SPACE_SCHEME << "_10.dat";
  of.open(ss.str()); ss.str(std::string());
  for ( auto i=0; i<ee.size() ; ++i ) {
    of << i*dt <<" " << ee[i] << "\n";
  }
  of.close();
  ss << FOLDER << "H_" << SCHEME << "_" << SPACE_SCHEME << ".dat";
  of.open(ss.str()); ss.str(std::string());
  for ( auto i=0; i<H.size() ; ++i ) {
    of << i*dt <<" " << (H[i]-H[0])/std::abs(H[0]) << "\n";
  }

  //double h = std::abs(*std::max_element( H.begin() , H.end() , [&](double a,double b){return ( std::abs((a-H[0])/std::abs(H[0])) < std::abs((b-H[0])/std::abs(H[0])) );} ));
  //std::cout << dt << " " << std::abs((h-H[0])/std::abs(H[0])) << "\n";

  of.close();
  ss << FOLDER << "Emax_" << SCHEME << "_" << SPACE_SCHEME << "_10.dat";
  of.open(ss.str()); ss.str(std::string());
  std::transform( Emax.begin() , Emax.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();

	return 0;
}
