#ifndef _RK_H_
#define _RK_H_

#include <array>
#include <fstream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array.hpp>

#include <boost/iterator/zip_iterator.hpp>

#include "array_view.h"
#include "field.h"
#include "fft.h"

using namespace boost::numeric;

namespace lawson {

// RK(3,3) ====================================================================
  template < typename Poisson_Solver >
  struct RK32 {
    std::size_t Nx;
    std::size_t Nv;
    Poisson_Solver poisson_solver;
    ublas::vector<double> E;
    ublas::vector<double> vk;
    ublas::vector<double> kx;
    field<double,1> f1;
    field<double,1> f2;
    field<double,1> Edvf;
    fft::spectrum_ hfn;
    fft::spectrum_ hf1;
    fft::spectrum_ hf2;
    fft::spectrum_ hEdvf;
    std::string label;

    std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L;

    RK32 ( std::size_t N_x , std::size_t N_v , double l , const std::size_t * shape ,
      ublas::vector<double> const& v_k , ublas::vector<double> const& k_x ,
      std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L_trp )
    :
      Nx(N_x) , Nv(N_v) , // size
      poisson_solver(N_x,l) , // Poisson solver (could be identity for rotation)
      E(N_x,0) , vk(v_k) , kx(k_x) ,
      f1(tools::array_view<const std::size_t>(shape,2)) ,
      f2(tools::array_view<const std::size_t>(shape,2)) ,
      Edvf(tools::array_view<const std::size_t>(shape,2)) ,
      hfn(N_x) , hf1(N_x) , hf2(N_x) , hEdvf(N_x) ,
      L(L_trp)
    { label = "RK32"; }

    field<double,1>
    operator () ( field<double,1> & fn , double dt )
    {
      static const std::complex<double> & I = std::complex<double>(0.,1.);
      
      // first stage
      E = poisson_solver( fn.density() );
      Edvf = L(fn,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i< Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf1[i] = std::exp(0.5*L*dt)*( hfn[i]-0.5*dt*hEdvf[i] );
        }

        hf1.ifft(&(f1[k][0]));
      }

      // second stage
      E = poisson_solver( f1.density() );
      Edvf = L(f1,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        //hf1.fft(&(f1[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf2[i] = std::exp(0.5*L*dt)*hfn[i] - 0.5*dt*hEdvf[i];
        }

        hf2.ifft(&(f2[k][0]));
      }

      // thrid stage
      E = poisson_solver( f2.density() );
      Edvf = L(f2,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        //hf1.fft(&(f1[k][0]));
        //hf2.fft(&(f2[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hfn[i] = std::exp(L*dt)*hfn[i] - dt*std::exp(0.5*L*dt)*hEdvf[i];
        }

        hfn.ifft(&(fn[k][0]));
      }

      return fn;
    }
  };

// RK(3,3) ====================================================================
  template < typename Poisson_Solver >
  struct RK33 {
    std::size_t Nx;
    std::size_t Nv;
    Poisson_Solver poisson_solver;
    ublas::vector<double> E;
    ublas::vector<double> vk;
    ublas::vector<double> kx;
    field<double,1> f1;
    field<double,1> f2;
    field<double,1> Edvf;
    fft::spectrum_ hfn;
    fft::spectrum_ hf1;
    fft::spectrum_ hf2;
    fft::spectrum_ hEdvf;
    std::string label;

    std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L;

    RK33 ( std::size_t N_x , std::size_t N_v , double l , const std::size_t * shape ,
      ublas::vector<double> const& v_k , ublas::vector<double> const& k_x ,
      std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L_trp )
    :
      Nx(N_x) , Nv(N_v) , // size
      poisson_solver(N_x,l) , // Poisson solver (could be identity for rotation)
      E(N_x,0) , vk(v_k) , kx(k_x) ,
      f1(tools::array_view<const std::size_t>(shape,2)) ,
      f2(tools::array_view<const std::size_t>(shape,2)) ,
      Edvf(tools::array_view<const std::size_t>(shape,2)) ,
      hfn(N_x) , hf1(N_x) , hf2(N_x) , hEdvf(N_x) ,
      L(L_trp)
    { label = "RK33"; }

    field<double,1>
    operator () ( field<double,1> & fn , double dt )
    {
      static const std::complex<double> & I = std::complex<double>(0.,1.);
      
      // first stage
      E = poisson_solver( fn.density() );
      Edvf = L(fn,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i< Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf1[i] = std::exp(L*dt)*( hfn[i]-dt*hEdvf[i] );
        }

        hf1.ifft(&(f1[k][0]));
      }

      // second stage
      E = poisson_solver( f1.density() );
      Edvf = L(f1,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf2[i] = 0.75*std::exp(0.5*L*dt)*hfn[i] + 0.25*std::exp(-0.5*L*dt)*( hf1[i]-dt*hEdvf[i] );
        }

        hf2.ifft(&(f2[k][0]));
      }

      // thrid stage
      E = poisson_solver( f2.density() );
      Edvf = L(f2,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hfn[i] = (1./3.)*std::exp(L*dt)*hfn[i] + (2./3.)*std::exp(0.5*L*dt)*( hf2[i]-dt*hEdvf[i] );
        }

        hfn.ifft(&(fn[k][0]));
      }

      return fn;
    }
  };

// RK(3,3) eq19 ===============================================================
  template < typename Poisson_Solver >
  struct RK33_eq19 {
    std::size_t Nx;
    std::size_t Nv;
    Poisson_Solver poisson_solver;
    ublas::vector<double> E;
    ublas::vector<double> vk;
    ublas::vector<double> kx;
    field<double,1> f1;
    field<double,1> f2;
    field<double,1> Edvf;
    fft::spectrum_ hfn;
    fft::spectrum_ hf1;
    fft::spectrum_ hf2;
    fft::spectrum_ hEdvf;
    std::string label;

    std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L;

    RK33_eq19 ( std::size_t N_x , std::size_t N_v , double l , const std::size_t * shape ,
      ublas::vector<double> const& v_k , ublas::vector<double> const& k_x ,
      std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L_trp )
    :
      Nx(N_x) , Nv(N_v) , // size
      poisson_solver(N_x,l) , // Poisson solver (could be identity for rotation)
      E(N_x,0) , vk(v_k) , kx(k_x) ,
      f1(tools::array_view<const std::size_t>(shape,2)) ,
      f2(tools::array_view<const std::size_t>(shape,2)) ,
      Edvf(tools::array_view<const std::size_t>(shape,2)) ,
      hfn(N_x) , hf1(N_x) , hf2(N_x) , hEdvf(N_x) ,
      L(L_trp)
    { label = "RK33_eq19"; }

    field<double,1>
    operator () ( field<double,1> & fn , double dt )
    {
      static const std::complex<double> & I = std::complex<double>(0.,1.);
      
      // first stage
      E = poisson_solver( fn.density() );
      Edvf = L(fn,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i< Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf1[i] = 0.5*std::exp((2./3.)*L*dt)*hfn[i] + 0.5*std::exp((2./3.)*dt*L)*( hfn[i] - (4./3.)*dt*hEdvf[i] );
        }

        hf1.ifft(&(f1[k][0]));
      }

      // second stage
      E = poisson_solver( f1.density() );
      Edvf = L(f1,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf2[i] = (2./3.)*std::exp((2./3.)*dt*L)*hfn[i] + (1./3.)*( hf1[i] - (4./3.)*dt*hEdvf[i] );
        }

        hf2.ifft(&(f2[k][0]));
      }

      // thrid stage
      E = poisson_solver( f2.density() );
      Edvf = L(f2,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hfn[i] = (59./128.)*std::exp(L*dt)*hfn[i] + (15./128.)*std::exp(L*dt)*( 2.*hf1[i]*std::exp(-(2./3.)*L*dt) - hfn[i] ) + (27./64.)*std::exp((1./3.)*dt*L)*( hf2[i] - (4./3.)*dt*hEdvf[i] );
        }

        hfn.ifft(&(fn[k][0]));
      }

      return fn;
    }
  };

// RK(4,4) ====================================================================
  template < typename Poisson_Solver >
  struct RK44 {
    std::size_t Nx;
    std::size_t Nv;
    Poisson_Solver poisson_solver;
    ublas::vector<double> E;
    ublas::vector<double> vk;
    ublas::vector<double> kx;
    field<double,1> f1;
    field<double,1> f2;
    field<double,1> f3;
    field<double,1> Edvf;
    fft::spectrum_ hfn;
    fft::spectrum_ hf1;
    fft::spectrum_ hf2;
    fft::spectrum_ hf3;
    fft::spectrum_ hEdvf;
    std::string label;

    std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L;

    RK44 ( std::size_t N_x , std::size_t N_v , double l , const std::size_t * shape ,
      ublas::vector<double> const& v_k , ublas::vector<double> const& k_x ,
      std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L_trp )
    :
      Nx(N_x) , Nv(N_v) , // size
      poisson_solver(N_x,l) , // Poisson solver (could be identity for rotation)
      E(N_x,0) , vk(v_k) , kx(k_x) ,
      f1(tools::array_view<const std::size_t>(shape,2)) ,
      f2(tools::array_view<const std::size_t>(shape,2)) ,
      f3(tools::array_view<const std::size_t>(shape,2)) ,
      Edvf(tools::array_view<const std::size_t>(shape,2)) ,
      hfn(N_x) , hf1(N_x) , hf2(N_x) , hf3(N_x) , hEdvf(N_x) ,
      L(L_trp)
    { label = "RK44"; }

    field<double,1>
    operator () ( field<double,1> & fn , double dt )
    {
      static const std::complex<double> & I = std::complex<double>(0.,1.);
      
      // first stage
      E = poisson_solver( fn.density() );
      Edvf = L(fn,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i< Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf1[i] = std::exp(0.5*L*dt)*( hfn[i] - 0.5*dt*hEdvf[i] );
        }

        hf1.ifft(&(f1[k][0]));
      }

      // second stage
      E = poisson_solver( f1.density() );
      Edvf = L(f1,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        //hf1.fft(&(f1[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf2[i] = std::exp(0.5*L*dt)*hfn[i] - 0.5*dt*hEdvf[i];
        }

        hf2.ifft(&(f2[k][0]));
      }

      // thrid stage
      E = poisson_solver( f2.density() );
      Edvf = L(f2,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        //hf1.fft(&(f1[k][0]));
        //hf2.fft(&(f2[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf3[i] = std::exp(L*dt)*hfn[i] - dt*std::exp(0.5*L*dt)*hEdvf[i];
        }

        hf3.ifft(&(f3[k][0]));
      }

      // fourth stage
      E = poisson_solver( f3.density() );
      Edvf = L(f3,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hf3.fft(&(f3[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hfn[i] = -(1./3.)*std::exp(L*dt)*hfn[i] + (1./3.)*std::exp(0.5*L*dt)*hf1[i] + (2./3.)*std::exp(0.5*L*dt)*hf2[i] + (1./3.)*hf3[i] - (1./6.)*dt*hEdvf[i];
        }

        hfn.ifft(&(fn[k][0]));
      }

      return fn;
    }
  };

// RK(4,4) 3/8 rule ===========================================================
  template < typename Poisson_Solver >
  struct RK44_38 {
    std::size_t Nx;
    std::size_t Nv;
    Poisson_Solver poisson_solver;
    ublas::vector<double> E;
    ublas::vector<double> vk;
    ublas::vector<double> kx;
    field<double,1> f1;
    field<double,1> f2;
    field<double,1> f3;
    field<double,1> Edvf;
    fft::spectrum_ hfn;
    fft::spectrum_ hf1;
    fft::spectrum_ hf2;
    fft::spectrum_ hf3;
    fft::spectrum_ hEdvf;
    std::string label;

    std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L;

    RK44_38 ( std::size_t N_x , std::size_t N_v , double l , const std::size_t * shape ,
      ublas::vector<double> const& v_k , ublas::vector<double> const& k_x ,
      std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L_trp )
    :
      Nx(N_x) , Nv(N_v) , // size
      poisson_solver(N_x,l) , // Poisson solver (could be identity for rotation)
      E(N_x,0) , vk(v_k) , kx(k_x) ,
      f1(tools::array_view<const std::size_t>(shape,2)) ,
      f2(tools::array_view<const std::size_t>(shape,2)) ,
      f3(tools::array_view<const std::size_t>(shape,2)) ,
      Edvf(tools::array_view<const std::size_t>(shape,2)) ,
      hfn(N_x) , hf1(N_x) , hf2(N_x) , hf3(N_x) , hEdvf(N_x) ,
      L(L_trp)
    { label = "RK44_38"; }

    field<double,1>
    operator () ( field<double,1> & fn , double dt )
    {
      static const std::complex<double> & I = std::complex<double>(0.,1.);
      
      // first stage
      E = poisson_solver( fn.density() );
      Edvf = L(fn,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i< Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf1[i] = std::exp((1./3.)*L*dt)*( hfn[i] - (1./3.)*dt*hEdvf[i] );
        }

        hf1.ifft(&(f1[k][0]));
      }

      // second stage
      E = poisson_solver( f1.density() );
      Edvf = L(f1,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf2[i] = 2.*std::exp((2./3.)*L*dt)*hfn[i] - std::exp((1./3.)*L*dt)*hf1[i] - dt*std::exp((1./3.)*L*dt)*hEdvf[i];
        }

        hf2.ifft(&(f2[k][0]));
      }

      // thrid stage
      E = poisson_solver( f2.density() );
      Edvf = L(f2,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        //hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf3[i] = 2.*std::exp((2./3.)*L*dt)*hf1[i] - std::exp((1./3.)*L*dt)*hf2[i] - dt*std::exp((1./3.)*L*dt)*hEdvf[i];
        }

        hf3.ifft(&(f3[k][0]));
      }

      // fourth stage
      E = poisson_solver( f3.density() );
      Edvf = L(f3,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        //hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hf3.fft(&(f3[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hfn[i] = -(1./8.)*std::exp(L*dt)*hfn[i] + 0.75*std::exp((1./3.)*L*dt)*hf2[i] + (3./8.)*hf3[i] - (1./8.)*dt*hEdvf[i];
        }

        hfn.ifft(&(fn[k][0]));
      }

      return fn;
    }
  };

// RK(5,3) ====================================================================
  template < typename Poisson_Solver >
  struct RK53 {
    std::size_t Nx;
    std::size_t Nv;
    Poisson_Solver poisson_solver;
    ublas::vector<double> E;
    ublas::vector<double> vk;
    ublas::vector<double> kx;
    field<double,1> f1;
    field<double,1> f2;
    field<double,1> f3;
    field<double,1> f4;
    field<double,1> Edvf;
    fft::spectrum_ hfn;
    fft::spectrum_ hf1;
    fft::spectrum_ hf2;
    fft::spectrum_ hf3;
    fft::spectrum_ hf4;
    fft::spectrum_ hEdvf;
    std::string label;

    std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L;

    RK53 ( std::size_t N_x , std::size_t N_v , double l , const std::size_t * shape ,
      ublas::vector<double> const& v_k , ublas::vector<double> const& k_x ,
      std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L_trp )
    :
      Nx(N_x) , Nv(N_v) , // size
      poisson_solver(N_x,l) , // Poisson solver (could be identity for rotation)
      E(N_x,0) , vk(v_k) , kx(k_x) ,
      f1(tools::array_view<const std::size_t>(shape,2)) ,
      f2(tools::array_view<const std::size_t>(shape,2)) ,
      f3(tools::array_view<const std::size_t>(shape,2)) ,
      f4(tools::array_view<const std::size_t>(shape,2)) ,
      Edvf(tools::array_view<const std::size_t>(shape,2)) ,
      hfn(N_x) , hf1(N_x) , hf2(N_x) , hf3(N_x) , hf4(N_x) , hEdvf(N_x) ,
      L(L_trp)
    { label = "RK53"; }

    field<double,1>
    operator () ( field<double,1> & fn , double dt )
    {
      static const std::complex<double> & I = std::complex<double>(0.,1.);
      
      // first stage
      E = poisson_solver( fn.density() );
      Edvf = L(fn,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i< Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf1[i] = std::exp((1./7.)*L*dt)*hfn[i] - (1./7.)*dt*std::exp((1./7.)*L*dt)*hEdvf[i];
        }

        hf1.ifft(&(f1[k][0]));
      }

      // second stage
      E = poisson_solver( f1.density() );
      Edvf = L(f1,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf2[i] = std::exp((3./16.)*L*dt)*hfn[i] - (3./16.)*dt*std::exp((5./112.)*L*dt)*hEdvf[i];
        }

        hf2.ifft(&(f2[k][0]));
      }

      // thrid stage
      E = poisson_solver( f2.density() );
      Edvf = L(f2,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf3[i] = std::exp((1./3.)*L*dt)*hfn[i] - (1./3.)*dt*std::exp((7./48.)*L*dt)*hEdvf[i];
        }

        hf3.ifft(&(f3[k][0]));
      }

      // fourth stage
      E = poisson_solver( f3.density() );
      Edvf = L(f3,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf4[i] = std::exp((2./3.)*L*dt)*hfn[i] - (2./3.)*dt*std::exp((1./3.)*L*dt)*hEdvf[i];
        }

        hf4.ifft(&(f4[k][0]));
      }

      // fifth stage
      E = poisson_solver( f4.density() );
      Edvf = L(f4,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hfn[i] = -0.75*std::exp(L*dt)*hfn[i] + 1.75*std::exp((6./7.)*L*dt)*hf1[i] - 0.75*dt*std::exp((1./3.)*L*dt)*hEdvf[i];
        }

        hfn.ifft(&(fn[k][0]));
      }
      return fn;
    }
  };

// DP5 ========================================================================
  template < typename Poisson_Solver >
  struct DP5 {
    std::size_t Nx;
    std::size_t Nv;
    Poisson_Solver poisson_solver;
    ublas::vector<double> E;
    ublas::vector<double> vk;
    ublas::vector<double> kx;
    field<double,1> f1;
    field<double,1> f2;
    field<double,1> f3;
    field<double,1> f4;
    field<double,1> f5;
    field<double,1> Edvf;
    fft::spectrum_ hfn;
    fft::spectrum_ hf1;
    fft::spectrum_ hf2;
    fft::spectrum_ hf3;
    fft::spectrum_ hf4;
    fft::spectrum_ hf5;
    fft::spectrum_ hEdvf;
    std::string label;

    std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L;

    DP5 ( std::size_t N_x , std::size_t N_v , double l , const std::size_t * shape ,
      ublas::vector<double> const& v_k , ublas::vector<double> const& k_x ,
      std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L_trp )
    :
      Nx(N_x) , Nv(N_v) , // size
      poisson_solver(N_x,l) , // Poisson solver (could be identity for rotation)
      E(N_x,0) , vk(v_k) , kx(k_x) ,
      f1(tools::array_view<const std::size_t>(shape,2)) ,
      f2(tools::array_view<const std::size_t>(shape,2)) ,
      f3(tools::array_view<const std::size_t>(shape,2)) ,
      f4(tools::array_view<const std::size_t>(shape,2)) ,
      f5(tools::array_view<const std::size_t>(shape,2)) ,
      Edvf(tools::array_view<const std::size_t>(shape,2)) ,
      hfn(N_x) , hf1(N_x) , hf2(N_x) , hf3(N_x) , hf4(N_x) , hf5(N_x) , hEdvf(N_x) ,
      L(L_trp)
    { label = "DP5"; }

    field<double,1>
    operator () ( field<double,1> & fn , double dt )
    {
      static const std::complex<double> & I = std::complex<double>(0.,1.);
      
      // first stage
      E = poisson_solver( fn.density() );
      Edvf = L(fn,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i< Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf1[i] = std::exp((1./5.)*L*dt)*hfn[i] - (1./5.)*dt*std::exp((1./5.)*L*dt)*hEdvf[i];
        }

        hf1.ifft(&(f1[k][0]));
      }

      // second stage
      E = poisson_solver( f1.density() );
      Edvf = L(f1,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf2[i] = (5./8.)*std::exp((3./10.)*L*dt)*hfn[i] + (3./8.)*std::exp((1./10.)*L*dt)*hf1[i] - (9./40.)*dt*std::exp((1./10.)*L*dt)*hEdvf[i];
        }

        hf2.ifft(&(f2[k][0]));
      }

      // thrid stage
      E = poisson_solver( f2.density() );
      Edvf = L(f2,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf3[i] = (175./27.)*std::exp((4./5.)*L*dt)*hfn[i] + (100./9.)*std::exp((3./5.)*L*dt)*hf1[i] - (448./27.)*std::exp(0.5*L*dt)*hf2[i] - (32./9.)*dt*std::exp(0.5*L*dt)*hEdvf[i];
        }

        hf3.ifft(&(f3[k][0]));
      }

      // fourth stage
      E = poisson_solver( f3.density() );
      Edvf = L(f3,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hf3.fft(&(f3[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf4[i] = (3551./6561.)*std::exp((8./9.)*L*dt)*hfn[i] + (7420./2187.)*std::exp((31./45.)*L*dt)*hf1[i] - (37376./6561.)*std::exp((53./90.)*L*dt)*hf2[i] + (2014./729.)*std::exp((4./45.)*L*dt)*hf3[i] + (212./729.)*dt*std::exp((4./45.)*L*dt)*hEdvf[i];
        }

        hf4.ifft(&(f4[k][0]));
      }

      // fifth stage
      E = poisson_solver( f3.density() );
      Edvf = L(f3,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hf3.fft(&(f3[k][0]));
        hf4.fft(&(f4[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf5[i] = (313397./335808.)*std::exp(L*dt)*hfn[i] + (424025./55968.)*std::exp((4./5.)*L*dt)*hf1[i] - (61400./5247.)*std::exp((7./10.)*L*dt)*hf2[i] + (96075./18656.)*std::exp((1./5.)*L*dt)*hf3[i] - (35721./37312.)*std::exp((1./9.)*L*dt)*hf4[i] + (5103./18656.)*dt*std::exp((1./9.)*L*dt)*hEdvf[i];
        }

        hf5.ifft(&(f5[k][0]));
      }

      // sixth stage
      E = poisson_solver( f3.density() );
      Edvf = L(f5,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hf3.fft(&(f3[k][0]));
        hf4.fft(&(f4[k][0]));
        hf5.fft(&(f5[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hfn[i] = -(563./3456.)*std::exp(L*dt)*hfn[i] - (575./252.)*std::exp((4./5.)*L*dt)*hf1[i] + (31400./10017.)*std::exp((7./10.)*L*dt)*hf2[i] + (325./1344.)*std::exp((1./5.)*L*dt)*hf3[i] - (7533./6784.)*std::exp((1./9.)*L*dt)*hf4[i] + (33./28.)*hf5[i] - (11./84.)*dt*hEdvf[i];
        }

        hfn.ifft(&(fn[k][0]));
      }

      return fn;
    }
  };

// RK(8,6) ====================================================================
  template < typename Poisson_Solver >
  struct RK86 {
    std::size_t Nx;
    std::size_t Nv;
    Poisson_Solver poisson_solver;
    ublas::vector<double> E;
    ublas::vector<double> vk;
    ublas::vector<double> kx;
    field<double,1> f1;
    field<double,1> f2;
    field<double,1> f3;
    field<double,1> f4;
    field<double,1> f5;
    field<double,1> f6;
    field<double,1> f7;
    field<double,1> Edvf;
    fft::spectrum_ hfn;
    fft::spectrum_ hf1;
    fft::spectrum_ hf2;
    fft::spectrum_ hf3;
    fft::spectrum_ hf4;
    fft::spectrum_ hf5;
    fft::spectrum_ hf6;
    fft::spectrum_ hf7;
    fft::spectrum_ hEdvf;
    std::string label;

    std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L;

    RK86 ( std::size_t N_x , std::size_t N_v , double l , const std::size_t * shape ,
      ublas::vector<double> const& v_k , ublas::vector<double> const& k_x ,
      std::function<field<double,1>(field<double,1>const&,ublas::vector<double>const&)> L_trp )
    :
      Nx(N_x) , Nv(N_v) , // size
      poisson_solver(N_x,l) , // Poisson solver (could be identity for rotation)
      E(N_x,0) , vk(v_k) , kx(k_x) ,
      f1(tools::array_view<const std::size_t>(shape,2)) ,
      f2(tools::array_view<const std::size_t>(shape,2)) ,
      f3(tools::array_view<const std::size_t>(shape,2)) ,
      f4(tools::array_view<const std::size_t>(shape,2)) ,
      f5(tools::array_view<const std::size_t>(shape,2)) ,
      f6(tools::array_view<const std::size_t>(shape,2)) ,
      f7(tools::array_view<const std::size_t>(shape,2)) ,
      Edvf(tools::array_view<const std::size_t>(shape,2)) ,
      hfn(N_x) , hf1(N_x) , hf2(N_x) , hf3(N_x) , hf4(N_x) , hf5(N_x) , hf6(N_x) , hf7(N_x) , hEdvf(N_x) ,
      L(L_trp)
    { label = "RK86"; }

    field<double,1>
    operator () ( field<double,1> & fn , double dt )
    {
      static const std::complex<double> & I = std::complex<double>(0.,1.);
      
      // first stage
      E = poisson_solver( fn.density() );
      Edvf = L(fn,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i< Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf1[i] = std::exp((1./9.)*L*dt)*hfn[i] - (1./9.)*dt*std::exp((1./9.)*dt*L)*hEdvf[i];
        }

        hf1.ifft(&(f1[k][0]));
      }

      // second stage
      E = poisson_solver( f1.density() );
      Edvf = L(f1,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf2[i] = (5./8.)*std::exp((1./6.)*L*dt)*hfn[i] + (3./8.)*std::exp((1./18.)*L*dt)*hf1[i] - (1./8.)*dt*std::exp((1./18.)*L*dt)*hEdvf[i];
        }

        hf2.ifft(&(f2[k][0]));
      }

      // thrid stage
      E = poisson_solver( f2.density() );
      Edvf = L(f2,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf3[i] = 2.*std::exp((1./3.)*L*dt)*hfn[i] + 3.*std::exp((2./9.)*L*dt)*hf1[i] -4.*std::exp((1./6.)*L*dt)*hf2[i] - (2./3.)*dt*std::exp((1./6.)*L*dt)*hEdvf[i];
        }

        hf3.ifft(&(f3[k][0]));
      }

      // fourth stage
      E = poisson_solver( f3.density() );
      Edvf = L(f3,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hf3.fft(&(f3[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf4[i] = (305./1268.)*std::exp(0.5*L*dt)*hfn[i] + (2817./1268.)*std::exp((7./18.)*L*dt)*hf1[i] - (927./317.)*std::exp((1./3.)*L*dt)*hf2[i] + (927./634.)*std::exp((1./6.)*L*dt)*hf3[i] - (321./1268.)*dt*std::exp((1./6.)*L*dt)*hEdvf[i];
        }

        hf4.ifft(&(f4[k][0]));
      }

      // fifth stage
      E = poisson_solver( f3.density() );
      Edvf = L(f3,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hf3.fft(&(f3[k][0]));
        hf4.fft(&(f4[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf5[i] = (3191./321.)*std::exp((2./3.)*L*dt)*hfn[i] - (2436./107.)*std::exp((5./9.)*L*dt)*hf1[i] - (2404./107.)*std::exp(0.5*L*dt)*hf2[i] + (12330./107.)*std::exp((1./3.)*L*dt)*hf3[i] - (25340./321.)*std::exp((1./6.)*L*dt)*hf4[i] - 8.*dt*std::exp((1./6.)*L*dt)*hEdvf[i];
        }

        hf5.ifft(&(f5[k][0]));
      }

      // sixth stage
      E = poisson_solver( f3.density() );
      Edvf = L(f5,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hf3.fft(&(f3[k][0]));
        hf4.fft(&(f4[k][0]));
        hf5.fft(&(f5[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf6[i] = -(15130159./6286464.)*std::exp((5./6.)*L*dt)*hfn[i] + (2014319./349248.)*std::exp((13./18.)*L*dt)*hf1[i] + (1194095./523872.)*std::exp((2./3.)*L*dt)*hf2[i] - (1471057./116416.)*std::exp(0.5*L*dt)*hf3[i] + (12601453./1571616.)*std::exp((1./3.)*L*dt)*hf4[i] - (433./19584.)*std::exp((1./6.)*L*dt)*hf5[i] - (33./1088.)*dt*std::exp((1./6.)*L*dt)*hEdvf[i];
        }

        hf6.ifft(&(f6[k][0]));
      }

      // seventh stage
      E = poisson_solver( f3.density() );
      Edvf = L(f5,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hf3.fft(&(f3[k][0]));
        hf4.fft(&(f4[k][0]));
        hf5.fft(&(f5[k][0]));
        hf6.fft(&(f6[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hf7[i] = (805187./78966.)*std::exp(L*dt)*hfn[i] - (2263766./48257.)*std::exp((8./9.)*L*dt)*hf1[i] + (2745422./144771.)*std::exp((5./6.)*L*dt)*hf2[i] + (2271108./48257.)*std::exp((2./3.)*L*dt)*hf3[i] - (13115270./434313.)*std::exp(0.5*L*dt)*hf4[i] - (227./2706.)*std::exp((1./3.)*L*dt)*hf5[i] + (888./451.)*std::exp((1./6.)*L*dt)*hf6[i] - (36./41.)*dt*std::exp((1./6.)*L*dt)*hEdvf[i];
        }

        hf7.ifft(&(f7[k][0]));
      }

      // eigth stage
      E = poisson_solver( f3.density() );
      Edvf = L(f5,E);
      for ( auto k=0 ; k<Nv ; ++k ) {
        hfn.fft(&(fn[k][0]));
        hf1.fft(&(f1[k][0]));
        hf2.fft(&(f2[k][0]));
        hf3.fft(&(f3[k][0]));
        hf4.fft(&(f4[k][0]));
        hf5.fft(&(f5[k][0]));
        hf6.fft(&(f6[k][0]));
        hf7.fft(&(f7[k][0]));
        hEdvf.fft(&(Edvf[k][0]));

        for ( auto i=0 ; i<Nx ; ++i ) {
          std::complex<double> L = -vk[k]*I*kx[i];
          hfn[i] = -(193999./179760.)*std::exp(L*dt)*hfn[i] + (2487363./329560.)*std::exp((8./9.)*L*dt)*hf1[i] - (847909./164780.)*std::exp((5./6.)*L*dt)*hf2[i] - (1600251./329560.)*std::exp((2./3.)*L*dt)*hf3[i] + (362713./98868.)*std::exp(0.5*L*dt)*hf4[i] + (109./1232.)*std::exp((1./3.)*L*dt)*hf5[i] + (186./385.)*std::exp((1./6.)*L*dt)*hf6[i] + (41./140.)*hf7[i] - (41./840.)*dt*hEdvf[i];
        }

        hfn.ifft(&(fn[k][0]));
      }

      return fn;
    }
  };

} // namespace lawson

#endif
