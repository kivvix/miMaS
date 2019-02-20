#ifndef MISC_TRP1D_BWENO_H_
#define MISC_TRP1D_BWENO_H_

namespace bweno
{

template <typename _T>
_T
beta ( _T ujm2 , _T ujm1 , _T uj , _T ujp1 , _T ujp2 )
{
  _T D14uj = 1./12.*( -ujp2 +  8.*ujp1          -  8.*ujm1 + ujm2 );
	_T D24uj = 1./12.*( -ujp2 + 16.*ujp1 - 30.*uj + 16.*ujm1 - ujm2 );
	_T D32uj =  1./2.*(  ujp2 -  2.*ujp1          +  2.*ujm1 - ujm2 );
	_T D42uj =        (  ujp2 -  4.*ujp1 +  6.*uj -  4.*ujm1 + ujm2 );

	return D14uj * ( D14uj +       D24uj +   1./3.*D32uj +    1./12.*D42uj )
	               + D24uj*( 4./3.*D24uj +   5./4.*D32uj +     2./5.*D42uj )
	                             + D32uj*( 83./60.*D32uj +   23./18.*D42uj )
	                                                     + 437./315.*D42uj*D42uj;
}

template <typename _T>
_T
ubjm12 ( _T ujm3 , _T ujm2 , _T ujm1 , _T uj , _T ujp1 , _T ujp2 , bool velocity_is_positive=true )
{
  _T betaL = beta(ujm3,ujm2,ujm1,uj,ujp1);
  _T betaR = beta(ujm2,ujm1,uj,ujp1,ujp2);

  _T aL = 0.5/(1e-6 + betaL);
  _T aR = 0.5/(1e-6 + betaR);

  _T wtL = aL/(aL+aR);
  _T wtR = aR/(aL+aR);

  _T wL = std::max(wtL,wtR);
  _T wR = std::min(wtL,wtR);
  if (!velocity_is_positive) { std::swap(wL,wR); }

  _T uL = 1./60.*( -3*ujp1 + 27.*uj   + 47.*ujm1 - 13*ujm2 + 2*ujm3 );
  _T uR = 1./60.*(  2.*ujp2 - 13.*ujp1 + 47.*uj   + 27*ujm1 - 3*ujm2 );

  return wL*uL + wR*uR;
}

template <typename Iterator>
inline auto
ubjm12 ( Iterator it , bool velocity_is_positive=true ) {
  return ubjm12(*(it-3),*(it-2),*(it-1),*(it),*(it+1),*(it+2),velocity_is_positive);
}

template <typename Container,typename _T>
Container
du ( Container const& u , _T const& v , _T dx )
{
  Container du(u.size());
  std::size_t N = u.size();

  for ( std::size_t i=0 ; i<3 ; ++i ) {
    du[i] = v*(ubjm12(u[(i-2+N)%N],u[(i-1+N)%N],u[i],u[i+1],u[i+2],u[i+3],v>0) - ubjm12(u[(i-3+N)%N],u[(i-2+N)%N],u[(i-1+N)%N],u[i],u[i+1],u[i+2],v>0))/dx;
  }
  for ( std::size_t i=3 ; i<u.size()-3 ; ++i ) {
    du[i] = v*( ubjm12( &u[i+1] , v>0 ) - ubjm12( &u[i] , v>0 ) )/dx;
  }
  for ( std::size_t i=u.size()-3 ; i<u.size() ; ++i ) {
    du[i] = v*( ubjm12(u[i-2],u[i-1],u[i],u[(i+1)%N],u[(i+2)%N],u[(i+3)%N],v>0) - ubjm12(u[i-3],u[i-2],u[i-1],u[i],u[(i+1)%N],u[(i+2)%N],v>0) )/dx; 
  }

  return du;
}


} // namespace bweno
#endif

