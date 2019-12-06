#ifndef _LAGRANGE5_H_
#define _LAGRANGE5_H_

namespace lagrange5 {

template < typename _T >
auto
generator ( _T y0 , _T y1 , _T y2 , _T y3 , _T y4 , _T y5 , _T dx , _T x3 )
{
  _T y0y1         = (-y0 +    y1)/dx;
  _T y0y1y2       = ( y0 - 2.*y1 +     y2)/(2.*dx*dx);
  _T y0y1y2y3     = (-y0 + 3.*y1 -  3.*y2 +     y3)/(6.*dx*dx*dx);
  _T y0y1y2y3y4   = ( y0 - 4.*y1 +  6.*y2 -  4.*y3 +    y4)/(24.*dx*dx*dx*dx);
  _T y0y1y2y3y4y5 = (-y0 + 5.*y1 - 10.*y2 + 10.*y3 - 5.*y4 + y5)/(120.*dx*dx*dx*dx*dx);

  return [dx,x3,y0,y0y1,y0y1y2,y0y1y2y3,y0y1y2y3y4,y0y1y2y3y4y5] ( double x ) -> double {
    /**
      i-3  i-2  i-1   i   i+1  i+2
     --|----|----|----|----|----|--
       0    1    2    3    4    5
       so :
       x1 = x3+dx
       x2 = x3+2*dx
       x3 = x3+3*dx
       x4 = x3+4*dx
       x5 = x3+5*dx
    **/
    double x0=x3-3.*dx , x1=x3-2.*dx , x2=x3-dx , x4=x3+dx , x5=x3+2.*dx;
    //std::cout << "\t" << x0 << " " << x1 << " " << x2 << " *" << x << "* " << x3 << " " << x4 << " " << x5 << "\n";
    double y = y0 + (x-x0)*( y0y1 + (x-x1)*( y0y1y2 + (x-x2)*( y0y1y2y3 + (x-x3)*( y0y1y2y3y4 + (x-x4)*y0y1y2y3y4y5 ) ) ) );

    return y;
  };
}

} // namespace lagrange5

#endif
