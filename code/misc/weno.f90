module weno
  use numeric
implicit none
contains

#define SQ(X) ((X)*(X))
function flux_WENO_p ( f ) result(fpip12)
  real(rp) , dimension(-2:2) , intent(in) :: f
  real(rp)                                :: fpip12

  real(rp) :: w0 , w1 , w2
  real(rp) :: sum_w
  real(rp) , parameter :: epsi = 1d-6

  w0 = 13._rp/12._rp*SQ( f(-2) - 2._rp*f(-1) + f(+0) ) &
            +0.25_rp*SQ( f(-2) - 4._rp*f(-1) + 3._rp*f(+0) )
  w1 = 13._rp/12._rp*SQ( f(-1) - 2._rp*f(+0) + f(+1) ) &
            +0.25_rp*SQ( f(-1) - f(+1) )
  w2 = 13._rp/12._rp*SQ( f(+0) - 2._rp*f(+1) + f(+2) ) &
            +0.25_rp*SQ( 3._rp*f(+0) - 4._rp*f(+1) + f(+2) )
  w0 = 0.1_rp /SQ(epsi+w0)
  w1 = 0.6_rp /SQ(epsi+w1)
  w2 = 0.3_rp /SQ(epsi+w2)

  sum_w = w0+w1+w2
  w0 = w0/sum_w
  w1 = w1/sum_w
  w2 = w2/sum_w

  fpip12 = w0*( (2._rp/6._rp)*f(-2) - (7._rp/6._rp)*f(-1) + (11._rp/6._rp)*f(+0) ) &
         + w1*(-(1._rp/6._rp)*f(-1) + (5._rp/6._rp)*f(+0) +  (2._rp/6._rp)*f(+1) ) &
         + w2*( (2._rp/6._rp)*f(+0) + (5._rp/6._rp)*f(+1) -  (1._rp/6._rp)*f(+2) )
end function flux_WENO_p

function flux_WENO_m ( f ) result(fmip12)
  real(rp) , dimension(-1:3) , intent(in) :: f
  real(rp)                                :: fmip12

  real(rp) :: w0 , w1 , w2
  real(rp) :: sum_w
  real(rp) , parameter :: epsi = 1d-6

  w0 = 13._rp/12._rp*SQ( f(+1) - 2._rp*f(+2) + f(+3) ) &
            +0.25_rp*SQ( 3._rp*f(+1) - 4._rp*f(+2) + f(+3) )
  w1 = 13._rp/12._rp*SQ( f(+0) - 2._rp*f(+1) + f(+2) ) &
            +0.25_rp*SQ( f(+0) - f(+2) )
  w2 = 13._rp/12._rp*SQ( f(-1) - 2._rp*f(+0) + f(+1) ) &
            +0.25_rp*SQ( f(-1) - 4._rp*f(+0) + 3._rp*f(+1) )
  w0 = 0.1_rp /SQ(epsi+w0)
  w1 = 0.6_rp /SQ(epsi+w1)
  w2 = 0.3_rp /SQ(epsi+w2)

  sum_w = w0+w1+w2
  w0 = w0/sum_w
  w1 = w1/sum_w
  w2 = w2/sum_w

  fmip12 = w2*(-(1._rp/6._rp)*f(-1) + (5._rp/6._rp)*f(+0) +  (2._rp/6._rp)*f(+1) ) &
         + w1*( (2._rp/6._rp)*f(+0) + (5._rp/6._rp)*f(+1) -  (1._rp/6._rp)*f(+2) ) &
         + w0*((11._rp/6._rp)*f(+1) - (7._rp/6._rp)*f(+2) +  (2._rp/6._rp)*f(+3) )
end function flux_WENO_m
#undef SQ

function WENO1d_p ( f , v_p ) result(fpip12)
  real(rp) , dimension(0:) :: f
  real(rp)                 :: v_p
  
  real(rp) , dimension(0:size(f)-1) :: fpip12

  integer :: i
  integer :: im2,im1,ip1,ip2

  do i=2,size(f)-3
    fpip12(i) = v_p*flux_WENO_p( f(i-2:i+2) )
  end do

  do i=0,1
    im2 = mod(i-2+size(f),size(f)) ; im1 = mod(i-1+size(f),size(f))
    fpip12(i) = v_p*flux_WENO_p( [f(im2),f(im1),f(i),f(i+1),f(i+2)] )
  end do

  do i=size(f)-2,size(f)-1
    ip1 = mod(i+1,size(f)) ; ip2 = mod(i+2,size(f))
    fpip12(i) = v_p*flux_WENO_p( [f(i-2),f(i-1),f(i),f(ip1),f(ip2)] )
  end do
end function WENO1d_p

function WENO1d_m ( f , v_m ) result(fmip12)
  real(rp) , dimension(0:) :: f
  real(rp)                 :: v_m
  
  real(rp) , dimension(0:size(f)-1) :: fmip12

  integer :: i
  integer :: im1,ip1,ip2,ip3

  do i=1,size(f)-4
    fmip12(i) = v_m*flux_WENO_m( f(i-1:i+3) )
  end do

  i=0
    im1 = size(f)-1 
    fmip12(i) = v_m*flux_WENO_m( [f(im1),f(i),f(i+1),f(i+2),f(i+3)] )

  do i=size(f)-3,size(f)-1
    ip1 = mod(i+1,size(f)) ; ip2 = mod(i+2,size(f)) ; ip3 = mod(i+3,size(f))
    fmip12(i) = v_m*flux_WENO_m( [f(i-1),f(i),f(ip1),f(ip2),f(ip3)] )
  end do
end function WENO1d_m

function trp1D_w( f , v , dx ) result(trpf)
  real(rp) , dimension(0:) :: f
  real(rp)                 :: v,dx

  real(rp) , dimension(0:size(f)-1) :: trpf

  real(rp) :: vp,vm
  real(rp) , dimension(0:size(f)-1) :: fp,fm

  integer  :: i
  integer  :: im1
  real(rp) :: ddx
  ddx = 1._rp/dx

  vp = max(v,0._rp)
  vm = min(v,0._rp)

  fp = WENO1d_p( f , vp )
  fm = WENO1d_m( f , vm )

  do i=1,size(f)-1
    trpf(i) = ddx*( fp(i)-fp(i-1) + fm(i)-fm(i-1) )
  end do

  i=0
    im1 = size(f)-1
    trpf(i) = ddx*( fp(i)-fp(im1) + fm(i)-fm(im1) )

end function trp1D_w

end module weno
