module time
  use numeric
implicit none
contains

function erkssp33 ( u , L , dt ) result(un1)
  real(rp) , dimension(0:) :: u
    interface
      function L(ui) result(uo)
        import rp
        real(rp) , dimension(0:) :: ui
        real(rp) , dimension(0:size(ui)-1) :: uo
    end function L
  end interface
  real(rp) :: dt

  real(rp) , dimension(0:size(u)-1) :: un1
  real(rp) , dimension(0:size(u)-1) :: u1,u2

  u1  = u + dt*L(u)
  u2  = 0.75_rp*u + 0.25_rp*u1 + 0.25*dt*L(u1)
  un1 = 1._rp/3._rp*u + 2._rp/3._rp*u2 + 2._rp/3._rp*dt*L(u2)
end function erkssp33

function erkssp43 ( u , L , dt ) result(un1)
  real(rp) , dimension(0:) :: u
    interface
      function L(ui) result(uo)
        import rp
        real(rp) , dimension(0:) :: ui
        real(rp) , dimension(0:size(ui)-1) :: uo
    end function L
  end interface
  real(rp) :: dt

  real(rp) , dimension(0:size(u)-1) :: un1

  un1 = u + 0.5_rp*dt*L(u)
  un1 = un1 + 0.5_rp*dt*L(un1)
  un1 = 2._rp/3._rp*u + 1._rp/3._rp*un1 + 1._rp/6._rp*dt*L(un1)
  un1 = un1 + 0.5_rp*dt*L(un1)
end function erkssp43

function erk76 ( u , L , dt ) result(un1)
  real(rp) , dimension(0:) :: u
    interface
      function L(ui) result(uo)
        import rp
        real(rp) , dimension(0:) :: ui
        real(rp) , dimension(0:size(ui)-1) :: uo
    end function L
  end interface
  real(rp) :: dt

  real(rp) , dimension(0:size(u)-1) :: un1

  real(rp) , parameter :: nu = 0.33333_rp
  real(rp) , dimension(0:size(u)-1) :: k1,k2,k3,k4,k5,k6,k7

  k1 = dt*L( u )
  k2 = dt*L( u + nu*k1 )
  k3 = dt*L( u + ((4._rp*nu-1)*k1 + k2)/(8._rp*nu) )
  k4 = dt*L( u + ((10._rp*nu-2._rp)*k1 + 2._rp*k2 + 8._rp*nu*k3)/(27._rp*nu) )
  k5 = dt*L( u + ( -((77._rp*nu-56._rp)+(17._rp*nu-8._rp)*sqrt(21._rp))*k1           &
                  -8._rp*(7._rp+sqrt(21._rp))*k2 + 48._rp*(7._rp+sqrt(21._rp))*nu*k3 &
                  -3._rp*(21._rp+sqrt(21._rp))*nu*k4 )/(392._rp*nu) )
  k6 = dt*L( u + ( -5._rp*((287._rp*nu-56._rp) - (59._rp*nu-8._rp)*sqrt(21._rp))*k1 &
              - 40._rp*(7._rp-sqrt(21._rp))*k2 + 320._rp*sqrt(21._rp)*nu*k3 + 3._rp*(21-121*sqrt(21._rp))*nu*k4 &
              + 392._rp*(6._rp-sqrt(21._rp))*nu*k5 )/(1960._rp*nu) )
  k7 = dt*L( u + ( 15._rp*((30._rp*nu-8._rp)-(7._rp*nu*sqrt(21._rp)))*k1 + 120._rp*k2 &
                   -40._rp*(5._rp+7._rp*sqrt(21._rp))*nu*k3 + 63._rp*(2._rp+3._rp*sqrt(21._rp))*nu*k4 &
                   - 14._rp*(49._rp-9._rp*sqrt(21._rp))*nu*k5 + 70._rp*(7._rp+sqrt(21._rp))*nu*k6)/(180._rp*nu) )

  un1 = u + (9._rp*k1 + 64._rp*k3 + 49._rp*k5 + 49._rp*k6 + 9._rp*k7)/180._rp
end function erk76

function euler ( u , L , dt ) result(un1)
  real(rp) , dimension(0:) :: u
    interface
      function L(ui) result(uo)
        import rp
        real(rp) , dimension(0:) :: ui
        real(rp) , dimension(0:size(ui)-1) :: uo
    end function L
  end interface
  real(rp) :: dt

  real(rp) , dimension(0:size(u)-1) :: un1
  
  un1 = u + dt*L(u)
end function euler

function fft_time ( u ,  x0 , x1 , dt ) result(un1)
  real(rp) , dimension(:) :: u
  real(rp) :: x0,x1
  real(rp) :: dt
  real(rp) , dimension(size(u)) :: un1

  real(rp) , dimension(2*size(u)+15) :: coefd
  real(rp) :: tmp,im,re
  integer  :: i
  real(rp) :: k
  integer  :: sizeu
  sizeu = size(u)
  tmp = 0.5_rp*((x1-x0)/pi) !/real(sizeu,rp)

  ! initiatilize working array
  call dffti(sizeu,coefd)

  ! compute FFT of rho-1
  call dfftf(sizeu,u,coefd)
  un1 = u
  un1(1) = 0._rp ! E est de moyenne nulle donc le premier coeffcient de Fourrier est nul
 
  do i=1,(sizeu)/2-1
    re = u(2*i) ; im=u(2*i+1)
    k = real(i,rp)/tmp
    un1(2*i)   = cos(k*dt)*re + sin(k*dt)*im
    un1(2*i+1) = cos(k*dt)*im - sin(k*dt)*re
  end do

  if (mod(sizeu,2)==0) un1(sizeu-1) = 0._rp
  call dfftb(sizeu,un1,coefd) ! fft inverse
  un1 = un1/sizeu
end function fft_time

end module time
