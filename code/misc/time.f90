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

!function poisson1D_per(rho) result(E)
!  ! pour se simiplifier la vie avec la FFTPACK
!  ! les indices commencent à 1 uniquement dans cette fonction
!  ! Fortran est capable de réadapter les indices au programme principal
!  real(rp) , dimension(d%nx) :: rho
!  real(rp) , dimension(d%nx) :: E
!
!  real(rp) , dimension(2*(d%nx)+15) :: coefd
!  real(rp) :: tmp,im,re
!  integer  :: i
!  tmp = 0.5_rp*((d%x_max-d%x_0)/pi)/real(d%nx,rp)
!
!  ! initiatilize working array
!  call dffti(d%nx,coefd)
!
!  do i=1,d%nx
!    E(i) = rho(i) - 1._rp
!  end do
!
!  ! compute FFT of rho-1
!  call dfftf(d%nx,E,coefd)
!  E(1) = 0._rp ! E est de moyenne nulle donc le premier coeffcient de Fourrier est nul
! 
!  do i=1,(d%nx)/2-1
!    re = E(2*i) ; im=E(2*i+1)
!    E(2*i)   =  tmp/real(i,rp)*im ! Re(\hat{E}) = Im(\hat{rho})/k
!    E(2*i+1) = -tmp/real(i,rp)*re ! Im(\hat{E}) = -Re(\hat{rho})/k
!  end do
!
!  if (mod(d%nx,2)==0) E(d%nx) = 0._rp
!  call dfftb(d%nx,E,coefd) ! fft inverse
!end function poisson1D_per

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
