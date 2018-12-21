program translate
use numeric
use weno
use time
implicit none

#define X(i) (real(i,rp)*dx+x0)

integer :: N
real(rp) , dimension(:) , allocatable :: u,u0
real(rp)                    :: T_final
real(rp)                    :: k,x0,x1
!real(rp)                    :: v,dx,dt

integer :: i, i_time, z


T_final = 10000_rp*pi
v = +1._rp
x0 = 0._rp
x1 = 2._rp*pi

k = 2._rp
dx = 1.38337934_rp/(k*2._rp*pi)

N = int((x1-x0)/dx)
allocate( u(0:N-1) , u0(0:N-1) )

!dt = 1.433_rp*dx
!dt = 1.606_rp*dx
dt = 0.5_rp*dx

!! INIT
do i = 0,N-1
  u0(i) = cos(k*X(i))
  !if (i<0.3*N) then
  !  u0(i) = 0.5_rp*X(i)
  !elseif (i<0.6*N) then
  !  u0(i) = 1._rp
  !else
  !  u0(i) = 0._rp
  !end if
end do
u=u0

open(newunit=z,file="u_init.dat")
do i=0,N-1
  write(z,*) X(i),u(i)
end do
close(z)

!! LOOP
print *, "(x0,x1): ",x0,x1
print *, "(dx,dt): ",dx,dt
print *, "N: " ,N
print *, "Tf: ",T_final 

write (*,*) T_final/dt
i_time = 0
do while ( i_time*dt < T_final )
!do while ( i_time <= 200 )
  write (*,"(A1,I0.4)",advance='no') 13,i_time
  u = erk76( u , L , dt )
  !u = euler( u , L , dt )
  i_time = i_time+1
end do

write (*,"(A1,I0.4)") 13,i_time
open(newunit=z,file="u_vp.dat")
do i=0,N-1
  write(z,*) X(i),u(i),u0(i)
  !write(z,*) X(i),u(i),cos(k*X(i)-(2._rp*real(i_time,rp))*dt)
end do
close(z)

deallocate(u,u0)

contains

  function L(u)
    real(rp) , dimension(0:) :: u
    real(rp) , dimension(0:size(u)-1) :: L

    L = -trp1D_w( u , v , dx )
  end function L

end program translate

