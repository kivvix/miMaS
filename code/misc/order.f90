program order 
use numeric
use weno
use time
implicit none

#define X(i) (real(i,rp)*dx+x0)
#define SQ(X) ((X)*(X))

integer  :: N,u
real(rp) :: t1, t2
open(newunit=u,file="order.dat")

t1 = 5._rp*2._rp*pi/10*1d-7
!t2 = 5*t1
do N=10,150,10
  !call trp_error_1step(N,u)
  call trp_error_time(N,t1,t2,u)
end do

close(u)

contains

function errors( un , us ) result(e)
  real(rp) , dimension(0:)  :: un,us
  real(rp) , dimension(1:3) :: e

  integer :: i
  e = 0._rp
  do i=0,size(un)-1
    e(1) = e(1) + abs(un(i)-us(i))*dx
    e(2) = e(2) + SQ(un(i)-us(i))*dx*dx
    e(3) = max( e(3) , abs(un(i)-us(i)) )
  end do
  e(2) = sqrt(e(2))
end function errors

subroutine trp_error_1step ( N , u )
  integer :: N ! number of points
  integer :: u ! unit of the file

  real(rp) :: x0,x1
  
  real(rp) , dimension(1:3) :: e
  integer  :: i,z

  real(rp) , dimension(0:N-1) :: un , us

  v = 1._rp
  x0 = 0._rp
  x1 = 2._rp*pi
  
  dx = (x1-x0)/real(N,rp)
  dt = 1d-7*dx

  ! INIT
  do i = 0,N-1
    un(i) = cos(X(i))
    us(i) = cos(X(i)-dt)
  end do

  ! LOOP
  un = ssprk33(un,L,dt)
  !un = euler(un,L,dt)

  !open(newunit=z,file="test.dat")
  !do i=0,N-1
  !  write(z,*) X(i),cos(X(i)),un(i),us(i)
  !end do
  !close(z)

  ! ERROR
  e = errors(un,us)
  write(u,"(I9,E15.5,E15.5,E15.5,E15.5)") N,dx,e(1),e(2),e(3)
  call flush()
end subroutine trp_error_1step

subroutine trp_error_time ( N , t1 , t2 , u )
  integer  :: N     ! number of points
  real(rp) :: t1,t2 ! times where we compute error
  integer  :: u     ! unit of the file

  real(rp) :: x0,x1
  
  real(rp) , dimension(1:3) :: e1,e2
  integer  :: i,i_time

  real(rp) , dimension(0:N-1) :: un , us1 , us2

  v  = 1._rp
  x0 = 0._rp
  x1 = 2._rp*pi
  
  dx = (x1-x0)/real(N,rp)
  dt = 1d-6*dx

  write(u,"(I3,E15.5)",advance='no') N,dx
  call flush()

  ! INIT
  do i = 0,N-1
    un(i)  = cos(X(i))
  end do

  ! LOOP 1
  i_time = 0
  do while ( dt*i_time <= t1 )
    un = ssprk33(un,L,dt)
    i_time = i_time+1
  end do

  do i=0,N-1
    us1(i) = cos(X(i)-real(i_time,rp)*dt)
  end do
  ! ERROR 1
  e1 = errors(un,us1)
  
  !write(u,"(I9,E15.5,E15.5,E15.5)",advance='no') i_time,e1(1),e1(2),e1(3)
  write(u,"(I9,E15.5,E15.5,E15.5)") i_time,e1(1),e1(2),e1(3)
  call flush()

  !! LOOP 2
  !do while ( dt*i_time <= t2 )
  !  un = ssprk33(un,L,dt)
  !  i_time = i_time+1
  !end do
  !do i=0,N-1
  !  us2(i) = cos(X(i)-real(i_time,rp)*dt)
  !end do
  !! ERROR 2
  !e2 = errors(un,us2)
  !
  !write(u,"(I9,E15.5,E15.5,E15.5)") i_time,e2(1),e2(2),e2(3)
  call flush()
end subroutine trp_error_time

function L(u)
  real(rp) , dimension(0:) :: u
  real(rp) , dimension(0:size(u)-1) :: L

  L = -trp1D_w( u , v , dx )
end function L

end program order

