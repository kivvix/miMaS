module numeric
implicit none
  integer , parameter :: sp_ = kind(1.0  )
  integer , parameter :: dp_ = kind(1.0d0)
  integer , parameter :: rp = dp_

  real(rp) , parameter :: pi = 3.1415926535897932384626433832795_rp

  real(rp) :: v,dt,dx
end module numeric
