module time
  use numeric
implicit none
contains

function ssprk33 ( u , L , dt ) result(un1)
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
end function ssprk33

function ssprk34 ( u , L , dt ) result(un1)
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
end function ssprk34

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

end module time
