module conductivity_function
use kind, only : rk
implicit none

private

public :: conductivity_fun

contains

  pure elemental real(rk) function conductivity_fun(T)
    real(rk), intent(in) :: T ! [K] temperature
    real(rk), parameter :: k0 = 1.2_rk
    conductivity_fun = k0
  endfunction conductivity_fun

endmodule conductivity_function
