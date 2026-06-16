module source_function
use kind, only : rk
implicit none

private

public :: source_fun

contains

  pure real(rk) function source_fun(x)
    use constants, only : pi
    real(rk), intent(in) :: x
    real(rk), parameter :: q0 = 20.0_rk
    real(rk), parameter :: L = 10.0_rk
    source_fun = q0 * cos(pi * x / L)
  endfunction source_fun

endmodule source_function
