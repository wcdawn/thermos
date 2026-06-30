module finite_element
use kind, only : ik, rk
implicit none

private

public :: finite_element_solve

contains

  subroutine finite_element_solve(geometry, nx, xcenter, dx, &
      bctype_left, bctype_right, bcval_left, bcval_right, &
      max_iter, tol_temperature, init_temperature, temperature, TCL)
    character(*), intent(in) :: geometry
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: xcenter(:) ! (nx)
    real(rk), intent(in) :: dx(:) ! (nx)

    character(*), intent(in) :: bctype_left, bctype_right
    real(rk), intent(in) :: bcval_left, bcval_right

    integer(ik), intent(in) :: max_iter
    real(rk), intent(in) :: tol_temperature ! [K]
    real(rk), intent(in) :: init_temperature ! [K]
    real(rk), intent(out) :: temperature(:) ! (nx) [K]
    real(rk), intent(out) :: TCL ! [K] centerline temperature
  endsubroutine finite_element_solve

endmodule finite_element
