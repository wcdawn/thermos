module geometry
use kind, only : rk, ik
implicit none

private

public :: geometry_calculate_coordinates

contains

  subroutine geometry_calculate_coordinates(length, nx, xcenter, dx)
    real(rk), intent(in) :: length
    integer(ik), intent(in) :: nx
    real(rk), intent(out) :: xcenter(:) ! Center coordinate for cell [cm]
    real(rk), intent(out) :: dx(:) ! Cell-width [cm]
  endsubroutine geometry_calculate_coordinates

endmodule geometry
