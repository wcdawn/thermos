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
    
    real(rk) :: h
    integer(ik) :: i

    ! NOTE: For now, this calculates uniformly spaced coordinates regardless of
    ! the geometry. In the future, it would be nice to support equal-area
    ! (cylindrical) and equal-volume (spherical) coordiantes for other
    ! geometries.

    h = length / nx
    dx = h
    do i = 1,nx
      xcenter(i) = (i - 0.5_rk) * h
    enddo ! i = 1,nx

  endsubroutine geometry_calculate_coordinates

endmodule geometry
