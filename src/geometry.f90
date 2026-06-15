module geometry
use kind, only : rk, ik
implicit none

private

public :: geometry_calculate_coordinates, geometry_refine, geometry_summary

contains

  subroutine geometry_calculate_coordinates(length, nx, xcenter, dx)
    real(rk), intent(in) :: length
    integer(ik), intent(in) :: nx
    real(rk), intent(out) :: xcenter(:) ! Center coordinate for cell [cm]
    real(rk), intent(out) :: dx(:) ! Cell-width [cm]
    
    real(rk) :: h

    ! NOTE: For now, this calculates uniformly spaced coordinates regardless of
    ! the geometry. In the future, it would be nice to support equal-area
    ! (cylindrical) and equal-volume (spherical) coordiantes for other
    ! geometries.

    h = length / nx
    dx = h
    call geometry_dx2xcenter(nx, dx, xcenter)
  endsubroutine geometry_calculate_coordinates

  subroutine geometry_refine(nx, xcenter, dx)
    integer(ik), intent(inout) :: nx
    real(rk), allocatable, intent(inout) :: xcenter(:) ! (nx)
    real(rk), allocatable, intent(inout) :: dx(:) ! (nx)

    integer(ik) :: i, idx

    integer(ik) :: nx_old
    real(rk), allocatable :: dx_old(:) ! (nx)

    ! NOTE: This will perform uniform refinement and split each element in half.
    ! This may be alright, but will not be "consistent" with equal-area or
    ! equal-volume discretizations.

    nx_old = nx
    allocate(dx_old(nx_old))

    dx_old = dx

    nx = 2*nx_old
    deallocate(dx)
    allocate(dx(nx))

    do i = 1,nx_old
      idx = 2*(i-1) + 1
      dx(idx+0) = 0.5_rk * dx_old(i)
      dx(idx+1) = 0.5_rk * dx_old(i)
    enddo ! i = 1,nx_old

    deallocate(xcenter)
    allocate(xcenter(nx))
    call geometry_dx2xcenter(nx, dx, xcenter)

    deallocate(dx_old)
  endsubroutine geometry_refine

  subroutine geometry_dx2xcenter(nx, dx, xcenter, xstart)
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    real(rk), intent(out) :: xcenter(:) ! (nx)
    real(rk), intent(in), optional :: xstart

    integer(ik) :: i
    real(rk) :: x0

    if (present(xstart)) then
      x0 = xstart
    else
      x0 = 0.0_rk
    endif

    do i = 1,nx
      xcenter(i) = x0 + 0.5_rk*dx(i)
      x0 = x0 + dx(i)
    enddo ! i = 1,nx
  endsubroutine geometry_dx2xcenter

  subroutine geometry_summary(nx, dx)
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    write(*, '(a)') '=== GEOMETRY SUMMARY ==='
    write(*, '(a, i0)') 'Number of cells: ', nx
    write(*, '(a, es13.6)') 'Minimum DX: ', minval(dx)
    write(*, '(a, es13.6)') 'Maximum DX: ', maxval(dx)
    write(*, *)
  endsubroutine

endmodule geometry
