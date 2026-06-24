module geometry
use kind, only : rk, ik
implicit none

private

public :: geometry_calculate_coordinates, geometry_refine, geometry_summary

contains

  subroutine geometry_calculate_coordinates(mesh_spacing, length, nx, xcenter, dx)
    use output, only : output_write
    character(*), intent(in) :: mesh_spacing
    real(rk), intent(in) :: length
    integer(ik), intent(in) :: nx
    real(rk), intent(out) :: xcenter(:) ! Center coordinate for cell [cm]
    real(rk), intent(out) :: dx(:) ! Cell-width [cm]
    
    select case (mesh_spacing)
      case ('uniform')
        call geometry_uniform_mesh(length, nx, dx)
      case ('area')
        call geometry_area_mesh(length, nx, dx)
      case ('volume')
        call geometry_uniform_mesh(length, nx, dx) ! TODO
      case default
        call output_write('ERROR: unknown mesh spacing: ' // &
          trim(adjustl(mesh_spacing)))
        stop
    endselect

    call geometry_dx2xcenter(nx, dx, xcenter)
  endsubroutine geometry_calculate_coordinates

  subroutine geometry_uniform_mesh(length, nx, dx)
    real(rk), intent(in) :: length
    integer(ik), intent(in) :: nx
    real(rk), intent(out) :: dx(:) ! (nx)
    dx = length / nx
  endsubroutine geometry_uniform_mesh

  subroutine geometry_area_mesh(length, nx, dx)
    use constants, only : pi
    real(rk), intent(in) :: length
    integer(ik), intent(in) :: nx
    real(rk), intent(out) :: dx(:) ! (nx)

    integer(ik) :: i
    real(rk) :: aequal
    real(rk) :: rprev, rthis

    aequal = pi * length**2 / nx

    rprev = 0.0_rk
    do i = 1,nx
      rthis = sqrt(aequal / pi + rprev**2)
      dx(i) = rthis - rprev
      rprev = rthis
    enddo
  endsubroutine geometry_area_mesh

  subroutine geometry_volume_mesh(length, nx, dx)
    use constants, only : pi
    real(rk), intent(in) :: length
    integer(ik), intent(in) :: nx
    real(rk), intent(out) :: dx(:) ! (nx)

    integer(ik) :: i
    real(rk) :: vequal
    real(rk) :: rprev, rthis

    vequal = 4.0_rk * pi * length**3 / 3.0_rk

    rprev = 0.0_rk
    do i = 1,nx
      rthis = (0.75_rk / pi * vequal + rprev**3)**(1.0_rk/3.0_rk)
    enddo ! i = 1,nx
  endsubroutine geometry_volume_mesh

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
    use output, only : output_write
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    character(1024) :: line
    call output_write('=== GEOMETRY SUMMARY ===')
    write(line, '(a, i0)') 'Number of cells: ', nx
    call output_write(line)
    write(line, '(a, es13.6)') 'Minimum DX: ', minval(dx)
    call output_write(line)
    write(line, '(a, es13.6)') 'Maximum DX: ', maxval(dx)
    call output_write(line)
    call output_write('')
  endsubroutine

endmodule geometry
