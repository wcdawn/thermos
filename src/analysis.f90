module analysis
use kind, only : rk, ik
implicit none

private

public :: analysis_analyze

contains

  subroutine analysis_analyze(analysis_name, nx, xcenter, temperature)
    use output, only : output_write
    use linalg, only : norm
    character(*), intent(in) :: analysis_name
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: xcenter(:) ! (nx)
    real(rk), intent(in) :: temperature(:) ! (nx)

    character(1024) :: line

    real(rk), allocatable :: texact(:) ! (nx)

    allocate(texact(nx))

    select case (analysis_name)
      case ('slab_cos')
        call temperature_exact_slab_cos(nx, xcenter, texact)
      case ('cyl_lin')
        call temperature_exact_cyl_lin(nx, xcenter, texact)
      case ('sph_lin')
        call temperature_exact_sph_lin(nx, xcenter, texact)
      case default
        call output_write('ERROR: unknown analysis name: ' // &
          trim(adjustl(analysis_name)))
        stop
    endselect

    call output_write('=== ANALYSIS SUMMARY ===')
    write(line, '(a,es9.2)') 'Linf error = ', norm(-1, temperature - texact)
    call output_write(line)
    call output_write('')

    deallocate(texact)
  endsubroutine analysis_analyze

  subroutine temperature_exact_slab_cos(nx, x, texact)
    use constants, only : pi
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: x(:) ! (nx)
    real(rk), intent(out) :: texact(:) ! (nx)

    ! TODO Try to be more clever about these
    real(rk), parameter :: k = 1.2_rk
    real(rk), parameter :: q0 = 20.0_rk
    real(rk), parameter :: T0 = 600.0_rk
    real(rk), parameter :: TL = 300.0_rk
    real(rk), parameter :: L = 10.0_rk

    integer(ik) :: i

    do i = 1,nx
      texact(i) = T0 &
        + (x(i)/L) * (TL - T0 + q0/k * (2.0_rk*L/pi)**2) &
        + q0/k * (2.0_rk*L/pi)**2 * (cos(pi*x(i)*0.5_rk/L) - 1.0_rk)
    enddo ! i = 1,nx
  endsubroutine temperature_exact_slab_cos

  subroutine temperature_exact_cyl_lin(nx, x, texact)
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: x(:) ! (nx)
    real(rk), intent(out) :: texact(:) ! (nx)

    integer(ik) :: i

    ! TODO Try to be more clever about these
    real(rk), parameter :: k0 = 0.5_rk
    real(rk), parameter :: q0 = 500.0_rk
    real(rk), parameter :: TR = 600.0_rk
    real(rk), parameter :: R = 0.75_rk

    do i = 1,nx
      texact(i) = -q0/k0 * (x(i)**2*0.25_rk - x(i)**3/(9.0_rk*R)) &
        + TR + q0/k0 * 5.0_rk/36.0_rk * R**2
    enddo ! i = 1,nx
  endsubroutine temperature_exact_cyl_lin

  subroutine temperature_exact_sph_lin(nx, x, texact)
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: x(:) ! (nx)
    real(rk), intent(out) :: texact(:) ! (nx)

    integer(ik) :: i

    ! TODO Try to be more clever about these
    real(rk), parameter :: k0 = 0.75_rk
    real(rk), parameter :: q0 = 200.0_rk
    real(rk), parameter :: TR = 600.0_rk
    real(rk), parameter :: R = 3.0_rk

    do i = 1,nx
      texact(i) = -q0/k0 * (x(i)**2/6.0_rk - x(i)**3/(12.0_rk*R)) &
        + TR + q0/k0 * R**2/12.0_rk
    enddo ! i = 1,nx
  endsubroutine temperature_exact_sph_lin

endmodule analysis
