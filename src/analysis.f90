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
        ! do something
        call temperature_exact_slab_cos(nx, xcenter, texact)
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
    enddo
  endsubroutine temperature_exact_slab_cos

endmodule analysis
