module analysis
use kind, only : rk, ik
implicit none

private

public :: analysis_analyze

contains

  subroutine analysis_analyze(analysis_name, fname_out, nx, xcenter, temperature)
    use output, only : output_write
    use linalg, only : norm
    character(*), intent(in) :: analysis_name
    character(*), intent(in) :: fname_out
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
      case ('slab_cos_klin')
        call temperature_exact_slab_cos_klin(nx, xcenter, texact)
      case ('cyl_lin_klin')
        call temperature_exact_cyl_lin_klin(nx, xcenter, texact)
      case default
        call output_write('ERROR: unknown analysis name: ' // &
          trim(adjustl(analysis_name)))
        stop
    endselect

    call output_write('=== ANALYSIS SUMMARY ===')
    write(line, '(a,es9.2)') 'Linf error = ', norm(-1, temperature - texact)
    call output_write(line)
    call output_write('Writing analysis output on: ' // trim(adjustl(fname_out)))
    call analysis_output(fname_out, nx, xcenter, temperature, texact)
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

  subroutine temperature_exact_slab_cos_klin(nx, x, texact)
    use constants, only : pi
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: x(:) ! (nx)
    real(rk), intent(out) :: texact(:) ! (nx)

    integer(ik) :: i

    ! Based on the Kirchoff transform.
    ! TODO Try to be more clever about these
    real(rk), parameter :: q0 = 400.0_rk ! [W/cc]
    real(rk), parameter :: k0 = 10.0_rk  ! [W/cm/K]
    real(rk), parameter :: beta = 1d-2 ! [W/cm/K^2]
    real(rk), parameter :: T0 = 500.0 ! [K]
    real(rk), parameter :: TL = 300.0 ! [K]
    real(rk), parameter :: L = 5.0 ! [cm]

    real(rk), parameter :: c1 = (q0 * (2.0_rk * L / pi)**2 &
      + 0.5_rk * beta * (T0**2 - TL**2) + k0 * (TL - T0)) / L
    real(rk), parameter :: c2 = 0.5_rk * (k0**2 - (beta * T0 - k0)**2) / beta &
      - q0 * (2.0_rk * L / pi)**2

    do i = 1,nx
      texact(i) = &
        (k0 - sqrt(k0**2 - 2.0_rk * beta * (q0 * (2.0_rk * L / pi)**2 &
        * cos(pi * x(i) * 0.5_rk / L) + c1 * x(i) + c2))) / beta
    enddo ! i = 1,nx
  endsubroutine temperature_exact_slab_cos_klin

  subroutine temperature_exact_cyl_lin_klin(nx, x, texact)
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: x(:) ! (nx)
    real(rk), intent(out) :: texact(:) ! (nx)

    integer(ik) :: i

    ! Based on the Kirchoff transform.
    ! TODO Try to be more clever about these
    real(rk), parameter :: q0 = 2e3 ! [W/cc]
    real(rk), parameter :: k0 = 67.2_rk  ! [W/cm/K]
    real(rk), parameter :: beta = 0.1_rk ! [W/cm/K^2]
    real(rk), parameter :: TR = 550.0_rk ! [K]
    real(rk), parameter :: R = 0.96_rk ! [cm]

    real(rk), parameter :: c2 = 0.5_rk * (k0**2 - (k0 - beta * TR)**2) / beta &
      + q0 * 5.0_rk / 36.0_rk * R**2

    do i = 1,nx
      texact(i) = &
        (k0 - sqrt(k0**2 - 2.0_rk * beta * (q0 * (x(i)**3/(9.0_rk * R) &
        - 0.25_rk * x(i)**2) + c2))) / beta
    enddo ! i = 1,nx
  endsubroutine temperature_exact_cyl_lin_klin

  subroutine analysis_output(fname, nx, xc, t, texact)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: xc(:) ! (nx)
    real(rk), intent(in) :: t(:) ! (nx)
    real(rk), intent(in) :: texact(:) ! (nx)

    integer(ik) :: i

    integer, parameter :: iout = 12

    call fileio_open_write(fname, iout)

    write(iout, '(a)') 'x [cm] , T [K] , Texact [K]'
    do i = 1,nx
      write(iout, '(es23.16, " , ", es23.16, " , ", es23.16)') &
        xc(i), t(i), texact(i)
    enddo ! i = 1,nx

    close(iout)
  endsubroutine analysis_output

endmodule analysis
