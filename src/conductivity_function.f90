module conductivity_function
use kind, only : ik, rk
implicit none

private

public :: conductivity_function_init, conductivity_function_cleanup, &
  conductivity_output_csv

abstract interface
  pure function k(T)
    use kind, only : rk
    real(rk) :: k ! [W/cm/K]
    real(rk), intent(in) :: T ! [K]
  endfunction k
endinterface

procedure(k), pointer, public :: conductivity_fun => null()

real(rk), allocatable :: coeff(:)

contains

  subroutine conductivity_function_init(conductivity_name, coeff_in)
    use exception_handler, only : exception_fatal
    character(*), intent(in) :: conductivity_name
    real(rk), intent(in) :: coeff_in(:)

    allocate(coeff(size(coeff_in)))
    coeff = coeff_in

    select case (conductivity_name)
      case ('constant')
        conductivity_fun => conductivity_fun_const
      case ('linear')
        conductivity_fun => conductivity_fun_linear
      case ('rational')
        conductivity_fun => conductivity_fun_rational
      case default
        call exception_fatal('Unknown name of conductivity function: ' &
          // trim(adjustl(conductivity_name)))
    endselect
  endsubroutine conductivity_function_init

  pure real(rk) function conductivity_fun_const(T)
    real(rk), intent(in) :: T ! [K] temperature
    real(rk) :: k0
    k0 = coeff(1)
    conductivity_fun_const = k0
  endfunction conductivity_fun_const

  pure real(rk) function conductivity_fun_linear(T)
    real(rk), intent(in) :: T ! [K] temperature
    real(rk) :: k0, beta
    k0 = coeff(1)
    beta = coeff(2)
    conductivity_fun_linear = k0 - beta * T
  endfunction conductivity_fun_linear

  pure real(rk) function conductivity_fun_rational(T)
    real(rk), intent(in) :: T ! [K] temperature
    real(rk) :: alpha, beta
    alpha = coeff(1)
    beta = coeff(2)
    conductivity_fun_rational = 1.0_rk / (alpha + beta * T)
  endfunction conductivity_fun_rational

  subroutine conductivity_function_cleanup()
    deallocate(coeff)
  endsubroutine conductivity_function_cleanup

  subroutine conductivity_output_csv(fname, nx, xcenter, temperature)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: xcenter(:) ! (nx)
    real(rk), intent(in) :: temperature(:) ! (nx)

    integer, parameter :: iout = 11

    integer(ik) :: i

    call fileio_open_write(fname, iout)
    write(iout, '(a)') 'xcenter [cm] , conductivity [W/cm/K]'
    do i = 1,nx
      write(iout, '(es23.16," , ",es23.6)') &
        xcenter(i), conductivity_fun(temperature(i))
    enddo ! i = 1,nx
    close(iout)
  endsubroutine conductivity_output_csv

endmodule conductivity_function
