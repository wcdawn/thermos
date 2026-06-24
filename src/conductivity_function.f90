module conductivity_function
use kind, only : rk
implicit none

private

public :: conductivity_function_init, conductivity_function_cleanup

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
    character(*), intent(in) :: conductivity_name
    real(rk), intent(in) :: coeff_in(:)

    allocate(coeff(size(coeff_in)))
    coeff = coeff_in

    select case (conductivity_name)
      case ('constant')
        conductivity_fun => conductivity_fun_const
      case ('linear')
        conductivity_fun => conductivity_fun_linear
      case default
        write(*,*) 'ERROR: Unknown name of conductivity function: ' &
          // trim(adjustl(conductivity_name))
        stop
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

  subroutine conductivity_function_cleanup()
    deallocate(coeff)
  endsubroutine conductivity_function_cleanup

endmodule conductivity_function
