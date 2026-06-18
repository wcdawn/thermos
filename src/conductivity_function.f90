module conductivity_function
use kind, only : rk
implicit none

private

public :: conductivity_function_init

abstract interface
  pure function k(T)
    use kind, only : rk
    real(rk) :: k ! [W/cm/K]
    real(rk), intent(in) :: T ! [K]
  endfunction k
endinterface

procedure(k), pointer, public :: conductivity_fun => null()

contains

  subroutine conductivity_function_init(conductivity_name)
    use output, only : output_write
    character(*), intent(in) :: conductivity_name
    select case (conductivity_name)
      case ('constant')
        conductivity_fun => conductivity_fun_const
      case default
        call output_write('ERROR: Unknown name of conductivity function: ' &
          // trim(adjustl(conductivity_name)))
        stop
    endselect
  endsubroutine conductivity_function_init

  pure real(rk) function conductivity_fun_const(T)
    real(rk), intent(in) :: T ! [K] temperature
    real(rk), parameter :: k0 = 1.2_rk
    conductivity_fun_const = k0
  endfunction conductivity_fun_const

endmodule conductivity_function
