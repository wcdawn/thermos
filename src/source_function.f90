module source_function
use kind, only : rk
implicit none

private

public :: source_function_init

abstract interface
  pure function q(x)
    use kind, only : rk
    real(rk) :: q ! q(x)
    real(rk), intent(in) :: x ! [cm]
  endfunction q
endinterface

procedure(q), pointer, public :: source_fun => null()

contains

  subroutine source_function_init(source_name)
    use output, only : output_write
    character(*), intent(in) :: source_name
    select case (source_name)
      case ('cos')
        source_fun => source_fun_cos
      case default
        call output_write('ERROR: Unknown name of source function: ' &
          // trim(adjustl(source_name)))
        stop
    endselect
  endsubroutine source_function_init

  pure real(rk) function source_fun_cos(x)
    use constants, only : pi
    real(rk), intent(in) :: x ! [cm] position
    real(rk), parameter :: q0 = 20.0_rk
    real(rk), parameter :: L = 10.0_rk
    source_fun_cos = q0 * cos(pi * x / L)
  endfunction source_fun_cos

endmodule source_function
