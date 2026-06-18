module source_function
use kind, only : rk
implicit none

private

public :: source_function_init, source_function_cleanup

abstract interface
  pure function q(x)
    use kind, only : rk
    real(rk) :: q ! q(x)
    real(rk), intent(in) :: x ! [cm]
  endfunction q
endinterface

procedure(q), pointer, public :: source_fun => null()

real(rk), allocatable :: coeff(:)

contains

  subroutine source_function_init(source_name, coeff_in)
    use output, only : output_write
    character(*), intent(in) :: source_name
    real(rk), intent(in) :: coeff_in(:)

    allocate(coeff(size(coeff_in)))
    coeff = coeff_in

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
    real(rk) :: q0
    real(rk) :: L
    q0 = coeff(1)
    L  = coeff(2)
    source_fun_cos = q0 * cos(pi * x * 0.5_rk / L)
  endfunction source_fun_cos

  subroutine source_function_cleanup()
    deallocate(coeff)
  endsubroutine source_function_cleanup

endmodule source_function
