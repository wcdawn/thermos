module analysis
use kind, only : rk, ik
implicit none

private

public :: analysis_analyze

contains

  subroutine analysis_analyze(analysis_name, nx, temperature)
    use output, only : output_write
    character(*), intent(in) :: analysis_name
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: temperature(:) ! (nx)

    select case (analysis_name)
      case ('slab_cos')
        ! do something
      case default
        call output_write('ERROR: unknown analysis name: ' // &
          trim(adjustl(analysis_name)))
        stop
    endselect

    call output_write('=== ANALYSIS SUMMARY ===')
    call output_write('')
  endsubroutine analysis_analyze

endmodule analysis
