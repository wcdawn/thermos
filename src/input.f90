module input
implicit none

private

character(16) :: geometry = 'cartesian'

public :: geometry
public :: input_parse, input_summary

contains

  subroutine input_parse(fname)
    character(*), intent(in) :: fname
    integer, parameter :: iunit = 11
  endsubroutine input_parse

  subroutine input_summary()
  endsubroutine input_summary

endmodule input
