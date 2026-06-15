module output
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit
implicit none

private

public :: output_open_file, output_close_file, output_write

integer, parameter, private :: output_file_unit = 99
integer, parameter, private :: output_list(2) = [ stdout, output_file_unit ]

contains

  subroutine output_open_file(fname)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    call fileio_open_write(fname, output_file_unit)
  endsubroutine output_open_file

  subroutine output_close_file()
    close(output_file_unit)
  endsubroutine output_close_file

  subroutine output_write(str)
    character(*), intent(in) :: str
    integer :: i
    do i = 1,size(output_list)
      write(output_list(i), '(a)') trim(adjustl(str))
    enddo ! i = 1,size(output_list)
  endsubroutine output_write

endmodule output
