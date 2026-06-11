module fileio
implicit none

private

public :: fileio_open_read, fileio_open_write

contains

  subroutine fileio_open_read(fname, iounit, replace)
    character(*), intent(in) :: fname
    integer, intent(in) :: iounit
    logical, intent(in), optional :: replace
    integer :: ios
    logical :: repl
    repl = .true.
    if (present(replace)) then
      repl = replace
    endif
    if (repl) then
      open(file = fname, unit = iounit, status = 'replace', iostat = ios)
    else
      open(file = fname, unit = iounit, status = 'new', iostat = ios)
    endif
    if (ios /= 0) then
      call fileio_error(fname, iounit, ios)
    endif
  endsubroutine fileio_open_read

  subroutine fileio_open_write(fname, iounit)
    character(*), intent(in) :: fname
    integer, intent(in) :: iounit
  endsubroutine fileio_open_read

  subroutine fileio_error(fname, iounit, iostat)
    character(*), intent(in) :: fname
    integer, intent(in) :: iounit
    integer, intent(in) :: iostat
    character(1024) :: line
    write(*,*) 'ERROR: Failed to interact with file "' // trim(adjustl(fname)) // ' on unit
  endsubroutine fileio_error

endmodule fileio
