module input
use kind, only : rk, ik
implicit none

private

real(rk) :: length ! Total length of the calculation domain [cm]. Either outer radius or slab length.
integer(ik) :: nx ! Number of spatial cells.
character(16) :: geometry = 'cartesian' ! ('cartesian', 'cylindrical', 'spherical')

public :: geometry, length, nx
public :: input_parse, input_summary

character(*), parameter :: comment_char = '#'

contains

  subroutine input_parse(fname)
    use fileio, only : fileio_open_read
    character(*), intent(in) :: fname
    integer, parameter :: iunit = 11
    integer :: ios
    character(1024) :: line, card

    call fileio_open_read(fname, iunit)

    do

      read(iunit, '(a)', iostat=ios) line
      if (ios /= 0) then
        exit
      endif

      line = adjustl(line)
      if (line(1:1) == comment_char) then
        cycle
      endif

      read(line, *) card
      backspace(iunit)

      select case (card)
        case ('geometry')
          read(iunit, *) card, geometry
        case ('length')
          read(iunit, *) card, length
        case default
          write(*, '(a,a)') 'ERROR: Unknown input card. Troublesome line follows.'
          write(*, '(a)') trim(adjustl(line))
          stop
      endselect
    enddo

    close(iunit)
  endsubroutine input_parse

  subroutine input_summary()
    write(*, '(a)') '=== INPUT SUMMARY ==='
    write(*, '(a,a)') 'Geometry: ', trim(adjustl(geometry))
    write(*, *)
  endsubroutine input_summary

endmodule input
