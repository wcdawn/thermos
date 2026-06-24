module output
use kind, only : rk, ik
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit
implicit none

private

public :: output_open_file, output_close_file, output_write, &
  output_temperature_csv, output_conductivity_csv, output_source_csv

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

  subroutine output_temperature_csv(fname, nx, xcenter, temperature)
    use fileio, only : fileio_open_write
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: xcenter(:) ! (nx)
    real(rk), intent(in) :: temperature(:) ! (nx)

    integer, parameter :: iout = 11

    integer(ik) :: i

    call fileio_open_write(fname, iout)
    write(iout, '(a)') 'xcenter [cm] , temperature [K]'
    do i = 1,nx
      write(iout, '(es23.16," , ",es23.6)') xcenter(i), temperature(i)
    enddo ! i = 1,nx
    close(iout)
  endsubroutine output_temperature_csv

  subroutine output_conductivity_csv(fname, nx, xcenter, temperature)
    use fileio, only : fileio_open_write
    use conductivity_function, only : conductivity_fun
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
  endsubroutine output_conductivity_csv

  subroutine output_source_csv(fname, nx, xcenter)
    use fileio, only : fileio_open_write
    use source_function, only : source_fun
    character(*), intent(in) :: fname
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: xcenter(:) ! (nx)

    integer, parameter :: iout = 11

    integer(ik) :: i

    call fileio_open_write(fname, iout)
    write(iout, '(a)') 'xcenter [cm] , source [W/cm^3]'
    do i = 1,nx
      write(iout, '(es23.16," , ",es23.6)') &
        xcenter(i), source_fun(xcenter(i))
    enddo ! i = 1,nx
    close(iout)
  endsubroutine output_source_csv

endmodule output
