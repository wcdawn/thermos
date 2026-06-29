module input
use kind, only : rk, ik
implicit none

private

real(rk) :: length ! Total length of the calculation domain [cm]. Either outer radius or slab length.
integer(ik) :: nx ! Number of spatial cells.
integer(ik) :: refine = 0 ! number of times to uniformly perform spatial refinement
character(16) :: geometry = 'cartesian' ! ('cartesian', 'cylindrical', 'spherical')
character(64) :: solver = 'finite_difference' ! ('finite_difference', 'finite_element')
integer(ik) :: max_iter = 10 ! Maximum number of Picard iterations
real(rk) :: tol_temperature = 0.5_rk ! [K] Tolerance for Picard iterations (max.  abs. diff.)
real(rk) :: init_temperature = 300.0_rk ! [K] Initial estimate for temperature (uniform)

! NOTE: For cylindrical and spherical geometry, bctype_left will always be
! treated as 'insulated' regardless of user input.
character(16) :: bctype_left = 'fixed' ! ('fixed', 'insulated')
character(16) :: bctype_right= 'fixed' ! ('fixed', 'insulated')
real(rk) :: bcval_left = 600.0_rk
real(rk) :: bcval_right = 300.0_rk

character(16) :: source_function_name = 'cos' ! ('cos', 'sin', 'linear')
real(rk) :: source_coeff(4) = [ 20.0, 10.0, 0.0, 0.0 ]
character(16) :: conductivity_function_name = 'constant' ! ('constant', 'linear', 'rational')
real(rk) :: conductivity_coeff(4) = [ 1.2, 0.0, 0.0, 0.0 ]

character(16) :: mesh_spacing = 'uniform' ! ('uniform', 'area', 'volume')

character(16) :: analysis_name = 'none'

! variables
public :: geometry, length, nx, refine, mesh_spacing
public :: solver, max_iter, tol_temperature, init_temperature
public :: bctype_left, bctype_right, bcval_left, bcval_right
public :: conductivity_function_name, conductivity_coeff
public :: source_function_name, source_coeff
public :: analysis_name

! subroutines
public :: input_parse, input_summary

character(*), parameter :: comment_char = '#'

contains

  subroutine input_parse(fname)
    use fileio, only : fileio_open_read
    use exception_handler, only : exception_fatal
    character(*), intent(in) :: fname
    integer, parameter :: iunit = 11
    integer :: ios
    character(1024) :: line, card, msg

    call fileio_open_read(fname, iunit)

    do

      read(iunit, '(a)', iostat=ios) line
      if (ios /= 0) then
        exit
      endif

      line = adjustl(line)

      if (line == '') then
        cycle
      endif

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
        case ('nx')
          read(iunit, *) card, nx
        case ('refine')
          read(iunit, *) card, refine

        case ('max_iter')
          read(iunit, *) card, max_iter
        case ('tol_temperature')
          read(iunit, *) card, tol_temperature
        case ('init_temperature')
          read(iunit, *) card, init_temperature

        case ('bctype_left')
          read(iunit, *) card, bctype_left
        case ('bctype_right')
          read(iunit, *) card, bctype_right
        case ('bcval_left')
          read(iunit, *) card, bcval_left
        case ('bcval_right')
          read(iunit, *) card, bcval_right
        case ('source_function')
          read(iunit, *) card, source_function_name
        case ('source_coeff')
          read(iunit, *) card, source_coeff ! NOTE: fortran array input
        case ('conductivity_function')
          read(iunit, *) card, conductivity_function_name
        case ('conductivity_coeff')
          read(iunit, *) card, conductivity_coeff ! NOTE: fortran array input
        case ('mesh_spacing')
          read(iunit, *) card, mesh_spacing
        case ('analysis_name')
          read(iunit, *) card, analysis_name
        case default
          write(msg, '(a)') 'Unknown input card. ' &
            // 'Troublesome line follows.' // new_line('a') // trim(adjustl(line))
          call exception_fatal(msg)
      endselect
    enddo

    if ((geometry == 'spherical') .or. (geometry == 'cylindrical')) then
      bctype_left = 'insulated'
    endif

    close(iunit)
  endsubroutine input_parse

  subroutine input_summary()
    use output, only : output_write
    character(1024) :: line
    integer :: i
    call output_write('=== INPUT SUMMARY ===')
    write(line, '(a,a)') 'Geometry: ', trim(adjustl(geometry))
    call output_write(line)
    write(line, '(a,es13.6)') 'Length: ', length
    call output_write(line)
    write(line, '(a,i0)') 'NX: ', nx
    call output_write(line)
    write(line, '(a,a)') 'Mesh spacing: ', trim(adjustl(mesh_spacing))
    call output_write(line)
    write(line, '(a,i0)') 'Uniform Refinement: ', refine
    call output_write(line)
    write(line, '(a,a)') 'Solver: ', trim(adjustl(solver))
    call output_write(line)
    write(line, '(a,i0)') 'Maximum Picard iterations: ', max_iter
    call output_write(line)
    write(line, '(a,es9.2)') 'Temperature Tolerance [K]: ', tol_temperature
    call output_write(line)
    write(line, '(a,es9.2)') 'Initial Temperature [K]: ', init_temperature
    call output_write(line)
    call output_write('Boundary conditions.' &
      // ' Left: ' // trim(adjustl(bctype_left)) &
      // ' Right: ' // trim(adjustl(bctype_right)))
    write(line, '(a,a,es9.2,a,es9.2)') 'Boundary values.', &
      ' Left: ', bcval_left, ' Right: ',bcval_right
    call output_write(line)

    write(line, '(a,a)') 'Source function q(x): ', &
      trim(adjustl(source_function_name))
    call output_write(line)
    write(line, '(a,es9.2)') 'Source coefficients: [ ', source_coeff(1)
    do i = 2,size(source_coeff)
      write(line, '(a,es9.2)') trim(adjustl(line)) // ' , ', source_coeff(i)
    enddo ! i = 2,size(source_coeff)
    line = trim(adjustl(line)) // ' ]'
    call output_write(line)

    write(line, '(a,a)') 'Conductivity function k(T): ', &
      trim(adjustl(conductivity_function_name))
    call output_write(line)
    write(line, '(a,es9.2)') 'Conductivity coefficients: [ ', conductivity_coeff(1)
    do i = 2,size(conductivity_coeff)
      write(line, '(a,es9.2)') trim(adjustl(line)) // ' , ', conductivity_coeff(i)
    enddo ! i = 2,size(conductivity_coeff)
    line = trim(adjustl(line)) // ' ]'
    call output_write(line)

    write(line, '(a,a)') 'Analytic comparison name: ', trim(adjustl(analysis_name))
    call output_write(line)

    call output_write('')
  endsubroutine input_summary

endmodule input
