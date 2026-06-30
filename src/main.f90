program thermos
use kind, only : rk, ik
use input, only : input_parse, input_summary, &
  geometry, length, nx, refine, mesh_spacing, &
  bctype_left, bctype_right, bcval_left, bcval_right, &
  solver, max_iter, tol_temperature, init_temperature, &
  conductivity_function_name, conductivity_coeff, &
  source_function_name, source_coeff, &
  analysis_name
use geometry, only : geometry_calculate_coordinates, geometry_refine, geometry_summary
use output, only : output_open_file, output_close_file, output_write, &
  output_temperature_csv
use finite_difference, only : finite_difference_solve
use finite_element, only : finite_element_solve
use source_function, only : source_function_init, source_function_cleanup, &
  source_output_csv
use conductivity_function, only : conductivity_function_init, conductivity_function_cleanup, &
  conductivity_output_csv
use analysis, only : analysis_analyze
use exception_handler, only : exception_fatal, exception_summary
implicit none

character(1024) :: fname_input, fname_stub, fname_out, &
  fname_temperature, fname_analysis, fname_conductivity, fname_source
character(1024) :: line

integer(ik) :: i
real(rk), allocatable :: xcenter(:), dx(:)
real(rk), allocatable :: temperature(:)
real(rk) :: TCL

if (command_argument_count() == 0) then
  stop 'missing input filename'
endif

fname_input = ''
call get_command_argument(1, fname_input)

! find last dot
i = index(trim(adjustl(fname_input)), '.', back=.true.) - 1
fname_stub = fname_input(:i)
fname_out = trim(adjustl(fname_stub)) // '.out'
fname_temperature = trim(adjustl(fname_stub)) // '_temperature.csv'
fname_analysis = trim(adjustl(fname_stub)) // '_analysis.csv'
fname_conductivity = trim(adjustl(fname_stub)) // '_conductivity.csv'
fname_source = trim(adjustl(fname_stub)) // '_source.csv'

call output_open_file(fname_out)

call output_write('Begin Thermos')
call output_write('Input file: ' // trim(adjustl(fname_input)))
call output_write('')

call input_parse(fname_input)
call input_summary()

allocate(xcenter(nx))
allocate(dx(nx))
call geometry_calculate_coordinates(mesh_spacing, length, nx, xcenter, dx)

call output_write('(before refinement)')
call geometry_summary(nx, dx)

do i = 1,refine
  call geometry_refine(nx, xcenter, dx)
enddo ! i = 1,refine

call output_write('(after refinement)')
call geometry_summary(nx, dx)

call source_function_init(source_function_name, source_coeff)
call conductivity_function_init(conductivity_function_name, conductivity_coeff)

allocate(temperature(nx))

select case (solver)
  case ('finite_difference')
    call finite_difference_solve(geometry, nx, xcenter, dx, &
      bctype_left, bctype_right, bcval_left, bcval_right, &
      max_iter, tol_temperature, init_temperature, temperature, TCL)
  case ('finite_element')
    call finite_element_solve(geometry, nx, xcenter, dx, &
      bctype_left, bctype_right, bcval_left, bcval_right, &
      max_iter, tol_temperature, init_temperature, temperature, TCL)
  case default
    call exception_fatal('unknown solver selection: ' // trim(adjustl(solver)))
endselect

if (TCL > 0.0_rk) then
  write(line, '(a,es13.6,a)') 'Tcenterline = ', TCL, ' [K]'
  call output_write(line)
  call output_write('')
endif

if (analysis_name /= 'none') then
  call analysis_analyze(analysis_name, fname_analysis, nx, xcenter, temperature, TCL)
endif

call output_write('Writing temperautre output on: ' // &
  trim(adjustl(fname_temperature)))
call output_temperature_csv(fname_temperature, nx, xcenter, temperature)

call output_write('Writing thermal conductivity output on: ' // &
  trim(adjustl(fname_conductivity)))
call conductivity_output_csv(fname_conductivity, nx, xcenter, temperature)

call output_write('Writing source output on: ' // &
  trim(adjustl(fname_source)))
call source_output_csv(fname_source, nx, xcenter)

call output_write('')

call exception_summary()

call output_write('Normal Termination :)')
call output_write('End Thermos')

call source_function_cleanup()
call conductivity_function_cleanup()
call output_close_file()
deallocate(xcenter, dx)
deallocate(temperature)

endprogram thermos
