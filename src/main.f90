program thermos
use kind, only : rk, ik
use input, only : input_parse, input_summary, &
  geometry, length, nx, refine, &
  bctype_left, bctype_right, bcval_left, bcval_right, &
  solver, max_iter, tol_temperature, init_temperature
use geometry, only : geometry_calculate_coordinates, geometry_refine, geometry_summary
use output, only : output_open_file, output_close_file, output_write, &
  output_temperature_csv
use finite_difference, only : finite_difference_solve_cartesian
use source_function, only : source_function_init
use conductivity_function, only : conductivity_function_init
implicit none

character(1024) :: fname_input, fname_stub, fname_out, fname_temperature

integer(ik) :: i
real(rk), allocatable :: xcenter(:), dx(:)
real(rk), allocatable :: temperature(:)

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

call output_open_file(fname_out)

call output_write('Begin Thermos')
call output_write('Input file: ' // trim(adjustl(fname_input)))
call output_write('')

call input_parse(fname_input)
call input_summary()

allocate(xcenter(nx))
allocate(dx(nx))
call geometry_calculate_coordinates(length, nx, xcenter, dx)

call output_write('(before refinement)')
call geometry_summary(nx, dx)

do i = 1,refine
  call geometry_refine(nx, xcenter, dx)
enddo ! i = 1,refine

call output_write('(after refinement)')
call geometry_summary(nx, dx)

call source_function_init('cos')
call conductivity_function_init('constant')

allocate(temperature(nx))

select case (trim(adjustl(solver)) // '_' // trim(adjustl(geometry)))
  case ('finite_difference_cartesian')
    call finite_difference_solve_cartesian(nx, xcenter, dx, &
      bctype_left, bctype_right, bcval_left, bcval_right, &
      max_iter, tol_temperature, init_temperature, temperature)
  case default
    call output_write('ERROR: unknown solver selection: ' // trim(adjustl(solver)))
    stop
endselect

call output_write('Writing temperautre output on: ' // &
  trim(adjustl(fname_temperature)))
call output_temperature_csv(fname_temperature, nx, xcenter, temperature)
call output_write('')

call output_write('End Thermos')

call output_close_file()
deallocate(xcenter, dx)
deallocate(temperature)

endprogram thermos
