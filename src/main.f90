program thermos
use kind, only : rk, ik
use input, only : input_parse, input_summary, &
  length, nx, refine
use geometry, only : geometry_calculate_coordinates, geometry_refine, geometry_summary
use output, only : output_open_file, output_close_file, output_write
implicit none

character(1024) :: fname_input, fname_stub, fname_out

integer(ik) :: i
real(rk), allocatable :: xcenter(:), dx(:)

if (command_argument_count() == 0) then
  stop 'missing input filename'
endif

fname_input = ''
call get_command_argument(1, fname_input)

! find last dot
i = index(trim(adjustl(fname_input)), '.', back=.true.) - 1
fname_stub = fname_input(:i)
fname_out = trim(adjustl(fname_stub)) // '.out'

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

call output_write('End Thermos')

call output_close_file()
deallocate(xcenter, dx)

endprogram thermos
