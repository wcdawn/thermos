program thermos
use kind, only : rk
use input, only : input_parse, input_summary, &
  length, nx
use geometry, only : geometry_calculate_coordinates
implicit none

character(1024) :: fname_input

real(rk), allocatable :: xcenter(:), dx(:)

if (command_argument_count() == 0) then
  stop 'missing input filename'
endif

write(*, '(a)') 'Begin Thermos'
write(*, *)

fname_input = ''
call get_command_argument(1, fname_input)

call input_parse(fname_input)
call input_summary()

allocate(xcenter(nx))
allocate(dx(nx))
call geometry_calculate_coordinates(length, nx, xcenter, dx)

deallocate(xcenter, dx)

write(*,'(a)') 'End Thermos'

endprogram thermos
