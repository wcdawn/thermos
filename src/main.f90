program thermos
use kind, only : rk, ik
use input, only : input_parse, input_summary, &
  length, nx, refine
use geometry, only : geometry_calculate_coordinates, geometry_refine, geometry_summary
implicit none

character(1024) :: fname_input

integer(ik) :: i
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

write(*, '(a)') '(before refinement)'
call geometry_summary(nx, dx)

do i = 1,refine
  call geometry_refine(nx, xcenter, dx)
enddo ! i = 1,refine

write(*, '(a)') '(after refinement)'
call geometry_summary(nx, dx)

deallocate(xcenter, dx)

write(*,'(a)') 'End Thermos'

endprogram thermos
