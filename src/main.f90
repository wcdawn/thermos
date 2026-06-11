program thermos
use input, only : input_parse, input_summary
implicit none

character(1024) :: fname_input

if (command_argument_count() == 0) then
  stop 'missing input filename'
endif

write(*, '(a)') 'Begin Thermos'
write(*, *)

fname_input = ''
call get_command_argument(1, fname_input)

call input_parse(fname_input)
call input_summary()

write(*,'(a)') 'End Thermos'

endprogram thermos
