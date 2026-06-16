module finite_difference
use kind, only : rk, ik
implicit none

private

public :: finite_difference_solve_cartesian

contains

  subroutine finite_difference_build_matrix(nx, dx, sub, dia, sup)
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    real(rk), intent(out) :: sub(:) ! (nx-1)
    real(rk), intent(out) :: dia(:) ! (nx)
    real(rk), intent(out) :: sup(:) ! (nx-1)
  endsubroutine finite_difference_build_matrix

  subroutine finite_difference_build_source(nx, dx, src)
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    real(rk), intent(out) :: src(:) ! (nx)
  endsubroutine finite_difference_build_source

  subroutine finite_difference_solve_cartesian(nx, dx, &
      max_iter, tol_temperature, init_temperature, temperature)
    use linalg, only : trid, norm
    use output, only : output_write
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: max_iter
    real(rk), intent(in) :: tol_temperature ! [K]
    real(rk), intent(in) :: init_temperature ! [K]
    real(rk), intent(out) :: temperature(:) ! (nx)

    real(rk), allocatable :: sub(:) ! (nx-1)
    real(rk), allocatable :: dia(:) ! (nx)
    real(rk), allocatable :: sup(:) ! (nx-1)
    real(rk), allocatable :: q(:) ! (nx)
    real(rk), allocatable :: qcpy(:) ! (nx)
    real(rk), allocatable :: temperature_old(:) ! (nx)

    integer(ik) :: iter
    real(rk) :: conv

    character(1024) :: line

    allocate(sub(nx-1))
    allocate(dia(nx))
    allocate(sup(nx-1))
    allocate(q(nx), qcpy(nx))
    allocate(temperature_old(nx))

    call finite_difference_build_source(nx, dx, qcpy)

    temperature = init_temperature

    ! Picard iteration
    do iter = 1,max_iter

      temperature_old = temperature

      ! must rebuild matrix since thermal conductivity may change on each iteration
      ! must copy the source since it is used as scratch space by trid
      call finite_difference_build_matrix(nx, dx, sub, dia, sup)
      q = qcpy

      call trid(nx, sub, dia, sup, q, temperature)

      conv = norm(-1, temperature - temperature_old) ! max. abs. diff.

      write(line, '(a,i0,1x,a,es9.2)') 'iter=', iter, 'conv=', conv

      if (conv < tol_temperature) then
        call output_write('CONVERGENCE!!!')
        exit
      endif

    enddo ! iter = 1,max_iter

    deallocate(temperature_old)
    deallocate(q, qcpy)
    deallocate(sub, dia, sup)
  endsubroutine finite_difference_solve_cartesian

endmodule finite_difference
