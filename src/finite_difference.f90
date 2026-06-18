module finite_difference
use kind, only : rk, ik
implicit none

private

public :: finite_difference_solve

contains

  subroutine finite_difference_build_matrix(nx, dx, temperature, &
      bctype_left, bctype_right, sub, dia, sup)
    use conductivity_function, only : conductivity_fun
    use output, only : output_write
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    real(rk), intent(in) :: temperature(:) ! (nx)
    character(*), intent(in) :: bctype_left, bctype_right
    real(rk), intent(out) :: sub(:) ! (nx-1)
    real(rk), intent(out) :: dia(:) ! (nx)
    real(rk), intent(out) :: sup(:) ! (nx-1)

    integer(ik) :: i
    real(rk) :: kprev, kthis, knext

    ! BC at x=0, i=1
    select case (bctype_left)
      case ('fixed')
        kthis = conductivity_fun(temperature(1))
        knext = conductivity_fun(temperature(2))
        dia(1) = -2.0_rk * ((kthis/dx(1)) * (knext/dx(2)) &
          / (kthis/dx(1) + knext/dx(2)) &
          + kthis/dx(1))
        sup(1) = 2.0_rk * (kthis/dx(1)) * (knext/dx(2)) &
          / (kthis/dx(1) + knext/dx(2))
      case ('insulated')
        kthis = conductivity_fun(temperature(1))
        knext = conductivity_fun(temperature(2))
        dia(1) = -2.0_rk * (kthis/dx(1)) * (knext/dx(2)) &
          / (kthis/dx(1) + knext/dx(2))
        sup(1) = 2.0_rk * (kthis/dx(1)) * (knext/dx(2)) &
          / (kthis/dx(1) + knext/dx(2))
      case default
        call output_write('ERROR: unknown value of bctype_left in build_matrix: ' &
          // trim(adjustl(bctype_left)))
        stop
    endselect

    do i = 2,nx-1

      kprev = conductivity_fun(temperature(i-1))
      kthis = conductivity_fun(temperature(i))
      knext = conductivity_fun(temperature(i+1))

      sub(i-1) = 2.0_rk * (kthis/dx(i)) * (kprev/dx(i-1)) &
        / (kthis/dx(i) + kprev/dx(i-1))
      dia(i) = -2.0_rk* (kthis/dx(i) * (kprev/dx(i-1)) &
        / (kthis/dx(i) + kprev/dx(i-1)) &
        + kthis/dx(i) * knext/dx(i+1) & 
        / (kthis/dx(i) + knext/dx(i+1)))
      sup(i) = 2.0_rk * (kthis/dx(i)) * (knext/dx(i+1)) &
        / (kthis/dx(i) + knext/dx(i+1))

    enddo ! i = 2,nx-1

    select case (bctype_right)
      case ('fixed')
        kprev = conductivity_fun(temperature(nx-1))
        kthis = conductivity_fun(temperature(nx))
        sub(nx-1) = 2.0_rk * (kthis/dx(nx)) * (kprev/dx(nx-1)) &
          / (kthis/dx(nx) + kprev/dx(nx-1))
        dia(nx) = -2.0_rk * ((kthis/dx(nx)) * (kprev/dx(nx-1)) &
          / (kthis/dx(nx) + kprev/dx(nx-1)) &
          +  kthis/dx(1))
      case ('insulated')
        kprev = conductivity_fun(temperature(nx-1))
        kthis = conductivity_fun(temperature(nx))
        sub(nx-1) = 2.0_rk * (kthis/dx(nx)) * (kprev/dx(nx-1)) &
          / (kthis/dx(nx) + kprev/dx(nx-1))
        dia(nx) = -2.0_rk * (kthis/dx(nx)) * (kprev/dx(nx-1)) &
          / (kthis/dx(nx) + kprev/dx(nx-1))
      case default
        call output_write('ERROR: unknown value of bctype_right in build_matrix: ' &
          // trim(adjustl(bctype_right)))
        stop
    endselect

  endsubroutine finite_difference_build_matrix

  subroutine finite_difference_build_source(nx, xcenter, dx, &
      bctype_left, bctype_right, Tleft, Tright, src)
    use source_function, only : source_fun
    use conductivity_function, only : conductivity_fun
    use output, only : output_write
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: xcenter(:) ! (nx)
    real(rk), intent(in) :: dx(:) ! (nx)
    character(*), intent(in) :: bctype_left, bctype_right
    real(rk), intent(in) :: Tleft, Tright
    real(rk), intent(out) :: src(:) ! (nx)
    
    integer(ik) :: i

    select case(bctype_left)
      case ('fixed')
        src(1) = dx(1) * source_fun(xcenter(1)) &
          + 2.0_rk * conductivity_fun(Tleft) / dx(1) * Tleft
      case ('insulated')
        src(1) = dx(1) * source_fun(xcenter(1))
      case default
        call output_write('ERROR: unknown value of bctype_left in build_source: ' &
          // trim(adjustl(bctype_left)))
        stop
    endselect

    do i = 2,nx-1
      src(i) = dx(i) * source_fun(xcenter(i))
    enddo ! i = 2,nx-1

    select case(bctype_right)
      case ('fixed')
        src(nx) = dx(nx) * source_fun(xcenter(nx)) &
          + 2.0_rk * conductivity_fun(Tright) / dx(nx) * Tright
      case ('insulated')
        src(nx) = dx(nx) * source_fun(xcenter(nx))
      case default
        call output_write('ERROR: unknown value of bctype_right in build_source: ' &
          // trim(adjustl(bctype_right)))
        stop
    endselect
    src = -src

  endsubroutine finite_difference_build_source

  subroutine finite_difference_solve(nx, xcenter, dx, &
      bctype_left, bctype_right, bcval_left, bcval_right, &
      max_iter, tol_temperature, init_temperature, temperature)
    use linalg, only : trid, norm
    use output, only : output_write
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: xcenter(:) ! (nx)
    real(rk), intent(in) :: dx(:) ! (nx)

    character(*), intent(in) :: bctype_left, bctype_right
    real(rk), intent(in) :: bcval_left, bcval_right

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

    call finite_difference_build_source(nx, xcenter, dx, &
      bctype_left, bctype_right, bcval_left, bcval_right, qcpy)

    temperature = init_temperature

    ! Picard iteration
    do iter = 1,max_iter

      temperature_old = temperature

      ! must rebuild matrix since thermal conductivity may change on each iteration
      ! must copy the source since it is used as scratch space by trid
      call finite_difference_build_matrix(nx, dx, temperature, &
        bctype_left, bctype_right, sub, dia, sup)
      q = qcpy

      call trid(nx, sub, dia, sup, q, temperature)

      conv = norm(-1, temperature - temperature_old) ! max. abs. diff.

      write(line, '(a,i3,1x,a,es9.2)') 'iter=', iter, 'conv=', conv
      call output_write(line)

      if (conv < tol_temperature) then
        call output_write('CONVERGENCE!!!')
        exit
      endif

    enddo ! iter = 1,max_iter

    call output_write('')

    deallocate(temperature_old)
    deallocate(q, qcpy)
    deallocate(sub, dia, sup)
  endsubroutine finite_difference_solve

endmodule finite_difference
