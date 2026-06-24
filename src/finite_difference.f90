module finite_difference
use kind, only : rk, ik
implicit none

private

public :: finite_difference_solve

contains

  subroutine finite_difference_build_matrix_cartesian(nx, xcenter, dx, temperature, &
      bctype_left, bctype_right, sub, dia, sup)
    use conductivity_function, only : conductivity_fun
    use output, only : output_write
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: xcenter(:) ! (nx)
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
        call output_write('ERROR: unknown value of bctype_left in build_matrix_cartesian: ' &
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
          +  kthis/dx(nx))
      case ('insulated')
        kprev = conductivity_fun(temperature(nx-1))
        kthis = conductivity_fun(temperature(nx))
        sub(nx-1) = 2.0_rk * (kthis/dx(nx)) * (kprev/dx(nx-1)) &
          / (kthis/dx(nx) + kprev/dx(nx-1))
        dia(nx) = -2.0_rk * (kthis/dx(nx)) * (kprev/dx(nx-1)) &
          / (kthis/dx(nx) + kprev/dx(nx-1))
      case default
        call output_write('ERROR: unknown value of bctype_right in build_matrix_cartesian: ' &
          // trim(adjustl(bctype_right)))
        stop
    endselect
  endsubroutine finite_difference_build_matrix_cartesian

  subroutine finite_difference_build_matrix_cylindrical(nx, xcenter, dx, temperature, &
      bctype_left, bctype_right, sub, dia, sup)
    use conductivity_function, only : conductivity_fun
    use output, only : output_write
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: xcenter(:) ! (nx)
    real(rk), intent(in) :: dx(:) ! (nx)
    real(rk), intent(in) :: temperature(:) ! (nx)
    character(*), intent(in) :: bctype_left, bctype_right
    real(rk), intent(out) :: sub(:) ! (nx-1)
    real(rk), intent(out) :: dia(:) ! (nx)
    real(rk), intent(out) :: sup(:) ! (nx-1)

    integer(ik) :: i
    real(rk) :: kprev, kthis, knext
    real(rk) :: xprev, xnext

    ! BC at x=0, i=1
    select case (bctype_left)
      case ('fixed')
        kthis = conductivity_fun(temperature(1))
        knext = conductivity_fun(temperature(2))
        xnext = xcenter(1) + 0.5_rk * dx(1)
        xprev = xcenter(1) - 0.5_rk * dx(1)
        dia(1) = -2.0_rk * xnext * ((kthis/dx(1)) * (knext/dx(2)) &
          / (kthis/dx(1) + knext/dx(2))) &
          -2.0_rk * xprev *  kthis/dx(1)
        sup(1) = 2.0_rk * xnext * (kthis/dx(1)) * (knext/dx(2)) &
          / (kthis/dx(1) + xnext*knext/dx(2))
      case ('insulated')
        kthis = conductivity_fun(temperature(1))
        knext = conductivity_fun(temperature(2))
        xnext = xcenter(1) + 0.5_rk * dx(1)
        dia(1) = -2.0_rk * xnext * (kthis/dx(1)) * (knext/dx(2)) &
          / (kthis/dx(1) + knext/dx(2))
        sup(1) = 2.0_rk * xnext * (kthis/dx(1)) * (knext/dx(2)) &
          / (kthis/dx(1) + knext/dx(2))
      case default
        call output_write('ERROR: unknown value of bctype_left in build_matrix_cylindrical: ' &
          // trim(adjustl(bctype_left)))
        stop
    endselect

    do i = 2,nx-1

      kprev = conductivity_fun(temperature(i-1))
      kthis = conductivity_fun(temperature(i))
      knext = conductivity_fun(temperature(i+1))

      xprev = xcenter(i) - 0.5_rk * dx(i)
      xnext = xcenter(i) + 0.5_rk * dx(i)

      sub(i-1) = 2.0_rk * xprev * (kthis/dx(i)) * (kprev/dx(i-1)) &
        / (kthis/dx(i) + kprev/dx(i-1))
      dia(i) = -2.0_rk* (xprev * kthis/dx(i) * kprev/dx(i-1) &
        / (kthis/dx(i) + kprev/dx(i-1)) &
        + xnext * kthis/dx(i) * knext/dx(i+1) & 
        / (kthis/dx(i) + knext/dx(i+1)))
      sup(i) = 2.0_rk * xnext * (kthis/dx(i)) * (knext/dx(i+1)) &
        / (kthis/dx(i) + knext/dx(i+1))

    enddo ! i = 2,nx-1

    select case (bctype_right)
      case ('fixed')
        kprev = conductivity_fun(temperature(nx-1))
        kthis = conductivity_fun(temperature(nx))
        xprev = xcenter(nx) - 0.5_rk * dx(nx)
        xnext = xcenter(nx) + 0.5_rk * dx(nx)
        sub(nx-1) = 2.0_rk * xprev * kthis/dx(nx) * kprev/dx(nx-1) &
          / (kthis/dx(nx) + kprev/dx(nx-1))
        dia(nx) = -2.0_rk * xprev * kthis/dx(nx) * kprev/dx(nx-1) &
          / (kthis/dx(nx) + kprev/dx(nx-1)) &
          -  2.0_rk*xnext*kthis/dx(nx)
      case ('insulated')
        kprev = conductivity_fun(temperature(nx-1))
        kthis = conductivity_fun(temperature(nx))
        xprev = xcenter(nx) - 0.5_rk * dx(nx)
        sub(nx-1) = 2.0_rk * xprev * (kthis/dx(nx)) * (kprev/dx(nx-1)) &
          / (kthis/dx(nx) + kprev/dx(nx-1))
        dia(nx) = -2.0_rk * xprev * (kthis/dx(nx)) * (kprev/dx(nx-1)) &
          / (kthis/dx(nx) + kprev/dx(nx-1))
      case default
        call output_write('ERROR: unknown value of bctype_right in build_matrix_cylindrical: ' &
          // trim(adjustl(bctype_right)))
        stop
    endselect
  endsubroutine finite_difference_build_matrix_cylindrical

  subroutine finite_difference_build_matrix_spherical(nx, xcenter, dx, temperature, &
      bctype_left, bctype_right, sub, dia, sup)
    use conductivity_function, only : conductivity_fun
    use output, only : output_write
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: xcenter(:) ! (nx)
    real(rk), intent(in) :: dx(:) ! (nx)
    real(rk), intent(in) :: temperature(:) ! (nx)
    character(*), intent(in) :: bctype_left, bctype_right
    real(rk), intent(out) :: sub(:) ! (nx-1)
    real(rk), intent(out) :: dia(:) ! (nx)
    real(rk), intent(out) :: sup(:) ! (nx-1)

    integer(ik) :: i
    real(rk) :: kprev, kthis, knext
    real(rk) :: xprev, xnext

    ! BC at x=0, i=1
    select case (bctype_left)
      case ('fixed')
        kthis = conductivity_fun(temperature(1))
        knext = conductivity_fun(temperature(2))
        xprev = (xcenter(1) - 0.5_rk * dx(1))**2
        xnext = (xcenter(1) + 0.5_rk * dx(1))**2
        dia(1) = -2.0_rk * xnext * ((kthis/dx(1)) * (knext/dx(2)) &
          / (kthis/dx(1) + knext/dx(2))) &
          - 2.0_rk*xprev*kthis/dx(1)
        sup(1) = 2.0_rk * (kthis/dx(1)) * (knext/dx(2)) &
          / (kthis/dx(1) + knext/dx(2))
      case ('insulated')
        kthis = conductivity_fun(temperature(1))
        knext = conductivity_fun(temperature(2))
        xnext = (xcenter(1) + 0.5_rk*dx(1))**2
        dia(1) = -2.0_rk * xnext * (kthis/dx(1)) * (knext/dx(2)) &
          / (kthis/dx(1) + knext/dx(2))
        sup(1) = 2.0_rk * xnext * (kthis/dx(1)) * (knext/dx(2)) &
          / (kthis/dx(1) + knext/dx(2))
      case default
        call output_write('ERROR: unknown value of bctype_left in build_matrix_spherical: ' &
          // trim(adjustl(bctype_left)))
        stop
    endselect

    do i = 2,nx-1

      kprev = conductivity_fun(temperature(i-1))
      kthis = conductivity_fun(temperature(i))
      knext = conductivity_fun(temperature(i+1))

      xprev = (xcenter(i) - 0.5_rk*dx(i))**2
      xnext = (xcenter(i) + 0.5_rk*dx(i))**2

      sub(i-1) = 2.0_rk * xprev * (kthis/dx(i)) * (kprev/dx(i-1)) &
        / (kthis/dx(i) + kprev/dx(i-1))
      dia(i) = -2.0_rk* (xprev * kthis/dx(i) * kprev/dx(i-1) &
        / (kthis/dx(i) + kprev/dx(i-1)) &
        + xnext * kthis/dx(i) * knext/dx(i+1) & 
        / (kthis/dx(i) + knext/dx(i+1)))
      sup(i) = 2.0_rk * xnext * (kthis/dx(i)) * knext/dx(i+1) &
        / (kthis/dx(i) + knext/dx(i+1))

    enddo ! i = 2,nx-1

    select case (bctype_right)
      case ('fixed')
        kprev = conductivity_fun(temperature(nx-1))
        kthis = conductivity_fun(temperature(nx))
        xprev = (xcenter(nx) - 0.5_rk*dx(nx))**2
        xnext = (xcenter(nx) + 0.5_rk*dx(nx))**2
        sub(nx-1) = 2.0_rk * xprev * kthis/dx(nx) * kprev/dx(nx-1) &
          / (kthis/dx(nx) + kprev/dx(nx-1))
        dia(nx) = -2.0_rk * xprev * kthis/dx(nx) * kprev/dx(nx-1) &
          / (kthis/dx(nx) + kprev/dx(nx-1)) &
          - 2.0_rk * xnext * kthis/dx(nx)
      case ('insulated')
        kprev = conductivity_fun(temperature(nx-1))
        kthis = conductivity_fun(temperature(nx))
        xprev = (xcenter(nx) - 0.5_rk*dx(nx))**2
        sub(nx-1) = 2.0_rk * xprev * kthis/dx(nx) * kprev/dx(nx-1) &
          / (kthis/dx(nx) + kprev/dx(nx-1))
        dia(nx) = -2.0_rk * xprev * kthis/dx(nx) * kprev/dx(nx-1) &
          / (kthis/dx(nx) + kprev/dx(nx-1))
      case default
        call output_write('ERROR: unknown value of bctype_right in build_matrix_spherical: ' &
          // trim(adjustl(bctype_right)))
        stop
    endselect
  endsubroutine finite_difference_build_matrix_spherical

  subroutine finite_difference_build_source_cartesian(nx, xcenter, dx, &
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
        call output_write('ERROR: unknown value of bctype_left in build_source_cartesian: ' &
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
        call output_write('ERROR: unknown value of bctype_right in build_source_cartesian: ' &
          // trim(adjustl(bctype_right)))
        stop
    endselect
    src = -src
  endsubroutine finite_difference_build_source_cartesian

  subroutine finite_difference_build_source_cylindrical(nx, xcenter, dx, &
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
    real(rk) :: x_plus_half, x_minus_half

    select case(bctype_left)
      case ('fixed')
        x_plus_half = xcenter(1) + 0.5_rk * dx(1)
        x_minus_half = 0.0_rk
        src(1) = source_fun(xcenter(1)) * 0.5_rk * (x_plus_half**2 - x_minus_half**2) &
          + 2.0_rk * x_minus_half * conductivity_fun(Tleft) / dx(1) * Tleft
      case ('insulated')
        x_plus_half = xcenter(1) + 0.5_rk * dx(1)
        x_minus_half = 0.0_rk
        src(1) = source_fun(xcenter(1)) * 0.5_rk * (x_plus_half**2 - x_minus_half**2)
      case default
        call output_write('ERROR: unknown value of bctype_left in build_source_cylindrical: ' &
          // trim(adjustl(bctype_left)))
        stop
    endselect

    do i = 2,nx-1
      x_plus_half = xcenter(i) + 0.5_rk * dx(i)
      x_minus_half = xcenter(i) - 0.5_rk * dx(i)
      src(i) = source_fun(xcenter(i)) * 0.5_rk * (x_plus_half**2 - x_minus_half**2)
    enddo ! i = 2,nx-1

    select case(bctype_right)
      case ('fixed')
        x_plus_half = xcenter(nx) + 0.5_rk * dx(nx)
        x_minus_half = xcenter(nx) - 0.5_rk * dx(nx)
        src(nx) = source_fun(xcenter(nx)) * 0.5_rk * (x_plus_half**2 - x_minus_half**2) &
          + 2.0_rk * x_plus_half * conductivity_fun(Tright) / dx(nx) * Tright
      case ('insulated')
        x_plus_half = xcenter(nx) + 0.5_rk * dx(nx)
        x_minus_half = xcenter(nx) - 0.5_rk * dx(nx)
        src(nx) = source_fun(xcenter(nx)) * 0.5_rk * (x_plus_half**2 - x_minus_half**2)
      case default
        call output_write('ERROR: unknown value of bctype_right in build_source_cylindrical: ' &
          // trim(adjustl(bctype_right)))
        stop
    endselect
    src = -src
  endsubroutine finite_difference_build_source_cylindrical

  subroutine finite_difference_build_source_spherical(nx, xcenter, dx, &
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
    real(rk) :: x_plus_half, x_minus_half

    select case(bctype_left)
      case ('fixed')
        x_plus_half = xcenter(1) + 0.5_rk * dx(1)
        x_minus_half = 0.0_rk
        src(1) = source_fun(xcenter(1)) * (x_plus_half**3 - x_minus_half**3) / 3.0_rk &
          + 2.0_rk * x_minus_half * conductivity_fun(Tleft) / dx(1) * Tleft
      case ('insulated')
        x_plus_half = xcenter(1) + 0.5_rk * dx(1)
        x_minus_half = 0.0_rk
        src(1) = source_fun(xcenter(1)) *  (x_plus_half**3 - x_minus_half**3) / 3.0_rk
      case default
        call output_write('ERROR: unknown value of bctype_left in build_source_spherical: ' &
          // trim(adjustl(bctype_left)))
        stop
    endselect

    do i = 2,nx-1
      x_plus_half = xcenter(i) + 0.5_rk * dx(i)
      x_minus_half = xcenter(i) - 0.5_rk * dx(i)
      src(i) = source_fun(xcenter(i)) * (x_plus_half**3 - x_minus_half**3) / 3.0_rk
    enddo ! i = 2,nx-1

    select case(bctype_right)
      case ('fixed')
        x_plus_half = xcenter(nx) + 0.5_rk * dx(nx)
        x_minus_half = xcenter(nx) - 0.5_rk * dx(nx)
        src(nx) = source_fun(xcenter(nx)) * (x_plus_half**3 - x_minus_half**3) / 3.0_rk &
          + 2.0_rk * x_plus_half**2 * conductivity_fun(Tright) / dx(nx) * Tright
      case ('insulated')
        x_plus_half = xcenter(nx) + 0.5_rk * dx(nx)
        x_minus_half = xcenter(nx) - 0.5_rk * dx(nx)
        src(nx) = source_fun(xcenter(nx)) * (x_plus_half**3 - x_minus_half**3) / 3.0_rk
      case default
        call output_write('ERROR: unknown value of bctype_right in build_source_spherical: ' &
          // trim(adjustl(bctype_right)))
        stop
    endselect
    src = -src
  endsubroutine finite_difference_build_source_spherical

  subroutine finite_difference_solve(geometry, nx, xcenter, dx, &
      bctype_left, bctype_right, bcval_left, bcval_right, &
      max_iter, tol_temperature, init_temperature, temperature)
    use linalg, only : trid, norm
    use output, only : output_write
    character(*), intent(in) :: geometry
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

    select case (geometry)
      case ('cartesian')
        call finite_difference_build_source_cartesian(nx, xcenter, dx, &
          bctype_left, bctype_right, bcval_left, bcval_right, qcpy)
      case ('cylindrical')
        call finite_difference_build_source_cylindrical(nx, xcenter, dx, &
          bctype_left, bctype_right, bcval_left, bcval_right, qcpy)
      case ('spherical')
        call finite_difference_build_source_spherical(nx, xcenter, dx, &
          bctype_left, bctype_right, bcval_left, bcval_right, qcpy)
      case default
        call output_write('ERROR1 : Unknown geometry in finite_difference_solve: ' // trim(adjustl(geometry)))
    endselect

    temperature = init_temperature

    ! Picard iteration
    do iter = 1,max_iter

      temperature_old = temperature

      ! must rebuild matrix since thermal conductivity may change on each iteration
      ! must copy the source since it is used as scratch space by trid
      select case (geometry)
        case ('cartesian')
          call finite_difference_build_matrix_cartesian(nx, xcenter, dx, temperature, &
            bctype_left, bctype_right, sub, dia, sup)
        case ('cylindrical')
          call finite_difference_build_matrix_cylindrical(nx, xcenter, dx, temperature, &
            bctype_left, bctype_right, sub, dia, sup)
        case ('spherical')
          call finite_difference_build_matrix_spherical(nx, xcenter, dx, temperature, &
            bctype_left, bctype_right, sub, dia, sup)
        case default
          call output_write('ERROR: unknown geometry in finite_difference_solve: ' // trim(adjustl(geometry)))
          stop
      endselect
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
