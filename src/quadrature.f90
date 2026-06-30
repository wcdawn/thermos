module quadrature
use kind, only : rk, ik
implicit none

private

public :: quadrature_gauss_legendre, QuadraturePoint

type QuadraturePoint
  real(rk) :: x ! abscissa
  real(rk) :: w ! weight
endtype QuadraturePoint

contains

  subroutine quadrature_gauss_legendre(n, q)
    use exception_handler, only : exception_fatal
    integer(ik), intent(in) :: n ! quadrature order
    type(QuadraturePoint), allocatable, intent(out) :: q(:) ! (n)

    character(1024) :: line

    ! NOTE: that for Gauss-Legendre quadratures, the order of the quadrature
    ! is the same as the number of points in the quadrature.

    allocate(q(n))

    select case (n)
      case (1)
        q(:)%x = [ 0.0_rk ]
        q(:)%w = [ 2.0_rk ]
      case (2)
        q(:)%x = [ -1.0_rk/sqrt(3.0_rk), +1.0_rk/sqrt(3.0_rk) ]
        q(:)%w = [ 1.0_rk, 1.0_rk ]
      case (3)
        q(:)%x = [ -sqrt(3.0_rk*0.2_rk), 0.0_rk, +sqrt(3.0_rk*0.2_rk) ]
        q(:)%w = [ 5.0_rk/9.0_rk, 8.0_rk/9.0_rk, 5.0_rk/9.0_rk ]
      case (4)
        q(:)%x = [ -sqrt(3.0_rk/7.0_rk + 2.0_rk/7.0_rk*sqrt(6.0_rk*0.2_rk)), &
                   -sqrt(3.0_rk/7.0_rk - 2.0_rk/7.0_rk*sqrt(6.0_rk*0.2_rk)), &
                   +sqrt(3.0_rk/7.0_rk - 2.0_rk/7.0_rk*sqrt(6.0_rk*0.2_rk)), &
                   +sqrt(3.0_rk/7.0_rk + 2.0_rk/7.0_rk*sqrt(6.0_rk*0.2_rk)) ]
        q(:)%w = [ (18.0_rk - sqrt(30.0_rk))/36.0_rk, &
                   (18.0_rk + sqrt(30.0_rk))/36.0_rk, &
                   (18.0_rk + sqrt(30.0_rk))/36.0_rk, &
                   (18.0_rk - sqrt(30.0_rk))/36.0_rk ]
      case (5)
        q(:)%x = [ -sqrt(5.0_rk + 2.0_rk*sqrt(10.0_rk/7.0_rk))/3.0_rk, &
                   -sqrt(5.0_rk - 2.0_rk*sqrt(10.0_rk/7.0_rk))/3.0_rk, &
                   0.0_rk, &
                   +sqrt(5.0_rk - 2.0_rk*sqrt(10.0_rk/7.0_rk))/3.0_rk, &
                   +sqrt(5.0_rk + 2.0_rk*sqrt(10.0_rk/7.0_rk))/3.0_rk ]
        q(:)%w = [ (322.0_rk - 13.0_rk*sqrt(70.0_rk))/900.0_rk, &
                   (322.0_rk + 13.0_rk*sqrt(70.0_rk))/900.0_rk, &
                   128.0_rk/225.0_rk                          , &
                   (322.0_rk + 13.0_rk*sqrt(70.0_rk))/900.0_rk, &
                   (322.0_rk - 13.0_rk*sqrt(70.0_rk))/900.0_rk  ]
      case default
        write(line, '(a,i0)') 'Unacceptable Gauss-Legendre quadrature order: ', n
        call exception_fatal(line)
    endselect
  endsubroutine quadrature_gauss_legendre

endmodule quadrature
