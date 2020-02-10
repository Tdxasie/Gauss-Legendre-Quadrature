! Quadrature Module
! DM1 - CSA 1
! Pierre Garnier | pierre.garnier1618@gmail.com
! pgarnier.dev
! 10/02/2020
module quadrature
    use precision
    implicit none
    private
    integer :: nquad
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), allocatable :: wk(:), xk(:)
    public :: Quad_GL, init_quad, clean_mem, show_points, show_weights
contains
    function Quad_GL(f, a, b) ! computes the Guauss Legendre quadrature
        real(dp), intent(in) :: a, b
        integer :: k
        real(dp) :: f, Quad_GL, rescale, shift
        real(dp) :: wk(nquad), xk(nquad)
        call GL_xw(xk, wk)
        shift = (a+b)/2
        rescale = (b-a)/2
        do k = 1, nquad
            Quad_GL = Quad_GL + wk(k)*f(rescale*xk(k)+shift)
        end do
        Quad_GL = Quad_GL * rescale
    end function Quad_GL
!
    subroutine GL_xw(xk, wk) ! generates points and weights
        ! if n odd:
        !   dpn structure: sym + piv + reversed_sym
        !   x   structure: sym + 0 + reversed_sym*-1
        ! if n even:
        !   dpn structure: sym + reversed_sym*-1
        !   x   structure: sym + reversed_sym*-1
        real(dp), intent(inout) :: wk(nquad), xk(nquad)
        integer :: k, i
        real(dp) :: dpn, dx, p0, p1, pn, x
        real(dp) :: dpnv(nquad)
        if (mod(nquad,2) == 0) then ! even case
            do k = 1, nquad/2
                x = cos(((4*k-1)*pi)/(4*nquad+2))
                do ! compute the Legendre polynomial
                    p0 = 1.0_dp
                    p1 = x
                    do i = 2, nquad
                        pn = ((2*i-1)*x*p1-(i-1)*p0)/i
                        p0 = p1
                        p1 = pn
                    end do
                    dpn = (nquad*p0-x*nquad*pn)/(1-x*x)
                    dx = pn/dpn
                    x = x - dx
                    if (abs(dx)<1.0e-16_dp) exit
                end do
                xk(k) = x
                dpnv(k) = dpn
            end do
            ! using symmetry
            dpnv(nquad/2+1:nquad) = dpnv(nquad/2:1:-1)*(-1)
            xk(nquad/2+1:nquad) = xk(nquad/2:1:-1)*(-1)
        else ! odd case
            do k = 1, nquad/2 + 1 
                x = cos(((4*k-1)*pi)/(4*nquad+2))
                do ! compute the Legendre polynomial
                    p0 = 1.0_dp
                    p1 = x
                    do i = 2, nquad
                        pn = ((2*i-1)*x*p1-(i-1)*p0)/i
                        p0 = p1
                        p1 = pn
                    end do
                    dpn = (nquad*p0-x*nquad*pn)/(1-x*x)
                    dx = pn/dpn
                    x = x - dx
                    if (abs(dx)<1.0e-16_dp) exit
                end do
                xk(k) = x
                dpnv(k) = dpn
            end do
            ! using symmetry
            dpnv(nquad/2+2:nquad) = dpnv(nquad/2:1:-1)
            xk(nquad/2+2:nquad) = xk(nquad/2:1:-1)*(-1)
        end if
        wk = 2.0_dp/((1.0_dp-xk**2)*dpnv**2)
    end subroutine GL_xw
!
    subroutine init_quad(filename) ! reads order of quadrature in "settings" file and allocates memory
        character (*), intent(in) :: filename
        namelist /settings/ nquad
        open(unit=7, file=filename)
        read(7, settings)
        allocate(xk(nquad), wk(nquad))
        close(7)
    end subroutine init_quad
!
    subroutine clean_mem() ! frees allocated memory
        deallocate(wk, xk)
    end subroutine clean_mem
!
    subroutine show_weights()
        write (6,*) wk
    end subroutine show_weights
!
    subroutine show_points()
        write (6,*) xk
    end subroutine show_points
end module quadrature