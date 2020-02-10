! Functions Module
! DM1 - CSA 1
! Pierre Garnier | pierre.garnier1618@gmail.com
! pgarnier.dev
! 10/02/2020
module functions
    use precision
    implicit none
    private
    integer :: npoly
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), allocatable :: cf(:)
    public :: clean_coef, f1, f2, horner, read_poly
contains  
    function f1(x)
        real(dp), intent(in) :: x
        real(dp) :: f1
        f1 = 2.0_dp / (1+x**2)
        return
    end function f1
!
    function f2(x)
        real(dp), intent(in) :: x
        real(dp) :: f2
        f2 = sin(x)**2
        return
    end function f2
!
    function horner(x)
        real(dp), intent(in) :: x
        integer :: i
        real(dp) :: horner
        horner = cf(npoly)*x
        do i=1, npoly-1
            horner = (horner+cf(i))*x
        end do
        horner = horner + cf(0)
        return
    end function horner
!
    subroutine read_poly(filename)
        character (*), intent(in) :: filename
        integer :: i
        namelist /poly/ npoly
        open(unit=7, file=filename)
        read(7, poly) ! coefficients are ordered in increasing power of x
        allocate(cf(0:npoly))
        read(7,*) (cf(i), i=0, npoly)
        close(7)
    end subroutine read_poly
!    
    subroutine clean_coef()
        deallocate(cf)
    end subroutine clean_coef
!
end module functions