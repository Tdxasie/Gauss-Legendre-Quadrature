! Quadrature Module
! DM1 - CSA 1
! Pierre Garnier | pierre.garnier1618@gmail.com
! pgarnier.dev
! 10/02/2020
program quad_mod
    use quadrature
    use functions
    use precision
    implicit none
    real(dp) :: a, b
    real(dp), parameter :: pi = acos(-1.0_dp)
    call init_quad("settings.conf")
    call read_poly("poly.dat")
    a = -1.0_dp
    b = 1.0_dp
    write(6, *) Quad_GL(f1, a, b) ! expected output: pi
    a = pi*(-1)
    b = pi
    write(6, *) Quad_GL(f2, a, b) ! expected output: pi
    a = -4.0_dp
    b = 0.0_dp
    write(6, *) Quad_GL(horner, a, b) ! expected output: 22.0_dp
    call clean_mem()
    call clean_coef()
end program quad_mod