module gto

implicit none

contains


function evaluate_contracted_gto(L, nprim, alpha, coeffs, r) result(res)

    implicit none
    integer, intent(in) :: L
    integer, intent(in) :: nprim
    real(8), dimension(:), intent(in) :: alpha
    real(8), dimension(:), intent(in) :: coeffs
    real(8), intent(in) :: r
    real(8) :: res
    integer :: i
    
    res = 0.0d0
    do i = 1, nprim
        res = res + coeffs(i) * evaluate_gaussian(L, alpha(i), r)
    end do

end function evaluate_contracted_gto



function integrate_contracted_gto(L, nprim, alpha, coeffs, rpow) result(res)

    implicit none
    integer, intent(in) :: L
    integer, intent(in) :: nprim
    real(8), dimension(:), intent(in) :: alpha
    real(8), dimension(:), intent(in) :: coeffs
    integer, intent(in) :: rpow
    real(8) :: res
    real(8) :: r
    real(8) :: right
    real(8) :: step
    real(8), parameter :: EPS = 1e-8
    real(8) :: ds
    real(8) :: f_r
    
    res = 0.0d0
    r = 0.0d0
    right = 50.0d0
    step = 1.0d-3
    do while (r <= right)
        f_r = evaluate_contracted_gto(L, nprim, alpha, coeffs, r)
        ds = step * r**2 * (r**rpow) * f_r**2
        res = res + ds
        r = r + step
    end do

end function integrate_contracted_gto



! calculate value of the normalized Gaussian function N*r^l*exp{-alpha*r^2} at point 'r'
function evaluate_gaussian(L, alpha, r) result(res)

    implicit none
    integer, intent(in) :: L
    real(8), intent(in) :: alpha
    real(8), intent(in) :: r
    real(8) :: res
    
    res = norm_factor(L,alpha) * (r**L) * exp(-alpha * r**2)

end function evaluate_gaussian



! calculate normalization factor for the Gaussian function r^l exp{-alpha*r^2}
! See Helgaker et al. 2000 p.234 for details
function norm_factor(L, alpha) result(res)

    implicit none
    integer, intent(in) :: L
    real(8), intent(in) :: alpha
    real(8) :: res
    real(8), parameter :: PI = 3.1415926535897932_8

    res = 2.0d0 * ((2.0d0*alpha)**0.75d0) / (PI**0.25d0) * &
          sqrt(2.0d0**L/double_factorial(2*L+1)) *         &
          (sqrt(2.0d0*alpha)**L)

end function norm_factor



function double_factorial(n) result(ans)

    implicit none
    integer, intent(in) :: n ! input
    integer :: ans ! output
    integer :: i
  
    if (mod(n,2) > 0) then
        ! even
        ans = product((/(i, i=2,n,2)/))
    else
        ! odd
        ans = product((/(i, i=1,n,2)/))
    end if
  
end function double_factorial


end module gto