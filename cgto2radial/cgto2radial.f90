!
! cgto2radial.f90
!
! Reads atomic orbitals expanded in Gaussian type orbitals.
! Transforms them to the rR(r) radial functions and write these
! functions to file.
!
! Input file syntax:
! <L=S,P,D,F,G,H,I,K,L> <nprim> <ncols>
! <exponent-1> <col1-coef-1> ... <col-ncols-coef-1>
!   . . .
! <exponent-nprim> <col1-coef-nprim> ... <col-ncols-coef-nprim>
!
! Here L denotes orbital angular momentum (S, P, ...),
! left column - exponents,
! right columns - contraction coefficients of contracted GTOs number 1..n
!
! Example of input file:
! S  30 1
!      18.027623573072357700       0.001987370146652646
!       9.070832728948866830      -0.038090372916132827
!       6.596969257417358180       0.176680531125234311
!          . . .
!       0.001216666666666667       0.000295599517296362
! Here: S orbital, 1 contracted GTO, 30 pairs 'exp - coef' total
!
! author: Alexander Oleynichenko, 4 Aug 2016
!         alexvoleynichenko@gmail.com
!


program cgto2radial

    use gto
    implicit none
    
    character(len=32) :: path
    integer :: L
    integer :: i
    integer :: num_prim
    integer :: num_funs
    real(8), dimension(:), allocatable :: alpha
    real(8), dimension(:,:), allocatable :: coeffs
    real(8), dimension(:), allocatable :: r
    real(8) :: expect_r
    real(8) :: f_square
    
    interface
    subroutine read_contracted_gtos(path, L, num_prim, num_funs, alpha, coeffs)
        character(len=*), intent(in) :: path
        integer, intent(out) :: L
        integer, intent(out) :: num_prim
        integer, intent(out) :: num_funs
        real(8), dimension(:), allocatable :: alpha
        real(8), dimension(:,:), allocatable :: coeffs
    end subroutine read_contracted_gtos
    subroutine calculate_radial(L, nprim, nfun, alpha, coeffs, left, right, step)
    integer, intent(in) :: L
    integer, intent(in) :: nprim
    integer, intent(in) :: nfun
    real(8), dimension(:), intent(in) :: alpha
    real(8), dimension(:,:), intent(in) :: coeffs
    real(8), intent(in) :: left
    real(8), intent(in) :: right
    real(8), intent(in) :: step
    end subroutine calculate_radial
    end interface
    
    ! extract input file name from the command-line arguments
    if (iargc() /= 1) then
        call print_usage
    end if
    call getarg(1, path)
    print *, 'Input file           = ', path
    
    ! read GTO expansions
    call read_contracted_gtos(path, L, num_prim, num_funs, alpha, coeffs)
    print *, 'Angular momentum     = ', L
    print *, 'Number of primitives = ', num_prim
    print *, 'Number of functions  = ', num_funs
    
    call calculate_radial(L, num_prim, num_funs, alpha, coeffs, 0.0d0, 5.0d0, 0.001d0)
    
    print '(a5,a12,a10)', 'n', '|R(r)|', '<r>'
    do i = 1, num_funs
        f_square = integrate_contracted_gto(L, num_prim, alpha, coeffs(:,i), 0)
        expect_r = integrate_contracted_gto(L, num_prim, alpha, coeffs(:,i), 1)
        print '(i5,f12.6,f10.4)', i, f_square, expect_r
    end do
    

end program cgto2radial



subroutine print_usage()

    print *, 'Usage: python plot_cgto.py <input-file>'
    print *, 'See source code for more details (input format, etc).'
    stop

end subroutine print_usage



subroutine errquit(message)
    
    implicit none

    character(len=*), intent(in) :: message
    
    print *, 'Error: ', message
    stop
    
end subroutine errquit



subroutine read_contracted_gtos(path, L, num_prim, num_funs, alpha, coeffs)

    implicit none
    
    character(len=*), intent(in) :: path
    integer, intent(out) :: L
    integer, intent(out) :: num_prim
    integer, intent(out) :: num_funs
    real(8), dimension(:), allocatable :: alpha
    real(8), dimension(:,:), allocatable :: coeffs

    integer :: descr
    character(len=4) :: angular_symbol
    integer :: i, j
    
    descr = 15
    open(unit=descr, file=path, form='formatted', err=13)
    read(descr, *) angular_symbol, num_prim, num_funs

    allocate(alpha(num_prim))
    allocate(coeffs(num_prim, num_funs))
    do i = 1, num_prim
        read(descr,*) alpha(i), (coeffs(i,j),j=1,num_funs)
    end do

    call angular_momentum_to_int(angular_symbol, L)
    
    close(unit=descr)
    return
    
    13 continue
    call errquit('input file not found')    

end subroutine read_contracted_gtos



! S,P,D,... -> 0,1,2,...
subroutine angular_momentum_to_int(angular_symbol, L)

    implicit none
    
    character(len=*), intent(in) :: angular_symbol
    character :: symbol
    integer, intent(out) :: L
    
    symbol = angular_symbol(1:1)
    
    if (symbol == 'S') then
        L = 0
    else if (symbol == 'P') then
        L = 1
    else if (symbol == 'D') then
        L = 2
    else if (symbol == 'F') then
        L = 3
    else if (symbol == 'G') then
        L = 4
    else if (symbol == 'H') then
        L = 5
    else if (symbol == 'I') then
        L = 6
    else if (symbol == 'K') then
        L = 7
    else if (symbol == 'L') then
        L = 8
    else
        call errquit('wrong angular momentum symbol')
    end if
    
end subroutine angular_momentum_to_int



! calculates values of 'nfun' contracted Gaussians N*r^l*exp{-alpha*r^2} on the
! radial grid. The rR(r) functions are written to the 'R.dat' file
subroutine calculate_radial(L, nprim, nfun, alpha, coeffs, left, right, step)

    use gto
    implicit none
    integer, intent(in) :: L
    integer, intent(in) :: nprim
    integer, intent(in) :: nfun
    real(8), dimension(:), intent(in) :: alpha
    real(8), dimension(:,:), intent(in) :: coeffs
    real(8), intent(in) :: left
    real(8), intent(in) :: right
    real(8), intent(in) :: step
    real(8) :: r
    real(8), dimension(nfun) :: f_r
    integer :: i
    integer :: descr
    
    open(unit=descr, file='R.dat', form='formatted')
    
    r = left
    do while (r <= right)
        do i = 1, nfun
            f_r(i) = evaluate_contracted_gto(L, nprim, alpha, coeffs(:,i), r)
        end do
        write (descr,'(f12.6,10e24.12)'), r, (r*f_r(i),i=1,nfun)
        r = r + step
    end do
    
    close(unit=descr)

end subroutine calculate_radial








