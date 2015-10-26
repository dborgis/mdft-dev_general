! Defines everything around the minimizer
module module_minimizer

    use precision_kinds, only: i2b , dp

    IMPLICIT NONE

    type lbfgsb_type
        !     Declare variables and parameters needed by the code.
        !       Note thar we wish to have output at every iteration.
        !          iprint=1
        !
        !       We also specify the tolerances in the stopping criteria.
        !          factr  = 1.0d+7, pgtol  = 1.0d-5
        !
        !       A description of all these variables is given at the beginning
        !       of  the driver
        !
        integer                :: n, m=5, iprint=1
        real(dp)               :: factr  = 1.0d+7, pgtol  = 1.0d-5
        !
        character(len=60)      :: task, csave
        logical                :: lsave(4)
        integer                :: isave(44)
        ! real(dp)               :: f
        real(dp)               :: dsave(29)
        integer,  allocatable  :: nbd(:), iwa(:)
        real(dp), allocatable  :: x(:), l(:), u(:), g(:), wa(:)
        integer :: itermax
    end type lbfgsb_type

    type (lbfgsb_type) :: lbfgsb

contains

    !===================================================================================================================================
    subroutine init
        ! this subroutine gets the informations in input file and then allocate, prepare, compute needed data

        use module_solvent, only: solvent
        use module_grid, only: grid
        use module_input, only: getinput

        IMPLICIT NONE

        integer :: i, icg, is, io, iz, iy, ix, itermax ! dummy

        itermax = getinput%int("maximum_iteration_nbr", defaultvalue=50, assert=">0")

        lbfgsb%n = grid%nx * grid%ny * grid%nz * grid%no * solvent(1)%nspec
        lbfgsb%m = 5
        lbfgsb%iprint = 1
        allocate ( lbfgsb%nbd(lbfgsb%n), lbfgsb%x(lbfgsb%n), lbfgsb%l(lbfgsb%n), lbfgsb%u(lbfgsb%n), lbfgsb%g(lbfgsb%n) )
        allocate ( lbfgsb%iwa(3*lbfgsb%n) )
        allocate ( lbfgsb%wa(2*lbfgsb%m*lbfgsb%n + 5*lbfgsb%n + 11*lbfgsb%m*lbfgsb%m + 8*lbfgsb%m) )

        !   bounds
        do i=1, lbfgsb%n
            lbfgsb%nbd(i) = 1
            lbfgsb%l(i)   = epsilon(1._dp)
            !lbfgsb%u(i)   = 1.0d2
        end do

        !   init vector x (the densities)
        icg=0
        do is=1,solvent(1)%nspec
            do io=1,grid%no
                do iz=1,grid%nz
                    do iy=1,grid%ny
                        do ix=1,grid%nx
                            icg=icg+1
                            lbfgsb%x(icg) = solvent(is)%density(ix,iy,iz,io)
                        end do
                    end do
                end do
            end do
        end do
        if (icg /= lbfgsb%n) then
            print*, "icg should be == lbfgsb%n in init of module_minimizer.f90 (l.174)"
            stop
        end if

    end subroutine init

end module module_minimizer
