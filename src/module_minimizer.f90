! Defines everything around the minimizer
module module_minimizer
    use precision_kinds, only: i2b , dp
    implicit none
    private
    type lbfgsb_type
        !     Declare variables and parameters needed by mylbfgsb.f90>setulb

        integer                :: n, m=5, iprint=1
        real(dp)               :: factr  = 1.0d+12, pgtol  = 1.0d-5
        character(len=60)      :: task, csave
        logical                :: lsave(4)
        integer                :: isave(44)
        real(dp)               :: dsave(29)
        integer,  allocatable  :: nbd(:), iwa(:)
        real(dp), allocatable  :: x(:), l(:), u(:), g(:), wa(:)
        integer :: itermax
    end type lbfgsb_type

    type (lbfgsb_type) :: lbfgsb

    public :: lbfgsb, init_lbfgsb

contains

    !===================================================================================================================================
    subroutine init_lbfgsb
        ! this subroutine gets the informations in input file and then allocate, prepare, compute needed data

        use module_solvent, only: solvent
        use module_grid, only: grid
        use module_input, only: getinput

        implicit none

        integer :: i, icg, is, io, iz, iy, ix, itermax
        real(dp), parameter :: epsdp=epsilon(1._dp)

        itermax = getinput%int("maximum_iteration_nbr", defaultvalue=50, assert=">0")

        lbfgsb%n = grid%nx * grid%ny * grid%nz * grid%no * solvent(1)%nspec
        allocate ( lbfgsb%nbd(lbfgsb%n), lbfgsb%x(lbfgsb%n), lbfgsb%l(lbfgsb%n), lbfgsb%u(lbfgsb%n), lbfgsb%g(lbfgsb%n) )
        allocate ( lbfgsb%iwa(3*lbfgsb%n) )
        allocate ( lbfgsb%wa(2*lbfgsb%m*lbfgsb%n + 5*lbfgsb%n + 11*lbfgsb%m*lbfgsb%m + 8*lbfgsb%m) )

        if (.not. allocated(solvent)) then
            print*, "probleme dans module_minimizer > init_lbfgsb"
            print*, "solvent is not allocated nor initiated"
            error stop
        end if

        if (.not. allocated(solvent(1)%density) ) then
            print*, "probleme dans module_minimizer > init_lbfgsb"
            print*, "solvent is allocated but not solvent%density"
            error stop
        end if

        !   init vector x (the densities)
        icg=0
        do is=1,solvent(1)%nspec
            do io=1,grid%no
                do iz=1,grid%nz
                    do iy=1,grid%ny
                        do ix=1,grid%nx
                            icg=icg+1
                            if (solvent(is)%vext(ix,iy,iz,io) < solvent(is)%vext_threeshold) then
                                lbfgsb%x(icg) = solvent(is)%density(ix,iy,iz,io)
                                lbfgsb%nbd(icg) = 1 ! lower bounded
                                lbfgsb%l(icg) = 0._dp ! the lower bound
                            else
                                ! if (abs(solvent(is)%density(ix,iy,iz,io))>epsdp) then
                                !     print*, "Dans module_minimizer, je croyais que qd vext>threeshold, density init à 0?!!!!"
                                !     print*, "On a solvent(is)%vext(ix,iy,iz,io) =" ,solvent(is)%vext(ix,iy,iz,io)
                                !     print*, "solvent(is)%density(ix,iy,iz,io) =", solvent(is)%density(ix,iy,iz,io)
                                !     print*, "ix,iy,iz,io=",ix,iy,iz,io
                                !     error stop
                                ! end if
                                lbfgsb%x(icg) = solvent(is)%density(ix,iy,iz,io)
                                lbfgsb%nbd(icg) = 2 ! lower and upper bounded
                                lbfgsb%l(icg) = 0._dp!solvent(is)%density(ix,iy,iz,io)
                                lbfgsb%u(icg) = 0._dp!solvent(is)%density(ix,iy,iz,io) ! lower and upper bounds the same, I hope lbfgsb understands this means don't touch to this?
                            end if
                        end do
                    end do
                end do
            end do
        end do

        if (icg /= lbfgsb%n) then
            print*, "icg should be == lbfgsb%n in init of module_minimizer.f90 (l.174)"
            stop
        end if

        !     We start the iteration by initializing task.
        lbfgsb%task = 'START'

    end subroutine init_lbfgsb

end module module_minimizer
