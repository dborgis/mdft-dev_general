! Defines everything around the minimizer
module module_minimizer
    use precision_kinds, only: i2b , dp
    implicit none
    private
    type lbfgsb_type
        !     Declare variables and parameters needed by mylbfgsb.f90>setulb

        integer                :: n, m=1, iprint=1
        real(dp)               :: factr  = 1.0d+3 ! The iteration will stop when
        ! real(dp)               :: factr  = 1.0d+12 ! The iteration will stop when
                                                   ! (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
                                                   ! where epsmch is the machine precision, which is automatically
                                                   ! generated by the code. Typical values for factr: 1.d+12 for
                                                   ! low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
                                                   ! high accuracy. On exit factr is unchanged.
        real(dp)               :: pgtol  = 1.0d-7 ! The iteration will stop when max{|proj g_i | i = 1, ..., n} <= pgtol
        ! real(dp)               :: pgtol  = 1.0d-3 ! The iteration will stop when max{|proj g_i | i = 1, ..., n} <= pgtol
                                                  ! where pg_i is the ith component of the projected gradient.
        character(len=60)      :: task, csave
        logical                :: lsave(4)
        integer                :: isave(44)
        real(dp)               :: dsave(29)
        ! integer,   allocatable :: nbd(:)
        integer,   allocatable :: iwa(:)
        ! real(dp),  allocatable :: l(:)
        real(dp),  allocatable :: wa(:)
        ! real(dp), allocatable :: x(:)
        ! real(dp),  allocatable :: u(:)
        ! real(dp), allocatable :: g(:)
        integer :: itermax=100
    contains
        procedure, nopass :: init => init_lbfgsb
    end type lbfgsb_type

    type (lbfgsb_type) :: lbfgsb

    public :: lbfgsb

contains

    !===================================================================================================================================
    subroutine init_lbfgsb
        ! this subroutine gets the informations in input file and then allocate, prepare, compute needed data

        use module_solvent, only: solvent
        use module_grid, only: grid
        use module_input, only: getinput

        implicit none
        integer :: n, m
        
        lbfgsb%itermax = getinput%int("maximum_iteration_nbr", defaultvalue=100, assert=">0")

        lbfgsb%n = grid%nx * grid%ny * grid%nz * grid%no * solvent(1)%nspec
        n=lbfgsb%n
        m=lbfgsb%m
        ! allocate ( lbfgsb%nbd(n) )
        ! allocate ( lbfgsb%l(n) )
        ! allocate ( lbfgsb%u(n) )
        ! allocate !, lbfgsb%x(n), lbfgsb%g(n)
        allocate ( lbfgsb%iwa(3*n) )
        allocate ( lbfgsb%wa(2*m*n + 5*n + 11*m*m + 8*m) )

        if (.not. allocated(solvent)) then
            print*, "probleme dans module_minimizer > init_lbfgsb"
            print*, "solvent is not allocated nor initiated"
            error stop
        end if

        if (.not. allocated(solvent(1)%xi) ) then
            print*, "probleme dans module_minimizer > init_lbfgsb"
            print*, "solvent is allocated but not solvent%xi"
            error stop
        end if

        !
        !   L-BFGS-B needs to know for each var to optimize if it is lower bounded (nbd=1), lower and upper bounded (nbd=2)
        !   Give also the bounds (l and u)
        !
        ! icg=0
        ! do is=1,solvent(1)%nspec
        !     do iz=1,grid%nz
        !         do iy=1,grid%ny
        !             do ix=1,grid%nx
        !                 do io=1,grid%no
        !                     icg=icg+1
        !                     lbfgsb%nbd(icg) = 0 ! non-bounded
        !                 end do
        !             end do
        !         end do
        !     end do
        ! end do

        ! if (icg /= lbfgsb%n) then
        !     print*, "icg should be == lbfgsb%n in init of module_minimizer.f90 (l.174)"
        !     stop
        ! end if

        !     We start the iteration by initializing task.
        lbfgsb%task = 'START'

    end subroutine init_lbfgsb

end module module_minimizer
