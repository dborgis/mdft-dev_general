! Defines everything around the minimizer
module module_minimizer

    use precision_kinds, only: i2b , dp

    IMPLICIT NONE

    character(10) :: minimizer_type ! bfgs or cg
    integer(i2b):: ncg ! nbr of variables
    real(dp):: FF ! value of the minimum found for the functional
    real(dp), allocatable :: cg_vect_new(:,:,:,:,:), df_new(:,:,:,:,:)
    real(dp):: epsg
    integer(i2b) :: iter ! iteration number of the minimizer. Goes from 1 to itermax during execution.
    integer(i2b) :: itermax ! maximum number of iterations

contains

    !===========================================================================================================================
    subroutine from_cgvect_get_rho

        use system      ,only: solvent
        use module_grid, only: grid
        implicit none
        integer(i2b) :: i, j, k, io, s, icg, nx, ny, nz, no

        nx = grid%nx
        ny = grid%ny
        nz = grid%nz
        no = grid%no

        do concurrent (s=1:solvent(1)%nspec)
        if (.not. allocated (solvent(s)%rho)) then
            allocate (solvent(s)%rho(nx,ny,nz,no), source=0._dp ,stat=i)
        else
            i = 0
            solvent(s)%rho = 0._dp
        end if
        end do
        if (i /= 0) stop "ERROR in from_cgvect_get_rho"

        do concurrent( s= 1: solvent(1)%nspec)
            solvent(s)%rho = cg_vect_new(:,:,:,:,s)**2 * solvent(s)%rho0
        end do
    end subroutine from_cgvect_get_rho

    !===========================================================================================================================
    subroutine from_rho_get_n

        use system      ,only: solvent
        use module_grid, only: grid
        implicit none
        integer(i2b) :: i,j,k,o,p,s,io
        integer(i2b) :: nfft(3)

        nfft = grid%n_nodes

        do concurrent (s=1:solvent(1)%nspec)
        if (.not. allocated(solvent(s)%n)) then
            allocate (solvent(s)%n(nfft(1),nfft(2),nfft(3)), source=0._dp ,stat=i)
        else
            i = 0
            solvent(s)%n = 0._dp
        end if
        end do
        if (i /= 0) stop "ERROR in allocate solvent%n in from_rho_get_n"

        do s =1, solvent(1)%nspec
            do i =1, grid%n_nodes(1)
                do j =1, grid%n_nodes(2)
                    do k =1, grid%n_nodes(3)
                        do io = 1, grid%no
                            solvent(s)%n(i,j,k) = solvent(s)%n(i,j,k) + solvent(s)%rho(i,j,k,io) * grid%w(io)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine from_rho_get_n


    !===========================================================================================================================
    subroutine deallocate_solvent_rho
        use system  ,only: solvent
        integer(i2b) :: s
        do concurrent (s=1:solvent(1)%nspec)  ! solvent%rho is now useless
        deallocate (solvent(s)%rho)
        end do
    end subroutine deallocate_solvent_rho

    !===================================================================================================================================
    subroutine prepare_minimizer
        ! this subroutine gets the informations in input file and then allocate, prepare, compute needed data

        use precision_kinds, ONLY: dp
        use module_input, ONLY: getinput
        use system , ONLY: solvent
        use module_grid, only: grid

        IMPLICIT NONE

        integer :: nx, ny, nz, no, ns, io

        nx = grid%nx
        ny = grid%ny
        nz = grid%nz
        no = grid%no
        ns = solvent(1)%nspec

        ncg = nx*ny*nz*no*ns

        allocate( cg_vect_new( nx, ny, nz, no, ns) ,source=0._dp)
        allocate(      df_new( nx, ny, nz, no, ns) ,source=0._dp)

        FF = 0._dp

    END subroutine prepare_minimizer

end module module_minimizer
