!===================================================================================================================================
module minimizer
!===================================================================================================================================
! Defines minimizer parameters

    USE precision_kinds,only: i2b , dp

    IMPLICIT NONE

    character(10) :: minimizer_type ! bfgs or cg
    integer(i2b):: ncg ! nbr of variables
    real(dp):: FF ! value of the minimum found for the functional
    real(dp), allocatable , dimension(:) :: CG_vect ! functional. It has to be a one dimensional array
    real(dp), allocatable , dimension(:) :: dF ! gradient of the functional
    integer(i2b):: ITERMAX ! maximum number of iterations ! a terme le remplacer dans bfgs par minimizer_iter
    real(dp):: epsg , epsmch , factr ! precision parameters
    integer(i2b):: minimizer_iter
    ! parametres pour gradients conjugues BFGS
    ! integer ( kind = i2b ), parameter :: mmax=4 ! doit etre remplacé par mcg
    integer(i2b):: mcg, iprint
    integer(i2b), allocatable, dimension(:) :: nbd
    integer(i2b), allocatable, dimension(:) :: iwa
    integer(i2b), dimension (1:44) :: isave
    real(dp):: pgtol
    real(dp), allocatable, dimension(:) :: ll, uu
    real(dp), dimension (1:29) :: dsave
    real(dp), allocatable, dimension(:) :: wa
    character*60 :: task, csave
    logical, dimension (1:4) :: lsave
    integer(i2b), allocatable , dimension (:,:,:,:) :: indicg
    integer(i2b) :: iter ! iteration number of the minimizer. Goes from 1 to itermax during execution.

    contains

        !===========================================================================================================================
        SUBROUTINE finalizeMinimizer
        !===========================================================================================================================
            IF (ALLOCATED ( dF ) ) deallocate ( dF )
            IF (ALLOCATED ( nbd ) ) deallocate ( nbd )
            IF (ALLOCATED ( iwa ) ) deallocate ( iwa )
            IF (ALLOCATED ( ll ) ) deallocate ( ll )
            IF (ALLOCATED ( uu ) ) deallocate ( uu )
            IF (ALLOCATED ( wa ) ) deallocate ( wa )
            IF (ALLOCATED ( indicg ) ) deallocate ( indicg )
        END SUBROUTINE finalizeMinimizer
        !===========================================================================================================================

        !===========================================================================================================================
        subroutine from_cgvect_get_rho
        !===========================================================================================================================
            use system      ,only: nb_species, spacegrid, solvent, nb_species
            use quadrature  ,only: anggrid, molrotgrid
            implicit none
            integer(i2b) :: i, j, k, o, p, s, icg
            integer(i2b) :: nfft(3)

            nfft = spacegrid%n_nodes

            do concurrent (s=1:nb_species)
                if (.not. allocated (solvent(s)%rho)) then
                    allocate (solvent(s)%rho(nfft(1),nfft(2),nfft(3),anggrid%n_angles,molrotgrid%n_angles), source=0._dp ,stat=i)
                else
                    i = 0
                    solvent(s)%rho = 0._dp
                end if
            end do
            if (i /= 0) stop "ERROR in from_cgvect_get_rho"

            icg = 0
            do s =1, nb_species
                do i =1, spacegrid%n_nodes(1)
                    do j =1, spacegrid%n_nodes(2)
                        do k =1, spacegrid%n_nodes(3)
                            do o =1, anggrid%n_angles
                                do p =1, molrotgrid%n_angles
                                    icg = icg+1
                                    solvent(s)%rho(i,j,k,o,p) = cg_vect(icg)**2 *solvent(s)%rho0
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end subroutine from_cgvect_get_rho
        !===========================================================================================================================

        !===========================================================================================================================
        subroutine from_rho_get_n
        !===========================================================================================================================
            use system      ,only: nb_species, spacegrid, solvent, nb_species
            use quadrature  ,only: anggrid, molrotgrid
            implicit none
            integer(i2b) :: i,j,k,o,p,s
            integer(i2b) :: nfft(3)

            nfft = spacegrid%n_nodes

            do concurrent (s=1:nb_species)
                if (.not. allocated(solvent(s)%n)) then
                    allocate (solvent(s)%n(nfft(1),nfft(2),nfft(3)), source=0._dp ,stat=i)
                else
                    i = 0
                    solvent(s)%n = 0._dp
                end if
            end do
            if (i /= 0) stop "ERROR in allocate solvent%n in from_rho_get_n"

            do s =1, nb_species
                do i =1, spacegrid%n_nodes(1)
                    do j =1, spacegrid%n_nodes(2)
                        do k =1, spacegrid%n_nodes(3)
                            do o =1, anggrid%n_angles
                                do p =1, molrotgrid%n_angles
                                    solvent(s)%n(i,j,k) = solvent(s)%n(i,j,k)   &
                                        + solvent(s)%rho(i,j,k,o,p) * angGrid%weight(o) * molRotGrid%weight(p)
                                end do
                            end do
                        end do
                    end do
                end do
            end do

!~             do concurrent (s=1:nb_species, i=1:spacegrid%n_angles(1), j=1:spacegrid%n_angles(2), k=1:spacegrid%n_angles(3), &
!~                            o=1:anggrid%n_angles, p=1:molrotgrid%n_angles)  ! TODO WHY DOESNT IT WORK??
        end subroutine from_rho_get_n
        !===========================================================================================================================


        !===========================================================================================================================
        subroutine deallocate_solvent_rho
        !===========================================================================================================================
            use system  ,only: solvent, nb_species
            integer(i2b) :: s
            do concurrent (s=1:nb_species)  ! solvent%rho is now useless
                deallocate (solvent(s)%rho)
            end do
        end subroutine deallocate_solvent_rho
        !===========================================================================================================================

!===================================================================================================================================
        SUBROUTINE prepare_minimizer
        ! this subroutine gets the informations in input file and then allocate, prepare, compute needed data

            USE precision_kinds, ONLY: dp
            USE input, ONLY: input_int, input_dp,input_char
            USE system , ONLY: spaceGrid, nb_species
            USE quadrature, ONLY: angGrid, molRotGrid

            IMPLICIT NONE

            ncg = PRODUCT(spaceGrid%n_nodes) * angGrid%n_angles *molRotGrid%n_angles* nb_species! total number of variables to optimize
            itermax=input_int('maximum_iteration_nbr')
            epsg=input_dp('epsg')
            pgtol=input_dp('pgtol')
            minimizer_type = input_char('minimizer') ! get from input the type of minimization one wants to do
            minimizer_iter = 0
            ALLOCATE ( cg_vect(ncg), SOURCE=0._dp )
            ALLOCATE ( dF(ncg), SOURCE=0._dp )
            FF = 0._dp

            IF ( minimizer_type(1:4) == 'bfgs') THEN ! if minimizer is bfgs then init what has to be init
                mcg = input_int('number_of_memorized_Steps')
                IF (input_char('constraint')=='no_constraint') THEN
                            ALLOCATE ( nbd(ncg), SOURCE=0)!nbd(i) is : 0 for no constaint, 1 for lower bounded by ll, 3 for upper bounded by uu, 2 for lower and upper bounded
                ELSE IF (input_char('constraint')=='lower_bounded' ) THEN
                            WRITE(*,*) '***********************************************************************************'
                            WRITE(*,*) 'WARNING YOU ARE USING A LOWER_BOUND CONSTRAINED MINIMIZATION '
                            WRITE(*,*) '***********************************************************************************'
                            ALLOCATE ( nbd(ncg), SOURCE=1)!nbd(i) is : 0 for no constaint, 1 for lower bounded by ll, 3 for upper bounded by uu, 2 for lower and upper bounded
                ELSE IF (input_char('constraint')=='upper_bounded') THEN
                            WRITE(*,*) '***********************************************************************************'
                            WRITE(*,*) 'WARNING YOU ARE USING A UPPER_BOUND CONSTRAINED MINIMIZATION '
                            WRITE(*,*) '***********************************************************************************'
                            ALLOCATE ( nbd(ncg), SOURCE=3)!nbd(i) is : 0 for no constaint, 1 for lower bounded by ll, 3 for upper bounded by uu, 2 for lower and upper bounded
                ELSE IF (input_char('constraint')=='both_boundeded') THEN
                            WRITE(*,*) '***********************************************************************************'
                            WRITE(*,*) 'WARNING YOU ARE USING A LOWER AND UPPER-BOUND CONSTRAINED MINIMIZATION '
                            WRITE(*,*) '***********************************************************************************'
                            ALLOCATE ( nbd(ncg), SOURCE=2)!nbd(i) is : 0 for no constaint, 1 for lower bounded by ll, 3 for upper bounded by uu, 2 for lower and upper bounded
                END IF

                IF( input_dp('lower_bound') < 0._dp .OR. input_dp('upper_bound') < 0._dp ) THEN
                    STOP "In prepare_minimizer.f90, I suppose that lower_bound and upper_bound are positive or zero. \&
                        You found a bug. Thank you!"
                ELSE
                    ALLOCATE ( ll(ncg), SOURCE=SQRT(input_dp('lower_bound'))) ! lower bound of cg_vect
                    ALLOCATE ( uu(ncg), SOURCE=SQRT(input_dp('upper_bound') )) ! uppder bound of cg_vect
                END IF

                epsmch = EPSILON(1.0_dp)  !  Precision de la machine
                factr = epsg / epsmch ! convergence criteria over energy
                iprint = -1
                ALLOCATE ( iwa ( 3 * ncg ) ,SOURCE=0)
                ALLOCATE ( wa ( 2 * mcg * ncg + 4 * ncg + 11 * mcg **2 + 8 * mcg ) ,SOURCE=0._dp)
            END IF

        END SUBROUTINE prepare_minimizer
!===================================================================================================================================

end module minimizer
!===================================================================================================================================
