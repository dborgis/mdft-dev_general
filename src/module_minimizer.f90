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
            if (i /= 0) stop "ERROR in allocate solvent%rho in energy_nn_cs.f90 > from_cgvect_get_rho"

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
            if (i /= 0) stop "ERROR in allocate solvent%n in energy_nn_cs.f90 > from_rho_get_n"

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
        
end module minimizer
