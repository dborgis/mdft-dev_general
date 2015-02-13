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
    real(dp):: epsg
    integer(i2b) :: iter ! iteration number of the minimizer. Goes from 1 to itermax during execution.
    integer(i2b) :: itermax ! maximum number of iterations

    contains

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
            ! pgtol=input_dp('pgtol')
            minimizer_type = input_char('minimizer') ! get from input the type of minimization one wants to do
            ALLOCATE ( cg_vect(ncg), SOURCE=0._dp )
            ALLOCATE ( dF(ncg), SOURCE=0._dp )
            FF = 0._dp

        END SUBROUTINE prepare_minimizer
!===================================================================================================================================

end module minimizer
!===================================================================================================================================
