! this SUBROUTINE init the density as a function of position and orientation
! it should be initiated to exp(-beta*Vext_total) but
! Vext_q is the electrostatic part of Vext_total, and is pathologic (it sometimes diverges)
! we thus init the density not using vext, but Vext_total - Vext_q
subroutine init_density

  use precision_kinds, only : dp, i2b
  use system, only: thermocond, nb_species
  use module_grid, only: grid
  use minimizer, only: cg_vect_new
  use external_potential, only: Vext_total, Vext_q
  use input, only: input_log
  use mathematica, only: chop

  implicit none

  real(dp) :: local_density0 !> @var local_density0 is the density at a space and angular grid point
  integer(i2b) :: i, j, k, icg, io
  real(dp) :: Vext_total_local , Vext_total_local_min_Vext_q ! dummy
  integer(i2b) :: species, ios, s
  logical :: exists
  integer :: nfft1, nfft2, nfft3

  nfft1 = grid%nx
  nfft2 = grid%ny
  nfft3 = grid%nz

  ! Read the density from a previous run
  if (input_log('reuse_density')) then
    INQUIRE (file='input/density.bin.in', EXIST=exists)
        IF ( .NOT. exists) STOP "input/density.bin.in not found"
    OPEN (10, file = 'input/density.bin.in' , form = 'unformatted' , iostat=ios, status='OLD' )
        IF ( ios /= 0 ) then
            print *, 'problem while opening input/density.bin.in. bug at init_density.f90'
            stop
        END IF
        READ ( 10, iostat=ios ) cg_vect_new
            IF ( ios<0 ) THEN
                STOP "input/density.bin.in is empty"
            ELSE IF ( ios>0 ) THEN
                STOP "problem while trying to read cg_vect_new in input/density.bin.in"
            END IF
        PRINT*, '*** RESTART ***'
    CLOSE (10)
    OPEN (10, FILE = 'output/density.bin.in.out', FORM = 'unformatted')
        WRITE ( 10 ) cg_vect_new
    CLOSE (10)
    RETURN
  end if


! We minimize with respect to sqrt(density) and not directly with respect to the density.
! This allows the variable (sqrt(density)) to be strictly positive.
! icg is the index of the cg_vect_new 1dimensional array. it obliges us to loop over x then y then z then omega

! only Vext_total is used in the functional. We may deallocate it now.
!if ( allocated ( Vext_q ) ) deallocate ( Vext_q )

do s = 1, nb_species
    do io = 1, GRID%no
        do k = 1, GRID%nz
            do j = 1, GRID%ny
                do i = 1, GRID%nx

                    Vext_total_local = Vext_total (i,j,k,io,s)
                    if ( allocated( Vext_q) ) then
                        Vext_total_local_min_Vext_q = Vext_total_local - Vext_q (i,j,k,io,s)
                    else
                        Vext_total_local_min_Vext_q = Vext_total_local
                    end if

                    if ( Vext_total_local >= 100.0_dp .or. Vext_total_local_min_Vext_q >= 100.0_dp ) then
                        local_density0 = 0.0_dp!tiny ( 0.0_dp ) ! Don't put 0 as it induces problems in the calculations of log(density) in ideal part of F.
                    ELSE
                        local_density0 = chop( EXP(- thermocond%beta * Vext_total_local_min_Vext_q ) )
                    END IF
                    ! put result in cg_vect_new. note that we do prefer to minimize sqrt(density) than the density in order to avoid sign problems
                    cg_vect_new(i,j,k,io,s) = chop( sqrt( local_density0) )

                end do
            end do
        end do
    end do
end do

END SUBROUTINE init_density
