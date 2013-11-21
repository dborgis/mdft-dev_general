SUBROUTINE energy_hard_sphere_fmt
USE precision_kinds,only : dp , i2b
use system,only : nfft1 , nfft2 , nfft3 , weight_function_0_k , weight_function_1_k , weight_function_2_k , &
                    weight_function_3_k , deltav , Fexc_0 , kBT , muexc_0 , n_0 , nb_species , n_0_multispec , &
                    muexc_0_multispec , Fexc_0_multispec , mole_fraction , rho_0_multispec
use quadrature,only : sym_order , angGrid, molRotGrid
USE minimizer, ONLY: cg_vect , FF , dF
use constants,only : pi , FourPi , twopi
use fft,only : fftw3
use input, only : input_line
IMPLICIT NONE
integer(i2b):: icg , i , j , k , o , p ! dummy
integer(i2b):: species ! dummy between 1 and nb_species
real(dp):: Nk ! Total number of k points = nfft1*nfft2*nfft3
real(dp):: Fint !> @var Internal part of free energy
real(dp):: local_density, psi !> Dummy , psi is sqrt(local_density)
real(dp), allocatable , dimension ( : ) :: nb_molecules ! temp before full implementation of multispecies
real(dp), allocatable , dimension ( : , : , : , : ) :: rho ! density per angle (recall : rho_0 = n_0 / 4pi ) ! x y z nb_species
complex(dp), allocatable , dimension ( : , : , : , :) :: rho_k !> @var Density in k space
real(dp):: time0 , time1 ! time stamps
real(dp):: w0 , w1 , w2 , w3 , omw3 ! local 1-w3 ! dummy 
real(dp):: F_HS ! Perkus Yevick (or Carnahan Starling) expression for excess energy
real(dp)   , allocatable , dimension ( : , : , : ) :: weighted_density_0 , weighted_density_1 , weighted_density_2 ,&
                                                                  weighted_density_3 ! weighted densities
real(dp)   , allocatable , dimension ( : , : , : ) :: one_min_weighted_density_3 ! dummy for 1 - weighted density 3
real(dp)   , allocatable , dimension ( : , : , : ) :: dFHS_0 , dFHS_1 , dFHS_2 , dFHS_3
complex(dp), allocatable , dimension ( : , : , : ) :: dFHS_0_k , dFHS_1_k , dFHS_2_k , dFHS_3_k
complex(dp), allocatable , dimension ( : , : , : , : ) :: dFex_k ! gradient in k space
real(dp)   , allocatable , dimension ( : , : , : , : ) :: dFex ! gradient in real space
 character ( len = 2 ) :: hs_functional ! hard sphere functional = PY for Percus-Yevick or CS for Carnahan-Starling
! parameters for easy reading and avoiding 'magic numbers'
real(dp), parameter :: inv8pi  = 1.0_dp / (  8.0_dp * pi ) ! dummy constant
real(dp), parameter :: inv12pi = 1.0_dp / ( 12.0_dp * pi ) ! dummy constant
real(dp), parameter :: inv18pi = 1.0_dp / ( 18.0_dp * pi ) ! dummy constant
real(dp), parameter :: inv24pi = 1.0_dp / ( 24.0_dp * pi ) ! dummy constant
real(dp), parameter :: inv36pi = 1.0_dp / ( 36.0_dp * pi ) ! dummy constant
! init timer
call cpu_time ( time0 )
! init 
Nk = real ( nfft1 * nfft2 * nfft3 , dp ) ! total number of k points needed for inverse fft normalization
! put result from last minimization step as density from which one computes energy and gradients
allocate ( rho ( nfft1 , nfft2 , nfft3 , nb_species ) )
! get rho ( \vec{r} )
allocate ( nb_molecules ( nb_species ) )
nb_molecules = 0.0_dp
icg = 0
do species = 1 , nb_species
  do i = 1 , nfft1
    do j = 1 , nfft2
      do k = 1 , nfft3
        local_density = 0.0_dp
        do o = 1, angGrid%n_angles
          do p=1, molRotGrid%n_angles
          icg = icg + 1
          local_density = local_density + angGrid%weight (o) * cg_vect (icg) ** 2*molRotGrid%weight(p)
           END DO          
        END DO
        ! correct by 8*pi²/n as the integral over all orientations o and psi is 4pi and 2pi/n
        local_density = local_density*real(sym_order, dp) / (fourpi*twopi)
        ! at the same time integrate rho in order to count the total number of implicit molecules. here we forget the integration factor = n_0 * deltav
        rho ( i , j , k , species ) = local_density
      END DO
    END DO
  END DO
END DO
! total number of molecules of each species
DO CONCURRENT ( species = 1 : nb_species )
  nb_molecules ( species ) = sum ( rho ( : , : , : , species ) ) * n_0_multispec ( species ) * mole_fraction ( species ) * deltav
END DO
! tell user about the number of molecule of each species in the supercell
do species = 1 , nb_species
  write (*,*) 'nb_molecule (' , species , ') = ' , nb_molecules ( species )
END DO
write (*,*) 'total number of molecules = ' , sum ( nb_molecules )
! fourier transform the density rho => rho_k
allocate ( rho_k ( nfft1 / 2 + 1 , nfft2 , nfft3 , nb_species ) )
do species = 1 , nb_species
  fftw3%in_forward = rho ( : , : , : , species )
  call dfftw_execute ( fftw3%plan_forward )
  rho_k ( : , : , : , species ) = fftw3%out_forward * n_0_multispec ( species ) * mole_fraction ( species )
END DO
deallocate ( rho )
! inverse fourier transform the weighted densities
! allocate the arrays for receiving the FFT-1
allocate ( weighted_density_0  ( nfft1 , nfft2 , nfft3 ) )
weighted_density_0 = 0.0_dp
allocate ( weighted_density_1  ( nfft1 , nfft2 , nfft3 ) )
weighted_density_1 = 0.0_dp
allocate ( weighted_density_2  ( nfft1 , nfft2 , nfft3 ) )
weighted_density_2 = 0.0_dp
allocate ( weighted_density_3  ( nfft1 , nfft2 , nfft3 ) )
weighted_density_3 = 0.0_dp
do species = 1 , nb_species
    fftw3%in_backward = weight_function_0_k ( : , : , : , species ) * rho_k ( : , : , : , species ) ! = weighted density in kspace _0
    call dfftw_execute ( fftw3%plan_backward )
    weighted_density_0 = weighted_density_0 + fftw3%out_backward / Nk
    fftw3%in_backward = weight_function_1_k ( : , : , : , species ) * rho_k ( : , : , : , species ) ! = weighted density kspace _1
    call dfftw_execute ( fftw3%plan_backward )
    weighted_density_1 = weighted_density_1 + fftw3%out_backward / Nk
    fftw3%in_backward = weight_function_2_k ( : , : , : , species ) * rho_k ( : , : , : , species ) ! = weighted density kspace _2
    call dfftw_execute ( fftw3%plan_backward )
    weighted_density_2 = weighted_density_2 + fftw3%out_backward / Nk
    fftw3%in_backward = weight_function_3_k ( : , : , : , species ) * rho_k ( : , : , : , species ) ! = weighted density kspace _3
    call dfftw_execute ( fftw3%plan_backward )
    weighted_density_3 = weighted_density_3 + fftw3%out_backward / Nk
END DO

! check if the hard sphere functional is Percus-Yevick or Carnahan-Starling
! Get the free energy functional that should be used. For now Percus Yevick and Carnahan Starling only. May be expanded.
do i = 1 , size ( input_line )
  j = len ( 'hs_functional' )
  if ( input_line (i) (1:j) == 'hs_functional' ) read ( input_line (i) (j+4:j+5) , * ) hs_functional
END DO
! compute free intrinsic energy
! init Fint
Fint = 0.0_dp
do k = 1 , nfft3 ! please pay attention to inner / outer loop.
  do j = 1 , nfft2
    do i = 1 , nfft1
      w0 = weighted_density_0 ( i , j , k )
      w1 = weighted_density_1 ( i , j , k )
      w2 = weighted_density_2 ( i , j , k )
      w3 = weighted_density_3 ( i , j , k )
      ! nowhere should w3 be lower or equal than 1
      omw3 = 1.0_dp - w3
      if ( omw3 <= 0.0_dp ) then
        call error_message_energy_hard_sphere_fmt ( i , j , k , w0 , w1 , w2 , w3 )
        call process_output ! process output so that we know where we are !
      END IF
      ! here is the Perkus Yevick excess free energy
      if ( hs_functional == 'PY' .or. hs_functional == 'py' ) then
        F_HS = - w0 * log ( omw3 ) + w1 * w2 / omw3 + inv24pi * w2 ** 3 / omw3 ** 2
      ELSE IF ( hs_functional == 'CS' .or. hs_functional == 'cs' ) then
        F_HS = ( inv36pi * w2 ** 3 / w3 ** 2 - w0 ) * log ( omw3 ) + w1 * w2 / omw3 &
                 + inv36pi * w2 ** 3 / ( omw3 ** 2 * w3 )
      END IF
  
      Fint = Fint + F_HS
    END DO
  END DO
END DO
Fint = Fint * kBT * DeltaV - sum ( muexc_0_multispec * nb_molecules ) - sum ( Fexc_0_multispec * mole_fraction )
FF = FF + Fint ! FF is used in BFGS algorithm with CGvect(icg) and dF(icg)
! gradients
! dFHS_i and weighted_density_j are arrays of dimension (nfft1,nfft2,nfft3)
! one_min_weighted_density_3 is dummy for speeding up and simplifying while using unnecessary memory.
allocate ( one_min_weighted_density_3 ( nfft1 , nfft2 , nfft3 ) )
one_min_weighted_density_3 = 1.0_dp - weighted_density_3
! excess free energy per atom derived wrt weighted density 0 (eq. A5 in Sears2003)
allocate ( dFHS_0 ( nfft1 , nfft2 , nfft3 ) )
allocate ( dFHS_1 ( nfft1 , nfft2 , nfft3 ) )
allocate ( dFHS_2 ( nfft1 , nfft2 , nfft3 ) )
allocate ( dFHS_3 ( nfft1 , nfft2 , nfft3 ) )
! Perkus Yevick
if ( hs_functional == 'PY' .or. hs_functional == 'py' ) then
  dFHS_0 = - log ( one_min_weighted_density_3 )
  dFHS_1 = weighted_density_2 / one_min_weighted_density_3
  dFHS_2 = weighted_density_1 / one_min_weighted_density_3 + inv8pi * dFHS_1 ** 2
  dFHS_3 = ( weighted_density_0 + weighted_density_1 * dFHS_1 ) / one_min_weighted_density_3 + inv12pi * dFHS_1 ** 3
! Carnahan Starling
ELSE IF ( hs_functional == 'CS' .or. hs_functional == 'cs' ) then
  dFHS_0 = - log ( one_min_weighted_density_3 )
  dFHS_1 = weighted_density_2 / one_min_weighted_density_3
  dFHS_2 = - inv12pi * ( weighted_density_2 / weighted_density_3 ) ** 2 * dFHS_0 &
           + weighted_density_1 / one_min_weighted_density_3 + inv12pi * dFHS_1 ** 2 / weighted_density_3
  dFHS_3 = inv18pi * dFHS_0 * ( weighted_density_2 / weighted_density_3 ) ** 3 &
           - ( inv36pi * weighted_density_2 ** 3 / weighted_density_3 ** 2 - weighted_density_0 ) / one_min_weighted_density_3 &
           + weighted_density_1 * dFHS_1 / one_min_weighted_density_3 + inv36pi * &
           weighted_density_2 ** 3 / weighted_density_3 ** 2 * &
           ( 3.0_dp * weighted_density_3 - 1.0_dp ) / one_min_weighted_density_3 ** 3
END IF
! deallocate weighted_densities
deallocate ( weighted_density_0 )
deallocate ( weighted_density_1 )
deallocate ( weighted_density_2 )
deallocate ( weighted_density_3 )
deallocate ( one_min_weighted_density_3 )
! compute gradients in k space
allocate ( dFHS_0_k ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
allocate ( dFHS_1_k ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
allocate ( dFHS_2_k ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
allocate ( dFHS_3_k ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
! FFT dFHS for computing convolution
fftw3%in_forward = dFHS_0
call dfftw_execute ( fftw3%plan_forward )
dFHS_0_k = fftw3%out_forward
fftw3%in_forward = dFHS_1
call dfftw_execute ( fftw3%plan_forward )
dFHS_1_k = fftw3%out_forward
fftw3%in_forward = dFHS_2
call dfftw_execute ( fftw3%plan_forward )
dFHS_2_k = fftw3%out_forward
fftw3%in_forward = dFHS_3
call dfftw_execute ( fftw3%plan_forward )
dFHS_3_k = fftw3%out_forward
! deallocate useless
deallocate ( dFHS_0 )
deallocate ( dFHS_1 )
deallocate ( dFHS_2 )
deallocate ( dFHS_3 )
! compute final gradient in k-space
allocate ( dFex_k ( nfft1 / 2 + 1 , nfft2 , nfft3 , nb_species ) ) ! FFT of dFex (in fact dFex_k is known from which FFT-1 gives dFex)
do species = 1 , nb_species
  dFex_k ( :,:,:, species ) = &
         dFHS_0_k * weight_function_0_k ( :,:,:, species ) &
       + dFHS_1_k * weight_function_1_k ( :,:,:, species ) &
       + dFHS_2_k * weight_function_2_k ( :,:,:, species ) &
       + dFHS_3_k * weight_function_3_k ( :,:,:, species )
END DO
!> Deallocate useless
deallocate ( dFHS_0_k )
deallocate ( dFHS_1_k )
deallocate ( dFHS_2_k )
deallocate ( dFHS_3_k )
! inverse fourier transform gradient
allocate ( dFex ( nfft1 , nfft2 , nfft3 , nb_species ) )
DO species = 1, nb_species
  fftw3%in_backward = dFex_k ( : , : , : , species )
  call dfftw_execute ( fftw3%plan_backward )
  dFex ( : , : , : , species ) = fftw3%out_backward / Nk
END DO
deallocate ( dFex_k )
! transfer in rank 1 vector dF
icg = 0
do species = 1 , nb_species
  do i = 1, nfft1
    do j = 1, nfft2
      do k = 1, nfft3
        do o = 1, angGrid%n_angles
          do p=1, molRotGrid%n_angles
          icg = icg + 1
          psi = cg_vect ( icg )
! AHAAAAAAAAAAAAATTENTION ICI LE 12 SEPTEMBRE 2011 JE ME RENDS COMPTE QU iL Y A PEUT ETRE UN FACTEUR WEIGHT(o) QUI MANQUE, CA DEPEND DE LA DEF DE N0 par RAPPORT à rho_0
          dF ( icg ) = dF ( icg ) &
+ 2.0_dp * psi * rho_0_multispec ( species ) * deltav * ( kBT * dFex ( i , j , k , species ) - muexc_0_multispec ( species ) )*&
angGrid%weight(o)*molRotGrid%weight(p)
 ! ATTENTION J'AI MULTIPLIE PAR WEIGHT(O) QD PASSAGE A angGrid%n_angles /=1 LE 18 JUILLET 2011
         END DO  !p
  
        END DO  !o
      END DO  !i
    END DO  !j
  END DO   !k
END DO   !species
! deallocate useless gradient
deallocate ( dFex )
 !stop timer
call cpu_time ( time1 )
! print info for user
write (*,*) 'Fexc fmt    = ' , Fint , 'computed in (sec)' , time1 - time0

CONTAINS

    ! this SUBROUTINE prints error message related to SUBROUTINE excess_cs_hard_sphere
    ! it may stop program execution depending on the error.
    SUBROUTINE error_message_energy_hard_sphere_fmt ( i , j , k , w0 , w1 , w2 , w3 )
        USE precision_kinds,only : i2b , dp
        IMPLICIT NONE
        integer(i2b), intent(in) :: i , j , k
        real(dp), intent(in) :: w0 , w1 , w2 , w3
        write (*,*) 'i , j , k = ' , i , j , k
        write (*,*) 'log (1-w3<=0) in energy_hard_sphere_fmt.f90. Critical stop'
        write (*,*) 'w0 = ' , w0
        write (*,*) 'w1 = ' , w1
        write (*,*) 'w2 = ' , w2
        write (*,*) 'w3 = ' , w3
        stop
    END SUBROUTINE error_message_energy_hard_sphere_fmt

END SUBROUTINE energy_hard_sphere_fmt
