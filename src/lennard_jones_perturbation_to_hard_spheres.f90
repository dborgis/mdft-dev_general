! this SUBROUTINE computes the lennard jones perturbative contribution to the hard sphere Helmotz energy
! 201109121332 creation by Maximilien Levesque
! 201109151545 added the calculation of the perturbation potential
SUBROUTINE lennard_jones_perturbation_to_hard_spheres
use system , only : nfft1 , nfft2 , nfft3 , Lx , Ly , Lz , n_0 , radius , sig_solv , eps_solv , v_perturbation_k
use quadrature , only : angGrid
USE minimizer, ONLY: cg_vect , dF , FF
USE precision_kinds , only : dp , i2b
use constants , only : fourpi , twopi
use fft , only : fftw3
IMPLICIT NONE
real(dp), allocatable , dimension ( : , : , : ) :: rho_n ! local density
complex(dp), allocatable , dimension ( : , : , : ) :: rho_k ! fourier transformed rho_n
real(dp):: local_density ! dummy for rho_n
integer(i2b):: icg , i , j , k , o , l , m , n , m1 , m2 , m3 ! dummy
integer(i2b):: nf1 , nf2 , nf3 ! dummy for nfft1 /2 , nfft2 /2, nfft3 /2
real(dp):: deltav ! integration factor = Lx*Ly*Lz/(nfft1*nfft2*nfft3)
real(dp):: Nk ! total number of k points
real(dp):: kx2 , ky2 , kz2 , norm_k ! related to number k
complex(dp), allocatable , dimension ( : , : , : ) :: vk ! fourier transform of the lennard jones perturbation (WCA)
real(dp), allocatable , dimension ( : , : , : ) :: v_perturbation_r ! vk in real space
real(dp):: Fperturbation ! what we want : the perturbative contribution of the lennard jones attractive potential tail
real(dp):: time0 , time1 ! timers
real(dp):: twopiolx , twopioly , twopiolz ! dummy for speeding up loops
real(dp):: potential ! dFp / drho at rho=rho_0 in order the grand potential to be zero at rho = rho_0
real(dp):: nb_molecule ! total number of hard spheres ie integral of density over all space
! init
call cpu_time ( time0 )
nf1 = nfft1 / 2
nf2 = nfft2 / 2
nf3 = nfft3 / 2
deltav = Lx * Ly * Lz / real ( nfft1 * nfft2 * nfft3 , dp ) ! integration factor
! put result from last minimization step as density from which one computes energy and gradients
allocate ( rho_n ( nfft1 , nfft2 , nfft3 ) )
icg = 0
nb_molecule = 0.0_dp
do i = 1 , nfft1
  do j = 1 , nfft2
    do k = 1 , nfft3
      local_density = 0.0_dp
      do o = 1, angGrid%n_angles ! angGrid%n_angles=1
        icg = icg + 1
        local_density = local_density + angGrid%weight (o) * cg_vect (icg) ** 2
      END DO
      ! correct by fourpi as the integral over all orientations o is 4pi
      local_density = local_density / fourpi
      ! at the same time integrate rho_n in order to count the total number of implicit molecules. here we forget the integration factor = n_0 * deltav
      rho_n ( i , j , k ) = local_density
      nb_molecule = nb_molecule + local_density
    END DO
  END DO
END DO
nb_molecule = nb_molecule * n_0 * deltav ! normalization
! total number of k points needed for inverse fft normalization
Nk = real ( nfft1 * nfft2 * nfft3 , dp )
! fourier transform the density rho_n => rho_k
fftw3%in_forward = rho_n
call dfftw_execute ( fftw3%plan_forward )
allocate ( rho_k ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
rho_k = fftw3%out_forward
! compute lennard jones perturbation in k space
! if v_perturbation_k doesn't exist, then compute it and put it in Vk. 
if ( .not. allocated ( v_perturbation_k ) ) then
  allocate ( v_perturbation_k ( nfft1 / 2 + 1 , nfft2 , nfft3 ) )
  twopiOLx = twopi / Lx ! dummy for speeding up the loops
  twopiOLy = twopi / Ly ! dummy for speeding up the loops
  twopiOLz = twopi / Lz ! dummy for speeding up the loops
  do l = 1 , nf1 + 1
    m1 = l - 1
    if ( l > nf1 ) m1 = l - 1 - nfft1
    kx2 = ( twopiOLx * real ( m1 , dp ) ) ** 2
    do m = 1 , nfft2
      m2 = m - 1
      if ( m > nf2 ) m2 = m - 1 - nfft2
      ky2 = ( twopiOLy * real ( m2 , dp ) ) ** 2
  
      do n = 1 , nfft3
        m3 = n - 1
        if ( n > nf3 ) m3 = n - 1 - nfft3
        kz2 = ( twopiOLz * real ( m3 , dp ) ) ** 2
        norm_k = sqrt ( kx2 + ky2 + kz2 )
        v_perturbation_k ( l , m , n ) = vlj_wca_k ( norm_k , sig_solv(1) , eps_solv(1) ) ! in kJ/mol
      END DO
 
    END DO
  END DO
END IF
! put the backup in what we use (perhaps redondant)
allocate ( Vk ( nfft1 / 2 + 1 , nfft2 , nfft3 ) ) ! Vk is the one we use in this routine. It may be saved in v_perturbation_k in order not to compute it each time.
vk = v_perturbation_k
! once equation written, dFp / drho at rho = rho_0 is shown to be equal to rho_0 * Vk(k=0)
potential = n_0 * real ( Vk ( 1 , 1 , 1 ) )
! FFT-1 of perturbation
fftw3%in_backward = Vk * rho_k
deallocate ( Vk )
deallocate ( rho_k )
call dfftw_execute (fftw3%plan_backward)
allocate ( v_perturbation_r ( nfft1 , nfft2 , nfft3 ) )
v_perturbation_r = fftw3%out_backward / Nk
! Compute the perturbative energy
Fperturbation = 0.0_dp
do i = 1 , nfft1
  do j = 1 , nfft2
    do k = 1 , nfft3
      do o = 1, angGrid%n_angles ! angGrid%n_angles=1
        Fperturbation = Fperturbation + rho_n ( i , j , k ) * v_perturbation_r ( i , j , k )
      END DO
    END DO
  END DO
END DO
! normalize
Fperturbation = Fperturbation * 0.5_dp * deltav ! * n_0
deallocate ( rho_n )
! add perturbation energy to total energy
write(*,*)'nb_molecule = ',nb_molecule
write(*,*)'mu_p = ', potential
write(*,*)'- nb_molecule * potential = ', - nb_molecule * potential
FF = FF + Fperturbation - nb_molecule * potential
! gradient
icg = 0
do i = 1 , nfft1
  do j = 1 , nfft2
    do k = 1 , nfft3
      do o = 1 , angGrid%n_angles
        icg = icg + 1
        dF ( icg ) = dF ( icg ) + 2.0_dp * cg_vect ( icg ) * DeltaV * v_perturbation_r ( i , j , k ) ! 2011 09 18 23h16 deleted *n_0
      END DO
    END DO
  END DO
END DO
deallocate ( v_perturbation_r )
!> Close timer
call cpu_time ( time1 )
! warn user
write (*,*) 'Fperturbati = ' , Fperturbation , 'computed in (sec)' , time1 - time0

contains
  ! Here we compute Uperturbation in kspace. it's an integration we do numericaly

  function vlj_wca_k ( k , sigma_lj , epsilon_lj )
  USE precision_kinds , only : dp , i2b
  use constants , only : fourpi
  IMPLICIT NONE
  complex(dp):: vlj_wca_k ! which computes the reciprocal value of the potential 'vk'
  real(dp), intent(in) :: k ! one gives the k point 'k' to eat to the routine
  real(dp), intent(in) :: sigma_lj , epsilon_lj ! lennard jones parameters in Angstroms and KJ/mol
  real(dp):: cutoff ! cut off under which U(r) = constant = -epsilon_lj and after which U(r)=Vlj(r)
  real(dp):: dx ! width of the integration step
  real(dp):: borne_sup , borne_inf !sup and inf limits of the integration
  integer(i2b), parameter :: nstep = 1000 ! nb of integration steps ! TODO CHECK EFFECT
  real(dp):: ri !integrants
  real(dp):: sigmaori6 !dummy for (sigma_lj/ri)**6
  integer(i2b):: i !dummy for loop
  !the integration is theoreticaly between 2^(1/6)*sigma and infinity
  !we do it between sigma_lj and N*sigma for a beginning
  !perhaps 10*sigma is more than enough
  !we use 10**3 steps for the integration.
  borne_sup = 10.0_dp * sigma_lj ! TODO 40 is better and converged. check later and do analyticaly
  borne_inf = 0.0_dp
  ! TODO here one integrates over the whole range of r. one could calculate the first part (U=-epsilon_lj) analyticaly and thus speed up everything by a factor of d_wca/borne_sup
  dx = ( borne_sup - borne_inf ) / real ( nstep , dp )
  ! compute the value after which U = Vlj, which is the value of x for which Vlj is minimum, thus 2^(1/6)sigma
  cutoff = 2.0_dp ** ( 1.0_dp / 6.0_dp ) * sigma_lj
  ! init vk
  vlj_wca_k = ( 0.0_dp , 0.0_dp )
  ! integrate v(r) over all r in R^3
  do i = 1 , nstep
    ri = borne_inf + real ( i - 1 , dp ) * dx
    if ( ri < cutoff ) then ! if under cutoff
      if ( k > 0.0_dp ) then
        vlj_wca_k = vlj_wca_k + cmplx ( - epsilon_lj * ri * sin ( k * ri ) / k , 0.0_dp )
      ELSE ! if k = 0 then lim sin(k*r)/k = r
        vlj_wca_k = vlj_wca_k + cmplx ( - epsilon_lj * ri ** 2 , 0.0_dp )
      END IF
    ELSE ! if after cutoff
      sigmaori6 = ( sigma_lj / ri ) ** 6
      if ( k > 0.0_dp ) then
        vlj_wca_k = vlj_wca_k + cmplx ( 4.0_dp * epsilon_lj * ri * ( sigmaori6 ** 2 - sigmaori6 ) * sin ( k * ri ) / k , 0.0_dp )
      ELSE ! if k = 0 then lim sin(k*r)/k = r 
        vlj_wca_k = vlj_wca_k + cmplx ( 4.0_dp * epsilon_lj * ri * ( sigmaori6 ** 2 - sigmaori6 ) * ri , 0.0_dp ) ! lim ( sin(kr) / k , k->0 ) = r
      END IF
    END IF
  END DO
  ! integration factors
  vlj_wca_k = vlj_wca_k * cmplx ( dx * fourpi , 0.0_dp )
  end function vlj_wca_k
  
END SUBROUTINE lennard_jones_perturbation_to_hard_spheres
