!> System variables.

module system

use precision_kinds , only:i2b,dp

implicit none

integer(kind=i2b) :: nb_species

real(kind=dp) :: Lx , Ly , Lz ! Taille de la boite en Ang

real(kind=dp) :: DeltaV ! volume elementaire=Lx*Ly*Lz/(nfft1*nfft2*nfft3)

real(kind=dp) :: deltax, deltay, deltaz ! Lx/nfft1, Ly/nfft2, Lz/nfft3

integer(kind=i2b) :: nb_solute_sites, nb_solvent_sites  ! nombre de site pour le solute et pour le solvent

real(kind=dp), allocatable, dimension(:) :: x_mol, y_mol, z_mol ! positions des sites du solute dans la boite (repere absolu)

real(kind=dp), allocatable, dimension(:) :: chg_mol, sig_mol, eps_mol  ! charge partielle et parametres LJ pour chaque site du solute (kJ/mol)

real(kind=dp), allocatable, dimension(:) :: x_solv, y_solv, z_solv  ! positions des sites du solvent dans le repere propre ede la molecule

real(kind=dp), allocatable, dimension(:) :: chg_solv, sig_solv, eps_solv  ! charge partielle et parametres LJ pour chaque site du solvant (kJ/mol)

integer(kind=i2b), allocatable, dimension(:) :: id_solv, id_mol ! atom type of each solvent site and solute site for instance id_solv(1)=1, id_solv(2)=2 and id_solv(3)=2 for OH2

real(kind=dp) :: Rc ! Taille de la charge (Ang)

real(kind=dp) :: TEMP  ! temperature du systeme lue dans dft.in 'temperature : XXXX'

real(kind=dp) :: kBT , beta

! Grille pour FFT

integer(kind=i2b) :: nfft ! nbr de point par unite de longueur. lu dans input/dft.in. nfft1, nfft2 et nfft3 sont deduits de nfft.

integer(kind=i2b) :: nfft1,nfft2,nfft3 ! nbr de points sur chaque dimensions de la grille. deduits de nfft lu dans dft.in

integer(kind=i2b) :: nb_legendre ! ordre pour integration Gauss-legendre

real(kind=dp) :: delta_r


! Grille angulaire

integer(kind=i2b) :: nb_psi ! nb angle psi

integer(kind=i2b) :: nb_omega


!

real(kind=dp) :: n_0 , rho_0   ! Densite du fluide homogene en part/A3 et incluant orientation

real ( kind = dp ) , allocatable , dimension ( : ) :: n_0_multispec , rho_0_multispec ! here are the equivalent of n_0 and rho_0 in multispecies case

real(kind=dp), allocatable, dimension(:) :: c_s, c_delta, c_d !> Direct correlation functions of various rotational invariants

real(kind=dp), allocatable, dimension(:) ::chi_L, chi_T

real(kind=dp) , allocatable , dimension (:) :: c_s_hs ! c(2)(k) of a hard sphere

complex(kind=dp),allocatable,dimension(:,:,:) :: Vk !>@var perturabtion in kspace

real(kind=dp) :: delta_k ! distance between two k points in cs.in, cdelta.in, cd.in

integer(kind=i2b) :: nb_k ! nb of k points in cs.in, cdelta.in, cd.in



! Electrostatics
real ( kind=dp ) , allocatable , dimension ( : , : , : , : , : ) :: wigma ! function of n1 , n2 , n3 , omega , psi

! charge factor & molecule polarization factor
complex ( kind=dp ) , allocatable , dimension ( : , : , : ,:,:,: ) :: sigma_k


real ( kind = dp ) , allocatable , dimension ( : , : , : ) :: rho_c ! charge density at each grid node

complex ( kind = dp ) , allocatable, dimension ( : , : , : ) :: rho_c_k, rho_c_k_myway


!Polarization

complex (kind=dp) , allocatable,dimension (:,:,:,:) :: pola_tot_x_k , pola_tot_y_k , pola_tot_z_k

complex ( kind=dp ) , allocatable , dimension ( : , : , : ,:,: , :) ::  molec_polarx_k,molec_polary_k,molec_polarz_k 


!> 3body parameters
!NEW
!!a=100.0_dp
!!gam=2.0_dp/3.0_dp
!!rmin1=1.5_dp
!!rsw1=2.0_dp
!!rmax1=0.5_dp*(sig_mol(1)+sig_solv(1))+1.9_dp
!!rmin2=2.25_dp
!!rsw2=2.5_dp
!!rmax2=4.75_dp
!!lambda_fs=100.0_dp

!OLD
!gamma  :  2.9
!rsw1  :  2.0
!rsw2  :  2.75
!rmin1  :  1.5
!rmax1  :  5.0
!rmin2  :  2.25
!rmax2  :  5.00
!lambda  :  125.00
!lambda_fs  :  0.0
!a3  :  0.0
!b3  :  3.0







real(kind=dp), allocatable, dimension(:,:,:) ::  V_int

!real ( kind = dp ) , allocatable , dimension ( : , : , : , : ) :: rho_n ! density per angle (=n/4pi) of each species at each node (nfft1,nfft2,nfft3,nb_species)
! for maximum efficiency, rho_species(nfft1,nfft2,nfft3,nb_species) should be used in loops where species is the outer loop (varies slowliest)

!do species = 1 , nb_species
!do z
!do y 
!do x
!density(x,y,z,species)




!> Hard spheres

! TODO this should be put in another module

logical(4) :: HS

complex ( kind = dp ) , allocatable , dimension ( : , : , : , : ) :: weight_function_1_k

complex ( kind = dp ) , allocatable , dimension ( : , : , : , : ) :: weight_function_2_k

complex ( kind = dp ) , allocatable , dimension ( : , : , : , : ) :: weight_function_3_k

complex ( kind = dp ) , allocatable , dimension ( : , : , : , : ) :: weight_function_0_k

real(kind=dp), allocatable, dimension(:) :: radius

real(kind=dp) :: eta, muexc_0, Fexc_0

!real ( kind = dp ) , allocatable , dimension ( : ) :: eta_multispec ! packing fraction = 4/3 * pi * n_0 * radius **3 = pi/6 * n_0 * diameter**3

real ( kind = dp ) , allocatable , dimension ( : ) :: Fexc_0_multispec ! 

real ( kind = dp ) , allocatable , dimension ( : ) :: muexc_0_multispec ! temp var during implementation of multispecies

complex ( kind = dp ) , allocatable , dimension ( : , : , : ) :: v_perturbation_k ! fourier transform of the lennard jones perturbation (WCA)

real ( kind = dp ) , allocatable , dimension ( : ) :: mole_fraction ! mole fraction of each species "x_i"



!> for xsf printing

integer(kind=i2b), allocatable, dimension (:) :: atomic_nbr

real(kind=dp), allocatable, dimension (:,:,:) :: V_coulomb ! nfft1 nfft2 nfft3 nb_omega

!Three body terms

real(kind = dp), allocatable, dimension(:) :: lambda1_mol, lambda2_mol 

contains

  subroutine deallocate_everything_system

  implicit none

  if ( allocated ( x_mol ) ) deallocate ( x_mol )

  if ( allocated ( y_mol ) ) deallocate ( y_mol )

  if ( allocated ( z_mol ) ) deallocate ( z_mol )

  if ( allocated ( chg_mol ) ) deallocate ( chg_mol )

  if ( allocated ( sig_mol ) ) deallocate ( sig_mol )

  if ( allocated ( eps_mol ) ) deallocate ( eps_mol )

  if ( allocated ( x_solv ) ) deallocate ( x_solv )

  if ( allocated ( y_solv ) ) deallocate ( y_solv )

  if ( allocated ( z_solv ) ) deallocate ( z_solv )

  if ( allocated ( chg_solv ) ) deallocate ( chg_solv )

  if ( allocated ( sig_solv ) ) deallocate ( sig_solv )

  if ( allocated ( eps_solv ) ) deallocate ( eps_solv )

  if ( allocated ( id_solv ) ) deallocate ( id_solv )

  if ( allocated ( id_mol ) ) deallocate ( id_mol )

  if ( allocated ( n_0_multispec ) ) deallocate ( n_0_multispec )

  if ( allocated ( rho_0_multispec ) ) deallocate ( rho_0_multispec )

  if ( allocated ( c_s ) ) deallocate ( c_s )

  if ( allocated ( c_delta ) ) deallocate ( c_delta )

  if ( allocated ( c_d ) ) deallocate ( c_d )

  if ( allocated ( id_mol ) ) deallocate ( id_mol )

  if ( allocated ( Vk ) ) deallocate ( Vk )

  if ( allocated ( V_int ) ) deallocate ( V_int )

  if ( allocated ( weight_function_1_k ) ) deallocate ( weight_function_1_k )

  if ( allocated ( weight_function_2_k ) ) deallocate ( weight_function_2_k )

  if ( allocated ( weight_function_3_k ) ) deallocate ( weight_function_3_k )

  if ( allocated ( weight_function_0_k ) ) deallocate ( weight_function_0_k )

  if ( allocated ( radius ) ) deallocate ( radius )

!  if ( allocated ( eta_multispec ) ) deallocate ( eta_multispec )

  if ( allocated ( Fexc_0_multispec ) ) deallocate ( Fexc_0_multispec )

  if ( allocated ( muexc_0_multispec ) ) deallocate ( muexc_0_multispec )

  if ( allocated ( v_perturbation_k ) ) deallocate ( v_perturbation_k )

  if ( allocated ( atomic_nbr ) ) deallocate ( atomic_nbr )

  if ( allocated ( V_coulomb ) ) deallocate ( V_coulomb )  

  end subroutine deallocate_everything_system


end module system
