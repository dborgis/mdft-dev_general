! Variable relative to the system being studied
MODULE system

    USE precision_kinds , ONLY: i2b,dp

    IMPLICIT NONE

    TYPE :: some_soluteSite
        CHARACTER(100) :: name
        INTEGER(i2b) :: type
        REAL(dp), DIMENSION(3) :: r
        REAL(dp) :: q, sig, eps, lambda1, lambda2
        INTEGER(i2b) :: Z ! atomic number
    END TYPE some_soluteSite
    
    TYPE (some_soluteSite), ALLOCATABLE, DIMENSION(:), TARGET :: soluteSite

    INTEGER(i2b) :: nb_species ! number of solvents in the species, e.g. 2 if the solvent is a mixture of water and acetone
    INTEGER(i2b) :: nb_solute_sites, nb_solvent_sites  ! nombre de site pour le solute et pour le solvent
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: x_mol, y_mol, z_mol ! positions des sites du solute dans la boite (repere absolu)
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: chg_mol, sig_mol, eps_mol  ! charge partielle et parametres LJ pour chaque site du solute (kJ/mol)
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: lambda1_mol, lambda2_mol     !Three body terms
    INTEGER(i2b), ALLOCATABLE, DIMENSION (:) :: atomic_nbr
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: x_solv, y_solv, z_solv  ! positions des sites du solvent dans le repere propre ede la molecule
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: chg_solv, sig_solv, eps_solv  ! charge partielle et parametres LJ pour chaque site du solvant (kJ/mol)
    INTEGER(i2b), ALLOCATABLE, DIMENSION(:) :: id_solv, id_mol ! atom type of each solvent site and solute site for instance id_solv(1)=1, id_solv(2)=2 and id_solv(3)=2 for OH2
    REAL(dp) :: Rc ! Taille de la charge (Ang)
    REAL(dp) :: temp  ! temperature du systeme lue dans dft.in 'temperature : XXXX'
    REAL(dp) :: kBT , beta
    
    ! Grille pour FFT
    TYPE :: spaceGridType
        INTEGER(i2b), DIMENSION(3) :: n_nodes ! number of grid nodes in direction x, y and z
        REAL(dp), DIMENSION(3) :: length ! total length in direction x, y and z
        REAL(dp), DIMENSION(3) :: dl ! elemental distance between two nodes in direction x, y and z
        REAL(dp) :: dv ! elemental volume
    END TYPE spaceGridType
    TYPE( spaceGridType ), TARGET :: spaceGrid
    REAL(dp), POINTER :: Lx => spaceGrid%length(1), Ly => spaceGrid%length(2), Lz => spaceGrid%length(3)
    REAL(dp), POINTER :: DeltaV => spaceGrid%dv
    REAL(dp), POINTER :: deltax => spaceGrid%dl(1), deltay => spaceGrid%dl(2), deltaz => spaceGrid%dl(3)
    INTEGER(i2b), POINTER :: nfft1 => spaceGrid%n_nodes(1), nfft2 => spaceGrid%n_nodes(2), nfft3 => spaceGrid%n_nodes(3) ! deprecated. Should be removed at some point

    REAL(dp) :: n_0 , rho_0   ! Densite du fluide homogene en part/A3 et incluant orientation
    REAL(dp), ALLOCATABLE, DIMENSION (:) :: n_0_multispec , rho_0_multispec ! here are the equivalent of n_0 and rho_0 in multispecies case
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: c_s, c_delta, c_d !> Direct correlation functions of various rotational invariants
    REAL(dp), ALLOCATABLE, DIMENSION(:) ::chi_L, chi_T
    REAL(dp), ALLOCATABLE, DIMENSION (:) :: c_s_hs ! c(2)(k) of a hard sphere
    COMPLEX(dp),ALLOCATABLE, DIMENSION(:,:,:) :: Vk !>@var perturabtion in kspace
    REAL(dp) :: delta_k ! distance between two k points in cs.in, cdelta.in, cd.in
    INTEGER(i2b) :: nb_k ! nb of k points in cs.in, cdelta.in, cd.in
    ! Electrostatics
    REAL ( dp ) , ALLOCATABLE , DIMENSION (:,:,:,:,:) :: wigma ! function of n1 , n2 , n3 , omega , psi
    ! charge factor & molecule polarization factor
    COMPLEX ( dp ) , ALLOCATABLE , DIMENSION (:,:,:,:,:,:) :: sigma_k
    REAL(dp), ALLOCATABLE , DIMENSION (:,:,:) :: rho_c ! charge density at each grid node
    COMPLEX(dp), ALLOCATABLE, DIMENSION (:,:,:) :: rho_c_k, rho_c_k_myway
    !Polarization
    COMPLEX (dp) , ALLOCATABLE,DIMENSION (:,:,:,:) :: pola_tot_x_k , pola_tot_y_k , pola_tot_z_k
    COMPLEX ( dp ) , ALLOCATABLE , DIMENSION (:,:,:,:,:,:) ::  molec_polarx_k,molec_polary_k,molec_polarz_k 
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
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::  V_int
    !REAL(dp), ALLOCATABLE , DIMENSION (:,:,:, : ) :: rho_n ! density per angle (=n/4pi) of each species at each node (nfft1,nfft2,nfft3,nb_species)
    ! for maximum efficiency, rho_species(nfft1,nfft2,nfft3,nb_species) should be used in loops where species is the outer loop (varies slowliest)
    !do species = 1 , nb_species
    !do z
    !do y 
    !do x
    !density(x,y,z,species)
    !> Hard spheres
    ! TODO this should be put in another module
    logical(4) :: HS
    COMPLEX(dp), ALLOCATABLE , DIMENSION (:,:,:,:) :: weight_function_1_k
    COMPLEX(dp), ALLOCATABLE , DIMENSION (:,:,:,:) :: weight_function_2_k
    COMPLEX(dp), ALLOCATABLE , DIMENSION (:,:,:,:) :: weight_function_3_k
    COMPLEX(dp), ALLOCATABLE , DIMENSION (:,:,:,:) :: weight_function_0_k
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: radius
    REAL(dp) :: eta, muexc_0, Fexc_0
    !REAL(dp), ALLOCATABLE , DIMENSION (:) :: eta_multispec ! packing fraction = 4/3 * pi * n_0 * radius **3 = pi/6 * n_0 * diameter**3
    REAL(dp), ALLOCATABLE , DIMENSION (:) :: Fexc_0_multispec ! 
    REAL(dp), ALLOCATABLE , DIMENSION (:) :: muexc_0_multispec ! temp var during implementation of multispecies
    COMPLEX(dp), ALLOCATABLE , DIMENSION (:,:,:) :: v_perturbation_k ! fourier transform of the lennard jones perturbation (WCA)
    REAL(dp), ALLOCATABLE , DIMENSION (:) :: mole_fraction ! mole fraction of each species "x_i"
    !> for xsf printing
    REAL(dp), ALLOCATABLE, DIMENSION (:,:,:) :: V_coulomb ! nfft1 nfft2 nfft3 angGrid%n_angles


END MODULE system
