!This routine compute : -The solvent molecular charge density, which can be used into Vcoul_from_solvent_charge_density.f90 to
!evaluate the electrostatic potential.
!                       -The solvent molecular polarization (from Ranieriet al : J. Chem. Phys. 98 (11) 1993) that can be used
!into energy_polarization_..._.f90 to compute the (multipolar) polarization Free energy.

subroutine chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin

    use precision_kinds     ,ONLY: dp
    use system              ,ONLY: solvent, grid => spacegrid, nb_species
    use fft                 ,ONLY: kx, ky, kz, k2
    use mathematica         ,ONLY: factorial

    implicit none

    integer :: nx, ny, nz, no, ns
    integer :: i, j, k, n, s, io
    real(dp)     :: r(3), kr, kvec(3)
    complex(dp)  :: fac, X
    complex(dp), parameter :: zeroc = cmplx(0._dp,0._dp), ic = cmplx(0._dp,1._dp)
    real(dp), parameter :: epsdp = epsilon(1._dp)
    type :: smoother_type
        real(dp) :: radius = 0.5_dp ! dramaticaly important
        real(dp) :: factor
    end type smoother_type
    type (smoother_type) :: smoother

    nx = grid%nx
    ny = grid%ny
    nz = grid%nz
    no = grid%no
    ns = nb_species

    ! sigma_k is the Fourier transformed charge density of a single solvent molecule in the reference frame defined by solvent.in
    ! molec_polar_k is the Fourier transformed molecular polarization
    do concurrent( s=1:size(solvent) , sum(abs(solvent(s)%site%q))>0) ! mask elimitates solvent molecules without point charges
        allocate( solvent(s)%sigma_k       (   nx/2+1, ny, nz, no), SOURCE=zeroC )
        allocate( solvent(s)%molec_polar_k (3, nx/2+1, ny, nz, no), SOURCE=zeroC )
    end do

    do concurrent ( i=1:nx/2+1, j=1:ny, k=1:nz, s=1:size(solvent)  , sum(abs(solvent(s)%site%q))>epsdp) ! mask elimitates solvent molecules without point charges)

        kvec = [ kx(i), ky(j), kz(k) ]
        smoother%factor =  exp(-smoother%radius**2 *k2(i,j,k)/2)

        do concurrent ( io=1:grid%no, n=1:SIZE(solvent(s)%site), abs(solvent(s)%site(n)%q)>epsdp )
            r(1) = dot_product(   [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)]  ,  solvent(s)%site(n)%r  )
            r(2) = dot_product(   [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)]  ,  solvent(s)%site(n)%r  )
            r(3) = dot_product(   [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)]  ,  solvent(s)%site(n)%r  )
            kr = dot_product( kvec, r )
            X = -iC*kr
            solvent(s)%sigma_k(i,j,k,io) = solvent(s)%sigma_k(i,j,k,io) + solvent(s)%site(n)%q *exp(X) *smoother%factor ! exact
            ! solvent(s)%sigma_k(i,j,k,o,p) = solvent(s)%sigma_k(i,j,k,o,p) +solvent(s)%site(n)%q* sum([(X**i/factorial(i), i=0,4)])&
            ! * smoother%factor ! Series expansion of exp(x) at 0 => multipole expansion of Vcoul(x). i=4 :: hexadecapole (16)
            if ( abs(kr)<=epsdp ) then
                solvent(s)%molec_polar_k(:,i,j,k,io) = solvent(s)%molec_polar_k(:,i,j,k,io) + solvent(s)%site(n)%q *r
            else
                fac = -iC*(exp(iC*kr)-1._dp)/kr *smoother%factor
                solvent(s)%molec_polar_k(:,i,j,k,io) = solvent(s)%molec_polar_k(:,i,j,k,io) + fac*solvent(s)%site(n)%q *r
            end if
        end do
    end do

    call removeTraceOfMolecularPolarizationTensor

contains

    subroutine removeTraceOfMolecularPolarizationTensor
        implicit none
        integer :: i,j,k,s,d
        do concurrent (i=1:nx/2+1, j=1:ny, k=1:nz, s=1:size(solvent), d=1:3)
            solvent(s)%molec_polar_k(d,i,j,k,:) = solvent(s)%molec_polar_k(d,i,j,k,:)  &
            -sum( solvent(s)%molec_polar_k(d,i,j,k,:) ) /real(grid%no,dp)
        end do
    end subroutine removeTraceOfMolecularPolarizationTensor

END SUBROUTINE chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin
