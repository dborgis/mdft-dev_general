! In this subroutine one calls the different parts of the total energy
! first is computed the radial_part, then ... blabla.
! this subroutine is the one called by the minimization stuff
! for computing the total energy and associated gradient

subroutine energy_and_gradient


use precision_kinds , only : i2b , dp

use input , only : input_line,input_log,input_char

use cg , only : FF , dF , cg_vect

use fft , only : in_forward , in_backward , out_forward , out_backward , plan_forward , plan_backward

use system, only : nb_omega, nb_psi, sigma_k, nfft1, nfft2, nfft3, Lx, Ly, Lz, deltax, rho_c, deltaV, nb_species,&
                   rho_0_multispec , rho_c_k_myway

use constants, only : twopi, fourpi, qfact

use quadrature, only : weight, weight_psi, sym_order

implicit none

integer ( kind = i2b ) :: i , j , k , l, icg ,  nf1, nf2, nf3
! dummy

logical :: hydrophobicity_tag

real (kind=dp) :: twopioLx , twopioLy , twopioLz, deltaVk


! init total energy and gradient to 0
! FF is the TOTAL ENERGY of the system, it is thus the functional of the density that is minimized by solver
! dF is the gradient of FF with respect to all coordinates. Remember it is of the kind dF ( number of variables over density (ie angles etc))

FF = 0.0_dp ! scalar

dF = 0.0_dp ! array

if (.not. allocated (rho_c_k_myway) ) allocate (rho_c_k_myway(nfft1/2+1, nfft2, nfft3))

rho_c_k_myway=(0.0_dp, 0.0_dp)

icg=0

twopioLx=twopi/Lx

twopioLy=twopi/Ly

twopioLz=twopi/Lz

nf1=nfft1/2

nf2=nfft2/2

nf3=nfft3/2

deltaVk=twopi**3/(Lx*Ly*Lz)




call energy_external



! compute ideal part of the total energy

call energy_ideal



! compute radial part of the excess free energy
! test if solvent fluid is of kind "water (tabulated cs) " or hard sphere
!do i = 1 , size ( input_line )
!  j = len ( 'hard_sphere_fluid' )
!  if ( input_line (i) (1:j) == 'hard_sphere_fluid' .and. input_line (i) (j+4:j+4) == 'T' ) then
!    ! fluid is composed of hard spheres  
!    call energy_hard_sphere_fmt ! => pure hard sphere contribution
!    exit ! it is useless to continue the loop over i if you've found the tag
!  end if
!end do

if (input_log('hard_sphere_fluid') )then
call energy_hard_sphere_fmt ! => pure hard sphere contribution
end if

! test if there is a LJ perturbation to hard spheres ! WCA model etc. to implement more intelligently

!do k = 1 , size ( input_line )
!  l = len ( 'lennard_jones_perturbation_to_hard_spheres' )
!  if ( input_line (k) (1:l) == 'lennard_jones_perturbation_to_hard_spheres' .and. input_line (k) (l+4:l+4) == 'T' ) then
!    call lennard_jones_perturbation_to_hard_spheres ! lennard jones perturbative contribution => Weeks-Chandler-Anderson
!    exit
!  end if
!end do

if (input_log('lennard_jones_perturbation_to_hard_spheres') ) then
call lennard_jones_perturbation_to_hard_spheres ! lennard jones perturbative contribution => Weeks-Chandler-Anderson
end if


if (input_log('read_ck_or_chi')) then
  if (input_log('hydrophobicity')) then
    if (trim(adjustl(input_char('treatment_of_hydro')))=='C')  then
      call cs_plus_hydro
    else if (trim(adjustl(input_char('treatment_of_hydro')))=='VdW')  then
      call energy_hydro
    end if
  else
  call cs_from_dcf
  end if
end if


! bridge calculation: F(FMT)-F(c2hs)+F(c2H2O)

!do i = 1 , size ( input_line )
!  j = len ( 'bridge_hard_sphere' )
!  if ( input_line (i) (1:j) == 'bridge_hard_sphere' .and. input_line (i) (j+4:j+4) == 'T' ) then
!    ! compute effect of direct correlation function of reference fluid
!    call energy_cs_hard_sphere
!    exit ! it is useless to continue the loop over i if you've found the tag
!  end if
!end do

if (input_log('bridge_hard_sphere')) then
  call energy_cs_hard_sphere
end if


! Dipolar polarization. What user wants (use it or not) is checked in subroutine for clearer code.

call energy_polarization

call energy_polarization_myway


! Threebody term that is needed to empiricaly force H-bonding in water. What user wants (use it or not) is checked in subroutine for clearer code




call energy_threebody
call energy_threebody_faster



! inform user

write(*,*) 'FREE ENERGY = ',FF

write(*,*)'-----------------------'



end subroutine energy_and_gradient
