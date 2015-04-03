! this SUBROUTINE compute the external part of the free energy functional
subroutine energy_external (Fext)

  use precision_kinds   ,only: dp , i2b
  use system            ,only: nb_species, spaceGrid, solvent
  use quadrature        ,only: angGrid, molRotGrid
  use minimizer         ,only: CG_vect , FF , dF
  use external_potential,only: Vext_total
  use input             ,only: input_dp

  implicit none

  real(dp), intent(out) :: Fext
  integer(i2b) :: icg , i , j , k , o , p,s, nfft1, nfft2, nfft3
  real(dp) :: psi, wdfve, imposedChemPot

  Fext = 0.0_dp
  nfft1 = spacegrid%n_nodes(1)
  nfft2 = spacegrid%n_nodes(2)
  nfft3 = spacegrid%n_nodes(3)

  ! Impose a chemical potential
  imposedChemPot = input_dp('imposed_chempot',0._dp) ! 0 means that you don't impose any chemical potential (delta_chempot is built this way)
  if ( nb_species/=1 .AND. imposedChemPot/=0._dp) STOP "Imposing a chemical potential is valid only for single-species solvent"

  ! F_{ext}[\rho(\vec{r},\vec{\Omega})]=\int d \vec{r} d \vec{\Omega} V_{ext}(\vec{r},\vec{\Omega})\rho(\vec{r},\vec{\Omega})

  icg = 0
  do s = 1 , nb_species
    do i = 1 , nfft1
      do j = 1 , nfft2
        do k = 1 , nfft3
          do o = 1 , angGrid%n_angles
            do p = 1 , molRotGrid%n_angles
              icg   = icg + 1
              psi   = cg_vect(icg)
              wdfve = angGrid%weight(o) * molRotGrid%weight(p) * spaceGrid%dv * solvent(s)%rho0 * (Vext_total(i,j,k,o,p,s) - imposedChemPot)
              Fext  = Fext + psi**2 * wdfve
              dF(icg) = dF(icg) + 2.0_dp*psi*wdfve
            end do
          end do
        end do
      end do
    end do
  end do

  FF = FF + Fext

end subroutine energy_external
