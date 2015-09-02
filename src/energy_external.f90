! this SUBROUTINE compute the external part of the free energy functional
subroutine energy_external (Fext)

  use precision_kinds    ,only: dp, i2b
  use system             ,only: nb_species, spacegrid, solvent
  use minimizer          ,only: cg_vect_new, FF, dF_new
  use external_potential ,only: Vext_total
  use input              ,only: input_dp

  implicit none

  real(dp), intent(out) :: Fext
  integer(i2b) :: i, j, k, o, p, s, nx, ny, nz, no, ns, io
  real(dp) :: psi, wdfve, imposedChemPot

  nx = spacegrid%nx
  ny = spacegrid%ny
  nz = spacegrid%nz
  no = spacegrid%no
  ns = nb_species

  ! Impose a chemical potential
  imposedChemPot = input_dp( 'imposed_chempot', defaultvalue=0._dp) ! 0 means that you don't impose any chemical potential (delta_chempot is built this way)
  if ( nb_species/=1 .AND. imposedChemPot/=0._dp) STOP "Imposing a chemical potential is valid only for single-species solvent"

  ! F_{ext}[\rho(\vec{r},\vec{\Omega})]=\int d \vec{r} d \vec{\Omega} V_{ext}(\vec{r},\vec{\Omega})\rho(\vec{r},\vec{\Omega})

  Fext = 0._dp
  do concurrent( i=1:nx, j=1:ny, k=1:nz, io=1:spacegrid%no, s=1:ns )
      psi = cg_vect_new(i,j,k,io,s)
      wdfve = spacegrid%w(io) * spacegrid%dv * solvent(s)%rho0 * (Vext_total(i,j,k,io,s) - imposedChemPot)
      Fext  = Fext + psi**2 * wdfve
      dF_new(i,j,k,io,s) = dF_new(i,j,k,io,s) + 2.0_dp*psi*wdfve
  end do
  FF = FF + Fext

end subroutine energy_external
