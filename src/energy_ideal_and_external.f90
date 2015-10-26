! This SUBROUTINE computes the ideal part of the free energy functional.
subroutine energy_ideal_and_external (fid, fext, df)

    use precision_kinds, only: dp
    use system, ONLY: thermocond
    use module_solvent, only: solvent
    use module_grid, only: grid
    use module_input, only: getinput

    implicit none

    real(dp), intent(out) :: fid, fext
    real(dp), intent(inout) :: df(:,:,:,:,:)
    integer :: is, nx, ny, nz, io, ns, ix, iy, iz
    real(dp) :: x, x0, vext, volume, dv, kT, mu, w
    real(dp), parameter :: zerodp = 0._dp


    ns = solvent(1)%nspec
    kT = thermocond%kbT
    dv = grid%dv
    volume = grid%v ! volume
    mu = getinput%dp( 'imposed_chempot', defaultvalue=0._dp)
    if ( ns/=1 .AND. mu/=0._dp) STOP "Imposing a chemical potential is valid only for single-species solvent"
    !
    ! fid = integrate over whole space of   x.log(x/x0)-x+x0 = Int[x.(log(x/x0)-1)] + Int[x0]
    !
    ! fext = Int[x*(vext-mu)]
    !
    fid = zerodp
    do is=1,ns
        x0 = solvent(is)%n0 ! bulk density
        do io=1,grid%no
            w = grid%w(io) ! weight of the orientation
            do iz=1,grid%nz
                do iy=1,grid%nz
                    do ix=1,grid%nx
                        x = solvent(is)%density(ix,iy,iz,io)
                        vext = solvent(is)%vext(ix,iy,iz,io) - mu
                        fid = fid + kT*(x*log(x/x0)-x+x0)*dv*w
                        fext = fext + x*vext*dv*w
                        if (vext<36._dp) then
                            df (ix,iy,iz,io,is) = df (ix,iy,iz,io,is) + w*(kT*log(x/x0)+vext)
                        else
                            df (ix,iy,iz,io,is) = df (ix,iy,iz,io,is) + w*(kT*log(x/x0)+vext)
                        end if
                    end do
                end do
            end do
        end do
    end do


!     Fideal = 0.0_dp! init Fideal to zero and its gradient
!     do concurrent( i=1:nx, j=1:ny, k=1:nz, io=1:grid%no, s=1:solvent(1)%nspec )
!       psi = cg_vect_new(i,j,k,io,s)
!       rho = psi**2
!       fid = fid + Fideal_local (io,s,rho)
!       dF_new(i,j,k,io,s) = dF_new(i,j,k,io,s) + dFideal_local (io,s,psi)
!     end do
!
!     Fideal = Fideal * thermocond%kbT * grid%dv ! integration factor
!     FF = FF + Fideal
!
!
!     CONTAINS
!
! !===================================================================================================================================
!
!     PURE FUNCTION dFideal_local (io,s,psi)
!         INTEGER, INTENT(IN) :: io, s
!         real(dp), INTENT(IN) :: psi
!         real(dp) :: dFideal_local
!         IF (abs(psi) > epsilon(1._dp)) THEN
!             dFideal_local = 2.0_dp * psi * prefactor(io,s) * grid%dv * ( thermocond%kbT*LOG(psi**2) )
!         ELSE
!             dFideal_local = 0._dp
!         END IF
!     END FUNCTION dFideal_local
!
! !===================================================================================================================================
!
!     PURE FUNCTION Fideal_local (io,s,rho)
!         INTEGER, INTENT(IN) :: io, s
!         real(dp), INTENT(IN) :: rho
!         real(dp) :: Fideal_local
!         IF (abs(rho) > epsilon(1.0_dp) ) THEN
!             Fideal_local = prefactor(io,s) * (rho*LOG(rho)-rho+1.0_dp)
!         ELSE
!             Fideal_local = prefactor(io,s)
!         END IF
!     END FUNCTION Fideal_local
!
! !===================================================================================================================================
!
!     PURE FUNCTION prefactor (io,s)
!         INTEGER, INTENT(IN) :: io,s
!         real(dp) :: prefactor
!         prefactor = grid%w(io) * solvent(s)%rho0
!     END FUNCTION
!
! !===================================================================================================================================

END SUBROUTINE energy_ideal_and_external
