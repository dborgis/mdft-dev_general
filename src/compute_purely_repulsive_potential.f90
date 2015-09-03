! ! This subroutine computes the external potential induced by purely repulsive r^-12 spheres at each solute sites. You may be interesed by Dzubiella and Hansen, J. Chem. Phys. 121 (2004)
! ! The potential has the form V(r)=kbT*(r-R0)^-(12)
!
SUBROUTINE compute_purely_repulsive_potential
!
!     use precision_kinds     ,ONLY: dp, i2b
!     use module_input               ,ONLY: getinput%dp, verbose
!     use system              ,ONLY: thermoCond, solute, solvent
!     use external_potential  ,ONLY: Vext_total
!     use constants           ,ONLY: zerodp=>zero
    use module_grid, only: grid
!
!     IMPLICIT NONE
!
!     INTEGER(i2b):: i,j,k,o,p,m,n,s
!     INTEGER(i2b), POINTER :: nfft1=>grid%n_nodes(1), nfft2=>grid%n_nodes(2), nfft3=>grid%n_nodes(3)
!     REAL(dp):: x_grid,y_grid,z_grid ! coordinates of grid nodes
!     REAL(dp):: x_m,y_m,z_m ! solvent sites coordinates
!     REAL(dp):: x_nm,y_nm,z_nm ! coordinate of vecteur solute-solvent
!     REAL(dp):: r_nm2 ! norm**2 of vector x_nm;y_nm;z_nm
!     REAL(dp):: dx,dy,dz,V_psi
!     REAL(dp):: time0,time1
!     REAL(dp), ALLOCATABLE :: xmod(:,:,:), ymod(:,:,:), zmod(:,:,:), Vrep(:,:,:,:,:)
!     REAL(dp):: radius_of_purely_repulsive_solute, radius_of_purely_repulsive_solute2
!     REAL(dp), PARAMETER :: Vmax=100._dp
!
!
!     CALL CPU_TIME(time0)    ! init timer
!
!     ! get the radius of the purely repulsive solute
!     ! the radius is defined such as in Dzubiella and Hansen, J Chem Phys 121 , 2011
!     ! look for tag 'purely_repulsive_solute_radius' in dft.in for hard wall thickness
!     radius_of_purely_repulsive_solute=getinput%dp('radius_of_purely_repulsive_solute')
!     radius_of_purely_repulsive_solute2 = radius_of_purely_repulsive_solute**2
!
!     ! tabulate coordinates of solvent sites for each omega and psi angles
!     ALLOCATE ( xmod (solvent(1)%nsite, grid%no) ,SOURCE=zerodp)
!     ALLOCATE ( ymod (solvent(1)%nsite, grid%no) ,SOURCE=zerodp)
!     ALLOCATE ( zmod (solvent(1)%nsite, grid%no) ,SOURCE=zerodp)
!
!     DO CONCURRENT ( m=1:solvent(1)%nsite, io=1:grid%no )
!         xmod (m,io) = DOT_PRODUCT( [grid%Rotxx(io), grid%Rotxy(io), grid%Rotxz(io)] , solvent(1)%site(m)%r )
!         ymod (m,io) = DOT_PRODUCT( [grid%Rotyx(io), grid%Rotyy(io), grid%Rotyz(io)] , solvent(1)%site(m)%r )
!         zmod (m,io) = DOT_PRODUCT( [grid%Rotzx(io), grid%Rotzy(io), grid%Rotzz(io)] , solvent(1)%site(m)%r )
!     END DO
!
!     dx=grid%dl(1)
!     dy=grid%dl(2)
!     dz=grid%dl(3)
!     ALLOCATE ( Vrep (nfft1,nfft2,nfft3,grid%no) ,SOURCE=Vmax)
!     DO s=1,solvent(1)%nspec
!         DO k=1,nfft3
!             z_grid = REAL(k-1,dp) * dz
!             DO j=1,nfft2
!                 y_grid = REAL(j-1,dp) * dy
!                 DO i=1,nfft1
!                     x_grid = REAL(i-1,dp) * dx
!                     do io=1,grid%no
!                         V_psi = 0.0_dp
!                    psiloop: DO m=1, 1 ! nb_solvent_sites => FOR DZUBIELLA HANSEN ONLY Oxygen atom is taken into account
!                                 x_m = x_grid + xmod (m,io)
!                                 y_m = y_grid + ymod (m,io)
!                                 z_m = z_grid + zmod (m,io)
!                                 DO n=1, solute%nsite
!                                     x_nm = x_m - solute%site(n)%r(1)
!                                     y_nm = y_m - solute%site(n)%r(2)
!                                     z_nm = z_m - solute%site(n)%r(3)
!                                     r_nm2 = x_nm**2+y_nm**2+z_nm**2
!                                     IF ( r_nm2   <= radius_of_purely_repulsive_solute2 ) THEN
!                                         V_psi = Vmax
!                                     ELSE
!                                         V_psi = V_psi + thermoCond%kbT/(Sqrt(r_nm2)-radius_of_purely_repulsive_solute)**12
!                                     END IF
!                                     IF (V_psi >= Vmax) THEN
!                                         V_psi = Vmax
!                                         EXIT psiloop
!                                     END IF
!                                 END DO
!                             END DO psiloop
!                             Vrep(i,j,k,o,p) = V_psi
!                         END DO
!                     END DO
!                 END DO
!             END DO
!         END DO
!         Vext_total (:,:,:,:,:,s) = Vext_total (:,:,:,:,:,s) + Vrep (:,:,:,:,:)
!     END DO
!     DEALLOCATE (xmod,ymod,zmod)
!
!     CALL CPU_TIME (time1)
!
!
!     IF (verbose) THEN
!         BLOCK
!             use constants, ONLY : fourpi
!             REAL(dp), ALLOCATABLE :: temparray(:,:,:)
!             CHARACTER(50) :: filename
!             ! warn user about vrep extrema for debugging
!             WRITE(*,*) 'minval(Vrep) = ' , MINVAL( Vrep )
!             WRITE(*,*) 'maxval(Vrep) = ' , MAXVAL( Vrep )
!             ALLOCATE ( temparray ( nfft1 , nfft2 , nfft3 ) ,SOURCE=0._dp)
!             CALL mean_over_orientations ( Vrep , temparray )! mean over orientations and print
!             temparray = temparray / fourpi
!             filename = 'output/Vrep.cube'
!             CALL write_to_cube_file ( temparray , filename )
!             DEALLOCATE ( temparray )
!             PRINT*,"time to compute purely repulsive : ",time1-time0
!         END BLOCK
!     END IF
!
!     DEALLOCATE (Vrep)
!
!
END SUBROUTINE compute_purely_repulsive_potential
