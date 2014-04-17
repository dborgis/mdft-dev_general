! This subroutine computes the external potential induced by purely repulsive r^-12 spheres at each solute sites. You may be interesed by Dzubiella and Hansen, J. Chem. Phys. 121 (2004)
! The potential has the form V(r)=kbT*(r-R0)^-(12)

SUBROUTINE compute_purely_repulsive_potential (Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)

    USE precision_kinds     ,ONLY: dp, i2b
    USE input               ,ONLY: input_dp, verbose
    USE system              ,ONLY: x_solv, y_solv, z_solv, x_mol, y_mol, z_mol, beta, nb_solute_sites, nb_solvent_sites, spaceGrid,&
                                   nb_species
    USE external_potential  ,ONLY: Vext_total
    USE quadrature          ,ONLY: angGrid, molRotGrid

    IMPLICIT NONE

    REAL(dp), DIMENSION(angGrid%n_angles,molRotGrid%n_angles), INTENT(IN) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz ! rotation matrix for sites
    INTEGER(i2b):: i,j,k,o,p,m,n,s
    INTEGER(i2b), POINTER :: nfft1=>spaceGrid%n_nodes(1), nfft2=>spaceGrid%n_nodes(2), nfft3=>spaceGrid%n_nodes(3)
    REAL(dp):: x_grid,y_grid,z_grid ! coordinates of grid nodes
    REAL(dp):: x_m,y_m,z_m ! solvent sites coordinates
    REAL(dp):: x_nm,y_nm,z_nm ! coordinate of vecteur solute-solvent
    REAL(dp):: r_nm2 ! norm**2 of vector x_nm;y_nm;z_nm
    REAL(dp):: dx,dy,dz,kbT,V_psi
    REAL(dp):: time0,time1
    REAL(dp), ALLOCATABLE :: xmod(:,:,:), ymod(:,:,:), zmod(:,:,:), Vrep(:,:,:,:,:)
    REAL(dp):: radius_of_purely_repulsive_solute, radius_of_purely_repulsive_solute2
    REAL(dp), PARAMETER :: Vmax=100._dp
    

    CALL CPU_TIME(time0)    ! init timer

    ! get the radius of the purely repulsive solute
    ! the radius is defined such as in Dzubiella and Hansen, J Chem Phys 121 , 2011
    ! look for tag 'purely_repulsive_solute_radius' in dft.in for hard wall thickness
    radius_of_purely_repulsive_solute=input_dp('radius_of_purely_repulsive_solute')
    radius_of_purely_repulsive_solute2 = radius_of_purely_repulsive_solute**2

    ! tabulate coordinates of solvent sites for each omega and psi angles
    ALLOCATE ( xmod (nb_solvent_sites,molRotGrid%n_angles,angGrid%n_angles) ,SOURCE=0._dp)
    ALLOCATE ( ymod (nb_solvent_sites,molRotGrid%n_angles,angGrid%n_angles) ,SOURCE=0._dp)
    ALLOCATE ( zmod (nb_solvent_sites,molRotGrid%n_angles,angGrid%n_angles) ,SOURCE=0._dp)
    DO CONCURRENT (m=1:nb_solvent_sites, p=1:molRotGrid%n_angles, o=1:angGrid%n_angles)
        xmod (m,p,o) = Rotxx (o,p) * x_solv (m) + Rotxy (o,p) * y_solv (m) + Rotxz (o,p) * z_solv (m)
        ymod (m,p,o) = Rotyx (o,p) * x_solv (m) + Rotyy (o,p) * y_solv (m) + Rotyz (o,p) * z_solv (m)
        zmod (m,p,o) = Rotzx (o,p) * x_solv (m) + Rotzy (o,p) * y_solv (m) + Rotzz (o,p) * z_solv (m)
    END DO

    kbT = 1.0_dp/beta
    dx=spaceGrid%dl(1)
    dy=spaceGrid%dl(2)
    dz=spaceGrid%dl(3)    
    ALLOCATE ( Vrep (nfft1,nfft2,nfft3,angGrid%n_angles,molRotGrid%n_angles) ,SOURCE=0._dp)
    DO CONCURRENT (s=1:nb_species)
        DO CONCURRENT (k=1:nfft3)
            z_grid = REAL(k-1,dp) * dz
            DO CONCURRENT (j=1:nfft2)
                y_grid = REAL(j-1,dp) * dy
                DO CONCURRENT (i=1:nfft1)
                    x_grid = REAL(i-1,dp) * dx
                    DO CONCURRENT (o=1:angGrid%n_angles)
                        DO CONCURRENT (p=1:molRotGrid%n_angles)
                            V_psi = 0.0_dp
                   psiloop: DO m=1,nb_solvent_sites
                                x_m = x_grid + xmod (m,p,o)
                                y_m = y_grid + ymod (m,p,o)
                                z_m = z_grid + zmod (m,p,o)
                                DO n=1,nb_solute_sites
                                    x_nm = x_m - x_mol(n)
                                    y_nm = y_m - y_mol(n)
                                    z_nm = z_m - z_mol(n)
                                    r_nm2 = x_nm**2+y_nm**2+z_nm**2
                                    IF ( r_nm2 == radius_of_purely_repulsive_solute2 ) THEN
                                        V_psi = Vmax
                                    ELSE
                                        V_psi = V_psi + kbT*1._dp/(SQRT(r_nm2)-radius_of_purely_repulsive_solute)**12
                                    END IF
                                    IF (V_psi >= Vmax) THEN
                                        V_psi = Vmax
                                        EXIT psiloop
                                    END IF
                                END DO
                            END DO psiloop
                            Vrep(i,j,k,o,p) = V_psi
                        END DO
                    END DO
                END DO
            END DO
        END DO
        Vext_total (:,:,:,:,:,s) = Vext_total (:,:,:,:,:,s) + Vrep (:,:,:,:,:)
    END DO
    DEALLOCATE (xmod,ymod,zmod)

    IF (verbose) THEN
        BLOCK
            USE constants, ONLY : fourpi
            REAL(dp), ALLOCATABLE :: temparray(:,:,:)
            CHARACTER(50) :: filename
            ! warn user about vrep extrema for debugging
            WRITE(*,*) 'minval(Vrep) = ' , MINVAL( Vrep )
            WRITE(*,*) 'maxval(Vrep) = ' , MAXVAL( Vrep )
            WRITE(*,*) 'time for compute_purely_repulsive_potential ' , time1 - time0
            ALLOCATE ( temparray ( nfft1 , nfft2 , nfft3 ) ,SOURCE=0._dp) 
            CALL mean_over_orientations ( Vrep , temparray )! mean over orientations and print
            temparray = temparray / fourpi
            filename = 'output/Vrep.cube'
            CALL write_to_cube_file ( temparray , filename )
            DEALLOCATE ( temparray )
        END BLOCK
    END IF
    
    DEALLOCATE (Vrep)
    CALL CPU_TIME (time1)
    PRINT*,"time to compute purely repulsive : ",time1-time0;STOP
END SUBROUTINE compute_purely_repulsive_potential
