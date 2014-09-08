! This subroutine computes the external potential induced by purely repulsive r^-12 spheres at each solute sites. You may be interesed by Dzubiella and Hansen, J. Chem. Phys. 121 (2004)
! The potential has the form V(r)=kbT*(r-R0)^-(12)

SUBROUTINE compute_purely_repulsive_potential

    USE precision_kinds     ,ONLY: dp, i2b
    USE input               ,ONLY: input_dp, verbose
    USE system              ,ONLY: thermoCond, nb_solute_sites, nb_solvent_sites, spaceGrid,&
                                   nb_species, soluteSite, solvent
    USE external_potential  ,ONLY: Vext_total
    USE quadrature          ,ONLY: Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz, angGrid, molRotGrid
    USE constants           ,ONLY: zerodp=>zero

    IMPLICIT NONE

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
    CALL test_allocated_Rotxx_etc
    ALLOCATE ( xmod (nb_solvent_sites,molRotGrid%n_angles,angGrid%n_angles) ,SOURCE=zerodp)
    ALLOCATE ( ymod (nb_solvent_sites,molRotGrid%n_angles,angGrid%n_angles) ,SOURCE=zerodp)
    ALLOCATE ( zmod (nb_solvent_sites,molRotGrid%n_angles,angGrid%n_angles) ,SOURCE=zerodp)
    DO CONCURRENT (m=1:nb_solvent_sites, p=1:molRotGrid%n_angles, o=1:angGrid%n_angles)
        xmod (m,p,o) = DOT_PRODUCT( [Rotxx(o,p), Rotxy(o,p), Rotxz(o,p)] , solvent(1)%site(m)%r )
        ymod (m,p,o) = DOT_PRODUCT( [Rotyx(o,p), Rotyy(o,p), Rotyz(o,p)] , solvent(1)%site(m)%r )
        zmod (m,p,o) = DOT_PRODUCT( [Rotzx(o,p), Rotzy(o,p), Rotzz(o,p)] , solvent(1)%site(m)%r )
    END DO

    dx=spaceGrid%dl(1)
    dy=spaceGrid%dl(2)
    dz=spaceGrid%dl(3)    
    ALLOCATE ( Vrep (nfft1,nfft2,nfft3,angGrid%n_angles,molRotGrid%n_angles) ,SOURCE=Vmax)
    DO s=1,nb_species
        DO k=1,nfft3
            z_grid = REAL(k-1,dp) * dz
            DO j=1,nfft2
                y_grid = REAL(j-1,dp) * dy
                DO i=1,nfft1
                    x_grid = REAL(i-1,dp) * dx
                    DO o=1,angGrid%n_angles
                        DO p=1,molRotGrid%n_angles
                            V_psi = 0.0_dp
                   psiloop: DO m=1, 1 ! nb_solvent_sites => FOR DZUBIELLA HANSEN ONLY Oxygen atom is taken into account
                                x_m = x_grid + xmod (m,p,o)
                                y_m = y_grid + ymod (m,p,o)
                                z_m = z_grid + zmod (m,p,o)
                                DO n=1, nb_solute_sites
                                    x_nm = x_m - soluteSite(n)%r(1)
                                    y_nm = y_m - soluteSite(n)%r(2)
                                    z_nm = z_m - soluteSite(n)%r(3)
                                    r_nm2 = x_nm**2+y_nm**2+z_nm**2
                                    IF ( r_nm2 == 0._dp ) THEN!  <= radius_of_purely_repulsive_solute2 ) THEN
                                        V_psi = Vmax
                                    ELSE
                                        V_psi = V_psi + thermoCond%kbT/(r_nm2**6)!)-radius_of_purely_repulsive_solute)**12
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

    CALL CPU_TIME (time1)


    IF (verbose) THEN
        BLOCK
            USE constants, ONLY : fourpi
            REAL(dp), ALLOCATABLE :: temparray(:,:,:)
            CHARACTER(50) :: filename
            ! warn user about vrep extrema for debugging
            WRITE(*,*) 'minval(Vrep) = ' , MINVAL( Vrep )
            WRITE(*,*) 'maxval(Vrep) = ' , MAXVAL( Vrep )
            ALLOCATE ( temparray ( nfft1 , nfft2 , nfft3 ) ,SOURCE=0._dp) 
            CALL mean_over_orientations ( Vrep , temparray )! mean over orientations and print
            temparray = temparray / fourpi
            filename = 'output/Vrep.cube'
            CALL write_to_cube_file ( temparray , filename )
            DEALLOCATE ( temparray )
            PRINT*,"time to compute purely repulsive : ",time1-time0
        END BLOCK
    END IF
    
    DEALLOCATE (Vrep)




    CONTAINS
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE test_allocated_Rotxx_etc
        IF (ALL(Rotxx==zerodp) .AND. &
            ALL(Rotxy==zerodp) .AND. &
            ALL(Rotxz==zerodp) .AND. &
            ALL(Rotyx==zerodp) .AND. &
            ALL(Rotyy==zerodp) .AND. &
            ALL(Rotyz==zerodp) .AND. &
            ALL(Rotzx==zerodp) .AND. &
            ALL(Rotzy==zerodp) .AND. &
            ALL(Rotzz==zerodp)) STOP "The rotation matrix in compute_purely_repulsive_potential should not be zero everywhere"
    END SUBROUTINE test_allocated_Rotxx_etc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
END SUBROUTINE compute_purely_repulsive_potential
