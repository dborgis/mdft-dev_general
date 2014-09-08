!===================================================================================================================================
SUBROUTINE compute_vcoul_as_sum_of_pointcharges
!===================================================================================================================================
! Returns the direct sum of all qi*qj/rij 

    USE precision_kinds     ,ONLY: dp, i2b
    USE system              ,ONLY: spaceGrid, nb_species, solvent, soluteSite
    USE constants           ,ONLY: fourpi, qfact
    USE external_potential  ,ONLY: Vext_q
    USE quadrature          ,ONLY: angGrid, molRotGrid, Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz
    USE input               ,ONLY: verbose

    IMPLICIT NONE

    INTEGER(i2b)    :: i,j,k,o,p,m,n,s
    REAL(dp)        :: xgrid(3), xv(3), xuv(3), xuv2, V_psi
    REAL, PARAMETER :: hardShellRadiusSQ=1._dp ! charge pseudo radius **2   == Rc**2
    REAL(dp), DIMENSION( SIZE(solvent(1)%site), molRotGrid%n_angles, angGrid%n_angles) :: xmod, ymod, zmod

    IF (.NOT. ALLOCATED(Vext_q)) STOP "Vext_q should be allocated in SUBROUTINE compute_vcoul_as_sum_of_pointcharges"
    IF ( ANY(Vext_q/=0.0_dp) ) STOP "Vext_q should be zero at the beginning of SUBROUTINE compute_vcoul_as_sum_of_pointcharges"

    ! test if all solutes have zero charge then don't waste your time : go to end of SUBROUTINE
    IF ( ALL(soluteSite%q == 0._dp) .OR. ALL(solvent(1)%site%q == 0._dp) ) THEN
        IF (verbose) PRINT*,"The electrostatic potential energy is zero."
        RETURN
    END IF
    
    ! precompute Rot_ij(omega,psi)*k_solv(a) for speeding up
    DO CONCURRENT (m=1:SIZE(solvent(1)%site), p=1:molRotGrid%n_angles, o=1:angGrid%n_angles)
        xmod(m,p,o)= DOT_PRODUCT( [Rotxx(o,p),Rotxy(o,p),Rotxz(o,p)] , solvent(1)%site(m)%r )
        ymod(m,p,o)= DOT_PRODUCT( [Rotyx(o,p),Rotyy(o,p),Rotyz(o,p)] , solvent(1)%site(m)%r )
        zmod(m,p,o)= DOT_PRODUCT( [Rotzx(o,p),Rotzy(o,p),Rotzz(o,p)] , solvent(1)%site(m)%r )
    END DO

    
    DO CONCURRENT (s=1:nb_species, i=1:spaceGrid%n_nodes(1), j=1:spaceGrid%n_nodes(2), k=1:spaceGrid%n_nodes(3), &
                   o=1:angGrid%n_angles, p=1:molRotGrid%n_angles)

        xgrid = REAL([i,j,k]-1,dp)*spaceGrid%dl
        V_psi = 0._dp

        DO CONCURRENT (m=1:SIZE(solvent(1)%site), n=1:SIZE(soluteSite), (solvent(1)%site(m)%q/=0._dp .AND. soluteSite(n)%q/=0._dp))
            xv = xgrid + [xmod(m,p,o), ymod(m,p,o), zmod(m,p,o)]
            xuv = xv - soluteSite(n)%r
            xuv2 = SUM(xuv**2)
            IF (xuv2 < hardShellRadiusSQ) THEN
                V_psi = HUGE(1.0_dp)
                CYCLE
            ELSE
                V_psi = V_psi + qfact*soluteSite(n)%q*solvent(1)%site(m)%q/SQRT(xuv2)
            END IF
        END DO

        IF (V_psi>100._dp) THEN; V_psi=100._dp; END IF!ELSE IF (V_psi<-100._dp) THEN; V_psi=-100._dp; END IF
        Vext_q(i,j,k,o,p,s ) = V_psi

    END DO
    
    IF (verbose) THEN
        BLOCK
            REAL(dp), ALLOCATABLE :: temparray(:,:,:)
            CHARACTER(50) :: filename
            PRINT*,"minval and maxval of Vext_q :", MINVAL( Vext_q (:,:,:,:,:,:) ), MAXVAL( Vext_q (:,:,:,:,:,:) )
            ! Get the external potential over orientations and print it
            ALLOCATE ( temparray ( spaceGrid%n_nodes(1),spaceGrid%n_nodes(2),spaceGrid%n_nodes(3) ),SOURCE=0._dp)
            CALL mean_over_orientations ( Vext_q (:,:,:,:,:,1) , temparray )
            temparray = temparray/fourpi
            filename = 'output/Vcoul.cube'
            CALL write_to_cube_file ( temparray , filename )
            filename = 'output/Vcoul_along-z.dat'
            CALL compute_z_density ( temparray , filename )
            DEALLOCATE(temparray)
        END BLOCK
    END IF

END SUBROUTINE
!===================================================================================================================================
