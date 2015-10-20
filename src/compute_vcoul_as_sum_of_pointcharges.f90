!===================================================================================================================================
SUBROUTINE compute_vcoul_as_sum_of_pointcharges
!===================================================================================================================================
! Returns the direct sum of all qi*qj/rij

    use precision_kinds     ,ONLY: dp, i2b
    use module_solute, only: solute
    use module_solvent, only: solvent
    use constants           ,ONLY: qfact
    use external_potential  ,ONLY: Vext_q
    use module_input               ,only: verbose
    use module_grid, only: grid

    IMPLICIT NONE

    INTEGER(i2b)    :: i,j,k,o,p,m,n,s,io,no
    REAL(dp)        :: xgrid(3), xv(3), xuv(3), xuv2, V_psi
    REAL, PARAMETER :: hardShellRadiusSQ=1._dp ! charge pseudo radius **2   == Rc**2
    real, parameter :: epsdp = epsilon(1._dp)
    REAL(dp), DIMENSION( SIZE(solvent(1)%site), grid%no) :: xmod, ymod, zmod

    IF (.NOT. ALLOCATED(Vext_q)) STOP "Vext_q should be allocated in SUBROUTINE compute_vcoul_as_sum_of_pointcharges"
    IF ( ANY(Vext_q/=0.0_dp) ) STOP "Vext_q should be zero at the beginning of SUBROUTINE compute_vcoul_as_sum_of_pointcharges"

    if (size(solvent)/=1) stop "CRITICAL. Compute_vcoul_as_pointcharges implemented for one solvent species only"

    !... return if no solute sites wear partial charges
    IF ( ALL(abs(solute%site%q) <= epsilon(1.0_dp)) ) return

    !... return if all sites of all solvent species are zero
    block
        logical :: ishouldreturn
        ishouldreturn = .true.
        do s=1,size(solvent)
            if ( any( abs(solvent(s)%site%q) >= epsilon(1.0_dp) ) ) ishouldreturn = .false.
        end do
        if ( ishouldreturn ) then
            if (verbose) print*, "The electrostatic potential is zero because no solvent site wear partial charges"
            print*,"haaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
            return
        end if
    end block


    ! precompute rot_ij(omega,psi)*k_solv(a) for speeding up
    DO CONCURRENT ( m=1:solvent(1)%nsite, io=1:grid%no )
        xmod (m,io) = DOT_PRODUCT( [grid%rotxx(io), grid%rotxy(io), grid%rotxz(io)] , solvent(1)%site(m)%r )
        ymod (m,io) = DOT_PRODUCT( [grid%rotyx(io), grid%rotyy(io), grid%rotyz(io)] , solvent(1)%site(m)%r )
        zmod (m,io) = DOT_PRODUCT( [grid%rotzx(io), grid%rotzy(io), grid%rotzz(io)] , solvent(1)%site(m)%r )
    END DO


    do s=1,solvent(1)%nspec
        do io=1,grid%no
            do k=1,grid%nz
                do j=1,grid%ny
                    do i=1,grid%nx

                        xgrid = real([i,j,k]-1,dp)*grid%dl
                        v_psi = 0._dp

                        do m=1,solvent(1)%nsite
                            if ( abs(solvent(s)%site(m)%q) < epsdp ) cycle
                            do n=1,solute%nsite
                                if( abs(solute%site(n)%q) < epsdp ) cycle
                                xv = xgrid + [xmod(m,io), ymod(m,io), zmod(m,io)]
                                xuv = xv - solute%site(n)%r
                                xuv2 = SUM(xuv**2)
                                IF (xuv2 < hardShellRadiusSQ) THEN
                                    V_psi = HUGE(1.0_dp)
                                    CYCLE
                                ELSE
                                    V_psi = V_psi + qfact *solute%site(n)%q *solvent(s)%site(m)%q/SQRT(xuv2)
                                END IF
                            end do
                        end do

                        Vext_q(i,j,k,io,s) = min( 100._dp, v_psi )

                    end do
                end do
            end do
        end do
    end do


    !
    ! DO CONCURRENT (s=1:solvent(1)%nspec, i=1:grid%n_nodes(1), j=1:grid%n_nodes(2), k=1:grid%n_nodes(3), &
    !                io=1:grid%no)
    !
    !     xgrid = REAL([i,j,k]-1,dp)*grid%dl
    !     V_psi = 0._dp
    !
    !     DO CONCURRENT (m=1:SIZE(solvent(s)%site), n=1:SIZE(solute%site),&
    !                             ( abs(solvent(s)%site(m)%q)>=epsilon(1.0_dp) .AND. abs(solute%site(n)%q)>=epsilon(1.0_dp)) )
    !         xv = xgrid + [xmod(m,io), ymod(m,io), zmod(m,io)]
    !         xuv = xv - solute%site(n)%r
    !         xuv2 = SUM(xuv**2)
    !         IF (xuv2 < hardShellRadiusSQ) THEN
    !             V_psi = HUGE(1.0_dp)
    !             CYCLE
    !         ELSE
    !             V_psi = V_psi + qfact *solute%site(n)%q *solvent(s)%site(m)%q/SQRT(xuv2)
    !         END IF
    !     END DO
    !
    !     IF (V_psi>100._dp) THEN; V_psi=100._dp; END IF!ELSE IF (V_psi<-100._dp) THEN; V_psi=-100._dp; END IF
    !     Vext_q(i,j,k,io,s ) = V_psi
    !
    ! END DO

    ! IF (verbose) THEN
    !     BLOCK
    !         REAL(dp), ALLOCATABLE :: temparray(:,:,:)
    !         CHARACTER(50) :: filename
    !         PRINT*,"minval and maxval of Vext_q :", MINVAL( Vext_q (:,:,:,:,:,:) ), MAXVAL( Vext_q (:,:,:,:,:,:) )
    !         ! Get the external potential over orientations and print it
    !         ALLOCATE ( temparray ( grid%n_nodes(1),grid%n_nodes(2),grid%n_nodes(3) ),SOURCE=0._dp)
    !         CALL mean_over_orientations ( Vext_q (:,:,:,:,:,1) , temparray )
    !         temparray = temparray/fourpi
    !         filename = 'output/Vcoul.cube'
    !         CALL write_to_cube_file ( temparray , filename )
    !         filename = 'output/Vcoul_along-z.dat'
    !         CALL compute_z_density ( temparray , filename )
    !         DEALLOCATE(temparray)
    !     END BLOCK
    ! END IF

END SUBROUTINE
!===================================================================================================================================
