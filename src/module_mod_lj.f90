! In module mod_lj, we compute the array Vext_lj that contains the lennard jones part of the external potential
! It is perhaps not the best idea to have a module for that, but is simpler for me (Max) to code, at least at the beginning.

MODULE mod_lj

    USE external_potential, ONly: Vext_lj
    USE precision_kinds,    ONly: dp, i2b
    USE system,             ONly: solute, solvent
    use module_grid, only: grid
    USE constants,          ONly: fourpi
    USE input,              ONly: input_line, verbose

    IMPLICIT NONE
    INTEGER(i2b), PRIVATE :: nx, ny, nz, no, ns
    REAL(dp),     PRIVATE :: lx, ly, lz

    CONTAINS

        SUBROUTINE init
            IMPLICIT NONE
            nx = GRID%nx
            ny = GRID%ny
            nz = GRID%nz
            lx = GRID%lx
            ly = GRID%ly
            lz = GRID%lz
            no = GRID%no
            ns = solvent(1)%nspec
            IF (.NOT. ALLOCATED(Vext_lj)) ALLOCATE( Vext_lj(nx,ny,nz,no,ns) ,SOURCE=0._dp)
            call calculate
        END SUBROUTINE

        SUBROUTINE calculate
            IMPLICIT NONE
            INTEGER(i2b) :: i,j,k,s,v,u,a,b,c
            REAL(dp) :: x_grid,y_grid,z_grid ! coordinates of grid nodes
            REAL(dp) :: V_node,dx,dy,dz,sigij,epsij
            LOGICAL :: fullpbc

            ! compute lennard jones potential at each position and for each orientation, for each species => Vext_lj ( i , j , k , omega , species )
            ! we impose the simplification that only the first site of the solvent sites has a lennard jones WATER ONly TODO
            ! test if this simplification is true and stop if not
            ! the test is done over the sigma lj. they're positive, so that the sum over all sigma is the same as the sigma of first site
            ! only if they're all zero but for first site
            IF ( SIZE(solvent(1)%site)>1 )THEN
                IF ( ANY( solvent(1)%site(2:)%sig/=0._dp )) STOP "The solvent must have only 1 LJ site for now."
            END IF

            CALL only_one_lj_site (solvent(1)%site(1)%r) ! verify the solvant wears only one lennard jones site

            ! Test if the supercell is big enough considering the LJ range (given by sigma).
            ! at 2.5*sigma, the lj potential is only 0.0163*epsilon. Almost zero.
            ! It would have no sense to have a box dimension < 2.5 sigma
            IF( MIN(lx,ly,lz)<= 2.5_dp*MAX(MAXVAL(solute%site%sig),MAXVAL(solvent(1)%site%sig))) THEN
                PRINT*,'The supercell is small. We replicate it in all directions. Vext calculation will be long.'
                fullpbc=.TRUE. ! much slower
!~                 STOP
            ELSE
                fullpbc=.FALSE. ! much faster
            END IF

            DO v=1,SIZE(solvent(1)%site)
                IF ( solvent(1)%site(v)%eps==0._dp ) CYCLE ! if solvent site wear no LJ

                DO CONCURRENT( i=1:nx, j=1:ny, k=1:nz, s=1:ns )
                    x_grid=REAL(i-1,dp)*GRID%dl(1)
                    y_grid=REAL(j-1,dp)*GRID%dl(2)
                    z_grid=REAL(k-1,dp)*GRID%dl(3)

                    V_node=0.0_dp
                    soluteloop: DO u=1,SIZE(solute%site)
                        IF (solute%site(u)%eps==0.0_dp) CYCLE
                        sigij = arithmetic_mean( solvent(1)%site(v)%sig, solute%site(u)%sig )
                        epsij = geometric_mean(  solvent(1)%site(v)%eps, solute%site(u)%eps )
                        IF (fullpbc) THEN
                            DO a=-1,1; DO b=-1,1; DO c=-1,1
                                dx =x_grid-solute%site(u)%r(1)+a*lx
                                dy =y_grid-solute%site(u)%r(2)+b*ly
                                dz =z_grid-solute%site(u)%r(3)+c*lz
                                V_node =V_node + vlj( epsij, sigij, NORM2([dx,dy,dz]))
                                IF (V_node >= 100.0_dp) THEN ! limit maximum value of Vlj to 100 TODO magic number
                                    V_node = 100.0_dp
                                    EXIT soluteloop
                                END IF
                            END DO; END DO; END DO
                        ELSE
                            dx =ABS(x_grid-solute%site(u)%r(1)); DO WHILE (dx>lx/2._dp); dx =ABS(dx-lx); END DO
                            dy =ABS(y_grid-solute%site(u)%r(2)); DO WHILE (dy>ly/2._dp); dy =ABS(dy-ly); END DO
                            dz =ABS(z_grid-solute%site(u)%r(3)); DO WHILE (dz>lz/2._dp); dz =ABS(dz-lz); END DO
                            V_node =V_node + vlj( epsij, sigij, SQRT(dx**2+dy**2+dz**2))
                            IF (V_node >= 100.0_dp) THEN ! limit maximum value of Vlj to 100 TODO magic number
                                V_node = 100.0_dp
                                EXIT soluteloop
                            END IF
                        END IF
                    END DO soluteloop

                    Vext_lj (i,j,k,:,s) = V_node ! all angles are treated in the same time as the oxygen atom is not sensitive to rotation around omega and psi
                END DO
            END DO

        END SUBROUTINE calculate


        PURE FUNCTION vlj(eps,sig,d) ! v_lj(d) = 4ε[(σ/d)^12-(σ/d)^6]
            REAL(dp) :: vlj
            REAL(dp), INTENT(IN) :: eps, sig, d ! ε,σ,distance
            REAL(dp) :: div
            IF (d <= EPSILON(1._dp)) THEN
                vlj = HUGE(1._dp)
            ELSE
                div = (sig/d)**6
                vlj = 4._dp*eps*div*(div-1._dp)
            END IF
        END FUNCTION vlj

        PURE FUNCTION arithmetic_mean( A, B)
            ! = sum_i^N a_i/N
            REAL(dp) :: arithmetic_mean
            REAL(dp), intent(in) :: A, B
            arithmetic_mean = (A+B)/2._dp
        END FUNCTION arithmetic_mean

        PURE FUNCTION geometric_mean( A, B)
            ! = (product_i^N a_i)^(1/N)
            REAL(dp) :: geometric_mean
            REAL(dp), intent(in) :: A, B
            geometric_mean = SQRT(A*B)
        END FUNCTION geometric_mean

        SUBROUTINE only_one_lj_site (coo)
            IMPLICIT NONE
            REAL(dp), DIMENSION(3), INTENT(IN) :: coo
            IF(ANY(coo/=0._dp)) THEN
                PRINT*,'For now, the solvent must have only one Lennard-Jones site.'
                PRINT*,'This site should have coordinates 0. 0. 0., i.e., it must be on top of a grid node.'
                PRINT*,'Sadly, its coordinates are now ',coo
                STOP
            END IF
        END SUBROUTINE

END MODULE mod_lj
