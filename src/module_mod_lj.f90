! In module mod_lj, we compute the array Vext_lj that contains the lennard jones part of the external potential
! It is perhaps not the best idea to have a module for that, but is simpler for me (Max) to code, at least at the beginning.

MODULE mod_lj

    USE external_potential, ONLY: Vext_lj
    USE precision_kinds,    ONLY: dp, i2b
    USE system,             ONLY: eps_mol,x_mol,y_mol,z_mol,&
                                  sig_solv,sig_mol,nb_species,id_mol,id_solv,chg_mol,chg_solv,eps_solv,&
                                  x_solv,y_solv,z_solv,spaceGrid,soluteSite,solventSite
    USE constants,          ONLY: fourpi
    USE input,              ONLY: input_line, verbose
    USE quadrature,         ONLY: angGrid, molRotGrid
    
    IMPLICIT NONE
    INTEGER(i2b), PRIVATE :: nb_id_mol, nb_id_solv ! number of different kinds of solvent sites or solute sites
    INTEGER(i2b), PRIVATE :: nfft1,nfft2,nfft3,nOm,nPsi
    REAL(dp),     PRIVATE :: Lx,Ly,Lz
    
    CONTAINS
    
        SUBROUTINE init
            IMPLICIT NONE
            nb_id_mol  = SIZE( chg_mol  ) ! total number of solute types
            nb_id_solv = SIZE( chg_solv ) ! total number of solvent types
            nfft1 =spaceGrid%n_nodes(1)
            nfft2 =spaceGrid%n_nodes(2)
            nfft3 =spaceGrid%n_nodes(3)
            Lx=spaceGrid%length(1)
            Ly=spacegrid%length(2)
            Lz=spaceGrid%length(3)
            nOm   =angGrid%n_angles
            nPsi  =molRotGrid%n_angles
            IF (.NOT. ALLOCATED(Vext_lj)) ALLOCATE( Vext_lj(nfft1,nfft2,nfft3,nOm,nPsi,nb_species) ,SOURCE=0._dp)
            CALL calculate
        END SUBROUTINE
        
        SUBROUTINE calculate
            IMPLICIT NONE
            INTEGER(i2b) :: i,j,k,s,v,u,a,b,c
            REAL(dp) :: x_grid,y_grid,z_grid ! coordinates of grid nodes
            REAL(dp) :: V_node,dx,dy,dz,sigij,epsij
            LOGICAL :: fullpdb

            ! compute lennard jones potential at each position and for each orientation, for each species => Vext_lj ( i , j , k , omega , species ) 
            ! we impose the simplification that only the first site of the solvent sites has a lennard jones WATER ONLY TODO
            ! test if this simplification is true and stop if not
            ! the test is done over the sigma lj. they're positive, so that the sum over all sigma is the same as the sigma of first site
            ! only if they're all zero but for first site
            IF( SUM( sig_solv(:) ) /= sig_solv(1) ) then
                PRINT*,'The solvent must have 1 Lennard Jones site only.'
                STOP
            END IF
            
            CALL only_one_lj_site ([x_solv(1),y_solv(1),z_solv(1)]) ! verify the solvant wears only one lennard jones site

            ! Test if the supercell is big enough considering the LJ range (given by sigma).
            ! at 2.5*sigma, the lj potential is only 0.0163*epsilon. Almost zero.
            ! It would have no sense to have a box dimension < 2.5 sigma
            IF( MIN(Lx,Ly,Lz)<= 2.5_dp*MAX(MAXVAL(sig_mol),MAXVAL(sig_solv))) THEN
                PRINT*,'The supercell is small. We replicate it in all directions. Vext calculation will be long.'
                fullpdb=.TRUE. ! much slower
!~                 STOP
            ELSE
                fullpdb=.FALSE. ! much faster
            END IF

            DO v=1,SIZE(solventSite)
                IF ( solventSite(v)%eps==0._dp ) CYCLE ! if solvent site wear no LJ
                
                DO CONCURRENT( i=1:nfft1, j=1:nfft2, k=1:nfft3, s=1:nb_species )
                    x_grid=REAL(i-1,dp)*spaceGrid%dl(1)
                    y_grid=REAL(j-1,dp)*spaceGrid%dl(2)
                    z_grid=REAL(k-1,dp)*spaceGrid%dl(3)

                    V_node=0.0_dp
                    solute: DO u=1,SIZE(soluteSite)
                        IF (soluteSite(u)%eps==0.0_dp) CYCLE
                        sigij = arithmetic_mean( solventSite(v)%sig, soluteSite(u)%sig )
                        epsij = geometric_mean(  solventSite(v)%eps, soluteSite(u)%eps )
                        IF (fullpdb) THEN
                            DO a=-1,1; DO b=-1,1; DO c=-1,1
                                dx =x_grid-soluteSite(u)%r(1)+a*Lx
                                dy =y_grid-soluteSite(u)%r(2)+b*Ly
                                dz =z_grid-soluteSite(u)%r(3)+c*Lz
                                V_node =V_node + vlj( epsij, sigij, SQRT(dx**2+dy**2+dz**2))
                                IF (V_node >= 100.0_dp) THEN ! limit maximum value of Vlj to 100 TODO magic number
                                    V_node = 100.0_dp
                                    EXIT solute
                                END IF
                            END DO; END DO; END DO
                        ELSE
                            dx =ABS(x_grid-soluteSite(u)%r(1)); DO WHILE (dx>Lx/2._dp); dx =ABS(dx-Lx); END DO
                            dy =ABS(y_grid-soluteSite(u)%r(2)); DO WHILE (dy>Ly/2._dp); dy =ABS(dy-Ly); END DO
                            dz =ABS(z_grid-soluteSite(u)%r(3)); DO WHILE (dz>Lz/2._dp); dz =ABS(dz-Lz); END DO
                            V_node =V_node + vlj( epsij, sigij, SQRT(dx**2+dy**2+dz**2))
                            IF (V_node >= 100.0_dp) THEN ! limit maximum value of Vlj to 100 TODO magic number
                                V_node = 100.0_dp
                                EXIT solute
                            END IF
                        END IF
                    END DO solute
              
                    Vext_lj (i,j,k,:,:,s) = V_node ! all angles are treated in the same time as the oxygen atom is not sensitive to rotation around omega and psi
                END DO
            END DO

        END SUBROUTINE calculate


        PURE FUNCTION vlj(eps,sig,d) ! v_lj(d) = 4ε[(σ/d)^12-(σ/d)^6]
            REAL(dp) :: vlj
            REAL(dp), INTENT(IN) :: eps, sig, d ! ε,σ,distance
            REAL(dp) :: div
            div = (sig/d)**6
            vlj = 4._dp*eps*div*(div-1._dp)
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
