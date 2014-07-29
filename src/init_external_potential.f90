! This SUBROUTINE computes the external potential. It is one of the most time consuming routine.
!Warning there are two ways of calculating the electrostatic potential (Poisson solver and point charge) you should always have one tag
!T and one tag F for the electrostatic potential, if thera are 2 tags T, it is the last evaluation which counts, i.e Poisson solver. 

SUBROUTINE init_external_potential

    USE precision_kinds     ,ONLY: dp , i2b
    USE input               ,ONLY: input_log, input_char
    USE system              ,ONLY: chg_solv, soluteSite, spaceGrid, nb_species
    USE external_potential  ,ONLY: Vext_total, Vext_q, vextdef0, vextdef1
    USE mod_lj              ,ONLY: init_lennardjones => init
    USE quadrature          ,ONLY: Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz, angGrid, molRotGrid
    USE constants           ,ONLY: zerodp=>zero
    
    IMPLICIT NONE
    
    INTEGER(i2b) :: nfft(3),i
    CHARACTER(180) :: j
    nfft = spaceGrid%n_nodes

    ALLOCATE( Vext_total(nfft(1),nfft(2),nfft(3),angGrid%n_angles,molRotGrid%n_angles,nb_species), SOURCE=zerodp ,STAT=i,ERRMSG=j)
        IF (i/=0) THEN; PRINT*,j; STOP "I can't allocate Vext_total in subroutine init_external_potential"; END IF
    
    CALL external_potential_hard_walls ! Hard walls

    CALL init_electrostatic_potential ! ELECTROSTATIC POTENTIAL

    CALL init_lennardjones ! LENNARD-JONES POTENTIAL

    ! r^-12 only
    IF (input_log('purely_repulsive_solute')) THEN
        CALL compute_purely_repulsive_potential
    END IF

    
    IF (input_log('hard_sphere_solute')) CALL compute_vext_hard_sphere     ! hard sphere
    IF (input_log('hard_cylinder_solute')) CALL compute_vext_hard_cylinder ! hard cylinder
    IF (input_log('personnal_vext')) CALL compute_vext_perso               ! personnal vext as implemented in personnal_vext.f90
    IF( input_char('other_predefined_vext')=='vextdef0') CALL vextdef0
    IF( input_char('other_predefined_vext')=='vextdef1') CALL vextdef1

    CALL vext_total_sum ! compute total Vext(i,j,k,omega), the one used in the free energy functional
    
    CALL prevent_numerical_catastrophes

STOP "OH MY GOD"
    
    
    CONTAINS

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE init_electrostatic_potential
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            IF ( input_log('direct_sum') .AND. input_log('poisson_solver')) THEN
                STOP 'You ask for two different methods for computing the electrostatic potential: direct_sum and poisson'
            END IF
            
            IF ( input_log('direct_sum') ) THEN ! Charges : treatment as point charges
                IF (.NOT. ALLOCATED(Vext_q)) THEN
                    BLOCK
                        INTEGER(i2b), DIMENSION(6) :: al
                        al(1) = nfft(1)
                        al(2) = nfft(2)
                        al(3) = nfft(3)
                        al(4) = angGrid%n_angles
                        al(5) = molRotGrid%n_angles
                        al(6) = nb_species
                        ALLOCATE ( vext_q (al(1),al(2),al(3),al(4),al(5),al(6)), SOURCE=zerodp )
                    END BLOCK
                END IF
                CALL compute_vcoul_as_sum_of_pointcharges( Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz )
            END IF

            ! Charges : Poisson solver
            IF (input_log('poisson_solver')) THEN
                BLOCK
                    REAL(dp), DIMENSION (nfft(1),nfft(2),nfft(3)) :: soluteChargeDensity
                    IF (.NOT. ALLOCATED(Vext_q) ) &
                        ALLOCATE ( Vext_q ( nfft(1),nfft(2),nfft(3),angGrid%n_angles,molRotGrid%n_angles,nb_species), SOURCE=zerodp)
                    CALL soluteChargeDensityFromSoluteChargeCoordinates (soluteChargeDensity)
                    CALL poissonSolver (soluteChargeDensity)
                    CALL vext_q_from_v_c (Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)
                END BLOCK
            END IF
        END SUBROUTINE
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE prevent_numerical_catastrophes
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            USE system    ,ONLY: spacegrid, solventsite
            USE constants ,ONLY: qfact
            IMPLICIT NONE
            INTEGER(i2b) :: i,j
            REAL(dp) :: d(3)
            REAL(dp), ALLOCATABLE :: dnn(:),epsnn(:),signn(:),qnn(:),rblock(:,:)
            LOGICAL :: iFoundTheNN
            
            i=SIZE(solutesite)
            j=SIZE(solventsite)
            ALLOCATE (dnn(i)) ; dnn=HUGE(1._dp)
            ALLOCATE (epsnn(i),signn(i),qnn(i),rblock(i,j) ,source=0._dp)

            ! liste tous les sites de soluté qui ont une charge mais pas de potentiel répulsif dessus:
            PRINT*,"---------- prevent numerical catastrophe"
            DO i=1, SIZE(soluteSite)
                IF( solutesite(i)%q/=0. .and. (solutesite(i)%eps==0. .or. solutesite(i)%sig==0.) ) THEN
                    iFoundTheNN=.false.
                    PRINT*,"-- found one charged site that is purely attractive"
                    PRINT*,"-- its q, eps and sig are ",solutesite(i)%q, solutesite(i)%eps, solutesite(i)%sig                    
                    ! trouve le site répulsif de soluté le plus proche, sa distance (dnn), eps et sig lj (epsnn et signn)
                    DO j=1, SIZE(solutesite)
                        IF (j==i) CYCLE
                        d= MODULO(    ABS(solutesite(i)%r(:) -solutesite(j)%r(:))    ,spacegrid%length(:)/2._dp)
                        IF (NORM2(d)<dnn(i)) THEN
                            dnn(i)= NORM2(d)
                            epsnn(i) = solutesite(j)%eps
                            signn(i) = solutesite(j)%sig
                            qnn(i) = solutesite(j)%q
                            iFoundTheNN=.true.
                        END IF
                    END DO
                    
                    IF (iFoundTheNN) THEN
                        PRINT*,"-- nearest repulsive site is at (Ang)",dnn(i)
                        PRINT*,"-- its q, eps and sig are:",qnn(i),epsnn(i),signn(i)
                    ELSE
                        STOP "I did not find any repulsive site!"
                    END IF


                    DO CONCURRENT (j=1:SIZE(solventsite), solutesite(i)%q*solventsite(j)%q<0._dp)
                    
                        PRINT*,"----- solutesite",i," va attirer solventsite n°",j
                        PRINT*,"----- q eps and sig of j are ",solventsite(j)%q,solventsite(j)%eps,solventsite(j)%sig
                    
                        BLOCK
                            INTEGER(i2b) :: ir
                            REAL(dp) :: r,vr,e,s,q
                            DO ir=1,50000
                                r=REAL(ir,dp)/10000._dp
                                q=solutesite(i)%q *solventsite(j)%q
                                e=SQRT(epsnn(i)   *solventsite(j)%eps)
                                s=0.5_dp*(signn(i)+solventsite(j)%sig)
                                vr=qfact*q/r + 4._dp*e*((s/(r+dnn(i)))**12-(s/(r+dnn(i)))**6)
                                IF (vr>0._dp) THEN
                                    rblock(i,j)=r
                                    EXIT
                                END IF
                            END DO
                            PRINT*,"-----rblock(i,j) =",rblock(i,j)
                        END BLOCK
                    
                    END DO
                    
                END IF
            END DO
        END SUBROUTINE prevent_numerical_catastrophes
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    
END SUBROUTINE init_external_potential
