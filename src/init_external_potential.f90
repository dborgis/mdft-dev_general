! This SUBROUTINE computes the external potential. It is one of the most time consuming routine.
!Warning there are two ways of calculating the electrostatic potential (Poisson solver and point charge) you should always have one tag
!T and one tag F for the electrostatic potential, if thera are 2 tags T, it is the last evaluation which counts, i.e Poisson solver.

SUBROUTINE init_external_potential

    USE precision_kinds     ,only: dp, i2b
    USE input               ,only: input_log, input_char
    USE system              ,only: solute, nb_species, solvent
    USE external_potential  ,only: Vext_total, Vext_q, vextdef0, vextdef1
    USE mod_lj              ,only: init_lennardjones => init
    USE fastPoissonSolver   ,only: start_fastPoissonSolver => init
    use module_grid, only: grid

    IMPLICIT NONE

    INTEGER(i2b) :: nx, ny, nz, no, nfft(3), i,j,k,o,s, ns
    type errType
        integer(i2b) :: i
        character(180) :: msg
    end type errType
    type (errType) :: er
    real(dp), parameter :: zerodp = 0._dp

    nfft = GRID%n_nodes
    i = GRID%n_nodes(1)
    j = GRID%n_nodes(2)
    k = GRID%n_nodes(3)
    o = GRID%no
    s = size(solvent)

    allocate ( Vext_total(i,j,k,o,s), SOURCE=zerodp ,STAT=er%i,ERRMSG=er%msg)
    if (er%i/=0) then
        print*,"Reported error is:", er%msg
        stop "I can't allocate Vext_total in init_external_potential.f90:init_external_potential"
    end if

    do s = 1, size(solvent)
        allocate ( solvent(s)%vext(i,j,k,o), stat=er%i, errmsg=er%msg )
        if (er%i/=0) then
            print*,"Reported error is:", er%msg
            stop "I can't allocate solvent(:)%vext in init_external_potential.f90:init_external_potential"
        end if
    end do

    CALL external_potential_hard_walls ! HARD WALLS

    CALL init_electrostatic_potential ! ELECTROSTATIC POTENTIAL

    CALL init_lennardjones ! LENNARD-JONES POTENTIAL

    IF (input_log('purely_repulsive_solute')) CALL compute_purely_repulsive_potential ! r^-12 only
    IF (input_log('hard_sphere_solute')) CALL compute_vext_hard_sphere     ! hard sphere
    IF (input_log('hard_cylinder_solute')) CALL compute_vext_hard_cylinder ! hard cylinder
    IF (input_log('personnal_vext')) CALL compute_vext_perso               ! personnal vext as implemented in personnal_vext.f90
    IF (input_char('other_predefined_vext')=='vextdef0') CALL vextdef0
    IF (input_char('other_predefined_vext')=='vextdef1') CALL vextdef1

    CALL vext_total_sum ! compute total Vext(i,j,k,omega,s), the one used in the free energy functional

    CALL prevent_numerical_catastrophes

!~ STOP "OH MY GOD"


    CONTAINS

!===================================================================================================================================

        SUBROUTINE init_electrostatic_potential

            IMPLICIT NONE
            INTEGER(i2b) :: i, nx, ny, nz, no, ns
            CHARACTER(180) :: j

            IF ( input_log('direct_sum') .AND. input_log('poisson_solver')) THEN
                STOP 'You ask for two different methods for computing the electrostatic potential: direct_sum and poisson'
            END IF

            ! ALLOCATE THE ELECTROSTATIC POTENTIAL
            nx = GRID%nx
            ny = GRID%ny
            nz = GRID%nz
            no = GRID%no
            ns = nb_species
            allocate ( vext_q (nx,ny,nz,no,ns) , SOURCE=zerodp ,STAT=i, ERRMSG=j)
            IF (i/=0) THEN
                PRINT j
                STOP "I cant allocate vext_q in subroutine init_electrostatic_potential."
            END IF

            ! DIRECT SUMMATION, pot = sum of qq'/r
            IF ( input_log('direct_sum') ) THEN
                CALL compute_vcoul_as_sum_of_pointcharges
            END IF

            ! FAST POISSON SOLVER, -laplacian(pot) = solute charge density
            IF (input_log('poisson_solver')) THEN
                CALL start_fastPoissonSolver
            END IF
        END SUBROUTINE

!===================================================================================================================================

        SUBROUTINE prevent_numerical_catastrophes

            USE system    ,only: solvent
            use module_grid, only: grid
            USE constants ,only: qfact
            IMPLICIT NONE
            INTEGER(i2b) :: i,j
            REAL(dp) :: d(3)
            REAL(dp), ALLOCATABLE :: dnn(:),epsnn(:),signn(:),qnn(:),rblock(:,:)
            LOGICAL :: iFoundTheNN

            i=SIZE(solute%site)
            j=SIZE(solvent(1)%site)
            ALLOCATE (dnn(i)) ; dnn=HUGE(1._dp)
            ALLOCATE (epsnn(i),signn(i),qnn(i),rblock(i,j) ,source=0._dp)

            ! liste tous les sites de soluté qui ont une charge mais pas de potentiel répulsif dessus:
            PRINT*,"---------- prevent numerical catastrophe"
            DO i=1, SIZE(solute%site)
                IF( solute%site(i)%q/=0. .and. (solute%site(i)%eps==0. .or. solute%site(i)%sig==0.) ) THEN
                    iFoundTheNN=.false.
                    PRINT*,"-- found one charged site that is purely attractive"
                    PRINT*,"-- its q, eps and sig are ",solute%site(i)%q, solute%site(i)%eps, solute%site(i)%sig
                    ! trouve le site répulsif de soluté le plus proche, sa distance (dnn), eps et sig lj (epsnn et signn)
                    DO j=1, SIZE(solute%site)
                        IF (j==i) CYCLE
                        d= MODULO(    ABS(solute%site(i)%r(:) -solute%site(j)%r(:))    ,GRID%length(:)/2._dp)
                        IF (NORM2(d)<dnn(i)) THEN
                            dnn(i)= NORM2(d)
                            epsnn(i) = solute%site(j)%eps
                            signn(i) = solute%site(j)%sig
                            qnn(i) = solute%site(j)%q
                            iFoundTheNN=.true.
                        END IF
                    END DO

                    IF (iFoundTheNN) THEN
                        PRINT*,"-- nearest repulsive site is at (Ang)",dnn(i)
                        PRINT*,"-- its q, eps and sig are:",qnn(i),epsnn(i),signn(i)
                    ELSE
                        STOP "I did not find any repulsive site!"
                    END IF


                    DO CONCURRENT (j=1:SIZE(solvent(1)%site), solute%site(i)%q*solvent(1)%site(j)%q<0._dp)
                        PRINT*,"----- solute%site",i," va attirer solvent(1)%site n°",j
                        PRINT*,"----- q eps and sig of j are ",solvent(1)%site(j)%q,solvent(1)%site(j)%eps,solvent(1)%site(j)%sig
                        BLOCK
                            INTEGER(i2b) :: ir
                            REAL(dp) :: r,vr,e,s,q
                            DO ir=1,50000
                                r=REAL(ir,dp)/10000._dp
                                q=solute%site(i)%q *solvent(1)%site(j)%q
                                e=SQRT(epsnn(i)   *solvent(1)%site(j)%eps)
                                s=0.5_dp*(signn(i)+solvent(1)%site(j)%sig)
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

!===================================================================================================================================

END SUBROUTINE init_external_potential
