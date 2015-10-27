!
! Module dedicated to computing the external potential.
! The external potential is often solvent-dependent, for instance through lennard-jones parameters.
!
module module_vext
    implicit none
    private
    public :: init_vext
contains
    ! This SUBROUTINE computes the external potential. It is one of the most time consuming routine.
    !Warning there are two ways of calculating the electrostatic potential (Poisson solver and point charge) you should always have one tag
    !T and one tag F for the electrostatic potential, if thera are 2 tags T, it is the last evaluation which counts, i.e Poisson solver.

    subroutine init_vext
        use precision_kinds, only: dp, i2b
        use module_input, only: getinput
        use module_solute, only: solute
        use module_solvent, only: solvent
        ! use external_potential, only: Vext_total, Vext_q, vextdef0
        use module_grid, only: grid

        IMPLICIT NONE

        integer :: nx, ny, nz, no, nfft(3), i,j,k,o,s, ns
        type errType
            integer :: i
            character(180) :: msg
        end type errType
        type (errType) :: er
        real(dp), parameter :: zerodp = 0._dp

        nfft = grid%n_nodes
        i = grid%nx
        j = grid%ny
        k = grid%nz
        o = grid%no
        ns = solvent(1)%nspec

        if (.not. allocated(solvent)) then
            print*, "Impossible de continuer. Dans module_external_potential > init_external_potential,"
            print*, "solvent type is not yet allocated. Ne devrait pas arriver"
            print*, "J'appelle l'initialisation de solvent% mais ce n'est pas la route usuelle."
            print*, "en fait je ferai ça plus tard. bisous"
            ! call solvent%init
            error stop
        end if

        do s=1,ns
            allocate ( solvent(s)%vext(nx,ny,nz,no), stat=er%i, errmsg=er%msg )
            if (er%i/=0) then
                print*,"Reported error is:", er%msg
                stop "I can't allocate solvent(:)%vext in init_external_potential.f90:init_external_potential"
            end if
        end do

        if (getinput%int('vext_hard_walls', defaultvalue=0, assert=">=0") /= 0) then
            do s=1,ns
                call vext_hard_walls
            end do
        end if

        CALL init_electrostatic_potential ! ELECTROSTATIC POTENTIAL
        stop "The problem is in init_electrostatic_potential"

        CALL vext_lennardjones ! LENNARD-JONES POTENTIAL

        IF (getinput%log('purely_repulsive_solute', defaultvalue=.false.)) CALL compute_purely_repulsive_potential ! r^-12 only
        IF (getinput%log('hard_sphere_solute', defaultvalue=.false.)) CALL compute_vext_hard_sphere     ! hard sphere
        IF (getinput%log('hard_cylinder_solute', defaultvalue=.false.)) CALL compute_vext_hard_cylinder ! hard cylinder
        ! IF (getinput%char('other_predefined_vext')=='vextdef0') CALL vextdef0

        CALL vext_total_sum ! compute total Vext(i,j,k,omega,s), the one used in the free energy functional
        CALL prevent_numerical_catastrophes
    end subroutine init_vext

    subroutine init_electrostatic_potential
        use precision_kinds, only: dp
        use module_input, only: getinput
        use module_solvent, only: solvent
        use module_grid, only: grid
        use module_fastpoissonsolver, only: init_fastpoissonsolver
        implicit none
        integer :: s, i, nx, ny, nz, no, ns
        character(180) :: myerrormsg
        logical :: is_direct, is_poisson
        is_direct = getinput%log('direct_sum', defaultvalue=.false.)
        is_poisson = getinput%log('poisson_solver', defaultvalue=.true.)

        if (is_poisson .and. is_direct) then
            STOP 'You ask for two different methods for computing the electrostatic potential: direct_sum and poisson'
        end if

        if (.not.allocated(solvent)) then
            print*, "solvent% is not allocated in module_vext init_electrostatic_potential"
            print*, "wierd idea"
            error stop
        end if

        !
        ! Definitions:
        ! The electrostatic potential, aka electric potential, is due to source charges only.
        ! The electrostatic energy is the application of the electric potential on target charges.
        ! The electrostatic energy density is the electrostatic energy per volume unit. It is the one we call vextq
        !
        do s=1,solvent(1)%nspec
            if (allocated(solvent(s)%vextq)) then
                print*, "I'm arriving in init_electrostatic_potential. There, solvent(s)%vextq is already allocated."
                print*, "That should not be the case."
                error stop
            end if
            allocate (solvent(s)%vextq(grid%nx,grid%ny,grid%nz,grid%no) ,source=0._dp, stat=i, errmsg=myerrormsg)
            if (i/=0) then
                print*, myerrormsg
                error stop "in init_electrostatic_potential"
            end if
        end do

        !
        ! Chose the Coulomb solver you want
        ! This may go elsewhere
        !

        ! DIRECT SUMMATION, pot = sum of qq'/r
        IF ( getinput%log('direct_sum', defaultvalue=.false.) ) call compute_vcoul_as_sum_of_pointcharges

        ! FAST POISSON SOLVER, -laplacian(pot) = solute charge density
        IF (getinput%log('fastpoissonsolver', defaultvalue=.true.)) call init_fastpoissonsolver
stop "apres l'appel à init_poisson_solver dans module_vext"
    end subroutine init_electrostatic_potential

    subroutine prevent_numerical_catastrophes
        use precision_kinds, only: dp
        use module_solvent, only: solvent
        use module_solute, only: solute
        use module_grid, only: grid
        use constants, only: qfact
        IMPLICIT NONE
        INTEGER :: i,j
        REAL(dp) :: d(3)
        REAL(dp), ALLOCATABLE :: dnn(:),epsnn(:),signn(:),qnn(:),rblock(:,:)
        LOGICAL :: iFoundTheNN

        i=SIZE(solute%site)
        j=SIZE(solvent(1)%site)
        ALLOCATE (dnn(i)) ; dnn=HUGE(1._dp)
        ALLOCATE (epsnn(i),signn(i),qnn(i),rblock(i,j) ,source=0._dp)

        ! liste tous les sites de soluté qui ont une charge mais pas de potentiel répulsif dessus:
        PRINT*,"===== Prevent numerical catastrophe ====="
        DO i=1, SIZE(solute%site)
            IF( solute%site(i)%q/=0. .and. (solute%site(i)%eps==0. .or. solute%site(i)%sig==0.) ) THEN
                iFoundTheNN=.false.
                PRINT*,"-- found one charged site that is purely attractive"
                PRINT*,"-- its q, eps and sig are ",solute%site(i)%q, solute%site(i)%eps, solute%site(i)%sig
                ! trouve le site répulsif de soluté le plus proche, sa distance (dnn), eps et sig lj (epsnn et signn)
                DO j=1, SIZE(solute%site)
                    IF (j==i) CYCLE
                    d= MODULO(    ABS(solute%site(i)%r(:) -solute%site(j)%r(:))    ,grid%length(:)/2._dp)
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
                        INTEGER :: ir
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
    end subroutine prevent_numerical_catastrophes

    subroutine vext_hard_walls
        ! this subroutine computes the external potential created by several (0 to infty) hard walls
        ! the hard walls "plans" have coordinates read in input/dft.in
        ! for a brief description of the coordinates of a plan, see wikipedia:
        ! step 1/ we read how many walls are to be found in the supercell.
        ! step 2/ we read their coordinates
        ! step 3/ we read their thinkness (2*radius) in dft.in
        ! step 4/ we compute the external potential created by the walls, which depends upon the fluid radius

        use precision_kinds, only: i2b, dp
        use module_input, only: input_line, getinput
        use module_solvent, only: solvent
        use hardspheres, only: hs
        use module_grid, only: grid

        implicit none

        integer :: i, j, k, w, s, iostatint
        integer :: nwall ! number of hard walls in the supercell
        real(dp) :: dplan ! local distance between point m (grid point) and plan defining wall
        real(dp) :: dist_hard
        real(dp), allocatable, dimension(:)   :: thickness ! thickness of each wall (thickness = 2radius. don't mix them up)
        real(dp), allocatable, dimension(:,:) :: normal_vec ! normal vector of each plan defining wall
        real(dp), allocatable, dimension(:)   :: norm2_normal_vec ! norm of normal_vec
        real(dp), allocatable, dimension(:,:) :: oa ! a is a point of coordinates oa(x,y,z) which is in the plan normal to normal_vec
        real(dp), allocatable, dimension(:)   :: dot_product_normal_vec_oa ! dummy
        real(dp) :: om(3) ! coordinates of grid points

        !
        ! Check everything we USE (from modules) is allocated and ready
        !
        if (.not. allocated(solvent)) then
            print*, "Arf. Dans module_vext, solvent n'est pas défini"
            error stop
        end if
        if(.not. allocated(hs)) then
            print*, "In module_vext I need hs(r) which is not allocated"
            error stop
        end if

        ! read how many hard walls are wanted in the supercell
        do i = 1, size(input_line)
            j = len('vext_hard_walls')
            if ( input_line (i) (1:j) == 'vext_hard_walls' ) then
                read ( input_line (i) (j+4:j+5) , * ) nwall
                if ( nwall == 0 ) return ! no hard wall
                allocate ( thickness (nwall) ,source=0._dp) ! thickness of each wall (thickness = 2radius. don't mix them up)
                allocate ( normal_vec (3,nwall) ,source=0._dp) ! normal vector of each plan defining wall
                allocate ( norm2_normal_vec (nwall) ,source=0._dp) ! norm of normal_vec. dummy
                allocate ( oa (3,nwall) ,source=0._dp) ! a is a point of coordinates oa(1,2,3) which is in the plan normal to normal_vec
                allocate ( dot_product_normal_vec_oa (nwall) ,source=0._dp) ! dummy
                do w = 1 , nwall
                    read ( input_line(i+w),*,iostat=iostatint) thickness(w) , normal_vec(1:3,w) , oa (1:3,w) ! for each wall, read thickness, normal vector coordinates and point in plan
                    if(iostatint/=0) stop "i could not read line containing thickness of the hard wall, its normal vector and position"
                    norm2_normal_vec ( w ) = norm2(normal_vec(:,w)) ! pretabulate norm of normal vec
                    if (norm2_normal_vec(w)==0._dp) stop "your hard wall is defined by a strange normal vector"
                    dot_product_normal_vec_oa ( w ) = dot_product(-normal_vec(:,w) , oa(:,w))
                end do
                exit
            end if
        end do

        do concurrent ( i=1:grid%nx, j=1:grid%ny, k=1:grid%nz, s=1:solvent(1)%nspec, w=1:nwall)
            om = real([i,j,k]-1,dp) *grid%dl
            dplan = abs( dot_product( normal_vec(:,w),om ) + dot_product_normal_vec_oa(w) )/ norm2_normal_vec(w) ! compute distance between grid point and plan
            if (.not. allocated(hs)) then
                dist_hard = 0.5_dp*thickness(w)
            else ! if allocated(hs)
                dist_hard = 0.5_dp*thickness(w) + hs(s)%r
            end if
            if (dplan <= dist_hard) then
                solvent(s)%vext(i,j,k,:) = huge(1.0_dp) ! no dependency over orientations
            endif
        end do

        deallocate (thickness, normal_vec, norm2_normal_vec, oa, dot_product_normal_vec_oa)
    end subroutine vext_hard_walls

    subroutine vext_lennardjones
        ! In module mod_lj, we compute the array Vext_lj that contains the lennard jones part of the external potential
        ! It is perhaps not the best idea to have a module for that, but is simpler for me (Max) to code, at least at the beginning.
        use precision_kinds, only: dp, i2b
        use module_solute, only: solute
        use module_solvent, only: solvent
        use module_grid, only: grid
        use module_input, only: input_line, verbose
        implicit none
        INTEGER :: nx, ny, nz, no, ns
        REAL(dp) :: lx, ly, lz
        real(dp), parameter :: fourpi=4._dp*acos(-1._dp)
        real(dp), parameter :: epsdp = epsilon(1._dp)
        INTEGER :: i,j,k,s,v,u,a,b,c,t, nsolutesite, nsolventsites_with_lj
        real(dp), allocatable :: x(:), y(:), z(:)
        REAL(dp) :: V_node, dx, dy, dz
        LOGICAL :: fullpbc
        integer, allocatable :: index_of_lj_site_in_solvent(:)
        real(dp), allocatable :: siguv(:), epsuv(:)

        if (.not. allocated(solvent)) then
            print*, "solvent should be allocated in vext_lennardjones"
            error stop
        end if
        if (.not. grid%isinitiated) then
            print*, "grid is not isinitiated in vext_lennardjones"
            error stop
        end if

        nx = grid%nx
        ny = grid%ny
        nz = grid%nz
        lx = grid%lx
        ly = grid%ly
        lz = grid%lz
        no = grid%no
        ns = solvent(1)%nspec

        ! compute lennard jones potential at each position and for each orientation, for each species => Vext_lj ( i , j , k , omega , species )
        ! we impose the simplification that only the first site of the solvent sites has a lennard jones WATER only TODO
        ! test if this simplification is true and stop if not
        ! the test is done over the sigma lj. they're positive, so that the sum over all sigma is the same as the sigma of first site
        ! only if they're all zero but for first site
        do s=1,ns
            nsolventsites_with_lj = count( abs(solvent(s)%site(:)%eps) > epsdp)
            if (nsolventsites_with_lj > 1) then
                print*, "in vext_lennardjones, I have only implemented the special case where the solvent wears only 1 LJ site"
                error stop
            end if
        end do



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

        if (solvent(1)%nspec /= 1) then
            print*, "Dans vext_lennardjones, je n'ai pas encore implemente le multi especes"
            error stop
        end if

        !
        ! For each solvent species, I want the index of the (unique) Lennard jones site
        !
        allocate (index_of_lj_site_in_solvent(solvent(1)%nspec), source=0)
        do s=1,solvent(1)%nspec
            do t=1,size(solvent(1)%site)
                if (solvent(1)%site(t)%eps > epsdp) then
                    index_of_lj_site_in_solvent(s) = t
                    cycle
                end if
            end do
        end do
        if (any(index_of_lj_site_in_solvent<=0)) then
            print*, "In module_vext, some site index is negative. Bisous"
            error stop
        end if

        !
        ! The unique LJ solvent site must be (see loop below) on top of a grid node
        !
        do s=1,solvent(1)%nspec
            t = index_of_lj_site_in_solvent(s)
            if (norm2(solvent(s)%site(t)%r)>epsdp) then
                print*, "Dans module_vext, vext_lennardjones, le calcul implique implicitement que l"
                print*, "unique site LJ du solvant est sur un grid node"
                print*, "on a ici comme coordonnées du site LJ du solvent:"
                print*, solvent(s)%site(t)%r
                print*, "qui est > epsdp"
                error stop
            end if
        end do



        nsolutesite = size(solute%site)
        allocate (siguv(nsolutesite), source=0._dp)
        allocate (epsuv(nsolutesite), source=0._dp)
        v = index_of_lj_site_in_solvent(1)
        ! pour tous les sites de soluté, tabule les eps et sig LJ avec l'unique site LJ du solvent
        siguv(1:nsolutesite) = [((solvent(1)%site(v)%sig + solute%site(u)%sig)/2._dp , u=1,nsolutesite)]
        epsuv(1:nsolutesite) = [(sqrt( solvent(1)%site(v)%eps * solute%site(u)%eps ) , u=1,nsolutesite)]
        ! do u=1,nsolutesite
        !     siguv(u) = (solvent(1)%site(v)%sig + solute%site(u)%sig )/2._dp
        !     epsuv(u) = sqrt( solvent(1)%site(v)%eps * solute%site(u)%eps )
        ! end do

        x(1:nx) = [(real(i-1,dp)*grid%dx ,i=1,nx)]
        y(1:ny) = [(real(j-1,dp)*grid%dy ,j=1,ny)]
        z(1:nz) = [(real(k-1,dp)*grid%dz, k=1,nz)]

        do concurrent (s=1:ns)
            do concurrent (i=1:nx, j=1:ny, k=1:nz)
                v_node=0._dp
                soluteloop: do u=1,nsolutesite
                    IF (abs(solute%site(u)%eps) <= epsdp) CYCLE
                    IF (fullpbc) THEN ! PARFAIT POUR OpenMPZ
                        DO a=-1,1; DO b=-1,1; DO c=-1,1
                            dx =x(i)-solute%site(u)%r(1)+a*lx ! we have implicitely the unique LJ site of the solvent that is on top of a grid node
                            dy =y(j)-solute%site(u)%r(2)+b*ly
                            dz =z(k)-solute%site(u)%r(3)+c*lz
                            V_node =V_node + vlj( epsuv(u), siguv(u), dx**2+dy**2+dz**2)
                            IF (V_node >= 100.0_dp) THEN ! limit maximum value of Vlj to 100 TODO magic number
                                V_node = 100.0_dp
                                EXIT soluteloop
                            END IF
                        END DO; END DO; END DO
                    ELSE
                        dx =ABS(x(i)-solute%site(u)%r(1)); DO WHILE (dx>lx/2._dp); dx =ABS(dx-lx); END DO
                        dy =ABS(y(j)-solute%site(u)%r(2)); DO WHILE (dy>ly/2._dp); dy =ABS(dy-ly); END DO
                        dz =ABS(z(k)-solute%site(u)%r(3)); DO WHILE (dz>lz/2._dp); dz =ABS(dz-lz); END DO
                        V_node =V_node + vlj( epsuv(u), siguv(u), dx**2+dy**2+dz**2 )
                        IF (V_node >= 100.0_dp) THEN ! limit maximum value of Vlj to 100 TODO magic number
                            V_node = 100.0_dp
                            EXIT soluteloop
                        END IF
                    END IF
                end do soluteloop
                solvent(s)%vext(i,j,k,:) = solvent(s)%vext(i,j,k,:) + V_node ! all angles are treated in the same time as the oxygen atom is not sensitive to rotation around omega and psi
            end do ! z
        end do ! solvent species
    contains
        pure function vlj(eps,sig,rsq) ! I hope this function is inlined by the compiler
            implicit none
            real(dp) :: vlj
            real(dp), intent(in) :: eps, sig, rsq
            real(dp) :: div
            real(dp), parameter :: epsdp=1._dp
            if (rsq<=epsdp) then
                vlj = huge(1._dp)
            else
                div = sig**6/rsq**3 ! rsq is a distance²
                vlj = 4._dp*eps*div*(div-1._dp)
            end if
        end function vlj
    end subroutine vext_lennardjones

    ! ! This subroutine computes the external potential induced by purely repulsive r^-12 spheres at each solute sites. You may be interesed by Dzubiella and Hansen, J. Chem. Phys. 121 (2004)
    ! ! The potential has the form V(r)=kbT*(r-R0)^-(12)
    !
    SUBROUTINE compute_purely_repulsive_potential
        !
        !     use precision_kinds     ,only: dp, i2b
        !     use module_input               ,only: getinput%dp, verbose
        !     use system              ,only: thermoCond
        !     use module_solute, only: solute
        !     use module_solvent, only: solvent
        !     use external_potential  ,only: Vext_total
        !     use constants           ,only: zerodp=>zero
        use module_grid, only: grid
        ! use module_quadrature, only: mean_over_orientations
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
        !                    psiloop: DO m=1, 1 ! nb_solvent_sites => FOR DZUBIELLA HANSEN only Oxygen atom is taken into account
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
        !             use constants, only : fourpi
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


    !
    ! Returns the direct sum of all qi*qj/rij
    ! That's very slow and does not consider periodic boundary conditions
    !
    SUBROUTINE compute_vcoul_as_sum_of_pointcharges
        use precision_kinds     ,only: dp, i2b
        use module_solute, only: solute
        use module_solvent, only: solvent
        use constants           ,only: qfact
        use module_input               ,only: verbose
        use module_grid, only: grid

        IMPLICIT NONE

        INTEGER(i2b)    :: i,j,k,o,p,m,n,s,io,no
        REAL(dp)        :: xgrid(3), xv(3), xuv(3), xuv2, V_psi
        REAL, PARAMETER :: hardShellRadiusSQ=1._dp ! charge pseudo radius **2   == Rc**2
        real, parameter :: epsdp = epsilon(1._dp)
        REAL(dp), DIMENSION( SIZE(solvent(1)%site), grid%no) :: xmod, ymod, zmod

        IF (.NOT. ALLOCATED(solvent(1)%vextq)) STOP "Vext_q should be allocated in SUBROUTINE compute_vcoul_as_sum_of_pointcharges"
        IF ( ANY(solvent(1)%vextq/=0.0_dp) ) STOP "Vext_q should be zero at the beginning of SUBROUTINE compute_vcoul_as_sum_of_pointcharges"

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

                            solvent(s)%vextq(i,j,k,io) = min( 100._dp, v_psi )

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


    ! this subroutine gets the partial vext ( hard sphere + hard wall + hard cylinder + purely repulsive)
    ! and adds to it Vext_lj and Vext_q.
    ! then it gives a upper value (100 kJ/mol) to vext_total.
    SUBROUTINE vext_total_sum

        use precision_kinds,    only: dp
        use module_solvent, only: solvent
        use module_grid, only: grid
        use module_input, only: verbose
        ! use module_quadrature, only: mean_over_orientations

        implicit none

        real(dp), parameter :: vmax = huge(1.0_dp)
        real(dp), parameter :: fourpi=4._dp*acos(-1._dp)
        real(dp), parameter :: zero=0._dp
        real(dp), parameter :: epsdp=epsilon(1._dp)
        integer :: nx, ny, nz, no, ns, s

        nx = grid%nx
        ny = grid%ny
        nz = grid%nz
        no = grid%no
        ns = solvent(1)%nspec

        ! vext_total = 0.0_dp
        ! vext is the sum over all external potentials
        ! note that purely repulsive and hard potentials are already included in vext
        where (solvent(1)%vext<vmax)
            solvent(1)%vext = solvent(1)%vext + solvent(1)%vextq
        end where

        if (all(abs(solvent(1)%vext)<=epsdp)) then
            print*, "WARNING: the external potential, vext, is uniform = 0 everywhere. OK?"
        end if

        ! caught nan and inf
        if ( any(solvent(1)%vext/=solvent(1)%vext) ) stop "there is a nan somewhere in vext_total."
        if ( any(abs(solvent(1)%vext)>huge(1.0_dp)) ) stop "there is an infinity somewhere in vext_total."

        where (solvent(1)%vext >vmax)
            solvent(1)%vext = vmax
        end where

        ! IF (verbose) THEN
        !     BLOCK
        !         CHARACTER(50) :: filename ! dummy
        !         REAL(dp), DIMENSION (nx,ny,nz) :: temparray ! dummy
        !         PRINT*, MINVAL ( Vext_total ) , ' < Vext_total < ' , MAXVAL(Vext_total)! give vext extrema to user for visual debugging
        !         PRINT*, MINVAL ( Vext_lj ) , ' < Vext_lj < ' , MAXVAL( Vext_lj )
        !         IF ( ALLOCATED (Vext_q) ) PRINT*, MINVAL(Vext_q), ' < Vext_q < ' , MAXVAL(Vext_q)
        !         IF ( ALLOCATED (Vext_hard_core) ) PRINT*, MINVAL(Vext_hard_core), ' < Vext_hard_core < ', MAXVAL(Vext_hard_core)
        !         ! mean over orientations and print
        !         CALL mean_over_orientations ( Vext_total ( : , : , : , : , : , 1 ) , temparray )
        !         temparray = temparray / (fourpi**2/(2.0_dp*molRotSymOrder))
        !         filename = 'output/Vext.cube' ! care when HS or multispec
        !         CALL write_to_cube_file ( temparray , filename )
        !         filename = 'output/z_Vext.dat' ! care when HS or multispec'
        !         CALL compute_z_density ( temparray , filename )
        !     END BLOCK
        ! END IF

    end subroutine vext_total_sum

end module module_vext
