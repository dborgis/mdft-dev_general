!===================================================================================================================================
SUBROUTINE external_potential_hard_walls
!===================================================================================================================================
! this SUBROUTINE computes the external potential created by several (0 to infty) hard walls
! The hard walls "plans" have coordinates read in input/dft.in
! For a brief description of the coordinates of a plan, see wikipedia:
! Step 1/ we read how many walls are to be found in the supercell.
! Step 2/ we read their coordinates
! Step 3/ we read their thinkness (2*radius) in dft.in
! Step 4/ we compute the external potential created by the walls, which depends upon the fluid radius

    use precision_kinds    ,ONLY: i2b, dp
    use module_input              ,ONLY: input_line, getinput
    use module_solvent, only: solvent
    use external_potential ,ONLY: Vext_total
    use hardspheres        ,ONLY: hs
    use module_grid, only: grid

    IMPLICIT NONE

    INTEGER(i2b) :: i,j,k,w,s,nfft1,nfft2,nfft3,iostatint
    INTEGER(i2b) :: species ! dummy between 1 and solvent(1)%nspec
    INTEGER(i2b) :: nwall ! number of hard walls in the supercell
    REAL(dp) :: dplan ! local distance between point M (grid point) and plan defining wall
    REAL(dp), ALLOCATABLE, DIMENSION(:)   :: thickness ! thickness of each wall (thickness = 2radius. Don't mix them up)
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: normal_vec ! normal vector of each plan defining wall
    REAL(dp), ALLOCATABLE, DIMENSION(:)   :: norm2_normal_vec ! norm of normal_vec
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: OA ! A is a point of coordinates OA(x,y,z) which is in the plan normal to normal_vec
    REAL(dp), ALLOCATABLE, DIMENSION(:)   :: dot_product_normal_vec_OA ! dummy
    REAL(dp) :: OM(3) ! coordinates of grid points

    nfft1 = grid%n_nodes(1)
    nfft2 = grid%n_nodes(2)
    nfft3 = grid%n_nodes(3)

! read how many hard walls are wanted in the supercell
    DO i = 1, SIZE(input_line)
        j = LEN('hard_wall_number')
        IF ( input_line (i) (1:j) == 'hard_wall_number' ) THEN
            READ ( input_line (i) (j+4:j+5) , * ) nwall
            IF ( nwall == 0 ) RETURN ! no hard wall
            ALLOCATE ( thickness (nwall) ,SOURCE=0._dp) ! thickness of each wall (thickness = 2radius. Don't mix them up)
            ALLOCATE ( normal_vec (3,nwall) ,SOURCE=0._dp) ! normal vector of each plan defining wall
            ALLOCATE ( norm2_normal_vec (nwall) ,SOURCE=0._dp) ! norm of normal_vec. dummy
            ALLOCATE ( OA (3,nwall) ,SOURCE=0._dp) ! A is a point of coordinates OA(1,2,3) which is in the plan normal to normal_vec
            ALLOCATE ( dot_product_normal_vec_OA (nwall) ,SOURCE=0._dp) ! dummy
            DO w = 1 , nwall
                READ ( input_line(i+w),*,IOSTAT=iostatint) thickness(w) , normal_vec(1:3,w) , OA (1:3,w) ! for each wall, read thickness, normal vector coordinates and point in plan
                IF(iostatint/=0) STOP "I could not read line containing thickness of the hard wall, its normal vector and position"
                norm2_normal_vec ( w ) = NORM2(normal_vec(:,w)) ! pretabulate norm of normal vec
                IF (norm2_normal_vec(w)==0._dp) STOP "Your hard wall is defined by a strange normal vector"
                dot_product_normal_vec_OA ( w ) = DOT_PRODUCT(-normal_vec(:,w) , OA(:,w))
            END DO
            EXIT
        END IF
    END DO

! be sure Vext_total is allocated
    IF (.NOT. ALLOCATED(Vext_total)) THEN
        ALLOCATE( Vext_total (nfft1,nfft2,nfft3,grid%no,solvent(1)%nspec) ,SOURCE=0._dp)
    END IF
! be sure radius(:) (solvent radius) is already computed
  !  IF(.NOT. ALLOCATED(hs)) STOP 'Molecular radius of solvent not allocated. critial stop in external_potential_hard_walls.f90'


    ! compute the potential potential
    DO CONCURRENT ( i=1:nfft1, j=1:nfft2, k=1:nfft3, s=1:solvent(1)%nspec, w=1:nwall)
        OM = [ REAL(i-1,dp)*grid%dl(1) , REAL(j-1,dp)*grid%dl(2) , REAL(k-1,dp)*grid%dl(3) ]
        dplan = ABS( DOT_PRODUCT( normal_vec(:,w),OM ) + dot_product_normal_vec_OA(w) )/ norm2_normal_vec(w) ! compute distance between grid point and plan
        IF (.NOT. ALLOCATED(hs)) THEN
          IF ( dplan <= 0.5_dp*thickness(w)  ) Vext_total(i,j,k,:,s) = HUGE(1.0_dp)
            IF (i*j*k*s*w==1) THEN
             PRINT*, 'WARNING : Radius of Solvent is not allocated in external_potential_hard_wall, it is set to zero'
            END IF
          ELSE
          IF ( dplan <= 0.5_dp*thickness(w) + hs(s)%R ) Vext_total(i,j,k,:,s) = HUGE(1.0_dp)
        ENDIF
    END DO

END SUBROUTINE external_potential_hard_walls
