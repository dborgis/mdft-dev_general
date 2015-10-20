module module_solute
    implicit none
    type :: solute_type
        character(130) :: name
        integer :: molrotsymorder
        integer :: nsite ! number of site of the solvent molecule
        integer :: nspec ! number of solvent species
        real(dp) :: monopole, dipole(3), quadrupole(3,3), octupole(3,3,3), hexadecapole(3,3,3,3)
        real(dp) :: diameter ! hard sphere diameter, for instance
        type (site_type), allocatable :: site(:)
    end type
    type (solute_type) :: solute
    private
    public :: solute_type, solute, read_solute, soluteChargeDensityFromSoluteChargeCoordinates

contains
    
    !> read solute atomic positions, charge, and lennard jones values in solute.in
    !! charge in electron units, sigma in Angstroms, epsilon in KJ/mol.
    subroutine read_solute

      use precision_kinds, ONLY: i2b,dp
      use module_input,    ONLY: input_line, getinput
      use module_grid, only: grid
    !  use periodic_table,  ONLY: init_periodic_table, ptable

      implicit none

      integer :: n,i,stat

    !  call init_periodic_table
      ! print *, ptable ( 1 ) % name
      ! open and test if input/solute.in is ok

      open (5, FILE='input/solute.in', STATUS='old', IOSTAT=stat)
      IF (stat /= 0) THEN
        PRINT*,'solute.in cannot be opened ! => STOP !'
        STOP
      END IF

      READ (5,*) ! comment line
      READ (5,*) solute%nsite ! total number of atom sites of the solute
      ALLOCATE(solute%site(solute%nsite))
      READ (5,*)
      DO n = 1, solute%nsite
        READ(5,*) i, solute%site(n)%q, solute%site(n)%sig, solute%site(n)%eps, &
                  solute%site(n)%lambda1, solute%site(n)%lambda2, solute%site(n)%r, solute%site(n)%Z
        END DO
      CLOSE (5)
      solute%site%q = solute%site%q * getinput%dp('solute_charges_scale_factor', defaultvalue=1._dp)

      block
        real(dp) :: solutexmin, solutexmax, soluteymin, soluteymax, solutezmin, solutezmax, &
                    solutediameterx, solutediametery, solutediameterz, solutesigmaljmax
        solutexmin = minval(solute%site%r(1))
        soluteymin = minval(solute%site%r(2))
        solutezmin = minval(solute%site%r(3))
        solutexmax = maxval(solute%site%r(1))
        soluteymax = maxval(solute%site%r(2))
        solutezmax = maxval(solute%site%r(3))
        solutediameterx = sqrt((solutexmax-solutexmin)**2)
        solutediametery = sqrt((soluteymax-soluteymin)**2)
        solutediameterz = sqrt((solutezmax-solutezmin)**2)
    !    solutesigmaljmax = maxval(solute%site%sig)
        grid%buffer_length = getinput%dp("buffer_length", defaultvalue=15._dp)
        grid%length(1) = solutediameterx +2.*grid%buffer_length !+2*solutesigmaljmax
        grid%length(2) = solutediametery +2.*grid%buffer_length !+2*solutesigmaljmax
        grid%length(3) = solutediameterz +2.*grid%buffer_length !+2*solutesigmaljmax
      end block

    end subroutine read_solute

    ! if user asks for it (tag 'translate_solute_to_center'), add Lx/2, Ly/2, Lz/2 to all solute coordinates
    subroutine mv_solute_to_center
        use module_input  ,only: getinput
        use module_grid, only: grid
        implicit none
        if( getinput%log( 'translate_solute_to_center', defaultvalue=.true. )) then
            solute%site%r(1) = solute%site%r(1) + grid%length(1)/2.0_dp
            solute%site%r(2) = solute%site%r(2) + grid%length(2)/2.0_dp
            solute%site%r(3) = solute%site%r(3) + grid%length(3)/2.0_dp
        end if
    end subroutine mv_solute_to_center

    subroutine assert_solute_inside
        use module_grid, only: grid
        implicit none
        integer :: i
        ! check if some positions are out of the supercell
        !j is a test tag. We loop over this test until every atom is in the box.
        ! This allows for instance, if a site is two boxes too far to still be ok.
        do concurrent( i=1:solute%nsite )
            solute%site(i)%r(1) = MODULO ( solute%site(i)%r(1) , grid%length(1) )
            solute%site(i)%r(2) = MODULO ( solute%site(i)%r(2) , grid%length(2) )
            solute%site(i)%r(3) = MODULO ( solute%site(i)%r(3) , grid%length(3) )
        end do
    end subroutine assert_solute_inside

    ! This subroutine uses the positions of the charges of the solute to extrapolate them on the grid nodes in order to
    ! get the solute charge density soluteChargeDensity (i,j,k).

    SUBROUTINE soluteChargeDensityFromSoluteChargeCoordinates (gridnode, gridlen, soluteChargeDensity)

        use precision_kinds, only: i2b,dp
        use module_input, only: verbose
        use mathematica, only: distToFloorNode, floorNode, ceilingNode

        IMPLICIT NONE

        INTEGER(i2b), INTENT(IN) :: gridnode(3)
        REAL(dp), INTENT(IN) :: gridlen(3)
        REAL(dp), DIMENSION(gridnode(1),gridnode(2),gridnode(3)), INTENT(OUT) :: soluteChargeDensity
        INTEGER(i2b)  :: s, m(3), p(3)
        REAL(dp)      :: volumElem ! elementary volume in Poisson Grid space
        REAL(dp)      :: r(3) ! coordinates of the charge in indicial coordinates
        REAL(dp)      :: wm(3), wp(3) ! weight associated to each index

        soluteChargeDensity = 0._dp

        ! extrapolate each solute point charge to grid nodes
        DO s = 1 , SIZE(solute%site)

            IF ( solute%site(s)%q == 0.0_dp ) CYCLE ! if the solute does not have charge, go to next solute

            r = distToFloorNode (gridnode,gridlen,solute%site(s)%r,.TRUE.)
            m = floorNode       (gridnode,gridlen,solute%site(s)%r,.TRUE.)
            p = ceilingNode     (gridnode,gridlen,solute%site(s)%r,.TRUE.)

            wp = r ! weights
            wm = 1._dp - r

            ! increase density accordingly
            soluteChargeDensity (m(1),m(2),m(3)) = soluteChargeDensity (m(1),m(2),m(3)) + solute%site(s)%q * wm(1) * wm(2) * wm(3)
            soluteChargeDensity (p(1),m(2),m(3)) = soluteChargeDensity (p(1),m(2),m(3)) + solute%site(s)%q * wp(1) * wm(2) * wm(3)
            soluteChargeDensity (m(1),p(2),m(3)) = soluteChargeDensity (m(1),p(2),m(3)) + solute%site(s)%q * wm(1) * wp(2) * wm(3)
            soluteChargeDensity (m(1),m(2),p(3)) = soluteChargeDensity (m(1),m(2),p(3)) + solute%site(s)%q * wm(1) * wm(2) * wp(3)
            soluteChargeDensity (p(1),p(2),m(3)) = soluteChargeDensity (p(1),p(2),m(3)) + solute%site(s)%q * wp(1) * wp(2) * wm(3)
            soluteChargeDensity (p(1),m(2),p(3)) = soluteChargeDensity (p(1),m(2),p(3)) + solute%site(s)%q * wp(1) * wm(2) * wp(3)
            soluteChargeDensity (m(1),p(2),p(3)) = soluteChargeDensity (m(1),p(2),p(3)) + solute%site(s)%q * wm(1) * wp(2) * wp(3)
            soluteChargeDensity (p(1),p(2),p(3)) = soluteChargeDensity (p(1),p(2),p(3)) + solute%site(s)%q * wp(1) * wp(2) * wp(3)
        END DO

        volumElem = PRODUCT(gridlen/REAL(gridnode,dp))
        soluteChargeDensity = soluteChargeDensity / volumElem ! charge density is in charge per unit volume


        IF (verbose) THEN
            BLOCK
                CHARACTER(50) :: filename
                filename='output/soluteChargeDensity.cube'
                CALL write_to_cube_file ( soluteChargeDensity, filename  )
            END BLOCK
        END IF

    END SUBROUTINE soluteChargeDensityFromSoluteChargeCoordinates

end module module_solute
