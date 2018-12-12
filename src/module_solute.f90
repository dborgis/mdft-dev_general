module module_solute
use precision_kinds, only: dp
    use system, only: site_type
    implicit none
    type :: solute_type
        character(130) :: name
        character(80) :: nature, solvent_coupling  ! 'classical' point charge model or QM
        integer :: molrotsymorder
        integer :: nsite ! number of site of the solvent molecule
        integer :: nspec ! number of solvent species
        real(dp) :: monopole, dipole(3), quadrupole(3,3), octupole(3,3,3), hexadecapole(3,3,3,3)
        real(dp) :: diameter ! hard sphere diameter, for instance
        complex(dp), allocatable :: sigma_k(:,:,:) ! charge factor
        type (site_type), allocatable :: site(:)
        complex(dp), allocatable :: electron_density_k(:,:,:)
        real(dp) :: electronic_charge, nuclei_charge
    end type
    type (solute_type), protected :: solute
    private
    public :: solute_type, solute, read_solute, soluteChargeDensityFromSoluteChargeCoordinates,getreciprocalsolutechargedensity

contains

    !> read solute atomic positions, charge, and lennard jones values in solute.in
    !! charge in electron units, sigma in Angstroms, epsilon in KJ/mol.
    subroutine read_solute

      use iso_c_binding
      use precision_kinds, only: dp, i4b
      use module_input,    ONLY: getinput
      use module_grid, only: grid
      use module_cubefiles


      implicit none

      integer :: n,i,j,k,stat
      real(dp) :: sum_charge
      character(3) :: QM_flag
      character(80) :: filename
      real(dp), allocatable :: solute_electron_density(:,:,:)
      character(130) ::  solute_electron_density_filename
      integer(i4b) :: plan_forward  !for fftw3

      include "fftw3.f03"

    !  call init_periodic_table
      ! print *, ptable ( 1 ) % name
      ! open and test if solute.in is ok

      open (5, FILE='solute.in', STATUS='old', IOSTAT=stat, action="read")
      IF (stat /= 0) THEN
        PRINT*,'solute.in cannot be opened ! => STOP !'
        STOP
      END IF
      ! copy solute.in in the output folder for further reference
      call execute_command_line ("cp solute.in output/.")

      READ (5,*) QM_flag
      if ( QM_flag == '#QM' ) solute%nature = 'QM'
      READ (5,*) solute%nsite ! total number of atom sites of the solute
      ALLOCATE(solute%site(solute%nsite))
      READ (5,*)
      DO n = 1, solute%nsite
        READ(5,*) i, solute%site(n)%q, solute%site(n)%sig, solute%site(n)%eps, solute%site(n)%r, solute%site(n)%Z
      END DO

      ! if solute is described quantum mechanically, one needs to specify the electron_densty_filename
      ! and upload the corresponding electron density
GET_QUANTUM_ELECTRON_DENSITY: Block
      if( solute%nature == 'QM') then
        write(*,*) 'the solute is described quantum-mechanically by its electron density'
        read(5,*)  ! comment line
        read(5,*)  solute_electron_density_filename
        write(*,*) 'The electron_density_filename : ', solute_electron_density_filename,' was read'
        read(5,*)  ! comment line
        read(5,*)  solute%solvent_coupling
        print*, solute%solvent_coupling

      allocate ( solute_electron_density( grid%nx, grid%ny, grid%nz ) )

      CALL Read_Gaussian_cube_file( solute_electron_density, solute_electron_density_filename, solute%electronic_charge )
      print*,'cube file read: DONE. Total electronic charge =',solute%electronic_charge

      allocate( solute%electron_density_k( grid%nx/2 + 1, grid%ny, grid%nz) )

      ! DO FFT forward
      select case(dp)
      case(c_double)
        call dfftw_plan_dft_r2c_3d (plan_forward, grid%nx, grid%ny, grid%nz, solute_electron_density, solute%electron_density_k, fftw_estimate)
        call dfftw_execute(plan_forward)
      case(c_float)
        call sfftw_plan_dft_r2c_3d (plan_forward, grid%nx, grid%ny, grid%nz, solute_electron_density, solute%electron_density_k, fftw_estimate)
        call sfftw_execute(plan_forward)
      end select

      solute%electron_density_k = solute%electron_density_k*grid%dv
!
      deallocate( solute_electron_density )

     end if !

END Block GET_QUANTUM_ELECTRON_DENSITY

      CLOSE (5)

      if( solute%nature /= 'QM') &
              solute%site%q = solute%site%q * getinput%dp('solute_charges_scale_factor', defaultvalue=1._dp)

      block
        real(dp) :: solutexmin, solutexmax, soluteymin, soluteymax, solutezmin, solutezmax, &
                    solutediameterx, solutediametery, solutediameterz!, solutesigmaljmax
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
        ! grid%buffer_length = getinput%dp("buffer_length", defaultvalue=15._dp)
        ! grid%length(1) = solutediameterx +2.*grid%buffer_length !+2*solutesigmaljmax
        ! grid%length(2) = solutediametery +2.*grid%buffer_length !+2*solutesigmaljmax
        ! grid%length(3) = solutediameterz +2.*grid%buffer_length !+2*solutesigmaljmax
      end block

      !call align_solute_to_longer_diagonal_of_supercell ! align the solute so that its largest distance is along the longest diagonal of the supercell
      CALL mv_solute_to_center ! if user wants all the sites to be translated to the center of the box, ie by Lx/2, Ly/2, Lz/2
      CALL print_solute_xsf ! Print periodic XSF file to be read by VMD or equivalent
    end subroutine read_solute


    subroutine align_solute_to_longer_diagonal_of_supercell
        ! We want to rotate the solute so that it is aligned with the largest diagonal of the supercell.
        ! Works only for cubic cells right now, while it is desirable to have it work for orthorhombic cells anytime soon.
        use module_input, only: getinput
        implicit none
        integer :: i, j, farthersite1, farthersite2
        integer :: nsite ! number of solute sitess
        real(dp) :: costheta, sintheta, x, y, z, newx, newy, newz, u, v, w, r(3), vector_product(3), distance_between_sites, biggest_distance
        nsite = size( solute%site ) ! number of solute sitess
        if( getinput%log( 'align_solute_with_diagonal', defaultvalue=.true. )) then
            ! start by finding the two sites of the solute that are the farthest from each other.
            ! This defines the "largest size" of our solute.
            ! These two sites are called fartersite1 and farthersite2.
            ! The largest size is called "biggest_distance"
            if( nsite > 1) then ! size(solute%site) returns the number of sites in the solute
                farthersite1 = 0
                farthersite2 = 0
                biggest_distance = 0._dp
                do i = 1, nsite -1
                    do j = i+1, nsite
                        distance_between_sites =  norm2( solute%site(i)%r - solute%site(j)%r )
                        if( distance_between_sites > biggest_distance ) then
                            farthersite1 = i
                            farthersite2 = j
                            biggest_distance = distance_between_sites
                        end if
                    end do
                end do
            end if
            ! Translate the solute so that the farthersite1 is at the origin
            r(1:3) = solute%site(farthersite1)%r
            do i = 1, nsite
                solute%site(i)%r = solute%site(i)%r - r
            end do
            ! The rotation angle is given by the dot product of the vector between coordinates of site 2 and 1, and the vector <1 1 1> (the longest diagonal)
            u = 1._dp
            v = 1._dp
            w = 1._dp
            r = solute%site(farthersite2)%r
            costheta = dot_product( solute%site(farthersite2)%r, [u,v,w] ) / norm2( solute%site(farthersite2)%r ) / norm2([u,v,w])
            sintheta = sqrt( 1._dp - costheta**2 )
            ! The rotation axis is the one given by the vector product between <1 1 1> and the coordinate vector
            ! see http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/   section 5.1    for the rotation matrix
            vector_product = [ r(2)*w - r(3)*v,&
                               r(3)*u - r(1)*w,&
                               r(1)*v - r(2)*u ]
            u = vector_product(1)
            v = vector_product(2)
            w = vector_product(3)
            u = u/norm2([u,v,w])
            v = v/norm2([u,v,w])
            w = w/norm2([u,v,w])
            do i = 1, nsite
                x = solute%site(i)%r(1)
                y = solute%site(i)%r(2)
                z = solute%site(i)%r(3)
                newx = ( u*(u*x+v*y+w*z)*(1._dp-costheta)+(u**2+v**2+w**2)*x*costheta+sqrt(u**2+v**2+w**2)*(-w*y+v*z)*sintheta )/sqrt(u**2+v**2+w**2)
                newy = ( v*(u*x+v*y+w*z)*(1._dp-costheta)+(u**2+v**2+w**2)*y*costheta+sqrt(u**2+v**2+w**2)*( w*x-u*z)*sintheta )/sqrt(u**2+v**2+w**2)
                newz = ( w*(u*x+v*y+w*z)*(1._dp-costheta)+(u**2+v**2+w**2)*z*costheta+sqrt(u**2+v**2+w**2)*(-v*x+u*y)*sintheta )/sqrt(u**2+v**2+w**2)
                solute%site(i)%r = [newx, newy, newz]
            end do
        end if
    end subroutine align_solute_to_longer_diagonal_of_supercell

    ! if user asks for it (tag 'translate_solute_to_center'), add Lx/2, Ly/2, Lz/2 to all solute coordinates
    subroutine mv_solute_to_center
        use precision_kinds, only: dp
        use module_input  ,only: getinput
        use module_grid, only: grid
        implicit none
        integer :: i, d
        real(dp) :: coo_midbox_x, coo_midbox_y, coo_midbox_z
        real(dp) :: solute_mean_x, solute_mean_y, solute_mean_z
        if (.not.grid%isinitiated) then
            print*, "The derived type grid is not allocated in mv_solute_to_center (in module_solute). It should"
            stop "in module_solute, line 87"
        end if
        if( getinput%log( 'translate_solute_to_center', defaultvalue=.true. )) then
            ! what are the coordinates of the middle of the simulation box ?
            coo_midbox_x = grid%lx/2._dp
            coo_midbox_y = grid%ly/2._dp
            coo_midbox_z = grid%lz/2._dp
            ! what are the coordinates of the center of mass of the solute?
            ! we don't now the mass of the sites as of mdft-dev 2016-07-20.
            ! we'll thus say all sites have the same mass.
            ! Thus, the coordinates of the center of mass is the mean coordinate
            solute_mean_x = sum(solute%site(:)%r(1)) / real(solute%nsite,dp)
            solute_mean_y = sum(solute%site(:)%r(2)) / real(solute%nsite,dp)
            solute_mean_z = sum(solute%site(:)%r(3)) / real(solute%nsite,dp)
            ! Now, we translate this center of mass to the center of the box
            ! by shifting all coordinates (and thus the center of mass).
            ! Removing solute_mean_x, y and z translates the center of mass to coordinate 0,0,0
            ! Then add coo_midbox_x, y and z to translate the center of mass to center of the box.
            solute%site%r(1) = solute%site%r(1) + coo_midbox_x - solute_mean_x
            solute%site%r(2) = solute%site%r(2) + coo_midbox_y - solute_mean_y
            solute%site%r(3) = solute%site%r(3) + coo_midbox_z - solute_mean_z
            ! solute%site%r(1) = solute%site%r(1) + grid%length(1)/2.0_dp
            ! solute%site%r(2) = solute%site%r(2) + grid%length(2)/2.0_dp
            ! solute%site%r(3) = solute%site%r(3) + grid%length(3)/2.0_dp
        end if
        ! check if some positions are out of the supercell
        !j is a test tag. We loop over this test until every atom is in the box.
        ! This allows for instance, if a site is two boxes too far to still be ok.
        do concurrent( i=1:solute%nsite, d=1:3 )
            solute%site(i)%r(d) = MODULO ( solute%site(i)%r(d) , grid%length(d) )
        end do
    end subroutine mv_solute_to_center


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
        real(dp), parameter :: epsdp = epsilon(1._dp)

        soluteChargeDensity = 0._dp

        ! extrapolate each solute point charge to grid nodes
        do concurrent (s=1:size(solute%site), abs(solute%site(s)%q) > epsdp)
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
        end do

        volumElem = PRODUCT(gridlen/REAL(gridnode,dp))
        soluteChargeDensity = soluteChargeDensity / volumElem ! charge density is in charge per unit volume


    END SUBROUTINE soluteChargeDensityFromSoluteChargeCoordinates





    !> Print an XSF file of the supercell for visualisation in VMD for instance.
    !! Type vmd --xsf output/solute.xsf to visualise it.
    SUBROUTINE print_solute_xsf
        use precision_kinds, ONLY: i2b
        use module_grid, only: grid
        IMPLICIT NONE
        integer(i2b) :: i

        open(5,file='output/solute.xsf')
        write(5,*)"# this is the specification file of the supercell"
        write(5,*)"# lines beginning with # are commented. There cannot be comment lines within the sections"
        write(5,*)"# XSF format specifications can be found on the XCrySDen website http://www.xcrysden.org/doc/XSF.html"
        write(5,*)"# I strongly recommends to read this documentation."
        write(5,*)
        write(5,*)"# for periodic structures one has to begin with word CRYSTAL"
        write(5,*)"CRYSTAL"
        write(5,*)
        write(5,*)"# Then one needs to specify the lattice vectors"
        write(5,*)"# specification of PRIMVEC (in ANGSTROMS) like:"
        write(5,*)"#         ax, ay, az    (first lattice vector)"
        write(5,*)"#         bx, by, bz    (second lattice vector)"
        write(5,*)"#         cx, cy, cz    (third lattice vector)"
        write(5,*)"# pay attention to vectors as they are written in horizontal way which is quite unusual"
        write(5,*)"# for now only orthorhombic structures are allowed (free norms of lattice vectors, all angles are 90 degrees)"
        write(5,*)"PRIMVEC"
        write(5,*) grid%length(1), 0., 0.
        write(5,*) 0., grid%length(2), 0.
        write(5,*) 0., 0., grid%length(3)
        write(5,*)
        write(5,*)"# Then one needs to specify the atoms belonging to the unit cell. "
        write(5,*)"# First number stands for number of atoms in the primitive cell (2 in this case)."
        write(5,*)"# The second number is always 1 for PRIMCOORD coordinates."
        write(5,*)"# in angstroms and cartesian coordinates"
        write(5,*)"PRIMCOORD"
        write(5,*) SIZE(solute%site), 1
        do i = 1, SIZE(solute%site)
            write(5,*) solute%site(i)%Z, solute%site(i)%r
        end do
        close(5)
    END SUBROUTINE print_solute_xsf



    !This subroutine computes analitically and directly in k-space the charge
    !distribution of the solute molecule. It is quite  a strong repeating of
    !soluteChargeDensityFromSoluteChargeCoordinates for solvent but I did not
    !find a smart way yet to avoid this redundancy
    SUBROUTINE GetReciprocalSoluteChargeDensity()
        use module_grid, only: grid
        implicit none
        integer :: nx, ny, nz, no
        integer :: i, j, k, n, s, io, d
        real(dp)     :: r(3), kr, kvec(3)
        complex(dp)  :: fac, X
        complex(dp), parameter :: zeroc = (0._dp,0._dp), ic = (0._dp,1._dp)
        real(dp), parameter :: epsdp = epsilon(1._dp)
        real(dp) :: smootherfactor
        real(dp) :: smootherradius = 0.5_dp ! dramaticaly important
        nx = grid%nx
        ny = grid%ny
        nz = grid%nz
        no = grid%no

        ! sigma_k is the Fourier transformed charge density of a single solvent molecule in the reference frame defined by solvent.in
        ! molec_polar_k is the Fourier transformed molecular polarization
        if( sum( abs( solute%site%q )) > epsdp ) then
          allocate( solute%sigma_k (   nx/2+1, ny, nz), SOURCE=zeroC )
        else
          print*, "you want to compute the charge distribution of a neutral solute,there is a problem somewhere"
          stop
        end if
        
           !$omp parallel private(i, j, k, kvec, smootherfactor, r, kr, X, fac)
           !$omp do
           do k = 1, nz
              do j = 1, ny
                 do i = 1, nx/2+1
                    kvec = [ grid%kx(i), grid%ky(j), grid%kz(k) ]
                    !smootherfactor =  exp(-smootherradius**2 * sum( kvec**2 )/2._dp)
                    smootherfactor =  1.0_dp
                    
                    do n = 1, SIZE(solute%site)
                        if ( abs(solute%site(n)%q) > epsdp ) then
                           kr = dot_product( kvec, solute%site(n)%r )
                           X = -iC*kr
                           solute%sigma_k(i,j,k) = solute%sigma_k(i,j,k) + solute%site(n)%q *exp(X) *smootherfactor ! exact
                        end if
                    end do

                 end do
              end do
           end do
           !$omp end do
           !$omp end parallel
    END SUBROUTINE
!

end module module_solute
