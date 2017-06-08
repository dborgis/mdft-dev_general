module module_postprocessing
    use precision_kinds, only: dp
    implicit none
    private
    public :: init_postprocessing

contains
!
!     !
!     !   Post-processing of MDFT, following user's requirements
!     !
    subroutine init_postprocessing
        use precision_kinds, only: dp
        use module_solvent, only: solvent
        use module_solute, only: solute
        use module_grid, only: grid
        use module_orientation_projection_transform, only: angl2proj
        implicit none
        character(len=80) :: filename
        real(dp), allocatable :: density(:,:,:)
        integer :: nx, ny, nz, ix, iy, iz, is, isite, no
        real(dp), parameter :: pi=acos(-1._dp)

        nx=grid%nx
        ny=grid%ny
        nz=grid%nz
        no=grid%no


        !
        ! print density (in fact, rho/rho0)
        !
        allocate (density(nx,ny,nz))
        call grid%integrate_over_orientations( solvent(1)%xi**2 * solvent(1)%rho0, density)
        filename = "output/density.cube"
        call write_to_cube_file (density/solvent(1)%rho0/(4*pi**2), filename)
        print*, "New file output/density.cube"


        !
        ! print binary file one can use as a restart point
        !
        open(10,file='output/density.bin',form='unformatted')
        write(10) size(solvent)
        write(10) grid%mmax
        write(10) grid%no
        write(10) grid%np
        write(10) grid%nx, grid%ny, grid%nz
        write(10) grid%dx, grid%dy, grid%dz
        write(10) size(solute%site)
        do isite=1,size(solute%site)
          write(10) solute%site(isite)
        enddo
        block
            use module_grid, only: grid
            complex(dp) :: xi_p(grid%np)
            do is=1,size(solvent)
                do iz=1,nz
                    do iy=1,ny
                        do ix=1,nx
                            call angl2proj( solvent(is)%xi(:,ix,iy,iz), xi_p)
                            write(10) xi_p
                        end do
                    end do
                end do
            end do
        end block
        close(10)
        print*, "New file output/density.bin"



        !
        ! print polarization in each direction
        !
        block
        use module_input, only: getinput
        real(dp), allocatable, dimension(:,:,:,:) :: px, py, pz
        logical :: write_polarization_to_disk
        write_polarization_to_disk = getinput%log( "write_polarization_to_disk", defaultValue = .false. )
        if( write_polarization_to_disk ) then
            allocate(px(nx,ny,nz,1), py(nx,ny,nz,1), pz(nx,ny,nz,1), source=0._dp)
            call get_final_polarization(px,py,pz)
            filename = "output/Px.cube"
            call write_to_cube_file(px,filename)
            print*, "New file output/Px.cube"
            filename = "output/Py.cube"
            call write_to_cube_file(py,filename)
            print*, "New file output/Py.cube"
            filename = "output/Pz.cube"
            call write_to_cube_file(pz,filename)
            print*, "New file output/Pz.cube"
            if( size(solute%site) < 10 ) then ! plotting site site radial distribution functions (of the polarization here) for large molecules is not usefull
                filename = 'output/pnorm.xvg'
                call output_rdf ( sqrt(  px(:,:,:,1)**2 +py(:,:,:,1)**2 +pz(:,:,:,1)**2  ) , filename ) ! Get radial distribution functions
                print*, "New output file ", trim(adjustl(filename))
            end if
        end if
        end block



!         use system,             ONLY: thermocond
!         use module_solvent, only: solvent
!         use module_input,              ONLY: verbose, getinput
!         ! use solute_geometry,    ONLY: soluteIsPlanar => isPlanar, soluteIsLinear => isLinear
!         use constants,          ONLY: zerodp
!         use hardspheres,        only: hs
!         use module_grid, only: grid
!
!         IMPLICIT NONE
!
!         CHARACTER(50):: filename
!         REAL(dp), ALLOCATABLE , DIMENSION (:,:,:,:) :: neq, Px, Py, Pz ! equilibrium density, ie rho(r), and Pi polarization(r)
!         INTEGER(i2b) :: nfft1, nfft2, nfft3
!         INTEGER(i2b) :: s
!
!         nfft1 = grid%nx
!         nfft2 = grid%ny
!         nfft3 = grid%nz
!
!         CALL print_cg_vect_new ! print output/density.bin that contains cg_vect_new
!
!         allocate ( neq (nfft1,nfft2,nfft3,solvent(1)%nspec) ,SOURCE=zerodp)
!         do s=1,size(solvent)
!             call get_final_density (neq,s)
!         end do
!
!         DO s=1,solvent(1)%nspec
!           write(*,'(A,F12.2)') "Solvent molecules in supercell", SUM(neq)*grid%dv *solvent(s)%n0
!         END DO
!
!
!         IF (verbose) THEN
!             filename = 'output/density.cube'
!             CALL write_to_cube_file (neq(:,:,:,1), filename) ! TODO for now only write for the first species
!             IF ( getinput%char("polarization", defaultvalue="no") /= "no" ) THEN
!                 ALLOCATE ( Px (nfft1,nfft2,nfft3,solvent(1)%nspec) ,SOURCE=zerodp)
!                 ALLOCATE ( Py (nfft1,nfft2,nfft3,solvent(1)%nspec) ,SOURCE=zerodp)
!                 ALLOCATE ( Pz (nfft1,nfft2,nfft3,solvent(1)%nspec) ,SOURCE=zerodp)
!                 CALL get_final_polarization (Px,Py,Pz)
!                 filename='output/normP.cube' ; CALL write_to_cube_file ( (SQRT(Px(:,:,:,1)**2+Py(:,:,:,1)**2+Pz(:,:,:,1)**2)), filename) ! TODO for now only write for the first species
!                 filename='output/Px.cube' ; CALL write_to_cube_file(Px(:,:,:,1), filename)
!                 filename='output/z_Px.dat'; CALL compute_z_density(Px(:,:,:,1) , filename)
!                 DEALLOCATE (Px)
!                 filename='output/Py.cube' ; CALL write_to_cube_file(Py(:,:,:,1), filename)
!                 filename='output/z_Py.dat'; CALL compute_z_density(Py(:,:,:,1) , filename)
!                 DEALLOCATE (Py)
!                 filename='output/Pz.cube' ; CALL write_to_cube_file(Pz(:,:,:,1), filename)
!                 filename='output/z_Pz.dat'; CALL compute_z_density(Pz(:,:,:,1) , filename)
!                 DEALLOCATE (Pz)
!             END IF
!
!
!             ! If calculation is for hard sphere fluid in presence of a hard wall compute profile perp wall
!             ! TODO: DONT HAVE TIME TO WRITE THE TEST TODAY
!             filename = 'output/z_density.out'
!             CALL compute_z_density ( neq (:,:,:,1) , filename ) ! TODO for now only write for the first species
!
!             ! IF (soluteIsLinear()) THEN
!             !     ! nothing for now
!             ! END IF
!             !
!             ! IF( soluteIsPlanar() ) THEN
!             !     PRINT*,'This solute has planar symetry'
!             !     filename = 'output/planardensity.out'
!             !     CALL compute_planar_density ( neq (:,:,:,1) , filename ) ! TODO for now only write for the first species
!             ! END IF
!
!             IF ( getinput%char('other_predefined_vext')=='vextdef0' ) THEN
!                 filename = 'output/molecular_density_in_xy_plane.out'
!                 PRINT*,"I am writing file ",filename
!                 OPEN(378,FILE=filename)
!                     BLOCK
!                         INTEGER(i2b) :: i,j
!                         DO i=1,SIZE(neq,1)
!                             DO j=1,SIZE(neq,2)
!                                 WRITE(378,*)[i,j]*grid%length(1:2)/grid%n_nodes(1:2),neq(i,j,1,1)
!                             END DO
!                             WRITE(378,*)
!                         END DO
!                     END BLOCK
!                 CLOSE(378)
!             END IF
!         END IF
!
!
        block
            use module_solvent, only: solvent
            use module_input, only: getinput
            if( (solvent(1)%nsite < 10 .and. size(solute%site) < 10) .or. getinput%log ('write_rdf', defaultvalue=.false.) ) then ! For solutes and solvents with more than a few sites, site-site radial distribution functions are no longer meaningful.
                density = density / solvent(1)%n0
                filename = 'output/rdf.xvg'
                call output_rdf ( density , filename ) ! Get radial distribution functions
                print*, "New output file ", trim(adjustl(filename))
                call output_gsitesite
                call output_gOfRandCosThetaAndPsi
            end if
        end block

        deallocate (density)

        block
            use module_pressure_correction, only: pressure_correction
            call pressure_correction()
        end block

!         if( allocated(hs) ) then
!           block
!             real(dp)::x
!             x=hs(1)%pf
!             write(*,'(A,F7.2)') "packing fraction η =",x
!             write(*,'(A,F7.2)') "βP/n PY by pressure route        = (1+2η+3η²)           /(1-η)² =",(1+2*x+3*x**2)/(1-x)**2
!             write(*,'(A,F7.2)') "βP/n PY by compressibility route = (1+ η+ η²)           /(1-η)³ =",(1+x+x**2)/(1-x)**3
!             write(*,'(A,F7.2)') "βP/n CS                          = (1+ η+ η²-η³)        /(1-η)³ =",(1+x+x**2-x**3)/(1-x)**3
!             write(*,'(A,F7.2)') "βP/n CSK                         = (1+ η+ η²-2(1+η)η³/3)/(1-η)³ =",&
!               (1+x+x**2-(2./3.)*(1+x)*x**3)/((1.-x)**3)
!           end block
!         end if
!
    end subroutine init_postprocessing
!
!     SUBROUTINE print_cg_vect_new
!         use module_minimizer, ONLY: cg_vect_new
!         if ( .not. allocated ( cg_vect_new ) ) then
!             print *, 'cg_vect_new is not allocated in SUBROUTINE print_cg_vect_new in process_output.f90. STOP.'
!             stop
!         END IF
!         OPEN (10, file = 'output/density.bin.out' , form = 'unformatted' )
!             write ( 10 ) cg_vect_new
!         CLOSE (10)
!     END SUBROUTINE print_cg_vect_new
!
!
!

!
!
!
!     pure subroutine get_final_density ( neq , s)
!         use precision_kinds, only: dp
!         use module_solvent, only: solvent
!         use module_quadrature, only: mean_over_orientations
!         implicit none
!         integer, intent(in) :: s ! the solvent species of which we want the number density
!         real(dp), intent(out) :: neq(:,:,:)
!         call mean_over_orientations( solvent(solventspecies)%rho , neq)
!     end subroutine get_final_density
!
!
!
!
!
    subroutine get_final_polarization ( Px , Py , Pz )

        use precision_kinds, only: dp, i2b
        use module_solvent, only: solvent
        use module_grid, only: grid

        IMPLICIT NONE
        INTEGER(i2b) :: i, j, k, io, s
        REAL(dp) :: x, local_Px, local_Py, local_Pz
        REAL(dp), dimension(:,:,:,:), intent(out) :: Px, Py, Pz ! equilibrium polarization(r)
        real(dp), parameter :: zerodp = 0._dp


        Px = zerodp
        Py = zerodp
        Pz = zerodp

        DO s =1,solvent(1)%nspec
            DO i =1,grid%nx
                DO j =1,grid%ny
                    DO k =1,grid%nz
                        local_Px = 0.0_dp
                        local_Py = 0.0_dp
                        local_Pz = 0.0_dp
                        DO io =1,grid%no
                            x = solvent(s)%xi(io,i,j,k)**2*solvent(s)%rho0
                            local_Px = local_Px + grid%omx(io) * grid%w(io) * x
                            local_Py = local_Py + grid%omy(io) * grid%w(io) * x
                            local_Pz = local_Pz + grid%omz(io) * grid%w(io) * x
                        END DO
                        Px(i,j,k,s) = local_Px
                        Py(i,j,k,s) = local_Py
                        Pz(i,j,k,s) = local_Pz
                    END DO
                END DO
            END DO
        END DO

    end subroutine get_final_polarization

end module module_postprocessing
