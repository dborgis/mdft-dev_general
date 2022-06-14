module module_postprocessing
    use precision_kinds, only: dp
    use module_cubefiles 
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
        use module_input
        implicit none
        character(len=80) :: filename
        real(dp), allocatable :: density(:,:,:), charge_density
        integer :: nx, ny, nz, ix, iy, iz, is, isite, no,io
        real(dp), parameter :: pi=acos(-1._dp)
        logical:: output_full_density,write_density,write_angular_density
        nx=grid%nx
        ny=grid%ny
        nz=grid%nz
        no=grid%no

        
        !
        ! print density (in fact, rho/rho0)
        !


WRITE_ANGULAR_DENSITIES: BLOCK
        ! print binary file one can use as a restart point
        write_angular_density=getinput%log('write_angular_density',defaultvalue=.false.)
        IF(write_angular_density) then
        open(10,file='output/density.bin',form='unformatted')
        output_full_density=getinput%log('write_full_density',defaultvalue=.false.)
        if(output_full_density) then
           write(10) -size(solvent)
        else
           write(10) size(solvent)
        end if
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
            if (.not. output_full_density) then
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
            else
              do is=1,size(solvent)
                  do iz=1,nz
                      do iy=1,ny
                          do ix=1,nx
                            do io=1,no
                              write(10) solvent(is)%xi(io,ix,iy,iz)
                            end do
                          end do
                      end do
                  end do
            end do
          end if
        end block
        close(10)
        print*, "New file output/density.bin"
        END IF
END BLOCK WRITE_ANGULAR_DENSITIES


allocate (density(nx,ny,nz))
call grid%integrate_over_orientations( solvent(1)%xi**2 * solvent(1)%rho0, density)


WRITE_NUMBER_DENSITIES: BLOCK
        write_density=getinput%log('write_density', defaultvalue=.false.)
        IF(write_density) then
        filename = "output/density.cube"
        call write_to_cube_file (density/solvent(1)%n0, filename)
        print*, "New file output/density.cube. Try$ vmd -cube output/density.cube"
        END IF
END BLOCK WRITE_NUMBER_DENSITIES


WRITE_NUMBER_DENSITY_RDFs: BLOCK
if( (solvent(1)%nsite < 50 .and. size(solute%site) < 50) .or. getinput%log ('write_rdf', defaultvalue=.false.) ) then ! For solutes and solvents with more than a few sites, site-site radial distribution functions are no longer meaningful.
density = density / solvent(1)%n0
filename = 'output/rdf'
call output_rdf ( density , filename ) ! Get radial distribution functions
print*, "New file ", trim(adjustl(filename))
end if

if( (solvent(1)%nsite > 5 .or. size(solute%site) > 10) .and. getinput%log ('write_site-site_rdf', defaultvalue=.false.) ) then
   print*, 'site-site rdfs only computed if (solvent(1)%nsite < 5 .and. size(solute%site) < 10)'
end if
if( (solvent(1)%nsite < 5 .and. size(solute%site) < 10) .and. getinput%log ('write_site-site_rdf', defaultvalue=.false.) ) then
   call output_gsitesite ! may be very time-consuming for large supercells / solutes
   print*, "New file:  output/g-sitesite.out"
end if

!if( getinput%log("write_angular_rdf", defaultValue=.false.)) then
!call output_gOfRandCosThetaAndPsi ! may also be very time-consuming
!end if

END BLOCK WRITE_NUMBER_DENSITY_RDFs
deallocate (density)

        ! print polarization in each direction
        !

WRITE_POLARIZATION: BLOCK
        real(dp), allocatable, dimension(:,:,:,:) :: px, py, pz ! last dimension accoutns for the solvent id.
        logical :: write_polarization_to_disk
        write_polarization_to_disk = getinput%log( "write_polarization", defaultValue = .false. )
        if( write_polarization_to_disk ) then
            allocate(px(nx,ny,nz,size(solvent)), py(nx,ny,nz,size(solvent)), pz(nx,ny,nz,size(solvent)), source=0._dp)
            call get_final_polarization(px,py,pz)
            filename = "output/polarization/Px.cube"
            call write_to_cube_file(px,filename)
            print*, "New file output/Px.cube. Try$ vmd -cube output/Px.cube"
            filename = "output/polarization/Py.cube"
            call write_to_cube_file(py,filename)
            print*, "New file output/polarization/Py.cube. Try$ vmd -cube output/Py.cube"
            filename = "output/Pz.cube"
            call write_to_cube_file(pz,filename)
            print*, "New file output/polarization/Pz.cube. Try$ vmd -cube output/Pz.cube"
            filename='output/polarization/z_Pz.dat'; CALL compute_z_density(Pz(:,:,:,1) , filename)
            filename = "output/polarization/Pnorm.cube"
            call write_to_cube_file( sqrt( px(:,:,:,1)**2 +py(:,:,:,1)**2 +pz(:,:,:,1)**2  ), filename ) 
            print*, "New file output/polarization/Pnorm.cube. Try$ vmd -cube output/Pnorm.cube"
            if( size(solute%site) < 50 ) then ! plotting site site radial distribution functions (of the polarization here) for large molecules is not usefull
                filename = 'output/polarization/pnorm'
                call output_rdf ( sqrt(  px(:,:,:,1)**2 +py(:,:,:,1)**2 +pz(:,:,:,1)**2  ) , filename ) ! Get radial distribution functions
                print*, "New output file ", trim(adjustl(filename)), ". Try$ xmgrace output/pnorm.xvg"
            end if
        deallocate( px, py, pz)
        end if
END BLOCK WRITE_POLARIZATION


PRESSURE_CORRECTIONS: BLOCK
use module_pressure_correction, only: pressure_correction
call pressure_correction()
END BLOCK PRESSURE_CORRECTIONS

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
