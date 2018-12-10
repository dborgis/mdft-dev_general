module module_postprocessing
    use precision_kinds, only: dp
    use module_cubefiles 
    implicit none
    private
    public :: init_postprocessing, get_Solvent_Molecule_Pseudo_Charge_Density

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
real(dp), allocatable :: density(:,:,:), charge_density(:,:,:), pseudo_charge_density(:,:,:)
        integer :: nx, ny, nz, ix, iy, iz, is, isite, no,io
        real(dp), parameter :: pi=acos(-1._dp)
        logical:: output_full_density
        character(180) :: solvent_pseudo_charge_density

        nx=grid%nx
        ny=grid%ny
        nz=grid%nz
        no=grid%no

        !
        ! print density (in fact, rho/rho0)

        allocate ( density(nx, ny, nz) , charge_density(nx, ny, nz ), pseudo_charge_density(nx, ny, nz) )


WRITE_DENSITY: BLOCK
        call grid%integrate_over_orientations( solvent(1)%xi**2 * solvent(1)%rho0, density)

        filename = "output/density.cube"
        call write_to_cube_file (density/solvent(1)%rho0/(4*pi**2), filename)
        filename = 'output/z_density.out'
        CALL compute_z_density ( density , filename ) ! TODO for now only write for the first species
        print*, "New file output/density.cube. Try$ vmd -cube output/density.cube"
        filename = 'output/z_density.out'
        CALL compute_z_density ( density(:,:,:) , filename ) ! TODO for now only write for the first species


        output_full_density=getinput%log('write_full_density', defaultvalue=.false.)
        ! print binary file one can use as a restart point
        !
        open(10,file='output/density.bin',form='unformatted')
        if (output_full_density) then
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
END BLOCK WRITE_DENSITY


WRITE_SOLVENT_PSEUDO_CHARGE_DENSITY: BLOCK
use module_input, only: getinput
real(dp) :: sum_charges

logical :: compute_solvent_pseudo_charge_density

compute_solvent_pseudo_charge_density = getinput%log( "compute_solvent_pseudo_charge_density", defaultValue = .false. )

if( compute_solvent_pseudo_charge_density ) then

call get_final_Solvent_Pseudo_Charge_Density ( charge_density )

filename = "output/solvent_pseudo_charge_density.cube"
call write_to_cube_file( charge_density, filename )
print*, "New file created:", trim(adjustl(filename))

!  Compute and print corresponding RDFs
if( (solvent(1)%nsite < 10 .and. size(solute%site) < 10) ) then
filename = 'output/charge_density_rdf.xvg'
charge_density = charge_density/solvent(1)%n0
call output_rdf ( charge_density , filename ) ! Get radial distribution functions
print*, "New file created:", trim(adjustl(filename))
end if

!
call get_Solvent_Molecule_Pseudo_Charge_Density( charge_density , 1 )
filename = "output/molecule_pseudo_charge_density.cube"
call write_to_cube_file( charge_density, filename )
print*, "New file created:", trim(adjustl(filename))
if( (solvent(1)%nsite < 10 .and. size(solute%site) < 10) ) then

filename = 'output/molecule_pseudo_charge_density_rdf.xvg'
call output_rdf ( charge_density , filename ) ! Get radial distribution functions
print*, "New file created:", trim(adjustl(filename))
end if

end if

END BLOCK WRITE_SOLVENT_PSEUDO_CHARGE_DENSITY

        !
        ! print polarization in each direction
        !
WRITE_POLARIZATION: BLOCK
        use module_input, only: getinput
        real(dp), allocatable, dimension(:,:,:,:) :: px, py, pz ! last dimension accoutns for the solvent id.
        logical :: write_polarization_to_disk
        write_polarization_to_disk = getinput%log( "write_polarization_to_disk", defaultValue = .false. )
        if( write_polarization_to_disk ) then
            allocate(px(nx,ny,nz,size(solvent)), py(nx,ny,nz,size(solvent)), pz(nx,ny,nz,size(solvent)), source=0._dp)
            call get_final_polarization(px,py,pz)
            filename = "output/Px.cube"
            call write_to_cube_file(px,filename)
            print*, "New file output/Px.cube. Try$ vmd -cube output/Px.cube"
            filename = "output/Py.cube"
            call write_to_cube_file(py,filename)
            print*, "New file output/Py.cube. Try$ vmd -cube output/Py.cube"
            filename = "output/Pz.cube"
            call write_to_cube_file(pz,filename)
            print*, "New file output/Pz.cube. Try$ vmd -cube output/Pz.cube"
            filename='output/z_Pz.dat'; CALL compute_z_density(Pz(:,:,:,1) , filename)
            filename = "output/Pnorm.cube"
            call write_to_cube_file( sqrt( px(:,:,:,1)**2 +py(:,:,:,1)**2 +pz(:,:,:,1)**2  ), filename ) 
            print*, "New file output/Pnorm.cube. Try$ vmd -cube output/Pnorm.cube"
            if( size(solute%site) < 50 ) then ! plotting site site radial distribution functions (of the polarization here) for large molecules is not usefull
                filename = 'output/pnorm.xvg'
                call output_rdf ( sqrt(  px(:,:,:,1)**2 +py(:,:,:,1)**2 +pz(:,:,:,1)**2  ) , filename ) ! Get radial distribution functions
                print*, "New output file ", trim(adjustl(filename)), ". Try$ xmgrace output/pnorm.xvg"
            end if
        end if
END BLOCK WRITE_POLARIZATION
!
!
WRITE_DENSITY_RDFs: BLOCK
            use module_solvent, only: solvent
            use module_input, only: getinput
            if( (solvent(1)%nsite < 10 .and. size(solute%site) < 10) .or. getinput%log ('write_rdf', defaultvalue=.false.) ) then ! For solutes and solvents with more than a few sites, site-site radial distribution functions are no longer meaningful.
                density = density / solvent(1)%n0
                filename = 'output/rdf.xvg'
                call output_rdf ( density , filename ) ! Get radial distribution functions
                print*, "New file ", trim(adjustl(filename))
                if( getinput%log("write_angular_rdf", defaultValue=.false.)) then
                    call output_gsitesite ! may be very time-consuming for large supercells / solutes
                    call output_gOfRandCosThetaAndPsi ! may also be very time-consuming
                end if

            end if
END BLOCK WRITE_DENSITY_RDFs

!        deallocate (density)

PRESSURE_CORRECTIONS: BLOCK
            use module_pressure_correction, only: pressure_correction
            call pressure_correction()
END BLOCK PRESSURE_CORRECTIONS

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
!
subroutine get_final_Solvent_Pseudo_Charge_Density( charge_density )

use iso_c_binding
use precision_kinds
use module_solvent, only: solvent
use module_grid, only: grid

IMPLICIT NONE
INTEGER(i2b) :: i, j, k, io, s, ix, iy, iz
integer(i2b) :: nx, ny, nz, no, ns

REAL(dp), allocatable, intent(out) :: charge_density(:,:,:) ! solvent pseudo-density
complex(dp), allocatable :: charge_density_k(:,:,:)  !  same in k-space

real(dp), allocatable :: fftw3inforward(:,:,:), fftw3outbackward(:,:,:)
complex(dp), allocatable :: fftw3outforward(:,:,:)
integer(i4b) :: plan_forward, plan_backward ! fftw3 plan forward or backward
complex(dp), parameter :: zeroc=(0._dp,0._dp)
real(dp), parameter :: zero = 0._dp
real(dp), parameter :: twopi  =2._dp*acos(-1._dp)
real(dp), parameter :: fourpi =2._dp*twopi
real(dp) :: sum_charges
complex(dp) :: sum_charges_c

include "fftw3.f03"

nx = grid%nx
ny = grid%ny
nz = grid%nz
no = grid%no
ns = solvent(1)%nspec

print*, 'passed in get_final_Solvent_Pseudo_Charge_Density'

if( ns /= 1 ) stop 'subroutine get_final_Solvent_Pseudo_Density only for a single solvent component'

! allocate the arrays needed as input for fft (in_forward) or output for fft (out_forward)
! or needed as input for inverse fft (in_backward) etc.
allocate ( fftw3inforward   (nx    , ny,nz), fftw3outbackward(nx, ny, nz) )
allocate ( fftw3outforward  (nx/2+1, ny,nz))

allocate( charge_density (nx, ny, nz) ,source=zero )
allocate( charge_density_k (nx/2+1, ny, nz) ,source=zeroc)


! prepare plans needed by fftw3
select case(dp)
case(c_double)
call dfftw_plan_dft_r2c_3d (plan_forward, nx, ny, nz, fftw3inforward, fftw3outforward, fftw_estimate)
call dfftw_plan_dft_c2r_3d (plan_backward, nx, ny, nz, charge_density_k, fftw3outbackward, fftw_estimate )
case(c_float)
call sfftw_plan_dft_r2c_3d (plan_forward, nx, ny, nz, fftw3inforward, fftw3outforward, fftw_estimate)
call sfftw_plan_dft_c2r_3d (plan_backward, nx, ny, nz, charge_density_k, fftw3outbackward, fftw_estimate )
end select

 do io= 1, no

   fftw3inforward(:,:,:) = solvent(1)%xi(io,:,:,:)**2*solvent(1)%rho0

! forward FFT
     select case(dp)
     case(c_double)
     call dfftw_execute(plan_forward)
     case(c_float)
     call sfftw_execute(plan_forward)
     end select

    charge_density_k(:,:,:) =  charge_density_k(:,:,:) +  &
          fftw3OutForward(:,:,:)*grid%w(io)*solvent(1)%pseudo_charge_density_k(:,:,:,io)

  end do !end loop over angles

! FFT backward
select case(dp)
case(c_double)
call dfftw_execute(plan_backward)
case(c_float)
call sfftw_execute(plan_backward)
end select

charge_density = fftw3outbackward/real(nx*ny*nz,dp)

select case(dp)
case(c_double)
call dfftw_destroy_plan (plan_forward)
call dfftw_destroy_plan (plan_backward)
case(c_float)
call sfftw_destroy_plan (plan_forward)
call sfftw_destroy_plan (plan_backward)
end select

deallocate( charge_density_k, fftw3InForward, fftw3OutForward, fftw3outbackward )

end subroutine get_final_Solvent_Pseudo_Charge_Density


subroutine get_Solvent_Molecule_Pseudo_Charge_Density( charge_density , io)

use iso_c_binding
use precision_kinds
use module_solvent, only: solvent
use module_grid, only: grid

IMPLICIT NONE
INTEGER(i2b) :: i, j, k, s
integer(i2b), intent(in) :: io

REAL(dp), allocatable, intent(out) :: charge_density(:,:,:) ! solvent pseudo-density
complex(dp), allocatable :: charge_density_k(:,:,:)  !  same in k-space

integer(i4b) :: plan_forward, plan_backward ! fftw3 plan forward or backward
complex(dp), parameter :: zeroc=(0._dp,0._dp), ic=(0._dp,1._dp)
real(dp), parameter :: zero = 0._dp
real(dp), parameter :: twopi  =2._dp*acos(-1._dp)
real(dp), parameter :: fourpi =2._dp*twopi
real(dp), parameter :: espdp=epsilon(1._dp)
integer :: nx, ny, nz, no, ns
real(dp) :: sum_charges, r(3), kr, kvec(3)

include "fftw3.f03"

nx = grid%nx
ny = grid%ny
nz = grid%nz
no = grid%no
ns = solvent(1)%nspec

!print*, 'passed in get_Solvent_Molecule_Pseudo_Charge_Density'

if( ns /= 1 ) stop 'subroutine get_Solvent_Molecule_Pseudo_Density only for a single solvent component'
if(io > no) stop 'in get_Solvent_Molecule_Pseudo_Charge_Density: io cannot be smaller than grid%no '

! allocate the arrays needed as input for fft (in_forward) or output for fft (out_forward)
! or needed as input for inverse fft (in_backward) etc.

allocate( charge_density (nx, ny, nz) ,source=zero )
allocate( charge_density_k (nx/2+1, ny, nz) ,source=zeroc)

! prepare plans needed by fftw3
select case(dp)
case(c_double)
call dfftw_plan_dft_c2r_3d (plan_backward, nx, ny, nz, charge_density_k, charge_density, fftw_estimate )
case(c_float)
call sfftw_plan_dft_c2r_3d (plan_backward, nx, ny, nz, charge_density_k, charge_density, fftw_estimate )
end select

charge_density_k(:,:,:) =  solvent(1)%pseudo_charge_density_k(:,:,:,io)

! shift to center of the box
do k = 1, nz
  do j = 1, ny
    do i = 1, nx/2+1

     kvec = [ grid%kx(i), grid%ky(j), grid%kz(k) ]

     r(1) = - grid%length(1)/2._dp
     r(2) = - grid%length(2)/2._dp
     r(3) = - grid%length(3)/2._dp
     kr =  dot_product( kvec, r )

     charge_density_k(i, j, k) = charge_density_k(i, j, k)  * exp(-ic*kr)

    end do !loop nx
  end do ! loop ny
end do !loop nz


! FFT backward
select case(dp)
case(c_double)
call dfftw_execute(plan_backward)
case(c_float)
call sfftw_execute(plan_backward)
end select

charge_density = charge_density*grid%dv

!check molecule charge density
   sum_charges = zero
   do k = 1, nz
     do j = 1, ny
       do i = 1, nx
          sum_charges = sum_charges + charge_density(i, j, k)
       end do
     end do
   end do
!'print*, '***************'
print*, 'check solvent charge : io= ', io,'total charge= ',sum_charges
!'print*, '***************'

select case(dp)
case(c_double)
call dfftw_destroy_plan (plan_backward)
case(c_float)
call sfftw_destroy_plan (plan_backward)
end select

deallocate( charge_density_k )

end subroutine get_Solvent_Molecule_Pseudo_Charge_Density
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



!  COMMENTED PIECES OF CODES
!
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

