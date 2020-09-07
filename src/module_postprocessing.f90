module module_postprocessing
    use precision_kinds, only: dp
    use module_cubefiles 
    implicit none
    private
    public :: init_postprocessing ! , get_Solvent_Molecule_Pseudo_Charge_Density

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
        real(dp), allocatable :: density(:,:,:), charge_density(:,:,:), pseudo_charge_density(:,:,:), electron_density(:,:,:)
        real(dp), allocatable :: solvent_electrostatic_potential(:,:,:)
        complex(dp), allocatable :: molecule_charge_density_k(:,:,:,:)
        integer :: nx, ny, nz, ix, iy, iz, is, isite, no,io
        real(dp), parameter :: pi=acos(-1._dp)
        complex(dp), parameter :: zeroc=(0._dp,0._dp)
        real(dp), parameter :: angtobohr = 1.889725989_dp
        real(dp) :: charge_smoothing_radius
        logical:: output_full_density,write_density,write_angular_density
  !      character(180) :: solvent_pseudo_charge_density
  !      character(80) :: charge_densities_directory_path

        nx=grid%nx
        ny=grid%ny
        nz=grid%nz
        no=grid%no
        charge_smoothing_radius = grid%dl(1)   ! TAKE CARE for non-cubic boxes

        !
        ! print density (in fact, rho/rho0)

        allocate ( density(nx, ny, nz) , charge_density(nx, ny, nz ), pseudo_charge_density(nx, ny, nz), electron_density(nx, ny, nz ) )


WRITE_NUMBER_DENSITIES: BLOCK
 !   use module_input
  !  logical:: output_full_density,output_density,write_angular_density

    write_density=getinput%log('write_density', defaultvalue=.false.)
    if(write_density) then
        call grid%integrate_over_orientations( solvent(1)%xi**2 * solvent(1)%rho0, density)

        filename = "output/density.cube"
        call write_to_cube_file (density/solvent(1)%rho0/(4*pi**2), filename)
        filename = 'output/z_density.out'
        CALL compute_z_density ( density , filename ) ! TODO for now only write for the first species
        print*, "New file output/density.cube. Try$ vmd -cube output/density.cube"
        filename = 'output/z_density.out'
        CALL compute_z_density ( density(:,:,:) , filename ) ! TODO for now only write for the first species
    end if

    write_angular_density=getinput%log('write_angular_density', defaultvalue=.false.)

! print binary file one can use as a restart point
    if(write_angular_density) then
        open(10,file='output/density.bin',form='unformatted')
        output_full_density=getinput%log('write_full_density', defaultvalue=.false.)
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
    end if
END BLOCK WRITE_NUMBER_DENSITIES

allocate ( molecule_charge_density_k(nx/2 +1, ny, nz, no), source = zeroC )

WRITE_SOLVENT_PSEUDO_CHARGE_DENSITY: BLOCK
use module_input, only: getinput
real(dp) :: sum_charges

logical :: write_solvent_pseudo_charge_density

write_solvent_pseudo_charge_density = getinput%log( "write_solvent_pseudo_charge_density", defaultValue = .false. )

    if( write_solvent_pseudo_charge_density ) then

        call Get_molecule_reciprocal_pseudo_charge_density( molecule_charge_density_k, charge_smoothing_radius )
        call get_final_Solvent_Pseudo_Charge_Density ( molecule_charge_density_k, pseudo_charge_density )

        filename = "output/charge_density/solvent_pseudo_charge_density.cube"
        ! write charge density as punctual charges
        call write_to_cube_file( pseudo_charge_density/solvent(1)%n0, filename )
        print*, "New file created:", trim(adjustl(filename))

        !WRITE FOR GAUSSIAN FILES
        !filename = "output/charge_density/solvent_pseudo_charge_density_GAUSSIAN.cube"
        !call write_to_Gaussian_cube_file( pseudo_charge_density*grid%dv, filename )
        !print*, "New file created:", trim(adjustl(filename)),  '   for visualization  !!!!'

        !filename = "output/charge_density/solvent_pseudo_charge_density_for_GAUSSIAN.out"
        !call write_to_Gaussian_charge_file( pseudo_charge_density*grid%dv, filename )
        !print*, "New file created:", trim(adjustl(filename)),  '    as input to Gaussian  !!!!'

        !  Compute and print corresponding RDFs
        if( (solvent(1)%nsite < 50 .and. size(solute%site) < 50) ) then
            filename = 'output/charge_density/solvent_pseudo_charge_density_rdf'
            !pseudo_charge_density = pseudo_charge_density/solvent(1)%n0
            !pseudo_charge_density = pseudo_charge_density*grid%dv
            call output_rdf ( pseudo_charge_density/angtobohr**3 , filename ) ! Get radial distribution functions
            print*, "New file created:", trim(adjustl(filename))
        end if

    end if
END BLOCK WRITE_SOLVENT_PSEUDO_CHARGE_DENSITY


WRITE_SOLVENT_CHARGE_DENSITY: BLOCK
use module_input, only: getinput
real(dp) :: sum_charges

logical :: write_solvent_charge_density
character(80) :: solvent_electron_density_type

write_solvent_charge_density = getinput%log( "write_solvent_charge_density", defaultValue = .false. )
solvent_electron_density_type = getinput%char('solvent_electron_density_type', defaultvalue='point_charges')

    if( write_solvent_charge_density ) then

    call Get_SPC_water_molecule_reciprocal_charge_density( molecule_charge_density_k, solvent_electron_density_type )
    call get_final_Solvent_Pseudo_Charge_Density ( molecule_charge_density_k, charge_density )

    filename = "output/charge_density/solvent_charge_density.cube"
    ! write charge density as punctual charges
    call write_to_cube_file( charge_density/angtobohr**3, filename )
    print*, "New file created:", trim(adjustl(filename))
    print*,'minval, maxval of charge_density = ',minval(charge_density/angtobohr**3), maxval(charge_density/angtobohr**3)

    !  Compute and print corresponding RDFs
    if( (solvent(1)%nsite < 50 .and. size(solute%site) < 50) ) then
        filename = 'output/charge_density/solvent_charge_density_rdf'
        !charge_density = charge_density/solvent(1)%n0
        !pseudo_charge_density = pseudo_charge_density*grid%dv
        call output_rdf ( charge_density/solvent(1)%n0 , filename ) ! Get radial distribution functions
        print*, "New file created:", trim(adjustl(filename))
    end if

    end if
END BLOCK WRITE_SOLVENT_CHARGE_DENSITY


WRITE_SOLVENT_ELECTRON_DENSITY: BLOCK
use module_input, only: getinput
real(dp) :: sum_charges

logical :: write_solvent_electron_density
character(80) :: solvent_electron_density_type

write_solvent_electron_density = getinput%log( "write_solvent_electron_density", defaultValue = .false. )
solvent_electron_density_type = getinput%char('solvent_electron_density_type', defaultvalue='point_charges')

if( write_solvent_electron_density ) then

    call Get_water_molecule_reciprocal_electron_density( molecule_charge_density_k, solvent_electron_density_type )
    call get_final_Solvent_Pseudo_Charge_Density ( molecule_charge_density_k, electron_density )

    filename = "output/charge_density/solvent_electron_density.cube"
    ! write charge density as punctual charges
    call write_to_cube_file( electron_density/angtobohr**3, filename )
    print*, "New file created:", trim(adjustl(filename))
    print*,'minval, maxval of solvent_electron_density = ',minval(electron_density/angtobohr**3), maxval(electron_density/angtobohr**3)

    !  Compute and print corresponding RDFs
    if( (solvent(1)%nsite < 50 .and. size(solute%site) < 50) ) then
        filename = 'output/charge_density/solvent_electron_density_rdf'
        !electron_density = electron_density/solvent(1)%n0
        !pseudo_charge_density = pseudo_charge_density*grid%dv
        call output_rdf ( electron_density/solvent(1)%n0 , filename ) ! Get radial distribution functions
        print*, "New file created:", trim(adjustl(filename))
    end if

end if
END BLOCK WRITE_SOLVENT_ELECTRON_DENSITY

deallocate ( molecule_charge_density_k )

WRITE_SOLVENT_ELECTROSTATIC_POTENTIAL: BLOCK
use module_input, only: getinput

logical :: write_solvent_electrostatic_potential

write_solvent_electrostatic_potential = getinput%log( "write_solvent_electrostatic_potential", defaultValue = .false. )

if( write_solvent_electrostatic_potential ) then

allocate( solvent_electrostatic_potential(nx, ny, nz) )

call Get_solvent_electrostatic_potential( charge_density, solvent_electrostatic_potential )


filename = "output/charge_density/solvent_electrostatic_potential.cube"
! write charge density as punctual charges
call write_to_cube_file( solvent_electrostatic_potential/angtobohr, filename )
print*, "New file created:", trim(adjustl(filename))
print*,'minval, maxval of solvent_electrostatic_potential = ',minval(solvent_electrostatic_potential/angtobohr), maxval(solvent_electrostatic_potential/angtobohr)

!  Compute and print corresponding RDFs
if( (solvent(1)%nsite < 50 .and. size(solute%site) < 50) ) then
filename = 'output/charge_density/solvent_electrostatic_potential_rdf'
call output_rdf ( solvent_electrostatic_potential, filename ) ! Get radial distribution functions
print*, "New file created:", trim(adjustl(filename))
end if

end if
END BLOCK WRITE_SOLVENT_ELECTROSTATIC_POTENTIAL

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
        end if
END BLOCK WRITE_POLARIZATION
!
!
WRITE_DENSITY_RDFs: BLOCK
            use module_solvent, only: solvent
            use module_input, only: getinput
            if( (solvent(1)%nsite < 11 .and. size(solute%site) < 11) .or. getinput%log ('write_rdf', defaultvalue=.false.) ) then ! For solutes and solvents with more than a few sites, site-site radial distribution functions are no longer meaningful.
                density = density / solvent(1)%n0
                filename = 'output/rdf'
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
subroutine get_final_Solvent_Pseudo_Charge_Density(molecule_charge_density_k, charge_density )

use iso_c_binding
use precision_kinds
use module_solvent, only: solvent
use module_grid, only: grid

IMPLICIT NONE
INTEGER(i2b) :: i, j, k, io, s, ix, iy, iz
integer(i2b) :: nx, ny, nz, no, ns

REAL(dp), allocatable, intent(out) :: charge_density(:,:,:) ! solvent pseudo-density
complex(dp),  allocatable  :: molecule_charge_density_k(:,:,:,:)
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

!print*, 'passed in get_final_Solvent_Pseudo_Charge_Density'

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
          fftw3OutForward(:,:,:)*grid%w(io)*molecule_charge_density_k(:,:,:,io)

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


subroutine get_Centered_Solvent_Molecule_Pseudo_Charge_Density( charge_density , io)

! pseudo-charge density of a single molecule at the center of the box
! for illustration only: takes solvent(1)%pseudo_charge_density_k(:,:,:,io) and
! shift it to (Lx/2,Ly/2,Lz/2)

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

charge_density = charge_density/grid%dv/real(nx*ny*nz,dp)

!check molecule charge density
   sum_charges = zero
   do k = 1, nz
     do j = 1, ny
       do i = 1, nx
          sum_charges = sum_charges + charge_density(i, j, k)*grid%dv
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

end subroutine get_Centered_Solvent_Molecule_Pseudo_Charge_Density
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


subroutine Get_water_molecule_reciprocal_electron_density( water_molecule_electron_density_k , electron_density_type )
use iso_c_binding
use precision_kinds
use module_solvent, only: solvent
use module_grid, only: grid

implicit none
integer :: nx, ny, nz, no, ns
integer :: i, j, k, n, s, io, d, m
real(dp)     :: r(3), kr, kvec(3), q2, f_o, f_h
complex(dp), allocatable, intent(out) :: water_molecule_electron_density_k(:,:,:,:)
complex(dp)  :: fac, X
real(dp) :: GaussianFactor
complex(dp), parameter :: zeroc = (0._dp,0._dp), ic = (0._dp,1._dp)
real(dp), parameter :: fourpi = 4._dp*acos(-1._dp)
real(dp), parameter :: q_h = 0.4238_dp, q_o = -0.8476_dp !SPC/E
!real(dp), parameter :: q_h = 0.41_dp, q_o = -0.82_dp !SPC
real(dp), parameter, dimension(4) :: a_o = (/3.0485, 2.2868, 1.5463, 0.867 /) , b_o =(/ 13.2771, 5.7011, 0.3229, 32.9089 /)
real(dp), parameter, dimension(4)  :: a_h = (/ 0.48918, 0.262003, 0.196767, 0.049879 /), b_h(4) = (/ 20.6593, 7.74029, 49.5519, 2.20159 /)
real(dp), parameter :: c_o = 0.2508,  c_h = 0.001305
character(80) :: electron_density_type

nx = grid%nx
ny = grid%ny
nz = grid%nz
no = grid%no
ns = size(solvent) ! Count of solvent species

allocate ( water_molecule_electron_density_k(nx/2+1, ny, nz, no), source = zeroC )

!print*, 'passed in subroutine get_solvent_molecule_pseudo_charge_density, created by Daniel on 7-12-2018 for introducing QM/MM electron-water pseudopotential'

! At this stage: pseudo_charge_density is the Fourier transformed charge density of a single water molecule in the reference frame defined by solvent.in

if( ns /= 1) stop 'init_solvent_molecule_pseudo_charge_density only works for ns = 1'
if( solvent(1)%name /= 'spce' ) stop 'init_solvent_molecule_pseudo_charge_density only work for spc or spce'


if(electron_density_type == 'point_charges') then
!$omp parallel private(i, j, k, kvec, smootherfactor, r, kr, X, fac)
!$omp do

    do io = 1, no

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx/2+1
    kvec = [ grid%kx(i), grid%ky(j), grid%kz(k) ]
    q2 = (grid%kx(i)**2 + grid%ky(j)**2 + grid%kz(k)**2)/fourpi**2

    n= 1 ! oxygen site


    r(1) = dot_product(   [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)]  ,  solvent(1)%site(n)%r  )
    r(2) = dot_product(   [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)]  ,  solvent(1)%site(n)%r  )
    r(3) = dot_product(   [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)]  ,  solvent(1)%site(n)%r  )
    kr = dot_product( kvec, r )
    X = -iC*kr
    water_molecule_electron_density_k(i,j,k,io) = water_molecule_electron_density_k(i,j,k,io) + exp(X) *(8.0 - q_o)

    do n= 2, 3 ! sum of hydrogen sites site


    r(1) = dot_product(   [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)]  ,  solvent(1)%site(n)%r  )
    r(2) = dot_product(   [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)]  ,  solvent(1)%site(n)%r  )
    r(3) = dot_product(   [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)]  ,  solvent(1)%site(n)%r  )
    kr = dot_product( kvec, r )
    X = -iC*kr
    water_molecule_electron_density_k(i,j,k,io) =  water_molecule_electron_density_k(i,j,k,io)  + exp(X)*(1.0 - q_h)

    end do ! loop over hydrogen sites


    end do !loop nx
    end do ! loop ny
    end do !loop nz

    end do !loop no

!$omp end do
!$omp end parallel

elseif (electron_density_type == 'dressed') then

!$omp parallel private(i, j, k, kvec, smootherfactor, r, kr, X, fac)
!$omp do

    do io = 1, no

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx/2+1
    kvec = [ grid%kx(i), grid%ky(j), grid%kz(k) ]
    q2 = (grid%kx(i)**2 + grid%ky(j)**2 + grid%kz(k)**2)/fourpi**2

    f_o = 0._dp
       do m = 1, 4
         f_o = f_o + a_o(m)*exp(-b_o(m)*q2)
       end do
    f_o = f_o + c_o
    f_o = f_o/8._dp

    f_h = 0._dp
    do m = 1, 4
    f_h = f_h + a_h(m)*exp(-b_h(m)*q2)
    end do
    f_h = f_h + c_h


    n= 1 ! oxygen site


    r(1) = dot_product(   [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)]  ,  solvent(1)%site(n)%r  )
    r(2) = dot_product(   [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)]  ,  solvent(1)%site(n)%r  )
    r(3) = dot_product(   [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)]  ,  solvent(1)%site(n)%r  )
    kr = dot_product( kvec, r )
    X = -iC*kr
    water_molecule_electron_density_k(i,j,k,io) = water_molecule_electron_density_k(i,j,k,io) + exp(X) *(8.0 - q_o)*f_o

    do n= 2, 3 ! sum of hydrogen sites site


    r(1) = dot_product(   [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)]  ,  solvent(1)%site(n)%r  )
    r(2) = dot_product(   [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)]  ,  solvent(1)%site(n)%r  )
    r(3) = dot_product(   [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)]  ,  solvent(1)%site(n)%r  )
    kr = dot_product( kvec, r )
    X = -iC*kr
    water_molecule_electron_density_k(i,j,k,io) =  water_molecule_electron_density_k(i,j,k,io)  + exp(X)*(1.0 - q_h)*f_h

    end do ! loop over hydrogen sites


    end do !loop nx
    end do ! loop ny
    end do !loop nz

    end do !loop no

!$omp end do
!$omp end parallel

else
    STOP 'option erroe in subroutine Get_water_molecule_reciprocal_electron_density'

endif

end subroutine Get_water_molecule_reciprocal_electron_density



subroutine Get_SPC_water_molecule_reciprocal_charge_density( SPC_water_charge_density_k, electron_density_type )

use module_grid, only: grid
use module_solvent, only: solvent
implicit none
integer :: nx, ny, nz, no, ns
integer :: i, j, k, n, s, io, d, m
real(dp)     :: r(3), kr, kvec(3), q2, f_o, f_h
complex(dp)  :: fac, X
complex(dp), parameter :: zeroc = (0._dp,0._dp), ic = (0._dp,1._dp)
real(dp), parameter :: epsdp = epsilon(1._dp)
real(dp) :: smootherfactor
real(dp) :: smootherradius
complex(dp), allocatable :: SPC_water_charge_density_k(:,:,:,:)
real(dp) :: fourpi = 4.0*acos(-1._dp)
real(dp), parameter :: q_h = 0.4238_dp, q_o = -0.8476_dp !SPC/E
real(dp), parameter, dimension(4) :: a_o = (/3.0485, 2.2868, 1.5463, 0.867 /) , b_o =(/ 13.2771, 5.7011, 0.3229, 32.9089 /)
real(dp), parameter, dimension(4)  :: a_h = (/ 0.48918, 0.262003, 0.196767, 0.049879 /), b_h(4) = (/ 20.6593, 7.74029, 49.5519, 2.20159 /)
real(dp), parameter :: c_o = 0.2508,  c_h = 0.001305
character(80) :: electron_density_type


nx = grid%nx
ny = grid%ny
nz = grid%nz
no = grid%no
ns = size(solvent) ! Count of solvent species

if( ns /= 1) stop 'init_solvent_molecule_pseudo_charge_density only works for ns = 1'

!allocate( SPC_water_charge_density_k(nx/2+1, ny, nz, no), SOURCE=zeroC )
SPC_water_charge_density_k = zeroC

if(electron_density_type == 'point_charges') then


!$omp parallel private(i, j, k, kvec, smootherfactor, r, kr, X, fac)
!$omp do
    do io = 1, no

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx/2+1

    do n=1, solvent(1)%nsite

    kvec = [ grid%kx(i), grid%ky(j), grid%kz(k) ]
    smootherfactor =  exp(-smootherradius**2 * sum( kvec**2 )/2._dp)

    r(1) = dot_product(   [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)]  ,  solvent(1)%site(n)%r  )   !  - grid%length(1)/2._dp
    r(2) = dot_product(   [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)]  ,  solvent(1)%site(n)%r  )   !  - grid%length(2)/2._dp
    r(3) = dot_product(   [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)]  ,  solvent(1)%site(n)%r  )   !  - grid%length(3)/2._dp
    kr = dot_product( kvec, r )
    X = -iC*kr
    SPC_water_charge_density_k(i,j,k,io) = SPC_water_charge_density_k(i,j,k,io) + solvent(1)%site(n)%q*exp(X) ! exact
    end do ! loop over solvent sites


    end do !loop nx
    end do ! loop ny
    end do !loop nz

    end do !loop no
!$omp end do
!$omp end parallel

elseif(electron_density_type == 'dressed') then
    smootherradius =  grid%dl(1)  ! smmothed over a grid bin length;  TAKE CARE for non cubic grids !!


!$omp parallel private(i, j, k, kvec, smootherfactor, r, kr, X, fac)
!$omp do

    do io = 1, no

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx/2+1
    kvec = [ grid%kx(i), grid%ky(j), grid%kz(k) ]
    smootherfactor =  exp(-smootherradius**2 * sum( kvec**2 )/2._dp)
    q2 = (grid%kx(i)**2 + grid%ky(j)**2 + grid%kz(k)**2)/fourpi**2

    f_o = 0._dp
    do m = 1, 4
    f_o = f_o + a_o(m)*exp(-b_o(m)*q2)
    end do
    f_o = f_o + c_o
    f_o = f_o/8._dp

    f_h = 0._dp
    do m = 1, 4
    f_h = f_h + a_h(m)*exp(-b_h(m)*q2)
    end do
    f_h = f_h + c_h


    n= 1 ! oxygen site


        r(1) = dot_product(   [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)]  ,  solvent(1)%site(n)%r  )
        r(2) = dot_product(   [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)]  ,  solvent(1)%site(n)%r  )
        r(3) = dot_product(   [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)]  ,  solvent(1)%site(n)%r  )
        kr = dot_product( kvec, r )
        X = -iC*kr
        SPC_water_charge_density_k(i,j,k,io) = SPC_water_charge_density_k(i,j,k,io) + exp(X)*( -(8.0_dp - q_o)*f_o + 8.0_dp)

    do n= 2, 3 ! sum of hydrogen sites site


        r(1) = dot_product(   [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)]  ,  solvent(1)%site(n)%r  )
        r(2) = dot_product(   [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)]  ,  solvent(1)%site(n)%r  )
        r(3) = dot_product(   [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)]  ,  solvent(1)%site(n)%r  )
        kr = dot_product( kvec, r )
        X = -iC*kr
        SPC_water_charge_density_k(i,j,k,io) = SPC_water_charge_density_k(i,j,k,io) + exp(X)*( -(1.0_dp - q_h)*f_h + 1.0_dp )

    end do ! loop over hydrogen sites


    end do !loop nx
    end do ! loop ny
    end do !loop nz

    end do !loop no

!$omp end do
!$omp end parallel

else
    STOP 'option error in subroutine Get_SPC_water_molecule_reciprocal_charge_density'
end if

end subroutine Get_SPC_water_molecule_reciprocal_charge_density

subroutine Get_molecule_reciprocal_pseudo_charge_density( molecule_charge_density_k, smootherradius )

use module_grid, only: grid
use module_solvent, only: solvent
implicit none
integer :: nx, ny, nz, no, ns
integer :: i, j, k, n, s, io, d, m
real(dp)     :: r(3), kr, kvec(3), q2, f_o, f_h
complex(dp)  :: fac, X
complex(dp), parameter :: zeroc = (0._dp,0._dp), ic = (0._dp,1._dp)
real(dp), parameter :: epsdp = epsilon(1._dp)
real(dp) :: smootherfactor
real(dp) :: smootherradius
complex(dp), allocatable :: molecule_charge_density_k(:,:,:,:)

nx = grid%nx
ny = grid%ny
nz = grid%nz
no = grid%no
ns = size(solvent) ! Count of solvent species

if( ns /= 1) stop 'init_solvent_molecule_pseudo_charge_density only works for ns = 1'

!allocate( SPC_water_charge_density_k(nx/2+1, ny, nz, no), SOURCE=zeroC )
molecule_charge_density_k = zeroC

!$omp parallel private(i, j, k, kvec, smootherfactor, r, kr, X, fac)
!$omp do
do io = 1, no

do k = 1, nz
do j = 1, ny
do i = 1, nx/2+1

do n=1, solvent(1)%nsite

kvec = [ grid%kx(i), grid%ky(j), grid%kz(k) ]
smootherfactor =  exp(-smootherradius**2 * sum( kvec**2 )/2._dp)

r(1) = dot_product(   [grid%Rotxx(io),grid%Rotxy(io),grid%Rotxz(io)]  ,  solvent(1)%site(n)%r  )   !  - grid%length(1)/2._dp
r(2) = dot_product(   [grid%Rotyx(io),grid%Rotyy(io),grid%Rotyz(io)]  ,  solvent(1)%site(n)%r  )   !  - grid%length(2)/2._dp
r(3) = dot_product(   [grid%Rotzx(io),grid%Rotzy(io),grid%Rotzz(io)]  ,  solvent(1)%site(n)%r  )   !  - grid%length(3)/2._dp
kr = dot_product( kvec, r )
X = -iC*kr
molecule_charge_density_k(i,j,k,io) = molecule_charge_density_k(i,j,k,io) + solvent(1)%site(n)%q*exp(X)*smootherfactor ! exact
end do ! loop over solvent sites


end do !loop nx
end do ! loop ny
end do !loop nz

end do !loop no
!$omp end do
!$omp end parallel

end subroutine Get_molecule_reciprocal_pseudo_charge_density


subroutine Get_solvent_electrostatic_potential( charge_density, electrostatic_potential )

use iso_c_binding
use precision_kinds
use module_solvent, only: solvent
use module_solute, only: solute
use module_grid, only: grid

IMPLICIT NONE
INTEGER(i2b) :: i, j, k, m1, m2, m3
integer(i2b) :: nx, ny, nz, no, ns
real(dp) :: kx, ky, kz, k2
real(dp) :: solute_net_charge

REAL(dp), allocatable :: charge_density(:,:,:)
REAL(dp), allocatable, intent(out) ::  electrostatic_potential(:,:,:) ! solvent pseudo-density
real(dp), allocatable :: fftw3inforward(:,:,:), fftw3outbackward(:,:,:)
complex(dp), allocatable :: fftw3outforward(:,:,:), fftw3inbackward(:,:,:)
integer(i4b) :: plan_forward, plan_backward ! fftw3 plan forward or backward
complex(dp), parameter :: zeroc=(0._dp,0._dp)
real(dp), parameter :: zero = 0._dp
real(dp), parameter :: twopi  =2._dp*acos(-1._dp)
real(dp), parameter :: fourpi =2._dp*twopi

include "fftw3.f03"

nx = grid%nx
ny = grid%ny
nz = grid%nz
no = grid%no
ns = solvent(1)%nspec
if( ns /= 1 ) stop 'subroutine get_final_Solvent_Pseudo_Density only for a single solvent component'


! allocate the arrays needed as input for fft (in_forward) or output for fft (out_forward)
! or needed as input for inverse fft (in_backward) etc.
allocate ( fftw3inforward   (nx    , ny,nz))
allocate ( fftw3outforward  (nx/2+1, ny,nz))
allocate ( fftw3outbackward (nx    , ny,nz))
allocate ( fftw3inbackward  (nx/2+1, ny,nz))
allocate( electrostatic_potential (nx, ny, nz) )

! prepare plans needed by fftw3
select case(dp)
case(c_double)
call dfftw_plan_dft_r2c_3d (plan_forward, nx, ny, nz, fftw3inforward, fftw3outforward, fftw_estimate)
call dfftw_plan_dft_c2r_3d (plan_backward, nx, ny, nz, fftw3inbackward , fftw3outbackward, fftw_estimate )
case(c_float)
call sfftw_plan_dft_r2c_3d (plan_forward, nx, ny, nz, fftw3inforward, fftw3outforward, fftw_estimate)
call sfftw_plan_dft_c2r_3d (plan_backward, nx, ny, nz, fftw3inbackward, fftw3outbackward, fftw_estimate )
end select


fftw3inforward = charge_density
select case(dp)
case(c_double)
call dfftw_execute(plan_forward)
case(c_float)
call sfftw_execute(plan_forward)
end select

DO k = 1, nz
    DO j = 1, ny
        DO i = 1, nx/2+1

        IF ( i<=nx/2 ) THEN
            m1 = i-1
            ELSE
            m1 = i-1-nx
        END IF

            IF ( j<=ny/2 ) THEN
                m2 = j-1
                ELSE
                m2 = j-1-ny
            END IF

        IF ( k<=nz/2 ) THEN
            m3 = k-1
            ELSE
            m3 = k-1-nz
        END IF

        kx = twopi*m1/grid%length(1)
        ky = twopi*m2/grid%length(2)
        kz = twopi*m3/grid%length(3)
        k2 = kx**2 + ky**2 + kz**2

        IF ( abs(k2) > epsilon(1._dp) ) THEN
            fftw3inbackward(i,j,k) = fftw3outforward(i,j,k) * fourpi/k2 ! in electrostatic units : V=-4pi rho
            ELSE
            fftw3inbackward(i,j,k) = (0._dp,0._dp)
        END IF

        END DO
    END DO
END DO

select case(dp)
case(c_double)
call dfftw_execute(plan_backward)
case(c_float)
call sfftw_execute(plan_backward)
end select

electrostatic_potential = fftw3outbackward/real(nx*ny*nz,dp)

! Only for cubic grids !
!solute_net_charge = sum(solute%site%q)
!If (solute_net_charge /= 0._dp) then
!   electrostatic_potential = electrostatic_potential &
!               - (1.0 - 1.0/solvent(1)%relativePermittivity)*(solute_net_charge &
!                 *2.837/grid%length(3) + twopi*solvent(1)%n0*0.8476/3.0 )
!   write(*,*) 'TypeB and C corrections applied to the electrostatic potential'
!end if

deallocate( fftw3inforward, fftw3outforward, fftw3outbackward, fftw3inbackward )
end subroutine Get_solvent_electrostatic_potential



end module module_postprocessing



