module module_coarse_grained_bridge 

    use iso_c_binding, only: C_PTR, C_INT, C_INT32_T, C_INTPTR_T, C_DOUBLE_COMPLEX, C_DOUBLE, C_FUNPTR, C_SIZE_T, C_FLOAT, &
                             C_FLOAT_COMPLEX, C_CHAR
    use precision_kinds, only: dp
!    use module_pressure_correction
    use module_energy_cproj_mrso, only: c000

    implicit none
    private
   
    !
    ! FFTW3 header - modern (fortran 2003) version. Expects iso_c_binding
    !
    include 'fftw3.f03'

    logical :: is_init = .false.

    type :: fft_type
        type(c_ptr) :: plan3dp, plan3dm 
    end type fft_type
    
    real(dp), allocatable :: density(:,:,:), density_cg(:,:,:), dfb_cg(:,:,:), dfb(:,:,:)
    complex(dp), allocatable :: kernel_k(:,:,:), density_k(:,:,:), dfb_cg_k(:,:,:)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!                                                          !!!!!
    !!!!!                 BRIDGE PARAMETERS A AND B                !!!!!
    !!!!!                   FOR MORE DETAILS SEE                   !!!!!
    !!!!!                    GAGEAT ET AL 2017                     !!!!!
    !!!!!                                                          !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type :: cgb_type
         real(dp) :: A3
         real(dp) :: B6, sigma
        character(3) :: version
    end type cgb_type
    type(cgb_type) :: cgb


    type(fft_type), protected :: fftRho, fftDfb

    public :: coarse_grained_bridge, cgb
   
   
contains

    !!! INITIALIZE !!!

    subroutine init
        use module_grid, only: grid
        use module_solvent, only: solvent
        use module_thermo, only: thermo
        use module_input, only: getinput

        implicit none
        
        real(dp) ::  n0, kT, B6, sigma
        integer :: is, ns
        integer :: nx, ny, nz

        nx = grid%nx
        ny = grid%ny
        nz = grid%nz
        ns=solvent(1)%nspec
        if(ns > 1) STOP 'Coarse-grained bridge for one component only: better check it out !'
        n0 = solvent(1)%n0
        kT   = thermo%kbT

        cgb%B6 = getinput%dp('cgb_B6', defaultvalue=0._dp)
        cgb%sigma = getinput%dp('cgb_sigma', defaultvalue=1._dp)
        cgb%version = getinput%char('cgb_version', defaultvalue='sla')

!         write(*,*) 'cgb_version =', cgb%version
!        if(cgb%version == 'wda') then
!            write(*,*) 'cgb: new WDA version'
!        else if (cgb%version == 'sla') then
!            write(*,*) 'cgb: semi local approximation'
!        else
!            STOP 'wrong coarse-grained bridge option'
!        endif

        allocate ( density(nx,ny,nz) )
        allocate ( density_cg(nx,ny,nz) )
        allocate ( dfb_cg(nx,ny,nz) )
        allocate ( dfb(nx,ny,nz) )
        
        allocate ( kernel_k(nx/2+1,ny,nz) )
        allocate ( density_k(nx/2+1,ny,nz) )
        allocate ( dfb_cg_k(nx/2+1,ny,nz) )

        


        select case(dp)
          case(c_double)
          
            call dfftw_plan_dft_r2c_3d( fftRho%plan3dp, nx, ny, nz, density, density_k, FFTW_MEASURE)
            call dfftw_plan_dft_c2r_3d( fftRho%plan3dm, nx, ny, nz, density_k, dfb_cg, FFTW_MEASURE)

            call dfftw_plan_dft_r2c_3d( fftDfb%plan3dp, nx, ny, nz, dfb_cg, dfb_cg_k, FFTW_MEASURE)
            call dfftw_plan_dft_c2r_3d( fftDfb%plan3dm, nx, ny, nz, dfb_cg_k, dfb, FFTW_MEASURE)
          
          case(c_float)
          
            call sfftw_plan_dft_r2c_3d( fftRho%plan3dp, nx, ny, nz, density, density_k, FFTW_MEASURE)
            call sfftw_plan_dft_c2r_3d( fftRho%plan3dm, nx, ny, nz, density_k, dfb_cg, FFTW_MEASURE)

            call sfftw_plan_dft_r2c_3d( fftDfb%plan3dp, nx, ny, nz, dfb_cg, dfb_cg_k, FFTW_MEASURE)
            call sfftw_plan_dft_c2r_3d( fftDfb%plan3dm, nx, ny, nz, dfb_cg_k, dfb, FFTW_MEASURE)
        
        end select


        call fill_kernel

        
        is_init = .true.
        
    end subroutine init
        

    subroutine fill_kernel
        
        use module_grid, only: grid

        implicit none

        real(dp) :: dSquare
        integer :: nx, ny, nz
        integer :: ix, iy, iz
        real(dp), allocatable :: xPbcSquare(:), yPbcSquare(:), zPbcSquare(:)
        real(dp) :: xPbc, yPbc, zPbc
        real(dp), allocatable :: kernel(:,:,:)
        type(c_ptr) :: fftPlan3d
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!                                                          !!!!!
        !!!!!              BRIDGE PARAMETER SIGMA (KERNEL)             !!!!!
        !!!!!                   FOR MORE DETAILS SEE                   !!!!!
        !!!!!                    GAGEAT ET AL 2017                     !!!!!
        !!!!!                                                          !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(dp) :: pi, gaussianPrefactor, oneOverTwoSigmaSquare, dSquareCutoff, sigma
        real(dp), parameter :: kernelCutoff = 1.0e-10_dp ! If the kernel value is less than this one, we put 0.0_dp

        pi=acos(-1._dp)
        sigma = cgb%sigma
        gaussianPrefactor = 1._dp/(sigma*sqrt(2._dp*pi))**3
        oneOverTwoSigmaSquare = 1._dp/(2.0_dp*sigma**2)
        dSquareCutoff = -1.0_dp / oneOverTwoSigmaSquare * log( kernelCutoff / gaussianPrefactor )

        nx=grid%nx
        ny=grid%ny
        nz=grid%nz

        allocate( xPbcSquare(nx) )
        allocate( yPbcSquare(ny) )
        allocate( zPbcSquare(nz) )

        do ix=1,nx
          xPbc = grid%dx * min( ix-1, nx-(ix-1) )
          xPbcSquare(ix) = xPbc*xPbc
        end do

        do iy=1,ny
          yPbc = grid%dy * min( iy-1, ny-(iy-1) )
          yPbcSquare(iy) = yPbc*yPbc
        end do

        do iz=1,nz
          zPbc = grid%dz * min( iz-1, nz-(iz-1) )
          zPbcSquare(iz) = zPbc*zPbc
        end do

        allocate (kernel(nx,ny,nz), source=0.0_dp)


        do iz=1,nz
          do iy=1,ny
            do ix=1,nx
               
              dSquare= xPbcSquare(ix) + yPbcSquare(iy) + zPbcSquare(iz)
              if ( dSquare .lt. dSquareCutoff ) then
                kernel(ix, iy, iz) = gaussianPrefactor * exp( - dSquare * oneOverTwoSigmaSquare )
              end if

            end do !ix
          end do !iy
        end do !iz

        deallocate(xPbcSquare)
        deallocate(yPbcSquare)
        deallocate(zPbcSquare)
       
        select case(dp)
          case(c_double)
            call dfftw_plan_dft_r2c_3d( fftPlan3d, nx, ny, nz, kernel, kernel_k, FFTW_MEASURE )
            call dfftw_execute_dft_r2c( fftPlan3d, kernel, kernel_k )
            call dfftw_destroy_plan( fftPlan3d )
          case(c_float)
            call sfftw_plan_dft_r2c_3d( fftPlan3d, nx, ny, nz, kernel, kernel_k, FFTW_MEASURE )
            call sfftw_execute_dft_r2c( fftPlan3d, kernel, kernel_k )
            call sfftw_destroy_plan( fftPlan3d )
        end select

        deallocate(kernel)
                
        kernel_k(:,:,:) = kernel_k(:,:,:) / real(nx*ny*nz) * (grid%lx * grid%ly * grid%lz) ! normalization to be use in convolution
    
    end subroutine fill_kernel
    
    
    !!! COMPUTING !!!
    
    
    subroutine coarse_grained_bridge(fb, df)

        use precision_kinds, only: dp
        use module_thermo, only: thermo
        use module_solvent, only: solvent 
        use module_grid, only: grid

        implicit none

        real(dp), intent(out) :: fb
        real(dp), intent(inout), contiguous, optional :: df(:,:,:,:,:)
        integer :: is, ix, iy, iz
        integer :: ns, nx, ny, nz
        real(dp) :: dv, kT,n0, Phib, A3, B6, sigma

        real(dp) :: n_cg, DeltaN_cg, deltaN, facteur
        real(dp) :: vVoxel, c00_000
        real(dp), parameter :: zero = 0._dp, one = 1.0_dp
        real(dp), parameter :: pi=acos(-1._dp)


        ns = solvent(1)%nspec
        kT = thermo%kbT
        dv = grid%dv

        fb = zero

        nx=grid%nx
        ny=grid%ny
        nz=grid%nz

        vVoxel = grid%dx * grid%dy * grid%dz
        n0 = solvent(1)%n0
        facteur = kT*vVoxel

        if(is_init .eqv. .FALSE.) then
            call init
        end if


        B6 = cgb%B6
        sigma = cgb%sigma

        A3 =  1.0_dp - c000/2.0_dp
!          write(*,*) "A3 = ",A3, ' B6 =', B6, ' sigma =', sigma, 'cgb_version = ', cgb%version
        if( solvent(1)%name == "spce" .and. grid%mmax == 0)  c00_000 = -14.64872
        if( solvent(1)%name == "spce" .and. grid%mmax > 0)   c00_000 = -13.75
        if( solvent(1)%name == "tip3p" ) c00_000 = -11.9078
        A3 =  1.0_dp - c00_000/2.0_dp




! We work here with reduced densities rho/rho_0

          call grid%integrate_over_orientations( solvent(1)%xi**2 * solvent(1)%rho0, density)
          density = density/n0
          call compute_coarse_grained_density( density, density_cg )

 if (cgb%version == 'sla') then

      IF(present(df)) then
          do iz=1,nz
            do iy=1,ny
              do ix=1,nx
 
                n_cg = density_cg(ix,iy,iz)
                deltaN = density( ix, iy, iz) - one
                deltaN_cg = density_cg(ix,iy,iz) - one

                  if ( n_cg < one ) then
                     Phib = A3 * deltaN_cg**3 + B6 * n_cg**2 * deltaN_cg**4
                     dfb_cg(ix,iy,iz) = kT*A3 * 3.0_dp * deltaN_cg**2 &
                                      + kT*B6 *  ( 2.0_dp * n_cg * deltaN_cg**4 &
                                                 + 4.0_dp * n_cg**2 * deltaN_cg**3 )
                  else
                    Phib = A3 * deltaN_cg**3
                    dfb_cg(ix,iy,iz) = kT*A3 * 3.0_dp * deltaN_cg**2
                  ! Phib = zero
                  ! dfb_cg(ix,iy,iz) = zero
                  end if

                  fb = fb + facteur * n0 *  Phib

              end do !ix
            end do !iy
          end do !iz


        select case(dp)
        case(c_double)
          call dfftw_execute_dft_r2c( fftDfb%plan3dp, dfb_cg, dfb_cg_k)
        case(c_float)
          call sfftw_execute_dft_r2c( fftDfb%plan3dp, dfb_cg, dfb_cg_k)
        end select

        dfb_cg_k(:,:,:) = dfb_cg_k(:,:,:) * kernel_k(:,:,:)

        select case(dp)
        case(c_double)
          call dfftw_execute_dft_c2r( fftDfb%plan3dm, dfb_cg_k, dfb)
        case(c_float)
          call sfftw_execute_dft_c2r( fftDfb%plan3dm, dfb_cg_k, dfb)
        end select

        dfb = dfb / real(grid%nx*grid%ny*grid%nz)

        call update_gradient(dfb, df)

     ELSE

     do iz=1,nz
       do iy=1,ny
         do ix=1,nx

         n_cg = density_cg(ix,iy,iz)
         deltaN = density( ix, iy, iz) - one
         deltaN_cg = density_cg(ix,iy,iz) - one

         if ( n_cg < one ) then
           Phib = A3 * deltaN_cg**3 + B6 * n_cg**2 * deltaN_cg**4
        else
           Phib = A3 * deltaN_cg**3
        end if

          fb = fb + facteur * n0 *  Phib

        end do !ix
      end do !iy
    end do !iz

  END IF  ! present(df)


ELSE IF (cgb%version == 'wda') then ! new WDA version

    IF(present(df)) then
       do iz=1,nz
         do iy=1,ny
           do ix=1,nx

            n_cg = density_cg(ix,iy,iz)
            deltaN = density( ix, iy, iz) - one
            deltaN_cg = density_cg(ix,iy,iz) - one

            if ( n_cg < one ) then
                    dfb_cg(ix,iy,iz) = A3 * 2.0_dp * deltaN_cg &
                                     + B6 * ( 2.0_dp * n_cg * deltaN_cg**3 &
                                              + 3.0_dp * n_cg**2 * deltaN_cg**2 )
            else
                    dfb_cg(ix,iy,iz) = A3 * 2.0_dp * deltaN_cg
        !            dfb_cg(ix,iy,iz) =  zero
            end if

            end do !ix
          end do !iy
       end do !iz


        select case(dp)
        case(c_double)
        call dfftw_execute_dft_r2c( fftDfb%plan3dp, dfb_cg, dfb_cg_k)
        case(c_float)
        call sfftw_execute_dft_r2c( fftDfb%plan3dp, dfb_cg, dfb_cg_k)
        end select

        dfb_cg_k(:,:,:) = dfb_cg_k(:,:,:) * kernel_k(:,:,:)

        select case(dp)
        case(c_double)
        call dfftw_execute_dft_c2r( fftDfb%plan3dm, dfb_cg_k, dfb)
        case(c_float)
        call sfftw_execute_dft_c2r( fftDfb%plan3dm, dfb_cg_k, dfb)
        end select

        dfb = dfb / real(grid%nx*grid%ny*grid%nz)

     END IF   ! present(df)

     do iz=1,nz
      do iy=1,ny
        do ix=1,nx

         n_cg = density_cg(ix,iy,iz)
         deltaN = density( ix, iy, iz) - one
         deltaN_cg = n_cg - one


         if( n_cg < one )then
               Phib = A3 * deltaN_cg**2 + B6 * n_cg**2 * deltaN_cg**3
         else
               Phib = A3 * deltaN_cg**2
  !            Phib = zero
         end if
              fb = fb + kT * vVoxel * n0 * deltaN * Phib

             if (present(df)) dfb_cg(ix,iy,iz) =  kT * ( deltaN * dfb(ix,iy,iz) + Phib )

         end do !ix
       end do !iy
     end do !iz

     if (present(df)) call update_gradient(dfb_cg, df)

ENDIF ! version gageat or new WDA

    end subroutine coarse_grained_bridge


    subroutine compute_coarse_grained_density ( density,  density_cg )

      use module_grid, only: grid

      implicit none

      real(dp), intent(in) :: density(:,:,:)
      real(dp), intent(out) :: density_cg(:,:,:)

      select case(dp)
        case(c_double)
          call dfftw_execute_dft_r2c( fftRho%plan3dp, density, density_k)
        case(c_float)
          call sfftw_execute_dft_r2c( fftRho%plan3dp, density, density_k)
      end select

      density_k(:,:,:) = density_k(:,:,:) * kernel_k(:,:,:)

      select case(dp)
        case(c_double)
          call dfftw_execute_dft_c2r( fftRho%plan3dm, density_k, density_cg)
        case(c_float)
          call sfftw_execute_dft_c2r( fftRho%plan3dm, density_k, density_cg)
      end select
    
      density_cg = density_cg / real(grid%nx * grid%ny * grid%nz)
    
    end subroutine compute_coarse_grained_density

    
    
    subroutine update_gradient ( dfb, df )
        
      use module_solvent, only: solvent 
      use module_grid, only: grid

      implicit none
 
      real(dp), intent(in) :: dfb(:,:,:)
      real(dp), intent(out) :: df(:,:,:,:,:)
    
    
      integer :: ix, iy, iz
      integer :: nx, ny, nz
      real(dp) :: rho0

        rho0 = solvent(1)%rho0

        do iz=1,grid%nz
          do iy=1,grid%ny
            do ix=1,grid%nx

              df(:,ix,iy,iz,1) = df(:,ix,iy,iz,1) + 2.0_dp * solvent(1)%xi(:,ix,iy,iz) * rho0 * dfb(ix, iy, iz) * grid%w(:)
            
            end do
          end do
        end do
   
    end subroutine update_gradient

    
end module  module_coarse_grained_bridge
