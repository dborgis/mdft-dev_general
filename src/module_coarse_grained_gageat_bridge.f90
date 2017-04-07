module module_coarse_grained_bridge

    use iso_c_binding, only: C_PTR, C_INT, C_INT32_T, C_INTPTR_T, C_DOUBLE_COMPLEX, C_DOUBLE, C_FUNPTR, C_SIZE_T, C_FLOAT, &
                             C_FLOAT_COMPLEX, C_CHAR
    use precision_kinds, only: dp

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
    real(dp), allocatable :: A(:)
    real(dp), parameter :: B = 15e-8_dp

    type(fft_type), protected :: fftRho, fftDfb

    public :: coarse_grained_bridge, init
   
   
contains

    !!! INITIALIZE !!!

    subroutine init
        use module_grid, only: grid
        use module_solvent, only: solvent
        use module_thermo, only: thermo

        implicit none
        
        real(dp) :: c00000, kT, n0
        integer :: is, ns
        integer :: nx, ny, nz

        nx = grid%nx
        ny = grid%ny
        nz = grid%nz
        ns=solvent(1)%nspec

        allocate ( density(nx,ny,nz) )
        allocate ( density_cg(nx,ny,nz) )
        allocate ( dfb_cg(nx,ny,nz) )
        allocate ( dfb(nx,ny,nz) )
        
        allocate ( kernel_k(nx/2+1,ny,nz) )
        allocate ( density_k(nx/2+1,ny,nz) )
        allocate ( dfb_cg_k(nx/2+1,ny,nz) )
        
        allocate ( A(ns) )
        


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


        !! TODO: l'extraire du fichier c directement mais cela implique de trop grosses modifs pour les permiers tests :)
        kT   = thermo%kbT
        
        do is=1,ns
        
          n0 = solvent(is)%n0
          c00000 = -13.6652260 / n0
          A(is) = kT*( 1.0_dp/n0**2 - c00000/(2.0_dp*n0) )

        end do
        
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
        real(dp), parameter :: sigma = 0.935_dp
        real(dp), parameter :: pi=acos(-1._dp)
        real(dp), parameter :: gaussianPrefactor = 1._dp/(sigma*sqrt(2._dp*pi))**3
        real(dp), parameter :: oneOverTwoSigmaSquare = 1._dp/(2.0_dp*sigma**2)
        real(dp), parameter :: kernelCutoff = 1.0e-10_dp ! If the kernel value is less than this one, we put 0.0_dp
        real(dp), parameter :: dSquareCutoff = -1.0_dp / oneOverTwoSigmaSquare * log( kernelCutoff / gaussianPrefactor )
        
        

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
        real(dp) :: dv, kT,n0
        real(dp) :: n_cg, deltaN_cg
        real(dp) :: vVoxel
        real(dp), parameter :: zerodp = 0._dp
        real(dp), parameter :: pi=acos(-1._dp)


        ns = solvent(1)%nspec
        kT = thermo%kbT
        dv = grid%dv

        fb = zerodp

        nx=grid%nx
        ny=grid%ny
        nz=grid%nz

        vVoxel = grid%dx * grid%dy * grid%dz

        if(is_init .eqv. .FALSE.) then
            call init
        end if

        
        do is=1,ns

          call grid%integrate_over_orientations( solvent(is)%xi**2 * solvent(is)%rho0, density)
          call compute_coarse_grained_density( density, density_cg )

          n0 = solvent(is)%n0

          do iz=1,nz
            do iy=1,ny
              do ix=1,nx
 
                n_cg = density_cg(ix,iy,iz)
          
                deltaN_cg = n_cg - n0
                fb = fb + A(is) * deltaN_cg**3 * vVoxel &
                        + B     * n_cg**2 * deltaN_cg**4 * vVoxel
                
                if(present(df)) then
                  if ( n_cg .lt. n0 ) then
                    dfb_cg(ix,iy,iz) = A(is) * 3.0_dp * deltaN_cg**2 &
                                     + B     *  ( 2.0_dp * n_cg * deltaN_cg**4 &
                                                + 4.0_dp * n_cg**2 * deltaN_cg**3 )
                  else 
                    dfb_cg(ix,iy,iz) = 0.0_dp
                  end if
                end if
                

              end do !ix
            end do !iy
          end do !iz
        end do !is

        if(present(df)) then
          call update_gradient(dfb_cg, df)
        end if

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
    
    
    
    
    
    
    
    subroutine update_gradient ( dfb_cg, df )
        
      use module_solvent, only: solvent 
      use module_grid, only: grid

      implicit none
 
      real(dp), intent(in) :: dfb_cg(:,:,:)
      real(dp), intent(out) :: df(:,:,:,:,:)
    
    
      integer :: is, ix, iy, iz
      integer :: ns, nx, ny, nz
      real(dp) :: wTot, rho0

      ns = solvent(1)%nspec
 
      nx=grid%nx
      ny=grid%ny
      nz=grid%nz
      
      wTot = sum(grid%w(:))
      
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
    
      do is=1,ns
        rho0 = solvent(is)%rho0
        do iz=1,nz
          do iy=1,ny
            do ix=1,nx

              df(:,ix,iy,iz,is) = df(:,ix,iy,iz,is) &
                                + 2.0_dp * solvent(is)%xi(:,ix,iy,iz) * rho0 * dfb(ix, iy, iz) * grid%w(:)
            
            end do
          end do
        end do
      end do
   
    end subroutine update_gradient
    
    
end module module_coarse_grained_bridge
