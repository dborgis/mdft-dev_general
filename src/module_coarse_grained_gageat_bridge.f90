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
    complex(dp), allocatable :: kernel_k(:,:,:)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!                                                          !!!!!
    !!!!!                 BRIDGE PARAMETERS A AND B                !!!!!
    !!!!!                   FOR MORE DETAILS SEE                   !!!!!
    !!!!!                    GAGEAT ET AL 2017                     !!!!!
    !!!!!                                                          !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(dp), allocatable :: A(:)
    real(dp), parameter :: B = 15e-8_dp

    
    public :: coarse_grained_bridge
   
   
contains

    !!! INITIALIZE !!!

    subroutine init
        use module_grid, only: grid
        use module_solvent, only: solvent
        use module_thermo, only: thermo

        implicit none
        
        real(dp) :: c00000, kT, n0
        integer :: is, ns


        ns=solvent(1)%nspec

        allocate ( kernel_k(grid%nx/2+1,grid%ny,grid%nz), source=(0._dp,0._dp) )
        
        allocate ( A(ns) )
        

        call fill_kernel( kernel_k) ! must be already allocated


        !! TODO: l'extraire du fichier c directement mais cela implique de trop grosses modifs pour les permiers tests :)
        kT   = thermo%kbT
        
        do is=1,ns
        
          n0 = solvent(is)%n0
          c00000 = -13.6652260 / n0
          A(is) = kT*( 1.0_dp/n0**2 - c00000/(2.0_dp*n0) )

        end do
        
        is_init = .true.
        
    end subroutine init
        

    subroutine fill_kernel( kernel_k)
        
        use module_grid, only: grid

        implicit none

        complex(dp), intent(out) :: kernel_k(:,:,:)
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

        if( .not. allocated( kernel) ) allocate (kernel(nx,ny,nz), source=0.0_dp)


        do iz=1,nz
          do iy=1,ny
            do ix=1,nx
               
              dSquare= xPbcSquare(ix) + yPbcSquare(iy) + zPbcSquare(iz)
              if ( dSquare .lt. dSquareCutoff ) then
                kernel(ix, iy, iz) = gaussianPrefactor * exp( - dSquare * oneOverTwoSigmaSquare )
              else
                kernel(ix, iy, iz) = 0._dp
              end if

            end do !ix
          end do !iy
        end do !iz
       
        deallocate(xPbcSquare)
        deallocate(yPbcSquare)
        deallocate(zPbcSquare)
       
        select case(dp)
          case(c_double)
            call dfftw_plan_dft_r2c_3d( fftPlan3d, nx, ny, nz, kernel, kernel_k, FFTW_ESTIMATE )
            call dfftw_execute_dft_r2c( fftPlan3d, kernel, kernel_k )
            call dfftw_destroy_plan( fftPlan3d )
          case(c_float)
            call sfftw_plan_dft_r2c_3d( fftPlan3d, nx, ny, nz, kernel, kernel_k, FFTW_ESTIMATE )
            call sfftw_execute_dft_r2c( fftPlan3d, kernel, kernel_k )
            call sfftw_destroy_plan( fftPlan3d )
        end select

        if( allocated(kernel)) deallocate(kernel)
                
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
        real(dp), intent(inout), optional :: df(:,:,:,:,:)
        integer :: is, ix, iy, iz
        integer :: ns, nx, ny, nz
        real(dp) :: dv, kT,n0
        real(dp) :: n_cg, deltaN_cg
        real(dp) :: vVoxel
        real(dp), parameter :: zerodp = 0._dp
        real(dp), parameter :: pi=acos(-1._dp)
        real(dp), allocatable :: density(:,:,:)
        real(dp), allocatable :: density_cg(:,:,:), dfb_cg(:,:,:)
        
        
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

        allocate ( density(grid%nx,grid%ny,grid%nz), source=0.0_dp )        
        allocate ( density_cg(grid%nx,grid%ny,grid%nz), source=0.0_dp )
        allocate ( dfb_cg(grid%nx,grid%ny,grid%nz), source=0.0_dp )

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

        deallocate ( density_cg )
        deallocate ( density )        


        if(present(df)) then
          call update_gradient(dfb_cg, df)
        end if
        
        deallocate ( dfb_cg )

    end subroutine coarse_grained_bridge
    
    
    
    
    
    
    
    
    subroutine compute_coarse_grained_density ( density,  density_cg )

      use module_grid, only: grid

      implicit none
 
      real(dp), intent(in) :: density(:,:,:)
      real(dp), intent(inout) :: density_cg(:,:,:)
      complex(dp), allocatable :: density_k(:,:,:)
      type :: fft_type
        type(c_ptr) :: plan3dp, plan3dm 
      end type fft_type
      type(fft_type) :: fftRho

  
  
      allocate ( density_k(grid%nx/2+1,grid%ny,grid%nz), source=(0._dp,0._dp) )  
    
      select case(dp)
        case(c_double)
          call dfftw_plan_dft_r2c_3d( fftRho%plan3dp, grid%nx, grid%ny, grid%nz, density, density_k, FFTW_ESTIMATE)
          call dfftw_execute_dft_r2c( fftRho%plan3dp, density, density_k)
          call dfftw_destroy_plan( fftRho%plan3dp )
       case(c_float)
          call sfftw_plan_dft_r2c_3d( fftRho%plan3dp, grid%nx, grid%ny, grid%nz, density, density_k, FFTW_ESTIMATE)
          call sfftw_execute_dft_r2c( fftRho%plan3dp, density, density_k)
          call sfftw_destroy_plan( fftRho%plan3dp )
      end select
        
      density_k(:,:,:) = density_k(:,:,:) * kernel_k(:,:,:)
        
      select case(dp)
        case(c_double)
          call dfftw_plan_dft_c2r_3d( fftRho%plan3dm, grid%nx, grid%ny, grid%nz, density_k, density_cg, FFTW_ESTIMATE)
          call dfftw_execute_dft_c2r( fftRho%plan3dm, density_k, density_cg)
          call dfftw_destroy_plan( fftRho%plan3dm )
        case(c_float)
          call sfftw_plan_dft_c2r_3d( fftRho%plan3dm, grid%nx, grid%ny, grid%nz, density_k, density_cg, FFTW_ESTIMATE)
          call sfftw_execute_dft_c2r( fftRho%plan3dm, density_k, density_cg)
          call sfftw_destroy_plan( fftRho%plan3dm )
      end select
    
      deallocate ( density_k )  

      density_cg = density_cg / real(grid%nx * grid%ny * grid%nz)
    
    end subroutine compute_coarse_grained_density
    
    
    
    
    
    
    
    subroutine update_gradient ( dfb_cg, df )
        
      use module_solvent, only: solvent 
      use module_grid, only: grid

      implicit none
 
      real(dp), allocatable, intent(in) :: dfb_cg(:,:,:)
      real(dp), intent(inout) :: df(:,:,:,:,:)
      real(dp), allocatable :: dfb(:,:,:)
      complex(dp), allocatable :: dfb_cg_k(:,:,:)
    
      type :: fft_type
        type(c_ptr) :: plan3dp, plan3dm 
      end type fft_type
      type(fft_type) :: fftDfb

    
      integer :: is, ix, iy, iz
      integer :: ns, nx, ny, nz
      real(dp) :: wTot, rho0

      ns = solvent(1)%nspec
 
      nx=grid%nx
      ny=grid%ny
      nz=grid%nz
      
      wTot = sum(grid%w(:))
      
      allocate ( dfb_cg_k(grid%nx/2+1,grid%ny,grid%nz), source=(0._dp,0._dp) )
      
      select case(dp)
        case(c_double)
          call dfftw_plan_dft_r2c_3d( fftDfb%plan3dp, grid%nx, grid%ny, grid%nz, dfb_cg, dfb_cg_k, FFTW_ESTIMATE)
          call dfftw_execute_dft_r2c( fftDfb%plan3dp, dfb_cg, dfb_cg_k)
          call dfftw_destroy_plan( fftDfb%plan3dp )
        case(c_float)
          call sfftw_plan_dft_r2c_3d( fftDfb%plan3dp, grid%nx, grid%ny, grid%nz, dfb_cg, dfb_cg_k, FFTW_ESTIMATE)
          call sfftw_execute_dft_r2c( fftDfb%plan3dp, dfb_cg, dfb_cg_k)
          call sfftw_destroy_plan( fftDfb%plan3dp )
      end select
        
      dfb_cg_k(:,:,:) = dfb_cg_k(:,:,:) * kernel_k(:,:,:)

      allocate ( dfb(grid%nx,grid%ny,grid%nz), source=0._dp )
        
      select case(dp)
        case(c_double)
          call dfftw_plan_dft_c2r_3d( fftDfb%plan3dm, grid%nx, grid%ny, grid%nz, dfb_cg_k, dfb, FFTW_ESTIMATE)
          call dfftw_execute_dft_c2r( fftDfb%plan3dm, dfb_cg_k, dfb)
          call dfftw_destroy_plan( fftDfb%plan3dm )
        case(c_float)
          call sfftw_plan_dft_c2r_3d( fftDfb%plan3dm, grid%nx, grid%ny, grid%nz, dfb_cg_k, dfb, FFTW_ESTIMATE)
          call sfftw_execute_dft_c2r( fftDfb%plan3dm, dfb_cg_k, dfb)
          call sfftw_destroy_plan( fftDfb%plan3dm )
      end select

      deallocate ( dfb_cg_k )


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
   
      deallocate ( dfb )
   
    end subroutine update_gradient
    
    
end module module_coarse_grained_bridge
