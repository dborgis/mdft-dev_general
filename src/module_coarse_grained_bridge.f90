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
    
    real(dp), allocatable :: A(:)
    real(dp) :: B


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


        ns=solvent(1)%nspec

        allocate ( density(grid%nx,grid%ny,grid%nz) )
        allocate ( density_cg(grid%nx,grid%ny,grid%nz) )
        allocate ( dfb_cg(grid%nx,grid%ny,grid%nz) )
        allocate ( dfb(grid%nx,grid%ny,grid%nz) )
        
        allocate ( kernel_k(grid%nx,grid%ny,grid%nz/2+1) )
        allocate ( density_k(grid%nx,grid%ny,grid%nz/2+1) )
        allocate ( dfb_cg_k(grid%nx,grid%ny,grid%nz/2+1) )
        
        allocate ( A(ns) )
        


        select case(dp)
          case(c_double)
          
            call dfftw_plan_dft_r2c_3d( fftRho%plan3dp, grid%nx, grid%ny, grid%nz, density, density_k, FFTW_MEASURE)
            call dfftw_plan_dft_c2r_3d( fftRho%plan3dm, grid%nx, grid%ny, grid%nz, density_k, dfb_cg, FFTW_MEASURE)

            call dfftw_plan_dft_r2c_3d( fftDfb%plan3dp, grid%nx, grid%ny, grid%nz, dfb_cg, dfb_cg_k, FFTW_MEASURE)
            call dfftw_plan_dft_c2r_3d( fftDfb%plan3dm, grid%nx, grid%ny, grid%nz, dfb_cg_k, dfb, FFTW_MEASURE)
          
          case(c_float)
          
            call sfftw_plan_dft_r2c_3d( fftRho%plan3dp, grid%nx, grid%ny, grid%nz, density, density_k, FFTW_MEASURE)
            call sfftw_plan_dft_c2r_3d( fftRho%plan3dm, grid%nx, grid%ny, grid%nz, density_k, dfb_cg, FFTW_MEASURE)

            call sfftw_plan_dft_r2c_3d( fftDfb%plan3dp, grid%nx, grid%ny, grid%nz, dfb_cg, dfb_cg_k, FFTW_MEASURE)
            call sfftw_plan_dft_c2r_3d( fftDfb%plan3dm, grid%nx, grid%ny, grid%nz, dfb_cg_k, dfb, FFTW_MEASURE)
        
        end select


        call fill_kernel


        !! m_thermo.kBT/puissance(m_solvent.rhoZero,2) - m_thermo.kBT*cZero/m_solvent.rhoZero/2.0 
        !! TODO: l'extraire du fichier c directement mais cela implique de trop grosses modifs pour les permiers tests :)
        kT   = thermo%kbT
        
        do is=1,ns
        
          n0 = solvent(is)%n0 ! 3.3289100974798203E-002
          c00000 = -13.6652260 / n0   ! -410.50150296973749
          A(is) = kT*( 1.0_dp/n0**2 - c00000/(2.0_dp*n0) )    ! c0 sphérique

        end do
        
        B = 15e-8_dp

        is_init = .true.
        
    end subroutine init
        

    subroutine fill_kernel
        
        use module_grid, only: grid

        implicit none

        real(dp) :: sigma, d
        integer :: i
        integer :: nx, ny, nz
        integer :: ix, iy, iz
        real(dp), allocatable :: x(:), y(:), z(:)
        real(dp) :: xPbc, yPbc, zPbc
        real(dp), allocatable :: kernel(:,:,:)
        type(c_ptr) :: fftPlan3d
        real(dp), parameter :: pi=acos(-1._dp)

        nx=grid%nx
        ny=grid%ny
        nz=grid%nz

        sigma = 0.935_dp !2.36_dp
        
        allocate( x(nx) ,source= [( (real(i-1,dp)*grid%dx) ,i=1,nx)] )
        allocate( y(ny) ,source= [( (real(i-1,dp)*grid%dy) ,i=1,ny)] )
        allocate( z(nz) ,source= [( (real(i-1,dp)*grid%dz), i=1,nz)] )

        allocate (kernel(grid%nx,grid%ny,grid%nz), source=0.0_dp)


        do iz=1,grid%nz
          zPbc = min( z(iz), grid%lz-z(iz) )
          do iy=1,grid%ny
            yPbc = min( y(iy), grid%ly-y(iy) )
            do ix=1,grid%nx
              xPbc = min( x(ix), grid%lx-x(ix) )
                
                d=sqrt( xPbc**2 + yPbc**2 + zPbc**2 )
                
               !kernel(ix, iy, iz) = exp( -0.5*(d/sigma)**2 ) ! formula in fourier space n 1d
               kernel(ix, iy, iz) = 1./(sigma*sqrt(2.*pi))**3 * exp( -(d**2)/(2.0_dp*sigma**2) )

            end do !ix
          end do !iy
        end do !iz
       
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
        
        kernel_k(:,:,:) = kernel_k(:,:,:) / real(grid%nx*grid%ny*grid%nz) * (grid%lx * grid%ly * grid%lz) ! normalization to be use in convolution
        
        deallocate(x)
        deallocate(y)
        deallocate(z)
        deallocate(kernel)
    
    end subroutine fill_kernel
    
    
    !!! COMPUTING !!!
    
    
    subroutine coarse_grained_bridge(fb, df)

        use precision_kinds, only: dp
        use module_thermo, only: thermo
        use module_solvent, only: solvent 
        use module_grid, only: grid
        !use module_input, only: getinput

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
          !density => 3.32891010E-02 without solute
          call compute_coarse_grained_density( density, density_cg )
          !density_cg => 3.32891010E-02 without solute

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
    
    
    
    
    
    
    
    
    subroutine compute_coarse_grained_density ( density,  density_cg )  !! CHECK => OK

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
    
      call write_to_cube_file (density, "output/toto.cube                                               ")
      call write_to_cube_file (density_cg, "output/toto_cg.cube                                                  ")

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

              df(:,ix,iy,iz,is) = df(:,ix,iy,iz,is) + 2.0_dp * solvent(is)%xi(:,ix,iy,iz) * rho0 * dfb(ix, iy, iz) * grid%w(:) / wTot 
            
            end do
          end do
        end do
      end do
   
    end subroutine update_gradient
    
    
end module module_coarse_grained_bridge
