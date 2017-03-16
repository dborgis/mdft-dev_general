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
    
    real(dp), allocatable :: density(:,:,:), coarseGrainedDensity(:,:,:), coarseGrainedDfb(:,:,:), dfb(:,:,:)
    complex(dp), allocatable :: fourierKernel(:,:,:), fourierDensity(:,:,:), fourierCoarseGrainedDfb(:,:,:)
    
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
        
        real(dp) :: c00000, kT, rho0
        integer :: is, ns


        ns=solvent(1)%nspec

        allocate ( density(grid%nx,grid%ny,grid%nz) )
        allocate ( coarseGrainedDensity(grid%nx,grid%ny,grid%nz) )
        allocate ( coarseGrainedDfb(grid%nx,grid%ny,grid%nz) )
        allocate ( dfb(grid%nx,grid%ny,grid%nz) )
        
        allocate ( fourierKernel(grid%nx,grid%ny,grid%nz/2+1) )
        allocate ( fourierDensity(grid%nx,grid%ny,grid%nz/2+1) )
        allocate ( fourierCoarseGrainedDfb(grid%nx,grid%ny,grid%nz/2+1) )
        
        allocate ( A(ns) )
        


        select case(dp)
          case(c_double)
          
            call dfftw_plan_dft_r2c_3d( fftRho%plan3dp, grid%nx, grid%ny, grid%nz, density, fourierDensity, FFTW_MEASURE)
            call dfftw_plan_dft_c2r_3d( fftRho%plan3dm, grid%nx, grid%ny, grid%nz, fourierDensity, coarseGrainedDfb, FFTW_MEASURE)

            call dfftw_plan_dft_r2c_3d( fftDfb%plan3dp, grid%nx, grid%ny, grid%nz, coarseGrainedDfb, fourierCoarseGrainedDfb, FFTW_MEASURE)
            call dfftw_plan_dft_c2r_3d( fftDfb%plan3dm, grid%nx, grid%ny, grid%nz, fourierCoarseGrainedDfb, dfb, FFTW_MEASURE)
          
          case(c_float)
          
            call dfftw_plan_dft_r2c_3d( fftRho%plan3dp, grid%nx, grid%ny, grid%nz, density, fourierDensity, FFTW_MEASURE)
            call dfftw_plan_dft_c2r_3d( fftRho%plan3dm, grid%nx, grid%ny, grid%nz, fourierDensity, coarseGrainedDfb, FFTW_MEASURE)

            call dfftw_plan_dft_r2c_3d( fftDfb%plan3dp, grid%nx, grid%ny, grid%nz, coarseGrainedDfb, fourierCoarseGrainedDfb, FFTW_MEASURE)
            call dfftw_plan_dft_c2r_3d( fftDfb%plan3dm, grid%nx, grid%ny, grid%nz, fourierCoarseGrainedDfb, dfb, FFTW_MEASURE)
        
        end select


        call fill_kernel


        !! m_thermo.kBT/puissance(m_solvent.rhoZero,2) - m_thermo.kBT*cZero/m_solvent.rhoZero/2.0 
        !! TODO: l'extraire du fichier c directement mais cela implique de trop grosses modifs pour les permiers tests :)
        c00000 = -13.6652260 ! -440 sur le code 1D :s
        kT   = thermo%kbT
        
        do is=1,ns

          rho0 = solvent(is)%rho0
          A(is) = kT*( 1.0_dp/rho0**2 - c00000/(2.0_dp*rho0) )    ! c0 sphérique
        
        end do
        
        B = -15e-8_dp
        

        is_init = .true.
        
    end subroutine init
        

    subroutine fill_kernel
        
        use module_grid, only: grid

        implicit none

        real(dp) :: sigma, d
        integer :: i
        integer :: nx, ny, nz
        integer :: ix, iy, iz
        real(dp), allocatable :: xSqr(:), ySqr(:), zSqr(:)
        real(dp), allocatable :: kernel(:,:,:)
        real(dp) :: cutoff, cutoffDist
        type(c_ptr) :: fftPlan3d

        nx=grid%nx
        ny=grid%ny
        nz=grid%nz

        sigma = 2.36_dp
        
        cutoff = log(huge(1._dp)) ! To avoid numerical errors
        cutoffDist = sqrt(2.0_dp * log(huge(1._dp)) ) /  sigma - 0.3_dp ! minus 0.3 because of numerical error above (empirical)

        allocate( xSqr(nx) ,source= [( (real(i-1,dp)*grid%dx)**2 ,i=1,nx)] )
        allocate( ySqr(ny) ,source= [( (real(i-1,dp)*grid%dy)**2 ,i=1,ny)] )
        allocate( zSqr(nz) ,source= [( (real(i-1,dp)*grid%dz)**2, i=1,nz)] )

        allocate (kernel(grid%nx,grid%ny,grid%nz), source=0.0_dp)


        do iz=1,grid%nz
          do iy=1,grid%ny
            do ix=1,grid%nx

                d=sqrt(xSqr(ix)+ySqr(iy)+zSqr(iz))
                
                if( d .lt. cutoffDist) then
                    kernel(ix, iy, iz) = exp( -0.5*(d*sigma)**2 )
                end if

            end do !ix
          end do !iy
        end do !iz
        
       
        select case(dp)
          case(c_double)
            call dfftw_plan_dft_r2c_3d( fftPlan3d, nx, ny, nz, kernel, fourierKernel, FFTW_MEASURE )
            call dfftw_execute_dft_r2c( fftPlan3d, kernel, fourierKernel )
            call dfftw_destroy_plan( fftPlan3d )
          case(c_float)
            call sfftw_plan_dft_r2c_3d( fftPlan3d, nx, ny, nz, kernel, fourierKernel, FFTW_MEASURE )
            call sfftw_execute_dft_r2c( fftPlan3d, kernel, fourierKernel )
            call sfftw_destroy_plan( fftPlan3d )
          end select
        
        
        deallocate(xSqr)
        deallocate(ySqr)
        deallocate(zSqr)
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
        real(dp) :: dv, kT,rho0
        real(dp) :: coarseGrainedXi, coarseGrainedRho, coarseGrainedDeltaRho
        real(dp), parameter :: zerodp = 0._dp


        ns = solvent(1)%nspec
        kT = thermo%kbT
        dv = grid%dv

        fb = zerodp

        nx=grid%nx
        ny=grid%ny
        nz=grid%nz


        if(is_init .eqv. .FALSE.) then
            call init
        end if

        
        do is=1,ns

          call grid%integrate_over_orientations( solvent(is)%xi**2 * solvent(is)%rho0, density)
          call compute_coarse_grained_density( density, coarseGrainedDensity )
          rho0 = solvent(is)%rho0

          do iz=1,nz
            do iy=1,ny
              do ix=1,nx
 
                coarseGrainedXi = coarseGrainedDensity(ix,iy,iz)
                coarseGrainedRho = coarseGrainedXi**2 *rho0
          
                coarseGrainedDeltaRho = coarseGrainedRho - rho0
                fb = fb + A(is) * coarseGrainedDeltaRho**3 &
                        + B     * coarseGrainedRho**2 * coarseGrainedDeltaRho**4
                
                if(present(df)) then
                  if ( coarseGrainedRho .lt. rho0 ) then
                    coarseGrainedDfb(ix,iy,iz) = A(is) * 3.0_dp * coarseGrainedDeltaRho**2 &
                                               + B     *  ( 2.0_dp * coarseGrainedRho * coarseGrainedDeltaRho**4 &
                                                          + 4.0_dp * coarseGrainedRho**2 * coarseGrainedDeltaRho**3 )
                  else 
                    coarseGrainedDfb(ix,iy,iz) = 0.0_dp
                  end if
                end if
                

              end do !ix
            end do !iy
          end do !iz
        end do !is

        if(present(df)) then
          call update_gradient(coarseGrainedDfb, df)
        end if

    end subroutine coarse_grained_bridge
    
    
    
    
    
    
    subroutine compute_coarse_grained_density ( density,  coarseGrainedDensity )

      implicit none
 
      real(dp), intent(in) :: density(:,:,:)
      real(dp), intent(out) :: coarseGrainedDensity(:,:,:)
    
      select case(dp)
        case(c_double)
          call dfftw_execute_dft_r2c( fftRho%plan3dp, density, fourierDensity)
        case(c_float)
          call sfftw_execute_dft_r2c( fftRho%plan3dp, density, fourierDensity)
      end select
        
      fourierDensity(:,:,:) = fourierDensity(:,:,:) * fourierKernel(:,:,:)
        
      select case(dp)
        case(c_double)
          call dfftw_execute_dft_c2r( fftRho%plan3dm, fourierDensity, coarseGrainedDensity)
        case(c_float)
          call sfftw_execute_dft_c2r( fftRho%plan3dm, fourierDensity, coarseGrainedDensity)
      end select
    
    end subroutine compute_coarse_grained_density
    
    
    
    
    subroutine update_gradient ( coarseGrainedDfb, df )
        
      use module_solvent, only: solvent 
      use module_grid, only: grid

      implicit none
 
      real(dp), intent(in) :: coarseGrainedDfb(:,:,:)
      real(dp), intent(out) :: df(:,:,:,:,:)
    
    
      integer :: is, ix, iy, iz
      integer :: ns, nx, ny, nz
      real(dp) :: wTot

      ns = solvent(1)%nspec
 
      nx=grid%nx
      ny=grid%ny
      nz=grid%nz
    
      wTot = sum(grid%w(:))
    
      select case(dp)
        case(c_double)
          call dfftw_execute_dft_r2c( fftDfb%plan3dp, coarseGrainedDfb, fourierCoarseGrainedDfb)
        case(c_float)
          call sfftw_execute_dft_r2c( fftDfb%plan3dp, coarseGrainedDfb, fourierCoarseGrainedDfb)
      end select
        
      fourierCoarseGrainedDfb(:,:,:) = fourierCoarseGrainedDfb(:,:,:) * fourierKernel(:,:,:)
        
      select case(dp)
        case(c_double)
          call dfftw_execute_dft_c2r( fftDfb%plan3dm, fourierCoarseGrainedDfb, dfb)
        case(c_float)
          call sfftw_execute_dft_c2r( fftDfb%plan3dm, fourierCoarseGrainedDfb, dfb)
      end select

    
    
      do is=1,ns
        do iz=1,nz
          do iy=1,ny
            do ix=1,nx

              df(:,ix,iy,iz,is) = df(:,ix,iy,iz,is) + dfb(ix, iy, iz) * grid%w(:) / wTot
            
            end do
          end do
        end do
      end do
   
    end subroutine update_gradient
    
    
end module module_coarse_grained_bridge
