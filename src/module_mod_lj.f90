! In module mod_lj, we compute the array Vext_lj that contains the lennard jones part of the external potential
! It is perhaps not the best idea to have a module for that, but is simpler for me (Max) to code, at least at the beginning.
module mod_lj
    use external_potential, only: Vext_lj
    use precision_kinds, only: dp, i2b
    use system, only: nfft1, nfft2, nfft3, nb_omega,deltax,deltay,deltaz,nb_solute_sites,eps_mol,x_mol,y_mol,z_mol,Lx,Ly,Lz&
                   ,sig_solv,sig_mol , nb_species, id_mol, chg_mol , chg_solv, eps_solv, nb_psi, nb_omega
    use constants, only:fourpi
    use input , only : input_line
    implicit none
    integer(i2b), private, parameter :: nrgrid = 100000 ! nb of point of the radial grid for tabulation of the potentials TODO magic number here ...
    integer(i2b), private :: nb_id_mol, nb_id_solv ! number of different kinds of solvent sites or solute sites
    real(dp), private :: drgrid, drgrid2 ! width in cartesian units of a radial bin

    contains
    
        subroutine init
            nb_id_mol  = size ( chg_mol  ) ! total number of solute types
            nb_id_solv = size ( chg_solv ) ! total number of solvent types
            drgrid  =   3._dp* max ( Lx , Ly , Lz ) * sqrt ( 3.0_dp ) * 2.0_dp / real ( nrgrid , kind=dp) ! distance between two grid points for tabulation
            drgrid2 = ( 3._dp * max ( Lx , Ly , Lz ) * sqrt ( 3.0_dp ) ) ** 2 * 2.0_dp / real( nrgrid ,kind=dp) ! distance **2 between two grid points
            if (.not. allocated( Vext_lj )) allocate( Vext_lj(nfft1,nfft2,nfft3,nb_omega,nb_psi,nb_species), source=0._dp ) ! Vext_lj is a *very* big array that should now be allocated
            call calculate ! compute Vext_lj
            ! deallocate the tabulated values of Vext_lj(r) and Vext_lj(r**2)
            !if ( allocated ( tabulated_ljsq ) ) deallocate ( tabulated_ljsq )
        end subroutine
        
        subroutine calculate
        
            real(dp), allocatable, dimension(:,:,:) :: tabulated_ljsq
            integer(i2b) :: i,j,k,n ! dummy
            real(dp) :: x_grid,y_grid,z_grid ! coordinates of grid nodes
            real(dp) :: V_node ! sum of all solute contributions to a given grid node
            integer(i2b) :: idm, i_r ! id of the solute site that is being used
            real(dp) :: r_nm2 ! norm**2 of vector x_nm;y_nm;z_nm
            real(dp) :: time0,time1 ! timer start and end
            real(dp), dimension(:,:,:), allocatable :: temparray
            character(50):: filename
            integer ( kind = i2b ) :: species ! dummy for loops over species
            integer ( kind = i2b ) :: pot_type
            real ( kind = dp ) :: eps_Lorentz_Berthelot , sigma_Lorentz_Berthelot 
            integer ( kind = i2b ) ::  i_mol,i_solv
            real ( kind = dp ) , dimension ( 0 : nrgrid ) :: pot !tabulated potential this subroutines computes
            real ( kind = dp ) :: drgrid, dx, dy, dz ! distance between two points in radial grid (in Angstroms) =abs(rcut-rmin)/nrgrid

            ! compute lennard jones potential at each position and for each orientation, for each species => Vext_lj ( i , j , k , omega , species ) 
            ! we impose the simplification that only the first site of the solvent sites has a lennard jones WATER ONLY TODO
            ! test if this simplification is true and stop if not
            ! the test is done over the sigma lj. they're positive, so that the sum over all sigma is the same as the sigma of first site
            ! only if they're all zero but for first site
            if( sum( sig_solv(:) ) /= sig_solv(1) ) then
                print*,'problem in vext_generator.f90'
                print*,'compute_vlj_ijko_from_tabulated IS ONLY VALID FOR SOLUTES WITH 1 LJ SITE'
                print*,'stop'
                stop
            end if

            ! Lennard Jones = 4*sqrt(eps1*eps2)*{[.5*(sig1+sig2)/r]^12 - [.5*(sig1+sig2)/r]^6}
            pot_type = 1 ! lennard jones
            allocate ( tabulated_ljsq ( nb_id_solv , nb_id_mol , 0 : nrgrid ) )
            do i_solv = 1 , nb_id_solv
                ! if epsilon is 0 then any potential involving this solvent site is zero
                if ( eps_solv ( i_solv ) == 0.0_dp ) then
                    tabulated_ljsq ( i_solv , : , : ) = 0.0_dp
                    cycle
                end if
                do i_mol = 1 , nb_id_mol
                    ! if epsilon solute site is zero then any potential involving this solvent site is zero
                    if ( eps_mol ( i_mol ) == 0.0_dp ) then
                        tabulated_ljsq ( i_solv , i_mol , : ) = 0.0_dp
                        cycle
                    end if
                    ! calculate epsilon_ij and sigma_ij using Lorentz-Berthelot mixing rules
                    eps_Lorentz_Berthelot = sqrt ( eps_solv ( i_solv ) * eps_mol ( i_mol ) ) !epsilonlj
                    sigma_Lorentz_Berthelot = 0.5_dp * ( sig_solv ( i_solv ) + sig_mol ( i_mol ) ) !sigmalj
                    ! tabulate Vij(r)
                    call potgensq ( pot_type , eps_Lorentz_Berthelot , sigma_Lorentz_Berthelot , nrgrid , drgrid2 , pot )
                    tabulated_ljsq ( i_solv , i_mol , : ) = pot (:)
                end do
            end do

            !> initiate
            call cpu_time(time0)

            !> Test if the supercell is big enough considering the LJ range (given by sigma).
            ! at 2.5*sigma, the lj potential is only 0.0163*epsilon. Almost zero.
            ! It would have no sense to have a box dimension < 2.5 sigma
            if( min(Lx,Ly,Lz)<= 2.5_dp*max(maxval(sig_mol),maxval(sig_solv))) then
                print*,'max(maxval(sig_mol),maxval(sig_solv)) = ' , max(maxval(sig_mol),maxval(sig_solv))
                print*,'problem detected in compute_vlj_ijko_from_tabulated.f90'
                print*,'the supercell is too small to use the minimum image convention with such large sigma values'
                stop
            end if

            do species = 1 , nb_species
                do k=1,nfft3
                    z_grid=real(k-1,dp)*deltaz
                    do j=1,nfft2
                        y_grid=real(j-1,dp)*deltay
                        do i=1,nfft1
                            x_grid=real(i-1,dp)*deltax
                            
                            V_node=0.0_dp
                            do n=1,nb_solute_sites
                                idm=id_mol(n)
                                if (eps_mol(idm)==0.0_dp) cycle
                                dx= x_grid-x_mol(n)
                                if ( dx > 0.5_dp*Lx) dx = dx - Lx
                                if ( dx < -0.5_dp*Lx) dx= Lx + dx
                                dy= y_grid-y_mol(n)
                                if ( dy > 0.5_dp*Ly) dy= dy - Ly
                                if ( dy < -0.5_dp*Ly) dy= Ly + dy
                                dz= z_grid-z_mol(n)
                                if ( dz > 0.5_dp*Lz) dz= dz - Lz
                                if ( dz < -0.5_dp*Lz) dz= Lz + dz
                                !r_nm2=x_nm2+y_nm2+z_nm2
                                r_nm2 = dx**2+dy**2+dz**2
                                i_r = int(r_nm2/drgrid2+0.5_dp)
                                ! Care that it is only valid for a solvent which has a lennard jones only on the first site
                                V_node = V_node + tabulated_ljsq(1,idm,i_r) ! ids=1 ONLY VALID IF ONLY ONE LJ SITE IN SOLVENT. OK FOR WATER AND STOCKMAYER

                                ! limit maximum value of Vlj to 100
                                if (V_node >= 100.0_dp) then
                                    V_node = 100.0_dp
                                    exit
                                end if
                            end do ! solute

                            ! all omegas are treated in the same time as the oxygen atom is not sensitive to rotation around omega and psi
                            Vext_lj ( i , j , k , : , : , species ) = V_node

                        end do ! i
                    end do ! j
                end do ! k
            end do ! species

            !> Get the lennard jones potential over orientations and print it
            allocate ( temparray ( nfft1 , nfft2 , nfft3 ) )
            temparray=Vext_lj(:,:,:,1,1,1)
            filename='output/Vlj.cube'
            call write_to_cube_file(temparray,filename)
            filename='output/Vlj_along-z.dat'
            call compute_z_density(temparray,filename)
            open(11,file='Vlj_aumillieu.dat')
            do i=1,nfft3
                write(11,*) i*deltaz, Vext_lj(nfft1/2,nfft2/2,i,1,1,1)
            end do
            close(11)
            deallocate(temparray)
            call cpu_time(time1)
            print *, 'Vext_lj : min = ' , minval(Vext_lj) , ' ; max = ' , maxval(Vext_lj) , ' ; in (sec) ' , time1-time0

        end subroutine
   
end module
