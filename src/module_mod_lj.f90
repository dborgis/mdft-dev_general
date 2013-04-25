! In module mod_lj, we compute the array Vext_lj that contains the lennard jones part of the external potential
! It is perhaps not the best idea to have a module for that, but is simpler for me (Max) to code, at least at the beginning.
module mod_lj
    use external_potential, only: Vext_lj
    use precision_kinds, only: dp, i2b
    use system, only: nfft1, nfft2, nfft3, nb_omega,deltax,deltay,deltaz,nb_solute_sites,eps_mol,x_mol,y_mol,z_mol,Lx,Ly,Lz&
                   ,sig_solv,sig_mol , nb_species, id_mol, id_solv, chg_mol , chg_solv, eps_solv, nb_psi, nb_omega,&
                   x_solv,y_solv,z_solv, nb_solvent_sites
    use constants, only:fourpi
    use input , only : input_line
    implicit none
    integer(i2b), private :: nb_id_mol, nb_id_solv ! number of different kinds of solvent sites or solute sites
    contains
    
        subroutine init
            nb_id_mol  = size ( chg_mol  ) ! total number of solute types
            nb_id_solv = size ( chg_solv ) ! total number of solvent types
            if (.not. allocated( Vext_lj )) allocate( Vext_lj(nfft1,nfft2,nfft3,nb_omega,nb_psi,nb_species), source=0._dp ) ! Vext_lj is a *very* big array that should now be allocated
            call calculate ! compute Vext_lj
        end subroutine
        
        subroutine calculate
        
            integer(i2b) :: i,j,k,n ! dummy
            real(dp) :: x_grid,y_grid,z_grid ! coordinates of grid nodes
            real(dp) :: V_node ! sum of all solute contributions to a given grid node
            integer(i2b) :: idm, ids ! id of the solute and solvent site that is being used
            real(dp) :: time0, time1 ! timer start and end
            integer(i2b) :: species ! dummy for loops over species
            integer(i2b) ::  i_mol, i_solv
            real(dp):: dx, dy, dz ! distance between two points in radial grid (in Angstroms) =abs(rcut-rmin)/nrgrid
            integer(i2b) :: solute_site, solvent_site
            ! compute lennard jones potential at each position and for each orientation, for each species => Vext_lj ( i , j , k , omega , species ) 
            ! we impose the simplification that only the first site of the solvent sites has a lennard jones WATER ONLY TODO
            ! test if this simplification is true and stop if not
            ! the test is done over the sigma lj. they're positive, so that the sum over all sigma is the same as the sigma of first site
            ! only if they're all zero but for first site
            if( sum( sig_solv(:) ) /= sig_solv(1) ) then
                print*,'problem in module mod_lj'
                print*,'ONLY VALID FOR SOLUTES WITH 1 LJ SITE'
                print*,'stop'
                stop
            end if
            
            ! Also, for now, the solvent site that wear the LJ potential should be on a grid node, and more precisely have coordinates 0 0 0.
            block
                real(dp), dimension(3) :: coo
                coo = [ x_solv(1), y_solv(1), z_solv(1) ]
                if( any(coo/=0._dp) ) then
                    print*,'For now, the solvent can only have one Lennard-Jones site.'
                    print*,'This site should have coordinates 0. 0. 0. in order to be on a grid node.'
                    print*,'These coordinates are now defined as',coo
                    print*,'For this reason, MDFT stops now.'
                    stop
                end if
            end block
            !> initiate
            call cpu_time(time0)
            !> Test if the supercell is big enough considering the LJ range (given by sigma).
            ! at 2.5*sigma, the lj potential is only 0.0163*epsilon. Almost zero.
            ! It would have no sense to have a box dimension < 2.5 sigma
            if( min(Lx,Ly,Lz)<= 2.5_dp*max(maxval(sig_mol),maxval(sig_solv))) then
                print*,'max(maxval(sig_mol),maxval(sig_solv)) = ' , max(maxval(sig_mol),maxval(sig_solv))
                print*,'problem detected in the module dedicated to computing the Lennard Jones contribution to external potential'
                print*,'the supercell is too small to use the minimum image convention with such large sigma values'
                print*,'MDFT stops now.'
                stop
            end if
            do species = 1, nb_species
                do solvent_site = 1, nb_solvent_sites
                    ids = id_solv(solvent_site)
                    if( eps_solv(ids) == 0._dp ) cycle ! if solvent site wear no LJ
                    do k=1,nfft3
                        z_grid=real(k-1,dp)*deltaz
                        do j=1,nfft2
                            y_grid=real(j-1,dp)*deltay
                            do i=1,nfft1
                                x_grid=real(i-1,dp)*deltax
                                V_node=0.0_dp
                                do solute_site = 1, nb_solute_sites
                                    idm=id_mol(solute_site)
                                    if (eps_mol(idm)==0.0_dp) cycle
                                    dx= x_grid-x_mol(solute_site)
                                    if ( dx > 0.5_dp*Lx) dx = dx - Lx ! should be replaced by mod() or modulo(). cant' remember the difference between them in Fortran. Look for on internet.
                                    if ( dx < -0.5_dp*Lx) dx= Lx + dx ! TODO : look at http://stackoverflow.com/questions/15069838/need-help-optimizing-code-minimum-image-convention
                                    dy= y_grid-y_mol(solute_site)
                                    if ( dy > 0.5_dp*Ly) dy= dy - Ly
                                    if ( dy < -0.5_dp*Ly) dy= Ly + dy
                                    dz= z_grid-z_mol(solute_site)
                                    if ( dz > 0.5_dp*Lz) dz= dz - Lz
                                    if ( dz < -0.5_dp*Lz) dz= Lz + dz
                                    V_node = V_node + vlj( geometric_mean( eps_solv(ids), eps_mol(idm) ),&
                                                           arithmetic_mean( sig_solv(ids), sig_mol(idm) ),&
                                                           norm2([dx,dy,dz]) ) ! rules of Lorentz-Berthelot
                                    if (V_node >= 100.0_dp) then ! limit maximum value of Vlj to 100
                                        V_node = 100.0_dp ! TODO magic number
                                        exit
                                    end if
                                end do ! solute
    
                                ! all omegas are treated in the same time as the oxygen atom is not sensitive to rotation around omega and psi
                                Vext_lj ( i , j , k , : , : , species ) = V_node
    
                            end do ! i
                        end do ! j
                    end do ! k
                end do ! solvent sites
            end do ! species
            
            block
                real(dp), dimension(:,:,:), allocatable :: temparray
                character(50):: filename
                !> Get the lennard jones potential over orientations and print it
                allocate ( temparray ( nfft1 , nfft2 , nfft3 ) )
                temparray=Vext_lj(:,:,:,1,1,1)
                filename='output/Vlj.cube'
                call write_to_cube_file(temparray,filename)
                filename='output/Vlj_along-z.dat'
                call compute_z_density(temparray,filename)
                open(11,file='Vlj_aumilieu.dat')
                do i=1,nfft3
                    write(11,*) i*deltaz, Vext_lj(nfft1/2,nfft2/2,i,1,1,1)
                end do
                close(11)
                deallocate(temparray)
            end block
            
            call cpu_time(time1)
            print *, 'Vext_lj : min = ' , minval(Vext_lj) , ' ; max = ' , maxval(Vext_lj) , ' ; in (sec) ' , time1-time0
        end subroutine
        
        pure function vlj(eps,sig,d)
            ! v_lj(d) = 4ε[(σ/d)^12-(σ/d)^6]
            real(dp) :: vlj
            real(dp), intent(in) :: eps, sig ! ε,σ
            real(dp), intent(in) :: d ! distance
            real(dp) :: div
            div = (sig/d)**6
            vlj = 4._dp*eps*div*(div-1._dp)
        end function
        pure function arithmetic_mean( A, B)
            ! = sum_i^N a_i/N
            real(dp) :: arithmetic_mean
            real(dp), intent(in) :: A, B
            arithmetic_mean = (A+B)/2._dp
        end function
        pure function geometric_mean( A, B)
            ! = (product_i^N a_i)^(1/N)
            real(dp) :: geometric_mean
            real(dp), intent(in) :: A, B
            geometric_mean = sqrt(A*B)
        end function
   
end module
