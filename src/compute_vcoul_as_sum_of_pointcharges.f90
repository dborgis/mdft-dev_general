! Here we compute the coulomb potential using the simple sum over all i, j of q_i*q_j/rij. i and j run over all charges in supercell and 26 neighbor supercells.
! See issue in Github for a proposition to improve this.
! Note that NO TABULATION IS USED. This name is here for historical reasons and should be modified at some point.

subroutine compute_vcoul_as_sum_of_pointcharges( Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz )

    use precision_kinds, only: dp,i2b
    use system, only: nfft1,nfft2,nfft3,deltax,deltay,deltaz,id_solv,id_mol,x_solv,y_solv,z_solv,x_mol,y_mol,z_mol,&
                        beta,nb_solute_sites,nb_solvent_sites,chg_mol,chg_solv,Lx,Ly,Lz , nb_species, RC
    use constants , only : fourpi , qfact
    use external_potential , only : Vext_q
    use quadrature, only: angGrid, molRotGrid
    USE input, ONLY: verbose

    implicit none

    real(dp), dimension(angGrid%n_angles,molRotGrid%n_angles), intent(in) :: Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz
    integer(i2b) :: i,j,k,o,p,m,n ! dummy
    real(dp),dimension(nb_solvent_sites,molRotGrid%n_angles,angGrid%n_angles) :: xmod,ymod,zmod
    real(dp) :: x_grid,y_grid,z_grid ! coordinates of grid nodes
    real(dp) :: x_m,y_m,z_m ! solvent sites coordinates
    real(dp) :: x_nm,y_nm,z_nm ! coordinate of vecteur solute-solvent
    real(dp) :: r_nm2 ! norm**2 of vector x_nm;y_nm;z_nm
    real(dp) :: V_psi ! dummy
    integer(i2b) :: ids,idm ! id of the solvent site or the solute site that is being used
    real(dp) :: time0,time1 ! timer start and end
    real(dp), dimension(:,:,:), allocatable :: temparray
    character(50):: filename
    integer(i2b) :: px,py,pz
    real(dp) :: qfactcc ! qfact*chg_solv()*chg_mol()
    real(dp) :: rc2 ! charge pseudo radius **2   == Rc**2
    real(dp) :: tempVcoul ! temporary Vcoul(i,j,k,o)
    integer(i2b):: species ! dummy for loops over species

    print*,'!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!'
    print*,'This routine use Rc which is defined in dft.in know what you are doing'
    write(*,*)'Rc = ',Rc 
    print*,'!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!!!WARNING!!!!!'

    call cpu_time(time0)! initiate
    Vext_q = 0.0_dp
    px=0
    py=0
    pz=0
    Rc2=Rc**2

    ! test if all solutes have zero charge then don't waste your time : go to end of subroutine
    IF (MINVAL(chg_mol)==0.0_dp .and. MAXVAL(chg_mol)==0.0_dp) THEN   ! solute is not charged
        IF (verbose) WRITE(*,*)'solute has no charge'
        RETURN
    ELSE IF (MINVAL(chg_solv)==0.0_dp .and. MAXVAL(chg_solv)==0.0_dp) THEN   ! solvent is not charged
        IF (verbose) WRITE(*,*)'solvent has no charge'
        RETURN
    END IF
    
    ! precompute Rot_ij(omega,psi)*k_solv(a) for speeding up
    do o=1,angGrid%n_angles
        do p=1,molRotGrid%n_angles
            do m=1,nb_solvent_sites
                xmod(m,p,o)= Rotxx(o,p)*x_solv(m) + Rotxy(o,p)*y_solv(m) + Rotxz(o,p)*z_solv(m)
                ymod(m,p,o)= Rotyx(o,p)*x_solv(m) + Rotyy(o,p)*y_solv(m) + Rotyz(o,p)*z_solv(m)
                zmod(m,p,o)= Rotzx(o,p)*x_solv(m) + Rotzy(o,p)*y_solv(m) + Rotzz(o,p)*z_solv(m)
            end do
        end do
    end do

    do species = 1 , nb_species
        do k=1,nfft3
            z_grid=real(k-1,dp)*deltaz
            !  if(z_grid>12.0_dp .and. z_grid<30.0_dp) then
            !    Vcoul(:,:,k,:)=100.0_dp
            !    cycle
            !  end if
            do j=1,nfft2
                y_grid=real(j-1,dp)*deltay
                do i=1,nfft1
                    x_grid=real(i-1,dp)*deltax
                    do o=1,angGrid%n_angles
                        tempVcoul=0.0_dp
                        ploop:  do p=1,molRotGrid%n_angles
                            V_psi=0.0_dp
                            do m=1,nb_solvent_sites
                                ids=id_solv(m)
                                if (chg_solv(ids)==0.0_dp) cycle
                                x_m = x_grid + xmod(m,p,o)
                                y_m = y_grid + ymod(m,p,o)
                                z_m = z_grid + zmod(m,p,o)
                                do n=1,nb_solute_sites
                                    idm=id_mol(n)
                                    if (chg_mol(idm)==0.0_dp) cycle
                                    qfactcc=qfact*chg_solv(ids)*chg_mol(idm)
                                    !   write(*,*) qfactcc, qfact
                                    !       consider every period in every direction => 3*3*3=27 solutes
                                    do px=-1,1
                                        x_nm= x_m-(x_mol(n)+px*Lx)
                                        do py=-1,1
                                            y_nm= y_m-(y_mol(n)+py*Ly)
                                            do pz=-1,1
                                                z_nm= z_m-(z_mol(n)+pz*Lz)
                                                r_nm2 = x_nm**2+y_nm**2+z_nm**2
                                                if(r_nm2<Rc2) then
                                                    V_psi = huge(1.0_dp)
                                                    cycle ploop ! sous-entendu tempVcoul=tempVcoul+exp(-beta*huge)=tempVcoul+0
                                                else
                                                    V_psi = V_psi + qfactcc/sqrt(r_nm2)
                                                end if
                                            end do
                                        end do
                                    end do
                                end do ! solute
                            end do ! solvent
                        !tempVcoul=tempVcoul+exp(-beta*V_psi)
                        !if (tempVcoul<1.e-10_dp) then
                        Vext_q(i,j,k,o,p,species ) = V_psi
                        !else
                        !   Vext_q(i,j,k,o,1)=-log(tempVcoul/molRotGrid%n_angles)/beta
                        !end if
                        ! Limit maximum value
                        if (Vext_q(i,j,k,o,p,species)>100.0_dp) Vext_q(i,j,k,o,p,species)=100.0_dp
                        !!!!!        if (Vcoul(i,j,k,o)<=-10.0_dp) then
                        !!!!!          write(*,*)'problem in compute_vcoul_ijko_from_tabulated.f90'
                        !!!!!          write(*,*)'reduced coordinate of problematic grid point : ',x_grid/Lx,y_grid/Ly,z_grid/Lz
                        !!!!!          write(*,*)'Vcoul(i,j,k,o) = ',Vcoul(i,j,k,o)
                        !!!!!        else if (Vcoul(i,j,k,o)>100.0_dp) then
                        !!!!!          Vcoul(i,j,k,o)=100.0_dp
                        !!!!!        end if
                        end do ploop ! psi
                    end do ! omega
                end do ! i
            end do ! j
        end do ! k
    end do ! species
!$OMP END PARALLEL DO
    
    999 continue
    
    call cpu_time(time1)
    
    IF (verbose) THEN
        write(*,*) 'minval(Vcoul)= ' , minval( Vext_q (:,:,:,:,:,:) )
        write(*,*) 'maxval(Vcoul)= ' , maxval( Vext_q (:,:,:,:,:,:) )
        write(*,*) 'runtime for compute_vcoul_ijko_from_tabulated ' , time1 - time0
        ! Get the external potential over orientations and print it
        allocate ( temparray ( nfft1 , nfft2 , nfft3 ) )
        call mean_over_orientations ( Vext_q (:,:,:,:,:,1) , temparray )
        temparray = temparray / fourpi
        filename = 'output/Vcoul.cube'
        call write_to_cube_file ( temparray , filename )
        filename = 'output/Vcoul_along-z.dat'
        call compute_z_density ( temparray , filename )
        deallocate(temparray)
    END IF

END SUBROUTINE
