! Here we compute the coulomb potential using the simple sum over all i, j of q_i*q_j/rij. i and j run over all charges in supercell and 26 neighbor supercells.
! See issue in Github for a proposition to improve this.
! Note that NO TABULATION IS USED. This name is here for historical reasons and should be modified at some point.

SUBROUTINE compute_vcoul_as_sum_of_pointcharges

    USE precision_kinds, only: dp,i2b
    use system, only: nfft1,nfft2,nfft3,deltax,deltay,deltaz,id_solv,id_mol,x_mol,y_mol,z_mol,&
                        beta,nb_solute_sites,nb_solvent_sites,chg_mol,chg_solv,Lx,Ly,Lz , nb_species, solventSite, soluteSite
    use constants,only : fourpi , qfact
    use external_potential,only : Vext_q
    use quadrature, only: angGrid, molRotGrid, Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz
    USE input, ONLY: verbose

    IMPLICIT NONE

    integer(i2b) :: i,j,k,o,p,m,n,s
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
    real(dp) :: qfactcc ! qfact*chg_solv()*chg_mol()
    real(dp) :: tempVcoul ! temporary Vcoul(i,j,k,o)
    REAL, PARAMETER :: Rc=1.0_dp, Rc2=Rc**2 ! charge pseudo radius **2   == Rc**2


    call cpu_time(time0)! initiate
    Vext_q = 0.0_dp

    ! test if all solutes have zero charge then don't waste your time : go to end of SUBROUTINE
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
                xmod(m,p,o)= Rotxx(o,p)*solventSite(m)%r(1) + Rotxy(o,p)*solventSite(m)%r(2) + Rotxz(o,p)*solventSite(m)%r(3)
                ymod(m,p,o)= Rotyx(o,p)*solventSite(m)%r(1) + Rotyy(o,p)*solventSite(m)%r(2) + Rotyz(o,p)*solventSite(m)%r(3)
                zmod(m,p,o)= Rotzx(o,p)*solventSite(m)%r(1) + Rotzy(o,p)*solventSite(m)%r(2) + Rotzz(o,p)*solventSite(m)%r(3)
            END DO
        END DO
    END DO

    do s = 1 , nb_species
        do k=1,nfft3
            z_grid=real(k-1,dp)*deltaz
            do j=1,nfft2
                y_grid=real(j-1,dp)*deltay
                do i=1,nfft1
                    x_grid=real(i-1,dp)*deltax
                    do o=1,angGrid%n_angles
                        tempVcoul=0.0_dp
                        ploop:  do p=1,molRotGrid%n_angles
                            V_psi=0.0_dp
                            do m=1,SIZE(solventSite)
                                if (solventSite(m)%q==0.0_dp) cycle
                                x_m = x_grid + xmod(m,p,o)
                                y_m = y_grid + ymod(m,p,o)
                                z_m = z_grid + zmod(m,p,o)
                                do n=1,SIZE(soluteSite)
                                    if (soluteSite(n)%q==0.0_dp) cycle
                                    qfactcc=qfact*solventSite(m)%q*soluteSite(n)%q
                                    x_nm = x_m - x_mol(n)
                                    y_nm = y_m - y_mol(n)
                                    z_nm = z_m - z_mol(n)
                                    r_nm2 = x_nm**2+y_nm**2+z_nm**2
                                    if(r_nm2<Rc2) then
                                        V_psi = huge(1.0_dp)
                                        cycle ploop ! sous-entendu tempVcoul=tempVcoul+exp(-beta*huge)=tempVcoul+0
                                    ELSE
                                        V_psi = V_psi + qfactcc/sqrt(r_nm2)
                                    END IF
                                END DO ! solute
                            END DO ! solvent
                        !tempVcoul=tempVcoul+exp(-beta*V_psi)
                        !if (tempVcoul<1.e-10_dp) then
                        Vext_q(i,j,k,o,p,s ) = V_psi
                        !ELSE
                        !   Vext_q(i,j,k,o,1)=-log(tempVcoul/molRotGrid%n_angles)/beta
                        !END IF
                        ! Limit maximum value
                        if (Vext_q(i,j,k,o,p,s)>100.0_dp) Vext_q(i,j,k,o,p,s)=100.0_dp
                        !!!!!        if (Vcoul(i,j,k,o)<=-10.0_dp) then
                        !!!!!          write(*,*)'problem in compute_vcoul_ijko_from_tabulated.f90'
                        !!!!!          write(*,*)'reduced coordinate of problematic grid point : ',x_grid/Lx,y_grid/Ly,z_grid/Lz
                        !!!!!          write(*,*)'Vcoul(i,j,k,o) = ',Vcoul(i,j,k,o)
                        !!!!!        ELSE IF (Vcoul(i,j,k,o)>100.0_dp) then
                        !!!!!          Vcoul(i,j,k,o)=100.0_dp
                        !!!!!        END IF
                        END DO ploop ! psi
                    END DO ! omega
                END DO ! i
            END DO ! j
        END DO ! k
    END DO ! species
    
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
