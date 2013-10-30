! This subroutine computes the external potential. It is one of the most time consuming routine.
!Warning there are two ways of calculating the electrostatic potential (Poisson solver and point charge) you should always have one tag
!T and one tag F for the electrostatic potential, if thera are 2 tags T, it is the last evaluation which counts, i.e Poisson solver. 

subroutine init_external_potential

    use precision_kinds , only : dp , i2b
    use input, only: input_line, input_log, input_char
    use system , only: chg_mol, chg_solv, eps_solv, eps_mol, sig_solv, sig_mol, Lx, Ly, Lz, nb_psi, &
                    nfft1, nfft2, nfft3, nb_species
    use external_potential, only: Vext_total, Vext_q
    use mod_lj, only: initLJ => init
    use quadrature, only: Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz, angGrid
    
    implicit none
    
    integer(i2b):: nb_id_mol , nb_id_solv ! nb of types of sites of solute and solvent
    integer(i2b):: i, j

    if(.not. allocated(Vext_total)) allocate( Vext_total(nfft1,nfft2,nfft3,angGrid%n_angles,nb_psi,nb_species), source=0._dp )

    nb_id_mol  = size ( chg_mol  ) ! total number of solute types
    nb_id_solv = size ( chg_solv ) ! total number of solvent types

    !call get_charge_factor ( Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz )
    ! If we microscopic description of solvent compute charge density in k space and electrostatic potential generated by a such distribution
    if (trim(adjustl(input_char('evaluate_polarization'))) == 'multi') then
        call get_charge_density_k ( Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz )  
    end if

    ! Pseudopotential
    !allocate(Vpseudo(nfft1,nfft2,nfft3,angGrid%n_angles))
    !call compute_vpseudo_ijko

    ! Hard walls
    call external_potential_hard_walls

    ! Charges : treatment as point charges
    if (input_log('point_charge_electrostatic')) then
        ! Compute Vcoul(i,j,k,omega)
        if (.not. allocated(Vext_q)) allocate( Vext_q(nfft1,nfft2,nfft3,angGrid%n_angles,nb_psi,nb_species), source=0._dp )
        call compute_vcoul_as_sum_of_pointcharges( Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz )
    end if

    ! Charges : Poisson solver
    if (input_log('poisson_solver')) then
        if (.not. allocated(Vext_q) ) allocate ( Vext_q ( nfft1 , nfft2 , nfft3 , angGrid%n_angles, nb_psi , nb_species ) )
        Vext_q=0.0_dp
        call electrostatic_potential_from_charge_density ! Electrostatic potential using FFT of Poisson equation Laplacian(V(r))= - charge_density / Epsilon_0
        call vext_q_from_v_c (Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)
    end if

    ! Lennard-Jones
    call initLJ

    ! r^-12 only
    if (input_log('purely_repulsive_solute')) then
        call compute_purely_repulsive_potential ( Rotxx , Rotxy , Rotxz , Rotyx , Rotyy , Rotyz , Rotzx , Rotzy , Rotzz )
    end if

    ! hard spherical solute
    if (input_log('hard_sphere_solute')) then
        call compute_vext_hard_sphere
    end if
    
    ! potential created by hard wall(s) ! CALLED TWICE BUG
    call external_potential_hard_walls

    ! hard cylinder
    if (input_log('hard_cylinder_solute')) then
        call compute_vext_hard_cylinder
    end if
    
    ! hard square well
    !do i = 1 , size ( input_line )
    !  j = len ( 'vext_hard_square_well' )
    !  if ( input_line (i) (1:j) == 'vext_hard_square_well' .and. input_line (i) (j+4:j+4) == 'T' ) call vext_hard_square_well
    !end do

    ! personnal vext as implemented in personnal_vext.f90
    if (input_log('personnal_vext')) then
        call compute_vext_perso
    end if
    
    ! compute total Vext(i,j,k,omega), the one used in the free energy functional
    call vext_total_sum
    
end subroutine init_external_potential
