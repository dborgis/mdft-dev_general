! In this subroutine one calls the different parts of the total energy
! first is computed the radial_part, then ... blabla.
! this subroutine is the one called by the minimization stuff
! for computing the total energy and associated gradient

subroutine energy_and_gradient

    use precision_kinds , only : i2b , dp
    use input , only : input_line,input_log,input_char
    use cg , only : FF , dF , cg_vect
    use fft , only : in_forward , in_backward , out_forward , out_backward , plan_forward , plan_backward
    use system, only : nb_omega, nb_psi, sigma_k, nfft1, nfft2, nfft3, Lx, Ly, Lz, deltax, rho_c, deltaV, nb_species,&
                    rho_0_multispec , rho_c_k_myway
    use constants, only : twopi, fourpi, qfact
    use quadrature, only : weight, weight_psi, sym_order

    implicit none
    integer(i2b):: i , j , k , l, nf1, nf2, nf3
    logical :: hydrophobicity_tag


    ! FF is the TOTAL ENERGY of the system, it is thus the functional of the density that is minimized by solver
    ! dF is the gradient of FF with respect to all coordinates. Remember it is of the kind dF ( number of variables over density (ie angles etc))
    call init_functional_and_gradient_to_zero( FF, dF )

    if (.not. allocated (rho_c_k_myway)) allocate (rho_c_k_myway(nfft1/2+1, nfft2, nfft3))
    rho_c_k_myway=(0.0_dp, 0.0_dp)

    call energy_external ! compute ideal part of the total energy
    call energy_ideal

    ! compute radial part of the excess free energy
    if (input_log('hard_sphere_fluid') ) call energy_hard_sphere_fmt ! => pure hard sphere contribution

    ! test if there is a LJ perturbation to hard spheres ! WCA model etc. to implement more intelligently
    if (input_log('lennard_jones_perturbation_to_hard_spheres') ) call lennard_jones_perturbation_to_hard_spheres ! lennard jones perturbative contribution => Weeks-Chandler-Anderson

    if (input_log('read_ck_or_chi')) then
        if (input_log('hydrophobicity')) then
            if (trim(adjustl(input_char('treatment_of_hydro')))=='C')  then
                call cs_plus_hydro
            else if (trim(adjustl(input_char('treatment_of_hydro')))=='VdW')  then
                call energy_hydro
            end if
        else
            call cs_from_dcf
        end if
    end if

    ! bridge calculation: F(FMT)-F(c2hs)+F(c2H2O)
    if (input_log('bridge_hard_sphere')) call energy_cs_hard_sphere ! better name should be given

    ! Dipolar polarization. What user wants (use it or not) is checked in subroutine for clearer code.
    call energy_polarization
    call energy_polarization_myway

    ! Threebody term that is needed to empiricaly force H-bonding in water. What user wants (use it or not) is checked in subroutine for clearer code
    call energy_threebody
    call energy_threebody_faster

    write(*,*)'SOLVATION FREE ENERGY AT THIS STEP = ',FF
    write(*,*)'-----------------------'

    contains
    
    pure subroutine init_functional_and_gradient_to_zero (FF,dF)
        real(dp), intent(out) :: FF, dF(:)
        FF = 0.0_dp ; dF = 0.0_dp ! functional and gradient
    end subroutine init_functional_and_gradient_to_zero
    
end subroutine energy_and_gradient
