subroutine get_psi_integration_roots_and_weights

    use system, only : nb_psi 
    use constants , only : twopi, pi
    use precision_kinds , only : i2b , dp
    use input , only : input_line
    use quadrature, only : x_psi , weight_psi, sym_order

    implicit none
    integer (kind = i2b ) :: psi

    allocate (weight_psi (nb_psi))
    allocate (x_psi      (nb_psi))

    do concurrent (psi=1:nb_psi) ! simple equidistant repartition between 0 and 2Pi
        weight_psi(psi) = twopi/real(nb_psi*sym_order,dp) ! weights
        x_psi(psi) = real(psi-1,dp)*twopi/real(nb_psi*sym_order,dp) ! roots
    end do

end subroutine get_psi_integration_roots_and_weights
