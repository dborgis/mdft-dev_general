SUBROUTINE get_final_polarization ( Px , Py , Pz )

    USE precision_kinds,    ONLY: dp,i2b
    USE system,             ONLY: nb_species, spacegrid
    USE constants,          ONLY: fourpi,twopi
    USE minimizer,          ONLY: cg_vect_new
    USE input,              ONLY: input_int

    IMPLICIT NONE
    INTEGER(i2b) :: i,j,k,io,icg,s,p, molRotSymOrder
    REAL(dp) :: rho_toto,local_Px,local_Py,local_Pz
    REAL(dp), DIMENSION(spacegrid%nx,spacegrid%ny,spacegrid%nz,nb_species), INTENT(OUT) :: Px,Py,Pz ! equilibrium polarization(r)
    INTEGER(i2b) :: nx, ny, nz, no, ns
    real(dp), parameter :: zerodp = 0._dp
    nx= spacegrid%n_nodes(1)
    ny= spacegrid%n_nodes(2)
    nz= spacegrid%n_nodes(3)
    no = spacegrid%no
    ns = nb_species


    Px = zerodp
    Py = zerodp
    Pz = zerodp

    icg = 0
    DO s =1,ns
        DO i =1,nx
            DO j =1,ny
                DO k =1,nz
                    local_Px = 0.0_dp
                    local_Py = 0.0_dp
                    local_Pz = 0.0_dp
                    DO io =1,spacegrid%no
                        icg = icg+1
                        rho_toto = cg_vect_new(i,j,k,io,s)**2 /(twopi*fourpi/spacegrid%molRotSymOrder)
                        local_Px = local_Px + spacegrid%omx(io) * spacegrid%w(io) * rho_toto
                        local_Py = local_Py + spacegrid%omy(io) * spacegrid%w(io) * rho_toto
                        local_Pz = local_Pz + spacegrid%omz(io) * spacegrid%w(io) * rho_toto
                    END DO
                    Px(i,j,k,s) = local_Px
                    Py(i,j,k,s) = local_Py
                    Pz(i,j,k,s) = local_Pz
                END DO
            END DO
        END DO
    END DO

END SUBROUTINE get_final_polarization
