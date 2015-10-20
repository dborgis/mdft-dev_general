SUBROUTINE get_final_polarization ( Px , Py , Pz )

    use precision_kinds,    ONLY: dp,i2b
    use module_solvent, only: solvent
    use constants,          ONLY: fourpi,twopi
    use module_minimizer,          ONLY: cg_vect_new
    use module_input,              ONLY: getinput
    use module_grid, only: grid

    IMPLICIT NONE
    INTEGER(i2b) :: i,j,k,io,icg,s,p, molRotSymOrder
    REAL(dp) :: rho_toto,local_Px,local_Py,local_Pz
    REAL(dp), DIMENSION(grid%nx,grid%ny,grid%nz,solvent(1)%nspec), INTENT(OUT) :: Px,Py,Pz ! equilibrium polarization(r)
    INTEGER(i2b) :: nx, ny, nz, no, ns
    real(dp), parameter :: zerodp = 0._dp
    nx= grid%n_nodes(1)
    ny= grid%n_nodes(2)
    nz= grid%n_nodes(3)
    no = grid%no
    ns = solvent(1)%nspec


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
                    DO io =1,grid%no
                        icg = icg+1
                        rho_toto = cg_vect_new(i,j,k,io,s)**2 /(twopi*fourpi/grid%molRotSymOrder)
                        local_Px = local_Px + grid%omx(io) * grid%w(io) * rho_toto
                        local_Py = local_Py + grid%omy(io) * grid%w(io) * rho_toto
                        local_Pz = local_Pz + grid%omz(io) * grid%w(io) * rho_toto
                    END DO
                    Px(i,j,k,s) = local_Px
                    Py(i,j,k,s) = local_Py
                    Pz(i,j,k,s) = local_Pz
                END DO
            END DO
        END DO
    END DO

END SUBROUTINE get_final_polarization
