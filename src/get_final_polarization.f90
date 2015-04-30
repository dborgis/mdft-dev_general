SUBROUTINE get_final_polarization ( Px , Py , Pz )

    USE precision_kinds,    ONLY: dp,i2b
    USE system,             ONLY: nb_species, spaceGrid
    USE constants,          ONLY: fourpi,twopi
    USE minimizer,          ONLY: CG_vect
    USE quadrature,         ONLY: Omx,Omy,Omz,angGrid,molRotGrid
    USE input,              ONLY: input_int
    
    IMPLICIT NONE
    INTEGER(i2b) :: i,j,k,o,icg,s,p,molRotSymOrder
    REAL(dp) :: rho_toto,local_Px,local_Py,local_Pz
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: weight_omx,weight_omy,weight_omz ! weight(:)*omx(:) for speeding up
    REAL(dp), DIMENSION(spacegrid%n_nodes(1),spacegrid%n_nodes(2),spacegrid%n_nodes(3),nb_species), INTENT(OUT) :: Px,Py,Pz ! equilibrium polarization(r)
    INTEGER(i2b) :: nfft1, nfft2, nfft3
    nfft1= spaceGrid%n_nodes(1)
    nfft2= spaceGrid%n_nodes(2)
    nfft3= spaceGrid%n_nodes(3)


    molRotSymOrder = input_int('molRotSymOrder', defaultvalue=1)
    Px = 0.0_dp
    Py = 0.0_dp
    Pz = 0.0_dp
    ALLOCATE( weight_omx ( angGrid%n_angles ), SOURCE=(angGrid%weight*Omx) )
    ALLOCATE( weight_omy ( angGrid%n_angles ), SOURCE=(angGrid%weight*Omy) )
    ALLOCATE( weight_omz ( angGrid%n_angles ), SOURCE=(angGrid%weight*Omz) )

    icg = 0
    DO s =1,nb_species
        DO i =1,nfft1
            DO j =1,nfft2
                DO k =1,nfft3
                    local_Px = 0.0_dp
                    local_Py = 0.0_dp
                    local_Pz = 0.0_dp
                    DO o =1,angGrid%n_angles
                        DO p =1,molRotGrid%n_angles
                            icg = icg+1
                            rho_toto = cg_vect(icg)**2 /(twopi*fourpi/molRotSymOrder)
                            local_Px = local_Px + weight_Omx(o) * rho_toto * molRotGrid%weight(p)
                            local_Py = local_Py + weight_Omy(o) * rho_toto * molRotGrid%weight(p)
                            local_Pz = local_Pz + weight_Omz(o) * rho_toto * molRotGrid%weight(p)
                        END DO
                    END DO
                    Px(i,j,k,s) = local_Px
                    Py(i,j,k,s) = local_Py
                    Pz(i,j,k,s) = local_Pz
                END DO
            END DO
        END DO
    END DO
    
END SUBROUTINE get_final_polarization
