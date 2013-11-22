SUBROUTINE init_solvent_polarization

    USE quadrature, ONLY: Rotxx, Rotxy, Rotxz, Rotyx, Rotyy, Rotyz, Rotzx, Rotzy, Rotzz, angGrid, molRotGrid
    USE input, ONLY: input_char

    IMPLICIT NONE

    ! call get_charge_factor ( Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz )
    ! If we microscopic description of solvent compute charge density in k space and electrostatic potential generated by a such distribution
    IF (trim(adjustl(input_char('polarization_order'))) == 'multi') then
        CALL chargeDensityAndMolecularPolarizationOfASolventMoleculeAtOrigin (Rotxx,Rotxy,Rotxz,Rotyx,Rotyy,Rotyz,Rotzx,Rotzy,Rotzz)
    END IF
    
END SUBROUTINE init_solvent_polarization