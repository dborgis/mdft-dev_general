subroutine adhoc_corrections_to_gsolv
! ... Here, we print all the adhoc corrections one should take into account before comparing MDFT results to MD and/or experiments.
  
    use precision_kinds, only: dp, sp
    use system, only: solute
  
    implicit none
  
    real(dp) :: PMcorr ! We use P-scheme instead of M-scheme for the electrostatics in MDFT. See Kastenholz and Hunenberger, JCP 124, 124106 (2006)
  
    if (  abs(sum(solute%site%q)) >= epsilon(1.0_dp)  ) then
        PMcorr = -79.8_dp*sum(solute%site%q) ! in kJ/mol
        print*,"You should add",real(PMcorr,sp),"kJ/mol to FREE ENERGY because of the P-scheme electrostatics"
    end if

end subroutine adhoc_corrections_to_gsolv
