!> read solute atomic positions, charge, and lennard jones values in solute.in
!! charge in electron units, sigma in Angstroms, epsilon in KJ/mol.
subroutine read_solute

  use precision_kinds, ONLY: i2b,dp
  use system,          ONLY: nb_solute_sites, solute, nb_species, spacegrid
  use input,           ONLY: input_line, input_log, input_dp
!  use periodic_table,  ONLY: init_periodic_table, ptable

  implicit none

  integer :: n,i,stat

!  call init_periodic_table
  ! print *, ptable ( 1 ) % name
  ! open and test if input/solute.in is ok

  open (5, FILE='input/solute.in', STATUS='old', IOSTAT=stat)
  IF (stat /= 0) THEN
    PRINT*,'solute.in cannot be opened ! => STOP !'
    STOP
  END IF

  READ (5,*) ! comment line
  READ (5,*) nb_solute_sites ! total number of atom sites of the solute
  ALLOCATE(solute%site(nb_solute_sites))
  READ (5,*)
  DO n = 1, nb_solute_sites
    READ(5,*) i, solute%site(n)%q, solute%site(n)%sig, solute%site(n)%eps, &
              solute%site(n)%lambda1, solute%site(n)%lambda2, solute%site(n)%r, solute%site(n)%Z
    END DO
  CLOSE (5)
  solute%site%q = solute%site%q * input_dp('solute_charges_scale_factor', defaultvalue=1._dp)

  block
    real(dp) :: solutexmin, solutexmax, soluteymin, soluteymax, solutezmin, solutezmax, &
                solutediameterx, solutediametery, solutediameterz, solutesigmaljmax
    solutexmin = minval(solute%site%r(1))
    soluteymin = minval(solute%site%r(2))
    solutezmin = minval(solute%site%r(3))
    solutexmax = maxval(solute%site%r(1))
    soluteymax = maxval(solute%site%r(2))
    solutezmax = maxval(solute%site%r(3))
    solutediameterx = sqrt((solutexmax-solutexmin)**2)
    solutediametery = sqrt((soluteymax-soluteymin)**2)
    solutediameterz = sqrt((solutezmax-solutezmin)**2)
!    solutesigmaljmax = maxval(solute%site%sig)
    spacegrid%buffer_length = input_dp("buffer_length", defaultvalue=15._dp)
    spacegrid%length(1) = solutediameterx +2.*spacegrid%buffer_length !+2*solutesigmaljmax 
    spacegrid%length(2) = solutediametery +2.*spacegrid%buffer_length !+2*solutesigmaljmax 
    spacegrid%length(3) = solutediameterz +2.*spacegrid%buffer_length !+2*solutesigmaljmax 
  end block

end subroutine read_solute
