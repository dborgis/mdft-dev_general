! Subroutine to initiate simulation parameters
SUBROUTINE init_simu

    use quadrature  ,ONLY: prepare_quadrature => init
    use fft         ,ONLY: prepare_corefft_for_MDFT => init, prepare_fft_threads => init_threads
    use dcf         ,ONLY: init_dcf => init
    use module_input       ,ONLY: getinput
    use module_minimizer   ,only: init_minimizer => init
    use hardspheres ,only: compute_hard_spheres_parameters
    ! use grid_mod    ,only: grid

    IMPLICIT NONE

!    type(grid) :: gr
!    call grid%build
!    print*, "dpsi and dphi =",grid%dpsi, grid%dphi
!    print*, "mmax=",grid%mmax
!    print*, "theta=",grid%theta
!    print*, "weigth=",grid%wtheta
!    stop "ici et la"

    CALL print_header ! package name, day & time
    call read_solute ! AWFUL TODO trick to default lx
    CALL allocate_from_input ! TODO this should be removed and done only when needed ! read dft.in, solute.in and solvent.in
    CALL prepare_quadrature


    CALL prepare_fft_threads ! prepare le multithread de fftw3
    CALL prepare_corefft_for_MDFT

    CALL print_input_to_output_folder! Print input parameters found in input files to output folder
    CALL read_solvent ! Read solvent atomic positions, charge and Lennard-Jones param
    CALL init_dcf
    IF (getinput%log('hard_sphere_fluid', defaultvalue=.false.)) CALL compute_hard_spheres_parameters ! If calculation based on Fundamental Measure Theory read FMT parameters and compute weight functions etc
    CALL init_minimizer ! Prepare minimization ! allocate tables and computes precision and so on
    CALL init_external_potential
    CALL init_density

END SUBROUTINE init_simu
