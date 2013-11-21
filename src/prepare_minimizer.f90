! this SUBROUTINE gets the informations in input file and then allocate, prepare, compute needed data

SUBROUTINE prepare_minimizer
    
    USE precision_kinds, ONLY: dp
    USE input, ONLY: input_int, input_dp,input_char
    USE cg, ONLY: nbd , iwa , ll , uu , wa , dF , cg_vect , epsmch , factr , epsg , ncg , mcg , iprint , minimizer_type ,&
                    itermax , minimizer_iter , pgtol, FF
    USE system , ONLY: spaceGrid, soluteSite
    USE quadrature, ONLY: angGrid, molRotGrid
    
    IMPLICIT NONE
    
    ncg = PRODUCT(spaceGrid%n_nodes) * angGrid%n_angles *molRotGrid%n_angles* SIZE(soluteSite)! total number of variables to optimize
    itermax=input_int('maximum_iteration_nbr')
    epsg=input_dp('epsg')
    pgtol=input_dp('pgtol')
    minimizer_type = TRIM(ADJUSTL(input_char('minimizer'))) ! get from input the type of minimization one wants to do
    minimizer_iter = 0
    ALLOCATE ( cg_vect(ncg), SOURCE=0._dp )
    ALLOCATE ( dF(ncg), SOURCE=0._dp )
    FF = 0._dp
    
    IF ( minimizer_type(1:4) == 'bfgs') THEN ! if minimizer is bfgs then init what has to be init
        mcg = 4
        ALLOCATE ( ll(ncg), SOURCE=0._dp ) ! lower bound of cg_vect
        ALLOCATE ( uu(ncg), SOURCE=0._dp ) ! uppder bound of cg_vect
        ALLOCATE ( nbd(ncg), SOURCE=0)!nbd(i) is : 0 for no constaint, 1 for lower bounded by ll, 3 for upper bounded by uu, 2 for lower and upper bounded
        epsmch = EPSILON(1.0_dp)  !  Precision de la machine
        factr = epsg / epsmch ! convergence criteria over energy
        iprint = 1
        ALLOCATE ( iwa ( 3 * ncg ) )
        ALLOCATE ( wa ( 2 * mcg * ncg + 4 * ncg + 11 * mcg **2 + 8 * mcg ) )
    END IF

END SUBROUTINE prepare_minimizer
