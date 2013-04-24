!> Defines the parameters of minimizer.
!! Could be CG+ or L-BFGS parameters, depending on input file.
module cg


use precision_kinds , only: i2b , dp


implicit none


! globale parameters

character(10) :: minimizer_type ! bfgs or cg

integer ( kind = i2b ) :: ncg ! nbr of variables

real ( kind = dp ) :: FF ! value of the minimum found for the functional

real ( kind = dp ) , allocatable , dimension(:) :: CG_vect ! functional. It has to be a one dimensional array

real ( kind = dp ) , allocatable , dimension(:) :: dF ! gradient of the functional

integer ( kind = i2b ) :: ITERMAX ! maximum number of iterations ! a terme le remplacer dans bfgs par minimizer_iter

real ( kind = dp ) :: epsg , epsmch , factr ! precision parameters

integer ( kind = i2b ) :: minimizer_iter



! parametres pour gradients conjugues BFGS

! integer ( kind = i2b ), parameter :: mmax=4 ! doit etre remplacé par mcg

integer ( kind = i2b ) :: mcg, iprint

integer ( kind = i2b ) , allocatable, dimension(:) :: nbd

integer ( kind = i2b ) , allocatable, dimension(:) :: iwa

integer ( kind = i2b ) , dimension (1:44) :: isave

real ( kind = dp ) :: pgtol

real ( kind = dp ) , allocatable, dimension(:) :: ll, uu

real ( kind = dp ) , dimension (1:29) :: dsave

real ( kind = dp ) , allocatable, dimension(:) :: wa

character*60 :: task, csave

logical, dimension (1:4) :: lsave


integer ( kind = i2b ) , allocatable , dimension (:,:,:,:) :: indicg




contains

  subroutine deallocate_everything_cg

    implicit none

    if ( allocated ( cg_vect ) ) deallocate ( cg_vect )

    if ( allocated ( dF ) ) deallocate ( dF )

    if ( allocated ( nbd ) ) deallocate ( nbd )

    if ( allocated ( iwa ) ) deallocate ( iwa )

    if ( allocated ( ll ) ) deallocate ( ll )

    if ( allocated ( uu ) ) deallocate ( uu )

    if ( allocated ( wa ) ) deallocate ( wa )

    if ( allocated ( indicg ) ) deallocate ( indicg )

  end subroutine deallocate_everything_cg


end module cg
