! variable relative to the system being studied
module system
    use precision_kinds , only: i2b,dp
    implicit none
    private

    type :: site_type
        character(100) :: name
        integer(i2b) :: type, n_sites
        real(dp), dimension(3) :: r
        real(dp) :: q, sig, eps, lambda1, lambda2
        integer(i2b) :: z ! atomic number
    end type site_type

    public :: site_type, mole_fraction
    !
    ! complex(dp),allocatable, dimension(:,:,:) :: vk !>@var perturabtion in kspace
    !
    ! logical(4) :: hs
    ! complex(dp), allocatable , dimension (:,:,:) :: v_perturbation_k ! fourier transform of the lennard jones perturbation (wca)
    real(dp), allocatable , dimension (:) :: mole_fraction ! mole fraction of each species "x_i"

end module system
