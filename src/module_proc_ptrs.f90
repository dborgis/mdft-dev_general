module proc_ptrs

 abstract interface
   
    subroutine excess_energy_t(ff, df, print_timers)
        use omp_lib
        use precision_kinds, only: dp
        use module_grid, only: grid
        use module_thermo, only: thermo
        use module_orientation_projection_transform, only: angl2proj, proj2angl, init_module_orientation_projection_transform => init
        use module_wigner_d, only: wigner_small_d
        implicit none

        real(dp), intent(out) :: ff
        real(dp), contiguous, intent(inout), optional :: df(:,:,:,:,:) ! x y z o s
        logical, intent(in), optional :: print_timers
    end subroutine

  end interface

  procedure(excess_energy_t), pointer :: excess_energy

end module
