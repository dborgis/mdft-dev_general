module module_debug
    implicit none
    private
    logical, protected :: debugmode
    public :: debugmode, init_debug
contains
    subroutine init_debug
        use module_input, only: getinput
        implicit none
        debugmode = getinput%log( "debugmode", defaultvalue=.true.)
    end subroutine
end module
