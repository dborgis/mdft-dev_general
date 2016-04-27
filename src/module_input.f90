! This module file is dedicated to reading input files for .
! It contains an input parser and functions to get some input, on demand.
! It is common to laboetie and mdft
! You're welcome to use it in another code, but please give credits for it

module module_input

  ! since this module should be shared between several codes,
  ! at least laboetie and mdft,
  ! it should have as few dependancies as possible.
  use iso_c_binding, only: dp=>c_double

  implicit none

  ! everything is private by default
  private

  character(len=12), parameter :: inputfilename="input/dft.in"

  character(len=100), public, allocatable, protected :: input_line(:) ! array containing all input lines
  logical, public, protected :: verbose

  type, public :: file_type
    integer :: unit
    character :: name
    integer :: nline = -1
  contains
    procedure :: line_count => file_line_count
  end type

  public :: n_linesInFile, deltaAbscissa,&
  input_dp, input_int, input_log, input_char, input_dp2, input_dp3, input_int2, input_int3

  type, private :: getinput_type
  contains
    procedure, nopass :: log => input_log
    procedure, nopass :: int => input_int
    procedure, nopass :: int2 => input_int2
    procedure, nopass :: int3 => input_int3
    procedure, nopass :: char => input_char
    procedure, nopass :: dp => input_dp
    procedure, nopass :: dp2 => input_dp2
    procedure, nopass :: dp3 => input_dp3
  end type
  type(getinput_type), public :: getinput

contains

  function file_line_count(file) result (line_count)
    implicit none
    class (file_type), intent(in) :: file
    integer :: line_count
    if (file%nline < 0) then
      line_count = 1
    else
      line_count = file%nline
    end if
  end function file_line_count

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function input_dp (tag, defaultvalue, assert)
    IMPLICIT NONE
    REAL(DP) :: input_dp
    CHARACTER(*), INTENT(IN) :: tag
    REAL(DP), optional, intent(in) :: defaultvalue
    INTEGER :: i, j
    logical :: tag_is_found
    character(*), optional, intent(in) :: assert
    if (.not. allocated(input_line) ) call put_input_in_character_array
    tag_is_found = .false.
    j=LEN(tag)
    DO i = 1, SIZE( input_line)
      IF( input_line( i)( 1:j) == tag  .AND. input_line(i)(j+1:j+1)==' ' ) then
        READ(input_line(i)(j+4:j+50),*) input_dp
        tag_is_found = .true.
        exit
      end if
    END DO
    if (tag_is_found .eqv. .false. .and. present(defaultvalue)) input_dp = defaultvalue
    if (present(assert)) then
      select case (assert)
      case (">0")
        if (input_dp <= 0) then
          print*, tag,"=",input_dp,". Must be >0"
          stop
        end if
      case (">=0")
        if (input_dp < 0) then
          print*, tag,"=",input_dp,". Must be >=0"
          stop
        end if
      case ("<0")
        if (input_dp >= 0) then
          print*, tag,"=",input_dp,". Must be <0"
          stop
        end if
      case ("<=0")
        if (input_dp > 0) then
          print*, tag,"=",input_dp,". Must be <=0"
          stop
        end if
      case default
        print*, 'I dont understand your assert for tag:', tag
        print*, 'defaultvalue', defaultvalue
        print*, 'assert=', assert
        stop
      end select
    end if
  END FUNCTION input_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function input_dp2 (tag, defaultvalue, assert)
    IMPLICIT NONE
    REAL(DP) :: input_dp2(2)
    CHARACTER(*), INTENT(IN) :: tag
    character(*), optional, intent(in) :: assert
    REAL(DP), optional, intent(in) :: defaultvalue(2)
    INTEGER :: i, j
    logical :: tag_is_found
    if (.not. allocated(input_line) ) call put_input_in_character_array
    tag_is_found = .false.
    j=LEN(tag)
    DO i = 1, SIZE( input_line)
      IF( input_line( i)( 1:j) == tag  .AND. input_line(i)(j+1:j+1)==' ' ) then
        READ(input_line(i)(j+4:j+50),*) input_dp2
        tag_is_found = .true.
        exit
      end if
    END DO
    if (tag_is_found .eqv. .false. .and. present(defaultvalue)) input_dp2 = defaultvalue
    if (present(assert)) then
      select case (assert)
      case (">0")
        if (any(input_dp2 <= 0)) then
          print*, tag,"=",input_dp2,". Must be >0"
          stop
        end if
      case (">=0")
        if (any(input_dp2 < 0)) then
          print*, tag,"=",input_dp2,". Must be >=0"
          stop
        end if
      case ("<0")
        if (any(input_dp2 >= 0)) then
          print*, tag,"=",input_dp2,". Must be <0"
          stop
        end if
      case ("<=0")
        if (any(input_dp2 > 0)) then
          print*, tag,"=",input_dp2,". Must be <=0"
          stop
        end if
      case default
        print*, 'I dont understand your assert for tag:', tag
        print*, 'defaultvalue', defaultvalue
        print*, 'assert=', assert
        stop
      end select
    end if
  END FUNCTION input_dp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function input_dp3 (tag, defaultvalue, assert)
    IMPLICIT NONE
    REAL(DP) :: input_dp3(3)
    character(*), optional, intent(in) :: assert
    CHARACTER(*), INTENT(IN) :: tag
    REAL(DP), optional, intent(in) :: defaultvalue(3)
    INTEGER :: i, j
    logical :: tag_is_found
    if (.not. allocated(input_line) ) call put_input_in_character_array
    tag_is_found = .false.
    j=LEN(tag)
    DO i = 1, SIZE( input_line)
      IF( input_line( i)( 1:j) == tag  .AND. input_line(i)(j+1:j+1)==' ' ) then
        READ(input_line(i)(j+4:j+50),*) input_dp3
        tag_is_found = .true.
        exit
      end if
    END DO
    if (tag_is_found .eqv. .false. .and. present(defaultvalue)) input_dp3 = defaultvalue
    if (present(assert)) then
      select case (assert)
      case (">0")
        if (any(input_dp3 <= 0)) then
          print*, tag,"=",input_dp3,". Must be >0"
          stop
        end if
      case (">=0")
        if (any(input_dp3 < 0)) then
          print*, tag,"=",input_dp3,". Must be >=0"
          stop
        end if
      case ("<0")
        if (any(input_dp3 >= 0)) then
          print*, tag,"=",input_dp3,". Must be <0"
          stop
        end if
      case ("<=0")
        if (any(input_dp3 > 0)) then
          print*, tag,"=",input_dp3,". Must be <=0"
          stop
        end if
      case default
        print*, 'I dont understand your assert for tag:', tag
        print*, 'defaultvalue', defaultvalue
        print*, 'assert=', assert
        stop
      end select
    end if
  END FUNCTION input_dp3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function input_int (tag, defaultvalue, assert)
    IMPLICIT NONE
    INTEGER :: input_int
    CHARACTER(*), INTENT(IN) :: tag
    integer, optional, intent(in) :: defaultvalue
    character(*), optional, intent(in) :: assert
    INTEGER :: i, j
    logical :: tag_is_found
    if (.not. allocated(input_line) ) call put_input_in_character_array
    tag_is_found = .false.
    j=LEN(tag)
    DO i = 1, SIZE( input_line)
      IF( input_line( i)( 1:j) == tag  .AND. input_line(i)(j+1:j+1)==' ' ) then
        READ(input_line(i)(j+4:j+50),*) input_int
        tag_is_found = .true.
        exit
      end if
    END DO
    if (.not. tag_is_found) then
      if (present(defaultvalue)) then
        input_int = defaultvalue
      else
        print*, "looking for tag", tag, "but unable to find it."
        stop
      end if
    end if

    if (present(assert)) then
      select case (assert)
      case (">0")
        if (input_int <= 0) then
          print*, tag,"=",input_int,". Must be >0"
          stop
        end if
      case (">=0")
        if (input_int < 0) then
          print*, tag,"=",input_int,". Must be >=0"
          stop
        end if
      case ("<0")
        if (input_int >= 0) then
          print*, tag,"=",input_int,". Must be <0"
          stop
        end if
      case ("<=0")
        if (input_int > 0) then
          print*, tag,"=",input_int,". Must be <=0"
          stop
        end if
      case default
        print*, 'I dont understand your assert for tag:', tag
        print*, 'defaultvalue', defaultvalue
        print*, 'assert=', assert
        stop
      end select
    end if

  end function input_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function input_int2 (tag, defaultvalue, assert)
    IMPLICIT NONE
    INTEGER :: input_int2(2)
    CHARACTER(*), INTENT(IN) :: tag
    character(*), optional, intent(in) :: assert
    integer, optional, intent(in) :: defaultvalue(2)
    INTEGER :: i, j
    logical :: tag_is_found
    if (.not. allocated(input_line) ) call put_input_in_character_array
    tag_is_found = .false.
    j=LEN(tag)
    DO i = 1, SIZE( input_line)
      IF( input_line( i)( 1:j) == tag  .AND. input_line(i)(j+1:j+1)==' ' ) then
        READ(input_line(i)(j+4:j+50),*) input_int2
        tag_is_found = .true.
        exit
      end if
    END DO
    if (tag_is_found.eqv..false. .and. present(defaultvalue)) input_int2 = defaultvalue
    if (present(assert)) then
      select case (assert)
      case (">0")
        if (any(input_int2 <= 0)) then
          print*, tag,"=",input_int2,". Must be >0"
          stop
        end if
      case (">=0")
        if (any(input_int2 < 0)) then
          print*, tag,"=",input_int2,". Must be >=0"
          stop
        end if
      case ("<0")
        if (any(input_int2 >= 0)) then
          print*, tag,"=",input_int2,". Must be <0"
          stop
        end if
      case ("<=0")
        if (any(input_int2 > 0)) then
          print*, tag,"=",input_int2,". Must be <=0"
          stop
        end if
      case default
        print*, 'I dont understand your assert for tag:', tag
        print*, 'defaultvalue', defaultvalue
        print*, 'assert=', assert
        stop
      end select
    end if

  END FUNCTION input_int2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function input_int3 (tag, defaultvalue, assert)
    IMPLICIT NONE
    INTEGER :: input_int3(3)
    character(*), optional, intent(in) :: assert
    CHARACTER(*), INTENT(IN) :: tag
    integer, optional, intent(in) :: defaultvalue(3)
    INTEGER :: i, j
    logical :: tag_is_found
    if (.not. allocated(input_line) ) call put_input_in_character_array
    tag_is_found = .false.
    j=LEN(tag)
    DO i = 1, SIZE( input_line)
      IF( input_line( i)( 1:j) == tag  .AND. input_line(i)(j+1:j+1)==' ' ) then
        READ(input_line(i)(j+4:j+50),*) input_int3
        tag_is_found = .true.
        exit
      end if
    END DO
    if (tag_is_found.eqv..false. .and. present(defaultvalue)) input_int3(1:3) = defaultvalue(1:3)
    if (present(assert)) then
      select case (assert)
      case (">0")
        if (any(input_int3 <= 0)) then
          print*, tag,"=",input_int3,". Must be >0"
          stop
        end if
      case (">=0")
        if (any(input_int3 < 0)) then
          print*, tag,"=",input_int3,". Must be >=0"
          stop
        end if
      case ("<0")
        if (any(input_int3 >= 0)) then
          print*, tag,"=",input_int3,". Must be <0"
          stop
        end if
      case ("<=0")
        if (any(input_int3 > 0)) then
          print*, tag,"=",input_int3,". Must be <=0"
          stop
        end if
      case default
        print*, 'I dont understand your assert for tag:', tag
        print*, 'defaultvalue', defaultvalue
        print*, 'assert=', assert
        stop
      end select
    end if

  END FUNCTION input_int3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION input_log (tag, defaultvalue)
    IMPLICIT NONE
    logical :: input_log
    CHARACTER(*), INTENT(IN) :: tag
    logical, intent(in), optional :: defaultvalue
    CHARACTER :: text
    INTEGER :: i, j, lentag
    logical :: found
    if (.not. allocated(input_line) ) call put_input_in_character_array
    found = .false.
    IF (tag=='point_charge_electrostatic') THEN
      STOP 'The tag point_charge_electrostatic in inputfilename must be renamed direct_sum since July 27th, 2014'
    END IF
    lentag = LEN(tag)
    DO i =1, SIZE( input_line)
      IF( input_line(i)(1:lentag)==tag .AND. input_line(i)(lentag+1:lentag+1)==' ' ) then
        READ( input_line(i)(lentag+4:lentag+50) ,*) text
        found = .true.
        exit
      end if
    END DO
    if( .not.found .and. present(defaultvalue) ) then
      input_log = defaultvalue
      return
    else if( .not.found ) then
      print*, "I could not find keyword '", tag,"' in ",inputfilename
      print*, "It should have been there associated to logical T or F"
      stop
    end if
    j = 999 ! means error in reading
    IF( text(1:1) == 'T' ) j = 1 ! means true, 2 means false
    IF( text(1:1) == 't' ) j = 1
    IF( text(1:1) == 'F' ) j = 2
    IF( text(1:1) == 'f' ) j = 2
    IF( j == 999 ) THEN
      PRINT*, 'Problem with logical tag "', tag,' in ',inputfilename
      print*, "It is here but is not logical. I read from it: '",text
      STOP
    END IF
    IF( j == 1 ) input_log = .TRUE.
    IF( j == 2 ) input_log = .FALSE.
  END FUNCTION input_log

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION input_char (tag, defaultvalue)
    IMPLICIT NONE
    CHARACTER(50) :: input_char
    CHARACTER(*), INTENT(IN) :: tag
    character(*), intent(in), optional :: defaultvalue
    INTEGER :: i,lentag,imax,iostatint
    logical :: tag_is_found
    tag_is_found = .false.
    if (.not. allocated(input_line) ) call put_input_in_character_array
    lentag=LEN(tag)
    i=0
    imax=SIZE(input_line)
    DO i=1,imax+1
      IF (i==imax+1) exit
      IF (input_line(i)(1:lentag)==tag .AND. input_line(i)(lentag+1:lentag+1)==' ') THEN
        READ(input_line(i)(lentag+4:lentag+50),*,IOSTAT=iostatint) input_char
        tag_is_found = .true.
        IF (iostatint/=0) THEN
          PRINT*,"I have a problem in reading input line:"
          PRINT*,TRIM(ADJUSTL(input_line(i)))
          IF (TRIM(ADJUSTL(input_line(i)(lentag+4:lentag+50)))=='') PRINT*,"I found nothing after sign ="
          STOP
        END IF
        EXIT
      END IF
    END DO
    if (tag_is_found) then
      IF (LEN(TRIM(ADJUSTL(input_char)))==0) THEN
        print*, "Problem while looking for input tag: ", tag
        print*, "The first input character is a whitespace"
        print*, "input_line(:) reads: ", input_line
        print*, "input_char(:) reads: ", input_char
        stop
      else IF (input_char(1:1)==' ') THEN
        print*, "Problem while looking for input tag: ", tag
        print*, "input_line(:) reads: ", input_line
        print*, "input_char(:) reads: ", input_char
        print*, "I only read white spaces!"
        stop
      end if
    else
      if (present(defaultvalue)) then
        input_char = defaultvalue
      else
        print*, "Failed to find the keyword for input tag: ", tag
        print*, "No default values have been given by developers."
        stop
      end if
    end if
  end function input_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION n_linesInFile (filename)
    IMPLICIT NONE
    INTEGER :: n_linesInFile
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER :: ios
    OPEN (77, FILE=filename)
    n_linesInFile = 0
    DO WHILE (.true.)
      READ (77,*,IOSTAT=ios)
      IF (ios>0) THEN
        WRITE(*,*)'Error in file:',filename
        STOP
      ELSE IF (ios<0) THEN ! end of file reached
        EXIT
      ELSE
        n_linesInFile = n_linesInFile +1
      END IF
    END DO
    CLOSE (77)
  END FUNCTION n_linesInFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION deltaAbscissa (filename)
    REAL(dp) :: abscissa, previousAbscissa, ordonates, deltaAbscissa
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER :: i, ios, n_lines
    OPEN (10, FILE=filename, IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE(*,*)"Cant open file ",filename," in FUNCTION deltaAbscissa"
    END IF

    DO i= 1, 2
      READ(10,*,IOSTAT=ios) previousAbscissa, ordonates
      IF (ios/=0) then
        PRINT*, 'Something went wrong while reading ', TRIM(ADJUSTL(filename)), ' in module_input>deltaAbscissa'
        STOP
      END IF
      READ(10,*, IOSTAT=ios) abscissa, ordonates
      IF (ios/=0) then
        PRINT*, 'Something went wrong while reading ', TRIM(ADJUSTL(filename)), ' in module_input>deltaAbscissa'
        STOP
      END IF
    END DO
    deltaAbscissa = abscissa-previousAbscissa
    CLOSE(10)
    n_lines = n_linesInFile(filename)
    OPEN (10, FILE=filename, IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE(*,*)"Cant open file ",filename," in FUNCTION deltaAbscissa"
    END IF

    READ(10,*,IOSTAT=ios) abscissa, ordonates
    IF (ios/=0) then
      PRINT*, 'Something went wrong while reading ', TRIM(ADJUSTL(filename)), ' in module_input>deltaAbscissa'
      STOP
    END IF
    DO i=1, n_lines-1
      previousAbscissa = abscissa
      READ(10,*,IOSTAT=ios) abscissa, ordonates
      IF (ios>0) then
        PRINT*, 'Something went wrong while reading ', TRIM(ADJUSTL(filename)), ' in module_input>deltaAbscissa'
        STOP
      ELSE IF (ios<0) THEN
        EXIT
      ELSE
        IF ((abscissa-previousAbscissa-deltaAbscissa)/deltaAbscissa > 1E-5) THEN
          PRINT*, abscissa, previousAbscissa,abscissa-previousAbscissa ,deltaAbscissa
          PRINT*, 'STOP. Non uniform absissa in ', filename
          STOP
        END IF
      END IF
    END DO
    CLOSE(10)
  end function deltaAbscissa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine put_input_in_character_array
    implicit none
    integer :: i, j, k, n, n_lines!, linelen
    character(len=100) :: text ! temporary input line
    character(len=100), allocatable, dimension(:) :: arraytemp  ! Temporary array to stock data for resizing input_line

    n_lines = n_linesInFile(inputfilename)
    allocate( input_line(n_lines), stat=i)
    if (i /= 0) then
      print *, "input_line: Allocation request denied"
      error stop
    end if

    open(unit=11, file=inputfilename, iostat=i, status="old", action="read")
    if ( i /= 0 ) then
      print*, "Error opening file ",inputfilename
      error stop
    end if

    do i=1, n_lines
      read(11,'(a)') text
      input_line (i) = trim(adjustl(text))
    end do

    close(unit=11, iostat=i)
    if ( i /= 0 ) then
      print*, "In module_input > put_input_in_character_array: Error closing file ",inputfilename
      error stop
    end if

    !  clean up comments in the lines (for instance, "option = 3 # blabla" should become "option = 3")
    DO i = 1, n_lines
      DO j = 1, len(text)
        IF ( input_line (i) (j:j) == '#' ) THEN
          DO CONCURRENT ( k=j:LEN(text) )
            input_line(i)(k:k) = ' '
          end do
          EXIT
        end if
      end do
      input_line(i) = TRIM( ADJUSTL( input_line(i) ))
    end do

    !Delete blank lines and count the size of the smallest array containing initial data
    n = 0
    do i = 1 , n_lines
      if ( input_line (i) (1:1) /= ' ' )  then
        input_line (n+1) = input_line(i)
        n = n + 1
      endif
    end do

    !Resize input_line to the smallest size by using a temporary array
    allocate (arraytemp(n))
    do i=1,n
      arraytemp(i)(:) = input_line(i)(:)
    end do
    deallocate (input_line)
    !  for some reason gfortran 5.2 is not happy with deallocate followed by allocate with source
    allocate (input_line(1:n))
    input_line(1:n) = arraytemp(1:n)
    deallocate (arraytemp)

    ! print what has been considered as input by the parser, that is what is contained by input_line(), to output dir.
    call execute_command_line("mkdir -p output", WAIT=.TRUE.)  ! just create folder. If it already exists, nothing happens.
    open(10, FILE='output/inputfile.out' )
    block
      integer :: i
      do i = 1 , SIZE( input_line ) ! print each line of input_line()
        write(10,*) input_line (i)
      end do
    end block
    close(10)

    ! check verbosity level. That's not the best place to read this, but certainly better than in module_init>allocate_from_input.
    verbose = getinput%log( "verbose", defaultvalue=.false.)

  end subroutine put_input_in_character_array

end module module_input
