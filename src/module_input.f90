module input

  use precision_kinds ,only: i2b,dp
  implicit none
  character(len=100), allocatable, dimension(:) :: input_line ! array containing all input lines
  logical :: verbose
  private
  public :: verbose, input_line, input_dp, input_int, input_log, input_char, n_linesinfile, deltaabscissa

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function input_dp (that)
    implicit none
    real(dp) :: input_dp
    character(*), intent(in) :: that
    integer(i2b) :: i, j
    j=len(that)
    do i =1,size(input_line)
      if( input_line( i)( 1:j) == that  .and. input_line(i)(j+1:j+1)==' ' ) read(input_line(i)(j+4:j+50),*) input_dp
    end do
  end function input_dp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function input_int (that, def)
    implicit none
    integer(i2b) :: input_int
    character(*), intent(in) :: that
    integer(i2b), optional, intent(in) :: def
    integer(i2b) :: i, j
    logical :: ifoundtag
    ifoundtag = .false.
    j=len(that)
    do i = 1, size( input_line)
      if( input_line( i)( 1:j) == that  .and. input_line(i)(j+1:j+1)==' ' ) then
        read(input_line(i)(j+4:j+50),*) input_int
        ifoundtag = .true.
        exit
      end if
    end do
    if (ifoundtag.eqv..false. .and. present(def)) input_int = def
  end function input_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function input_log (that)
    implicit none
    logical :: input_log
    character(*), intent(in) :: that
    character :: text
    integer(i2b) :: i, j
    if (that=='point_charge_electrostatic') then
      stop 'the tag point_charge_electrostatic in dft.in must be renamed direct_sum since july 27th, 2014'
    end if
    j=len(that)
    do i =1,size( input_line)
      if( input_line(i)(1:j)==that .and. input_line(i)(j+1:j+1)==' ' ) read( input_line (i) (j+4:j+50) , * ) text
    end do
    j = 999 ! means error in reading
    if( text(1:1) == 't' ) j = 1 ! means true, 2 means false
    if( text(1:1) == 't' ) j = 1
    if( text(1:1) == 'f' ) j = 2
    if( text(1:1) == 'f' ) j = 2
    if( j == 999 ) then
      print*, 'i did not find the tag ', that,' in dft.in'
      stop
    end if
    if( j == 1 ) input_log = .true.
    if( j == 2 ) input_log = .false.
  end function input_log

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function input_char (that)
    implicit none
    character(50) :: input_char
    character(*), intent(in) :: that
    integer(i2b) :: i,j,imax,iostatint
    j=len(that)
    i=0
    imax=size(input_line)
    do i=1,imax+1
      if (i==imax+1) then
        print*,"i didnt find keyword '",that,"' in dft.in"
        stop
      end if
      if (input_line(i)(1:j)==that .and. input_line(i)(j+1:j+1)==' ') then
        read(input_line(i)(j+4:j+50),*,iostat=iostatint) input_char
        if (iostatint/=0) then
          print*,"i have a problem in reading input line:"
          print*,trim(adjustl(input_line(i)))
          if (trim(adjustl(input_line(i)(j+4:j+50)))=='') print*,"i found nothing after sign ="
          stop
        end if
        exit
      end if
    end do
    if (input_char(1:1)==' ') then
      print*,"first character of ",that," is a whitespace"
      stop
    end if
    if (len(trim(adjustl(input_char)))==0) then
      print*,"tag after ",that," is only whitespaces."
      stop
    end if
  end function input_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function n_linesinfile (filename)
    implicit none
    integer :: n_linesinfile
    character(*), intent(in) :: filename
    integer :: ios
    open (77, file=filename)
    n_linesinfile = 0
    do while (.true.)
      read (77,*,iostat=ios)
      if (ios>0) then
        write(*,*)'error in file:',filename
        stop
        else if (ios<0) then ! end of file reached
          exit
        else
          n_linesinfile = n_linesinfile +1
        end if
      end do
      close (77)
  end function

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function deltaabscissa (filename)
    implicit none
    real(dp) :: abscissa, previousabscissa, ordonates, deltaabscissa
    character(*), intent(in) :: filename
    integer(i2b) :: i, ios, n_lines
    open (10, file=filename, iostat=ios)
    if (ios /= 0) then
      write(*,*)"cant open file ",filename," in function deltaabscissa"
    end if

    do i= 1, 2
      read(10,*,iostat=ios) previousabscissa, ordonates
      if (ios/=0) then
        print*, 'something went wrong while reading ', trim(adjustl(filename)), ' in module_input>deltaabscissa'
        stop
      end if
      read(10,*, iostat=ios) abscissa, ordonates
      if (ios/=0) then
        print*, 'something went wrong while reading ', trim(adjustl(filename)), ' in module_input>deltaabscissa'
        stop
      end if
    end do
    deltaabscissa = abscissa-previousabscissa
    close(10)
    n_lines = n_linesinfile(filename)
    open (10, file=filename, iostat=ios)
    if (ios /= 0) then
      write(*,*)"cant open file ",filename," in function deltaabscissa"
    end if

    read(10,*,iostat=ios) abscissa, ordonates
    if (ios/=0) then
      print*, 'something went wrong while reading ', trim(adjustl(filename)), ' in module_input>deltaabscissa'
      stop
    end if
    do i=1, n_lines-1
      previousabscissa = abscissa
      read(10,*,iostat=ios) abscissa, ordonates
      if (ios>0) then
        print*, 'something went wrong while reading ', trim(adjustl(filename)), ' in module_input>deltaabscissa'
        stop
        else if (ios<0) then
          exit
        else
          if ((abscissa-previousabscissa-deltaabscissa)/deltaabscissa > 1e-5) then
            print*, abscissa, previousabscissa,abscissa-previousabscissa ,deltaabscissa
            print*, 'stop. non uniform absissa in ', filename
            stop
          end if
        end if
      end do
      close(10)
    end function deltaabscissa

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module input
