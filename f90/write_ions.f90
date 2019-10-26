program write_ions

  use module_domain
  use module_ramses
  use module_mesh
  use module_gas_composition
  use module_select

  implicit none

  real(kind=8) :: xmin,xmax,ymin,ymax,zmin,zmax
  integer(kind=4),dimension(:),allocatable :: cpu_list
  integer(kind=4) :: ncpu_read
  character(2000),allocatable :: ion_data_path(:)
  character(2000),allocatable :: ions(:)
  character(2000) :: parameter_file
  real(kind=8) :: start, finish
  integer :: narg

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [CreateDomDump] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)             :: repository = './'          ! ramses run directory (where all output_xxxxx dirs are).
  character(2000)             :: ion_parameter_file = './ions.dat' 
  !integer(kind=4)             :: ion_number = 1             ! number of ions we want to read
  integer(kind=4)             :: snapnum = 1                ! ramses output number to use
  character(2000)             :: write_path = './cells'
  integer(kind=4)             :: max_cells = 1000000
  ! --- domain  
  real(kind=8),dimension(3)   :: comput_dom_pos       = (/0.5,0.5,0.5/)  ! center of domain [code units]
  real(kind=8)                :: comput_dom_rsp       = 0.3              ! radius of sphere [code units]

  ! --- miscelaneous
  logical                     :: verbose = .false.
  logical                     :: read_excited_state = .false.
  ! --------------------------------------------------------------------------

  call cpu_time(start)

  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: CreateDomDump params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_write_ions_params(parameter_file)
  if (verbose) call print_write_ions_params
  ! ------------------------------------------------------------


  xmax = comput_dom_pos(1) + comput_dom_rsp
  xmin = comput_dom_pos(1) - comput_dom_rsp
  ymax = comput_dom_pos(2) + comput_dom_rsp
  ymin = comput_dom_pos(2) - comput_dom_rsp
  zmax = comput_dom_pos(3) + comput_dom_rsp
  zmin = comput_dom_pos(3) - comput_dom_rsp

  call get_cpu_list_periodic(repository, snapnum, xmin,xmax,ymin,ymax,zmin,zmax, ncpu_read, cpu_list)
  call write_ion_leaf(repository, snapnum, ion_number, ion_data_path, ions, max_cells, write_path, ncpu_read, cpu_list)

  call cpu_time(finish)
  print '(" --> Time = ",f12.3," seconds.")',finish-start
  print*,' '

contains

  subroutine read_write_ions_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! ALSO read parameter form used modules (mesh)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present
    logical         :: ndomain_present 

    section_present = .false.
    ndomain_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:12) == '[write_ions]') then
          section_present = .true.
          exit
       end if
    end do
    ! read section if present
    if (section_present) then 
       do
          read (10,'(a)',iostat=err) line
          if(err/=0) exit
          if (line(1:1) == '[') exit ! next section starting... -> leave
          i = scan(line,'=')
          if (i==0 .or. line(1:1)=='#' .or. line(1:1)=='!') cycle  ! skip blank or commented lines
          name=trim(adjustl(line(:i-1)))
          value=trim(adjustl(line(i+1:)))
          i = scan(value,'!')
          if (i /= 0) value = trim(adjustl(value(:i-1)))
          select case (trim(name))
          case ('comput_dom_pos')
             read(value,*) comput_dom_pos(1),comput_dom_pos(2),comput_dom_pos(3)
          case ('comput_dom_rsp')
             read(value,*) comput_dom_rsp
          case ('verbose')
             read(value,*) verbose
          case ('repository')
             write(repository,'(a)') trim(value)
          case ('ion_parameter_file')
             write(ion_parameter_file,'(a)') trim(value)
          case ('write_path')
             write(write_path,'(a)') trim(value)
          case ('snapnum')
             read(value,*) snapnum
          case ('max_cells')
             read(value,*) max_cells
          case ('read_excited_state')
             read(value,*) read_excited_state
          end select
       end do
    end if
    close(10)

    call read_ramses_params(pfile)
    call read_gas_composition_params(pfile)
    call read_ion_params_file
    
    return

  end subroutine read_write_ions_params


  subroutine read_ion_params_file

    implicit none

    integer(kind=4) :: i

    allocate(ions(ion_number),ion_data_path(ion_number))

    open(unit=10,file=trim(ion_parameter_file),action='read',form='formatted')
    do i=1,ion_number
       read(10,'(a)') ions(i)
       read(10,'(a)') ion_data_path(i)
    end do
    
  end subroutine read_ion_params_file


  subroutine print_write_ions_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit
    character(100) :: fmt

    if (present(unit)) then 
       write(unit,'(a,a,a)')         '[write_ions]'
       write(unit,'(a)')             '# input / output parameters'
       write(unit,'(a,a)')           '  repository = ',trim(repository)
       write(unit,'(a,i5)')          '  snapnum    = ',snapnum
       write(unit,'(a,a)')           '  write_path = ',trim(write_path)
       write(unit,'(a,i5)')          '  max_cells  = ',max_cells
       write(unit,'(a,3(ES10.3,1x))')'  comput_dom_pos       = ',comput_dom_pos(1),comput_dom_pos(2),comput_dom_pos(3)
       write(unit,'(a,ES10.3,1x)')   '  comput_dom_rsp       = ',comput_dom_rsp
       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,L1)')          '  verbose         = ',verbose
       call print_ramses_params(unit)
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ' '
       write(*,'(a,a,a)')         '[write_ions]'
       write(*,'(a)')             '# input / output parameters'
       write(*,'(a,a)')           '  repository = ',trim(repository)
       write(*,'(a,i5)')          '  snapnum    = ',snapnum
       write(*,'(a,a)')           '  write_path = ',trim(write_path)
       write(*,'(a,i5)')          '  max_cells  = ',max_cells
       write(*,'(a,3(ES10.3,1x))')'  comput_dom_pos       = ',comput_dom_pos(1),comput_dom_pos(2),comput_dom_pos(3)
       write(*,'(a,ES10.3,1x)')   '  comput_dom_rsp       = ',comput_dom_rsp
       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,L1)')          '  verbose         = ',verbose
       call print_ramses_params
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
    end if

    return

  end subroutine print_write_ions_params



end program write_ions
