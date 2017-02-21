program RaysFromSourceModel

  use module_ray
  use module_utils
  use module_random
  use module_constants

  implicit none

  type(ray_type),dimension(:),allocatable :: rays
  integer                                 :: iran, i, narg
  real(kind=8)                            :: nu, scalar
  character(2000)                         :: parameter_file

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [RaysFromSourceModel] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: outputfile = 'RaysIC.dat'   ! file to which outputs will be written
  ! --- emission properties  
  real(kind=8),dimension(3) :: source_pos  = (/0.5d0,0.5d0,0.5d0/)   ! position of the source [code units]
  real(kind=8),dimension(3) :: source_vel  = (/0.0d0,0.0d0,0.0d0/)    ! velocity of the source [km/s]
  real(kind=8)              :: nu_min      = clight/1221d-8          ! min frequency (in source's frame) [Hz]
  real(kind=8)              :: nu_max      = clight/1210d-8          ! max frequency (in source's frame) [Hz]
  ! --- miscelaneous
  integer                   :: nrays       = 10000                   ! number of rays to generate
  integer                   :: ranseed = 1234                        ! seed for random generator
  logical                   :: verbose = .true.
  ! --------------------------------------------------------------------------


  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: RaysFromSourceModel path/to/params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_RaysFromSourceModel_params(parameter_file)
  if (verbose) call print_RaysFromSourceModel_params
  ! ------------------------------------------------------------


  ! --------------------------------------------------------------------------------------
  ! define rays
  ! --------------------------------------------------------------------------------------
  allocate(rays(nrays))
  iran = -abs(ranseed)
  do i = 1,nrays
     rays(i)%ID     = i
     nu             = ran3(iran) * (nu_max-nu_min) + nu_min 
     call isotropic_direction(rays(i)%k_em,iran)
     scalar         = rays(i)%k_em(1)*source_vel(1) + rays(i)%k_em(2)*source_vel(2) + rays(i)%k_em(3)*source_vel(3)
     rays(i)%nu_ext = nu / (1d0 - scalar/clight)
     rays(i)%x_em   = source_pos
  enddo


  ! --------------------------------------------------------------------------------------
  ! write ICs
  ! --------------------------------------------------------------------------------------
  if (verbose) write(*,*) '--> writing file'
  open(unit=14, file=trim(outputfile), status='unknown', form='unformatted', action='write')
  write(14) nrays
  write(14) (rays(i)%ID,i=1,nrays)
  write(14) (rays(i)%nu_ext,i=1,nrays)
  write(14) (rays(i)%x_em(:),i=1,nrays)
  write(14) (rays(i)%k_em(:),i=1,nrays)
  close(14)
  ! --------------------------------------------------------------------------------------

  deallocate(rays)

contains

  subroutine read_RaysFromSourceModel_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present
    
    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:24) == '[RaysFromSourceModel]') then
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
          case ('outputfile')
             write(outputfile,'(a)') trim(value)
          case ('source_pos')
             read(value,*) source_pos(1),source_pos(2),source_pos(3)
          case ('source_vel')
             read(value,*) source_vel(1),source_vel(2),source_vel(3)
          case ('verbose')
             read(value,*) verbose
          case ('ranseed')
             read(value,*) ranseed
          case ('nu_min')
             read(value,*) nu_min
          case ('nu_max')
             read(value,*) nu_max
          case ('nrays')
             read(value,*) nrays
          end select
       end do
    end if
    close(10)
    return

  end subroutine read_RaysFromSourceModel_params


  subroutine print_RaysFromSourceModel_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)')         '[RaysFromSourceModel]'
       write(unit,'(a)')             '# input / output parameters'
       write(unit,'(a,a)')           '  outputfile      = ',trim(outputfile)
       write(unit,'(a)')             '# emission properties'
       write(unit,'(a,3(ES9.3,1x))') '  source_pos      = ',source_pos(1),source_pos(2),source_pos(3)
       write(unit,'(a,3(ES9.3,1x))') '  source_vel      = ',source_vel(1),source_vel(2),source_vel(3)
       write(unit,'(a,es9.3,a)')     '  nu_min          = ',nu_min, ' ! [Hz]'
       write(unit,'(a,es9.3,a)')     '  nu_max          = ',nu_max, ' ! [Hz]'
       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,i8)')          '  nrays           = ',nrays
       write(unit,'(a,i8)')          '  ranseed         = ',ranseed
       write(unit,'(a,L1)')          '  verbose         = ',verbose
       write(unit,'(a)')             ' '
    else
       write(*,'(a,a,a)')         '[RaysFromSourceModel]'
       write(*,'(a)')             '# input / output parameters'
       write(*,'(a,a)')           '  outputfile    = ',trim(outputfile)
       write(*,'(a)')             '# emission properties'
       write(*,'(a,3(ES9.3,1x))') '  source_pos      = ',source_pos(1),source_pos(2),source_pos(3)
       write(*,'(a,3(ES9.3,1x))') '  source_vel      = ',source_vel(1),source_vel(2),source_vel(3)
       write(*,'(a,es9.3,a)')     '  nu_min          = ',nu_min, ' ! [Hz]'
       write(*,'(a,es9.3,a)')     '  nu_max          = ',nu_max, ' ! [Hz]'
       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,i8)')          '  nrays           = ',nrays
       write(*,'(a,i8)')          '  ranseed         = ',ranseed
       write(*,'(a,L1)')          '  verbose         = ',verbose
       write(*,'(a)')             ' '
    end if

    return

  end subroutine print_RaysFromSourceModel_params


end program RaysFromSourceModel
