program PhotonsFromSourceModel

  use module_photon
  use module_utils
  use module_random
  use module_constants

  implicit none

  type(photon_init),dimension(:),allocatable :: photgrid
  integer                                    :: iran, i, narg
  real(kind=8)                               :: nu, r1, r2
  character(2000)                            :: parameter_file

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [PhotonsFromSourceModel] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: outputfile = 'ICs_photons_n5e6.dat'   ! file to which outputs will be written

  ! --- source type 
  character(10)             :: source_type = 'pointlike'             ! type of source model
  real(kind=8),dimension(3) :: source_pos  = (/0.5d0,0.5d0,0.5d0/)   ! position of the source [code units]
  integer                   :: nphot       = 5000000                 ! number of photons to generate

  ! --- how source shines
  character(30)             :: spec_type = 'monochromatic'           ! how to draw frequencies
  ! ------ spec_type == 'monochromatic' : all photons at the same frequency nu_0
  real(kind=8)              :: nu_0    = clight/1215.6701d-8         ! emission frequency [Hz] 
  ! ------ spec_type == 'flat_fnu' : photons have a flat distribution in nu, between nu_min and nu_max
  real(kind=8)              :: nu_min  = clight/1221d-8              ! min frequency [Hz]
  real(kind=8)              :: nu_max  = clight/1210d-8              ! max frequency [Hz]
  ! ------ spec_type == 'gauss' : photons have a gaussian distribution in nu.
  real(kind=8)              :: nu_cen   = clight / 1215.6701d-8      ! central frequency [Hz]
  real(kind=8)              :: velwidth = 10.0                       ! line width in velocity [km/s]  
  ! --- miscelaneous
  integer                   :: ranseed = 1234                        ! seed for random generator
  logical                   :: verbose = .true.
  ! --------------------------------------------------------------------------



  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: PhotonsFromSourceModel path/to/params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_PhotonsFromSourceModel_params(parameter_file)
  if (verbose) call print_PhotonsFromSourceModel_params
  ! ------------------------------------------------------------


  ! --------------------------------------------------------------------------------------
  ! make source shines
  ! --------------------------------------------------------------------------------------
  allocate(photgrid(nphot))
  iran = -abs(ranseed)

  do i=1,nphot
     photgrid(i)%ID    = i
     select case(trim(spec_type))
     case('monochromatic')
        nu = nu_0
     case('flat_fnu')
        nu = ran3(iran) * (nu_max-nu_min) + nu_min 
     case('gauss')
        r1 = ran3(iran)
        r2 = ran3(iran)
        nu = sqrt(-log(r1)) * cos(2.0d0*pi*r2)
        nu = (velwidth * 1d5 * nu_cen / clight) * nu + nu_cen
     case default
        print*,'ERROR: unknown spec_type :',trim(spec_type)
     end select
     photgrid(i)%nu_em = nu
     select case(trim(source_type))
     case('pointlike')
        photgrid(i)%x_em  = source_pos
     case default
        print*,'ERROR: unknown source_type :',trim(source_type)
     end select
     photgrid(i)%iran  = -i !!iran
     call isotropic_direction(photgrid(i)%k_em,iran)
  enddo


  ! --------------------------------------------------------------------------------------
  ! write ICs
  ! --------------------------------------------------------------------------------------
  if (verbose) write(*,*) '--> writing file'
  open(unit=14, file=trim(outputfile), status='unknown', form='unformatted', action='write')
  write(14) nphot
  write(14) ranseed
  write(14) (photgrid(i)%ID,i=1,nphot)
  write(14) (photgrid(i)%nu_em,i=1,nphot)
  write(14) (photgrid(i)%x_em(:),i=1,nphot)
  write(14) (photgrid(i)%k_em(:),i=1,nphot)
  write(14) (photgrid(i)%iran,i=1,nphot)
  close(14)
  ! --------------------------------------------------------------------------------------

  deallocate(photgrid)

contains

  subroutine read_PhotonsFromSourceModel_params(pfile)

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
       if (line(1:24) == '[PhotonsFromSourceModel]') then
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
          case ('source_type')
             write(source_type,'(a)') trim(value)
          case ('verbose')
             read(value,*) verbose
          case ('ranseed')
             read(value,*) ranseed
          case ('spec_type')
             write(spec_type,'(a)') trim(value)
          case ('nu_0')
             read(value,*) nu_0
          case ('nu_cen')
             read(value,*) nu_cen
          case ('velwidth')
             read(value,*) velwidth
          case ('nu_min')
             read(value,*) nu_min
          case ('nu_max')
             read(value,*) nu_max
          case ('nphot')
             read(value,*) nphot
          end select
       end do
    end if
    close(10)

    return

  end subroutine read_PhotonsFromSourceModel_params


  subroutine print_PhotonsFromSourceModel_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)')         '[PhotonsFromSourceModel]'
       write(unit,'(a)')             '# input / output parameters'
       write(unit,'(a,a)')           '  outputfile      = ',trim(outputfile)
       write(unit,'(a)')             '# source type parameters'
       write(unit,'(a,a)')           '  source_type     = ',trim(source_type)
       write(unit,'(a,3(ES9.3,1x))') '  source_pos      = ',source_pos(1),source_pos(2),source_pos(3)
       write(unit,'(a)')             '# how source shines'
       write(unit,'(a,i8)')          '  nphot           = ',nphot
       write(unit,'(a,a)')           '  spec_type       = ',trim(spec_type)
       select case(trim(spec_type))
       case('monochromatic')
          write(unit,'(a,es9.3,a)')     '  nu_0            = ',nu_0, ' ! [Hz]'
       case('flat_fnu')
          write(unit,'(a,es9.3,a)')     '  nu_min          = ',nu_min, ' ! [Hz]'
          write(unit,'(a,es9.3,a)')     '  nu_max          = ',nu_max, ' ! [Hz]'
       case('gauss')
          write(unit,'(a,es9.3,a)')     '  nu_cen          = ',nu_cen, ' ! [Hz]'
          write(unit,'(a,es9.3,a)')     '  velwidth        = ',velwidth, ' ! [km/s]'
       case default
          print*,'ERROR: unknown spec_type :',trim(spec_type)
       end select
       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,i8)')          '  ranseed         = ',ranseed
       write(unit,'(a,L1)')          '  verbose         = ',verbose
       write(unit,'(a)')             ' '
    else
       write(*,'(a,a,a)')         '[PhotonsFromSourceModel]'
       write(*,'(a)')             '# input / output parameters'
       write(*,'(a,a)')           '  outputfile    = ',trim(outputfile)
       write(*,'(a)')             '# source type parameters'
       write(*,'(a,a)')           '  source_type   = ',trim(source_type)
       write(*,'(a,3(ES9.3,1x))') '  source_pos    = ',source_pos(1),source_pos(2),source_pos(3)
       write(*,'(a)')             '# how source shines'
       write(*,'(a,i8)')          '  nphot         = ',nphot
       write(*,'(a,a)')           '  spec_type     = ',trim(spec_type)
       select case(trim(spec_type))
       case('monochromatic')
          write(*,'(a,es9.3,a)')     '  nu_0            = ',nu_0, ' ! [Hz]'
       case('flat_fnu')
          write(*,'(a,es9.3,a)')     '  nu_min          = ',nu_min, ' ! [Hz]'
          write(*,'(a,es9.3,a)')     '  nu_max          = ',nu_max, ' ! [Hz]'
       case('gauss')
          write(*,'(a,es9.3,a)')     '  nu_cen          = ',nu_cen, ' ! [Hz]'
          write(*,'(a,es9.3,a)')     '  velwidth        = ',velwidth, ' ! [km/s]'
       case default
          print*,'ERROR: unknown spec_type :',trim(spec_type)
       end select
       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,i8)')          '  ranseed     = ',ranseed
       write(*,'(a,L1)')          '  verbose    = ',verbose
       write(*,'(a)')             ' '       
    end if

    return

  end subroutine print_PhotonsFromSourceModel_params


end program PhotonsFromSourceModel
