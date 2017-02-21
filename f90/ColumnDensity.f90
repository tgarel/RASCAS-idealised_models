program main

  ! Use RASCAS to compute column densities along sight-lines. 

  use module_photon
  use module_mesh
  use module_domain
  use module_uparallel
  use module_constants

  implicit none

  type(photon_current),dimension(:),allocatable :: photgrid
  type(mesh)                                    :: meshdom
  type(domain)                                  :: compute_dom
  integer                                       :: nphot
  real(kind=8)                                  :: start, tmptime, finish

  character(2000) :: parameter_file, line, file_compute_dom
  character(2000),dimension(:),allocatable :: mesh_file_list 
  integer(kind=4) :: narg, i, j, ndomain
  
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [ColumnDensity] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- inputs 
  character(2000)           :: DataDir      = './'                      ! where input files below are 
  character(2000)           :: PhotonICFile = 'Photon_IC_file.dat'      ! the file containing photons to cast.
  character(2000)           :: DomDumpFile  = 'MCLya_domain_params.dat' ! the file describing the outputs of CreateDomDump.
  ! --- outputs
  character(2000)           :: fileout = 'photons_done.dat'   ! output file ... 
  ! --- miscelaneous
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------

  call cpu_time(start)

  
  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: ColumnDensity params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_serial_params(parameter_file)
  if(verbose) call print_serial_params
  ! ------------------------------------------------------------

  

  ! -------------------- read ICs photons --------------------
  if (verbose) print *,'--> reading ICs photons in file: ',trim(PhotonICFile)
  call init_photons_from_file(PhotonICFile,photgrid)
  nphot = size(photgrid)
  if (verbose) print *,'--> Nphoton =',nphot
  ! ------------------------------------------------------------


  
  ! -------------------- Get domain properties --------------------
  ! here, we parse the file written out by CreateDomDump, which contains all file names and nb of domains.
  if (verbose) print *,'--> reading domain and mesh...'
  open(unit=18,file=DomDumpFile,status='old',form='formatted')
  read(18,'(a)') line ; i = scan(line,'=') ; file_compute_dom = trim(DataDir)//trim(adjustl(line(i+1:)))
  read(18,'(a)') line ; i = scan(line,'=') ; read(line(i+1:),*) ndomain
  if(ndomain/=1)then
     print *,'ndomain /= 1 use the MPI version'
     stop
  endif
  allocate(mesh_file_list(ndomain))
  do j = 1, ndomain
     read(18,'(a)') line ! this is .dom files -> we want the .mesh next line 
     read(18,'(a)') line ; i = scan(line,'=') ; mesh_file_list(j) = trim(DataDir)//trim(adjustl(line(i+1:)))
  end do
  close(18)
  
  call domain_constructor_from_file(file_compute_dom,compute_dom)
  if (verbose) print*,'computational domain built'
  call mesh_from_file(mesh_file_list(1),meshdom)
  if (verbose) print*,'mesh read'
  
  if (verbose) then
     print *,'--> Ndomain =',ndomain
     print *,'    |_ ',trim(file_compute_dom)
     print *,'    |_ ',trim(mesh_file_list(1))
  endif
  ! ------------------------------------------------------------

  call cpu_time(tmptime)
  if (verbose) print '(" --> Time = ",f12.3," seconds.")',tmptime-start


  if (verbose) print *,'--> starting RT...'
  ! do the RT stuff
  call MCRT(nphot,photgrid,meshdom,compute_dom)

  if (verbose) print *,'--> RT done'

  if (verbose) print *,' '
  if (verbose) print*,'--> writing results in file: ',trim(fileout)
  call dump_photons(fileout,photgrid)

  call cpu_time(finish)
  if (verbose) print '(" --> Time = ",f12.3," seconds.")',finish-start

  

contains

  subroutine read_serial_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! ALSO read parameter form used modules which have parameters (mesh)
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
       if (line(1:15) == '[ColumnDensity]') then
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
          case ('verbose')
             read(value,*) verbose
          case ('DataDir')
             write(DataDir,'(a)') trim(value)
          case ('PhotonICFile')
             write(PhotonICFile,'(a)') trim(value)
          case ('DomDumpFile')
             write(DomDumpFile,'(a)') trim(value)
          case ('fileout')
             write(fileout,'(a)') trim(value)
          end select
       end do
    end if
    close(10)

    ! add path (datadir) to input files 
    PhotonICFile = trim(DataDir)//trim(PhotonICFile)
    DomDumpFile  = trim(DataDir)//trim(DomDumpFile)
    
    call read_mesh_params(pfile)
    
    return

  end subroutine read_serial_params

  
  subroutine print_serial_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a)')             '[ColumnDensity]'
       write(unit,'(a,a)')           '  DataDir        = ',trim(DataDir)
       write(unit,'(a,a)')           '  PhotonICFile   = ',trim(PhotonICFile)
       write(unit,'(a,a)')           '  DomDumpFile    = ',trim(DomDumpFile)
       write(unit,'(a,a)')           '  fileout        = ',trim(fileout)
       write(unit,'(a,L1)')          '  verbose        = ',verbose
       write(unit,'(a)')             ' '
       call print_mesh_params(unit)
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ''
       write(*,'(a)')             '[ColumnDensity]'
       write(*,'(a,a)')           '  DataDir        = ',trim(DataDir)
       write(*,'(a,a)')           '  PhotonICFile   = ',trim(PhotonICFile)
       write(*,'(a,a)')           '  DomDumpFile    = ',trim(DomDumpFile)
       write(*,'(a,a)')           '  fileout        = ',trim(fileout)
       write(*,'(a,L1)')          '  verbose        = ',verbose
       write(*,'(a)')             ' '       
       call print_mesh_params
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
    end if

    return

  end subroutine print_serial_params
  
end program main
