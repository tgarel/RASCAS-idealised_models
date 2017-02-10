program main

  use module_parallel_mpi
  use module_master
  use module_worker

  implicit none

  real(kind=8)                             :: start,finish
  character(2000)                          :: parameter_file, line, file_compute_dom
  character(2000),allocatable,dimension(:) :: mesh_file_list, domain_file_list
  integer(kind=4)                          :: narg, i, j, ndomain

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [MCLya] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- inputs 
  character(2000)           :: DataDir      = 'test/'                   ! where input files below are 
  character(2000)           :: PhotonICFile = 'Photon_IC_file.dat'      ! the file containing photons to cast.
  character(2000)           :: DomDumpFile  = 'MCLya_domain_params.dat' ! the file describing the outputs of CreateDomDump.
  ! --- outputs
  character(2000)           :: fileout = 'photons_done.dat'   ! output file ... 
  ! --- miscelaneous
  integer(kind=4)           :: nbuffer = 1000
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------

  
  call cpu_time(start)

  call initialisation_mpi
  call type_derive

  nslave=nb_cpus-1

  if (rank == 0) then
     if(verbose) print*,'--> Nworker =',nslave
  end if

  
  ! -------------------- read parameters -----------------------------------------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: MCLya params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_mclya_params(parameter_file)
  if(rank==0) call print_mclya_params
  ! ------------------------------------------------------------------------------------------  

  
  ! -------------------- Read domain list from CreateDomDump param file --------------------
  if (verbose .and. rank==0) print *,'--> reading domain and mesh...'
  open(unit=18,file=DomDumpFile,status='old',form='formatted')
  read(18,'(a)') line ; i = scan(line,'=') ; file_compute_dom = trim(DataDir)//trim(adjustl(line(i+1:)))
  read(18,'(a)') line ; i = scan(line,'=') ; read(line(i+1:),*) ndomain
  allocate(mesh_file_list(ndomain),domain_file_list(ndomain))
  do j = 1, ndomain
     read(18,'(a)') line ; i = scan(line,'=') ; domain_file_list(j) = trim(DataDir)//trim(adjustl(line(i+1:)))
     read(18,'(a)') line ; i = scan(line,'=') ; mesh_file_list(j) = trim(DataDir)//trim(adjustl(line(i+1:)))
  end do
  close(18)
  ! ------------------------------------------------------------------------------------------

  
  call MPI_BARRIER(MPI_COMM_WORLD,code)

  ! Master - Worker separation
  if (rank == 0) then
     ! Master section, will dispatch the jobs.
     call master(file_compute_dom, ndomain, domain_file_list, PhotonICFile, nbuffer, fileout)
  else
     ! Worker section, will mostly do radiative transfer (MCLya)
     call worker(file_compute_dom, ndomain, mesh_file_list, nbuffer)
  end if

  ! write results: this is done by master

  ! deallocations
  deallocate(mesh_file_list,domain_file_list)
  
  call finalisation_mpi
  call cpu_time(finish)
  if(rank==0)then
     print*,'--> work done, MPI finalized'
     print '(" --> Time = ",f12.3," seconds.")',finish-start
  endif

contains

    subroutine read_mclya_params(pfile)

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
       if (line(1:7) == '[MCLya]') then
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
             write(DataDir,'(a,"/")') trim(value)
          case ('PhotonICFile')
             write(PhotonICFile,'(a)') trim(value)
          case ('DomDumpFile')
             write(DomDumpFile,'(a)') trim(value)
          case ('fileout')
             write(fileout,'(a)') trim(value)
          case ('nbuffer')
             read(value,*) nbuffer
          end select
       end do
    end if
    close(10)

    ! add path (datadir) to input files 
    PhotonICFile = trim(DataDir)//trim(PhotonICFile)
    DomDumpFile  = trim(DataDir)//trim(DomDumpFile)
    
    call read_master_params(pfile)
    call read_worker_params(pfile)
    
    return

  end subroutine read_mclya_params

  
  subroutine print_mclya_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a)')             '[MCLya]'
       write(unit,'(a,a)')           '  DataDir        = ',trim(DataDir)
       write(unit,'(a,a)')           '  PhotonICFile   = ',trim(PhotonICFile)
       write(unit,'(a,a)')           '  DomDumpFile    = ',trim(DomDumpFile)
       write(unit,'(a,a)')           '  fileout        = ',trim(fileout)
       write(unit,'(a,i8)')          '  nbuffer        = ',nbuffer
       write(unit,'(a,L1)')          '  verbose        = ',verbose
       write(unit,'(a)')             ' '
       call print_master_params(unit)
       write(unit,'(a)')             ' '
       call print_worker_params(unit)
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ''
       write(*,'(a)')             '[MCLya]'
       write(*,'(a,a)')           '  DataDir        = ',trim(DataDir)
       write(*,'(a,a)')           '  PhotonICFile   = ',trim(PhotonICFile)
       write(*,'(a,a)')           '  DomDumpFile    = ',trim(DomDumpFile)
       write(*,'(a,a)')           '  fileout        = ',trim(fileout)
       write(*,'(a,i8)')          '  nbuffer        = ',nbuffer
       write(*,'(a,L1)')          '  verbose        = ',verbose
       write(*,'(a)')             ' '       
       call print_master_params
       write(*,'(a)')             ' '
       call print_worker_params
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
    end if

    return

  end subroutine print_mclya_params

  
end program main
