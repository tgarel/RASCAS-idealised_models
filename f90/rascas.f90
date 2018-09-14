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
  ! user-defined parameters - read from section [RASCAS] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- inputs 
  character(2000)           :: DomDumpDir   = 'test/'                   ! where outputs of CreateDomDump are
  character(2000)           :: PhotonICFile = 'Photon_IC_file.dat'      ! the file containing photons to cast (incl. full path)
  character(2000)           :: DomDumpFile  = 'domain_decomposition_params.dat' ! the file describing the outputs of CreateDomDump (in DomDumpDir)
  ! --- outputs
  character(2000)           :: fileout = 'photons_done.dat'   ! output file ... 
  ! --- miscelaneous
  integer(kind=4)           :: nbundle = 1000
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------

  
  call cpu_time(start)

  call start_mpi
  call define_mpi_type

  nworker=nb_cpus-1
  if(nworker==0)then
     print*,'rascas is a parallel code, you should run it with MPI'
     stop
  end if

  ! -------------------- read parameters -----------------------------------------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: rascas params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_rascas_params(parameter_file)
  if(verbose .and. rank==0) call print_rascas_params
  ! ------------------------------------------------------------------------------------------  

  if (rank == 0 .and. verbose) print*,'--> Nworker =',nworker

  ! -------------------- Read domain list from CreateDomDump param file --------------------
  if (verbose .and. rank==0) print *,'--> reading domain list'
  open(unit=18,file=DomDumpFile,status='old',form='formatted')
  read(18,'(a)') line ; i = scan(line,'=') ; file_compute_dom = trim(DomDumpDir)//trim(adjustl(line(i+1:)))
  read(18,'(a)') line ; i = scan(line,'=') ; read(line(i+1:),*) ndomain
  allocate(mesh_file_list(ndomain),domain_file_list(ndomain))
  do j = 1, ndomain
     read(18,'(a)') line ; i = scan(line,'=') ; domain_file_list(j) = trim(DomDumpDir)//trim(adjustl(line(i+1:)))
     read(18,'(a)') line ; i = scan(line,'=') ; mesh_file_list(j) = trim(DomDumpDir)//trim(adjustl(line(i+1:)))
  end do
  close(18)
  ! ------------------------------------------------------------------------------------------

  
  call MPI_BARRIER(MPI_COMM_WORLD,code)

  ! Master - Worker separation
  if (rank == 0) then
     ! Master section, will dispatch the jobs.
     call master(file_compute_dom, ndomain, domain_file_list, PhotonICFile, nbundle, fileout)
  else
     ! Worker section, will mostly do radiative transfer (MCRT)
     call worker(file_compute_dom, ndomain, mesh_file_list, nbundle)
  end if

  ! write results: this is done by master

  ! deallocations
  deallocate(mesh_file_list,domain_file_list)
  
  call finish_mpi
  call cpu_time(finish)
  if(verbose .and. rank==0)then
     print*,' '
     print*,'--> work done, MPI finalized'
     print '(" --> Time = ",f12.3," seconds.")',finish-start
     print*,' '
  endif

contains

    subroutine read_rascas_params(pfile)

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
       if ((line(1:8) == '[RASCAS]').or.(line(1:8) == '[rascas]')) then
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
          case ('DomDumpDir')
             write(DomDumpDir,'(a,"/")') trim(value)
          case ('PhotonICFile')
             write(PhotonICFile,'(a)') trim(value)
          case ('DomDumpFile')
             write(DomDumpFile,'(a)') trim(value)
          case ('fileout')
             write(fileout,'(a)') trim(value)
          case ('nbundle')
             read(value,*) nbundle
          end select
       end do
    end if
    close(10)

    ! add path (DomDumpDir) to files generated by CreateDomDump
    DomDumpFile  = trim(DomDumpDir)//trim(DomDumpFile)
    
    call read_master_params(pfile)
    call read_worker_params(pfile)
    
    return

  end subroutine read_rascas_params

  
  subroutine print_rascas_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a)')             '[RASCAS]'
       write(unit,'(a,a)')           '  DomDumpDir     = ',trim(DomDumpDir)
       write(unit,'(a,a)')           '  PhotonICFile   = ',trim(PhotonICFile)
       write(unit,'(a,a)')           '  DomDumpFile    = ',trim(DomDumpFile)
       write(unit,'(a,a)')           '  fileout        = ',trim(fileout)
       write(unit,'(a,i8)')          '  nbundle        = ',nbundle
       write(unit,'(a,L1)')          '  verbose        = ',verbose
       write(unit,'(a)')             ' '
       call print_master_params(unit)
       write(unit,'(a)')             ' '
       call print_worker_params(unit)
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ''
       write(*,'(a)')             '[RASCAS]'
       write(*,'(a,a)')           '  DomDumpDir     = ',trim(DomDumpDir)
       write(*,'(a,a)')           '  PhotonICFile   = ',trim(PhotonICFile)
       write(*,'(a,a)')           '  DomDumpFile    = ',trim(DomDumpFile)
       write(*,'(a,a)')           '  fileout        = ',trim(fileout)
       write(*,'(a,i8)')          '  nbundle        = ',nbundle
       write(*,'(a,L1)')          '  verbose        = ',verbose
       write(*,'(a)')             ' '       
       call print_master_params
       write(*,'(a)')             ' '
       call print_worker_params
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
    end if

    return

  end subroutine print_rascas_params

  
end program main
