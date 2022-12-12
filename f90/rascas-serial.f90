program main

  ! serial version of rascas 

  use module_photon
  use module_mesh
  use module_domain
  use module_constants
  use module_utils, only : print_rascas_header
  use module_mock, only : read_mock_params, print_mock_params, peeling_off

  implicit none

  type(photon_current),dimension(:),allocatable :: photgrid
  type(mesh)                                    :: meshdom
  type(domain)                                  :: compute_dom
  integer(kind=4)                               :: nphot
  real(kind=8)                                  :: start, tmptime, finish
  character(2000)                               :: parameter_file, line, file_compute_dom
  character(2000),dimension(:),allocatable      :: mesh_file_list 
  integer(kind=4)                               :: narg, i, j, ndomain
  character(2000)                               :: DomDumpFile = 'domain_decomposition_params.dat' ! the file describing the outputs of CreateDomDump (in DomDumpDir)
  
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [RASCAS-serial] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- inputs 
  character(2000)           :: DomDumpDir   = 'test/'                   ! where outputs of CreateDomDump are
  character(2000)           :: PhotonICFile = 'Photon_IC_file.dat'      ! the file containing photons to cast (incl. full path)
  ! --- outputs
  character(2000)           :: fileout = 'photons_done.dat'   ! output file ... 
  ! --- miscelaneous
  logical                   :: verbose = .true.
  ! --------------------------------------------------------------------------

  call cpu_time(start)

  
  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: rascas-serial params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_serial_params(parameter_file)
  if(verbose)then
     call print_rascas_header
     call print_serial_params
     print *,'--> Working with the serial version of RASCAS'
  endif
  ! ------------------------------------------------------------
  
  
  
  ! -------------------- read ICs photons --------------------
  if (verbose) print *,'--> reading ICs photons in file: ',trim(PhotonICFile)
  call init_photons_from_file(PhotonICFile,photgrid)
  nphot = size(photgrid)
  if (verbose) print *,'--> Nphoton =',nphot
  ! ------------------------------------------------------------
  
  
  
  ! -------------------- Get domain properties --------------------
  ! here, we parse the file written out by CreateDomDump, which contains all file names and nb of domains.
  if (verbose) print *,'--> reading domains'
  open(unit=18,file=DomDumpFile,status='old',form='formatted')
  read(18,'(a)') line ; i = scan(line,'=') ; file_compute_dom = trim(DomDumpDir)//trim(adjustl(line(i+1:)))
  read(18,'(a)') line ; i = scan(line,'=') ; read(line(i+1:),*) ndomain
  if(ndomain/=1)then
     print *,'ERROR: ndomain /= 1 -> use the MPI version'
     stop
  endif
  allocate(mesh_file_list(ndomain))
  do j = 1, ndomain
     read(18,'(a)') line ! this is .dom files -> we want the .mesh next line 
     read(18,'(a)') line ; i = scan(line,'=') ; mesh_file_list(j) = trim(DomDumpDir)//trim(adjustl(line(i+1:)))
  end do
  close(18)
  call domain_constructor_from_file(file_compute_dom,compute_dom)
  call mesh_from_file(mesh_file_list(1),meshdom)
  if (verbose) then
     print *,'--> Ndomain =',ndomain
     print *,'    |_ ',trim(file_compute_dom)
     print *,'    |_ ',trim(mesh_file_list(1))
  endif
  ! ------------------------------------------------------------
  
  call cpu_time(tmptime)
  if (verbose) then
     print '(" --> Time = ",f12.3," seconds.")',tmptime-start
     print *,' '
     print *,'--> start processing photon packets, RT in progress'
  endif
  
  ! do the Monte Carlo Radiative Transfer
   call MCRT(nphot,photgrid,meshdom,compute_dom)

  if (verbose) print *,'--> RT done'

  ! some checks & logs
  if(verbose)then
     print *,' '
     print *,'--> some quick checks on the results...'
     ! test status of photons
     do i=1,nphot
        if(photgrid(i)%status == 0)print*,'ERROR: oh oh problem with photon status...',i,photgrid(i)%status
     enddo
     ! Some stats on photon status
     print *,'photon status:'
     print *,'# of photons             =',size(photgrid(:)%status)
     print *,'# of status=1 (escaped)  =',count(mask=(photgrid(:)%status==1))
     print *,'# of status=2 (absorbed) =',count(mask=(photgrid(:)%status==2))
     print *,'some diagnostics:'
     print *,'min max status      =',minval(photgrid%status),maxval(photgrid%status)
     print *,'min max pos x       =',minval(photgrid%xcurr(1)),maxval(photgrid%xcurr(1))
     print *,'min max pos y       =',minval(photgrid%xcurr(2)),maxval(photgrid%xcurr(2))
     print *,'min max pos z       =',minval(photgrid%xcurr(3)),maxval(photgrid%xcurr(3))
     print *,'min max nb scatt    =',minval(photgrid%nb_abs),maxval(photgrid%nb_abs)
     print *,'min max nu          =',minval(photgrid%nu_ext),maxval(photgrid%nu_ext)
     print *,'min max lambda      =',clight/maxval(photgrid%nu_ext)*cmtoA,clight/minval(photgrid%nu_ext)*cmtoA
     print *,'min max travel time =',minval(photgrid%time),maxval(photgrid%time)
     print *,'last scattering:'
     print *,'min max pos x       =',minval(photgrid%xlast(1)),maxval(photgrid%xlast(1))
     print *,'min max pos y       =',minval(photgrid%xlast(2)),maxval(photgrid%xlast(2))
     print *,'min max pos z       =',minval(photgrid%xlast(3)),maxval(photgrid%xlast(3))
  endif

  ! write results
  if (verbose) print *,' '
  if (verbose) print*,'--> writing results in file: ',trim(fileout)
  call dump_photons(fileout,photgrid)

  !--PEEL--
  ! and dump mocks
  if (peeling_off) then
     print*,'--> writing mock to file'
     call dump_mocks
     print*,'mock statistics:'
     print*,'   peels_count     = ',peels_count
     print*,'   rays_count      = ',rays_count
     print*,'   detectors_count = ',(detectors_count(i), i=1,nDirections)
  end if
  !--LEEP--

  call cpu_time(finish)
  if (verbose) then
     print*,'--> work done'
     print '(" --> Time = ",f12.3," seconds.")',finish-start
     print*,' '
  endif
  

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
       if ((line(1:15) == '[RASCAS-serial]').or.(line(1:15) == '[rascas-serial]')) then
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
             write(DomDumpDir,'(a)') trim(value)
          case ('PhotonICFile')
             write(PhotonICFile,'(a)') trim(value)
          case ('fileout')
             write(fileout,'(a)') trim(value)
          end select
       end do
    end if
    close(10)

    ! add path (DomDumpDir) to input files generated by CreateDomDump
    DomDumpFile  = trim(DomDumpDir)//trim(DomDumpFile)
    
    call read_mesh_params(pfile)

    call read_mock_params(pfile)
    
    return

  end subroutine read_serial_params

  
  subroutine print_serial_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a)')             '[RASCAS-serial]'
       write(unit,'(a,a)')           '  DomDumpDir     = ',trim(DomDumpDir)
       write(unit,'(a,a)')           '  PhotonICFile   = ',trim(PhotonICFile)
       write(unit,'(a,a)')           '  fileout        = ',trim(fileout)
       write(unit,'(a,L1)')          '  verbose        = ',verbose
       write(unit,'(a)')             ' '
       call print_mesh_params(unit)
       write(unit,'(a)')             ' '
       if (peeling_off) call print_mock_params(unit)
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ''
       write(*,'(a)')             '[RASCAS-serial]'
       write(*,'(a,a)')           '  DomDumpDir     = ',trim(DomDumpDir)
       write(*,'(a,a)')           '  PhotonICFile   = ',trim(PhotonICFile)
       write(*,'(a,a)')           '  fileout        = ',trim(fileout)
       write(*,'(a,L1)')          '  verbose        = ',verbose
       write(*,'(a)')             ' '       
       call print_mesh_params
       write(*,'(a)')             ' '
       if (peeling_off) call print_mock_params
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ' '
    end if

    return

  end subroutine print_serial_params
  
end program main
