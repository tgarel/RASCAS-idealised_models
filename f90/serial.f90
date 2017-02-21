program main

  ! serial version of rascas 

  use module_photon
  use module_mesh
  use module_domain
  use module_uparallel
  use module_constants

  implicit none

  type(photon_current),dimension(:),allocatable :: photgrid
  type(mesh)                                    :: meshdom
  type(domain)                                  :: compute_dom
  integer(kind=4)                               :: nphot
  real(kind=8)                                  :: start, tmptime, finish
  character(2000)                               :: parameter_file, line, file_compute_dom
  character(2000),dimension(:),allocatable      :: mesh_file_list 
  integer(kind=4)                               :: narg, i, j, ndomain
  
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [serial] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- inputs 
  character(2000)           :: DataDir      = 'test/'                   ! where input files below are 
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
     write(*,*)'You should type: serial params.dat'
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

  

#ifdef DEBUG
  print *,'--> check mesh dom'
  print *,meshdom%domain
  print *,meshdom%nCoarse,meshdom%nOct,meshdom%nLeaf,meshdom%nCell
  print *,minval(meshdom%xoct(:,:)),maxval(meshdom%xoct(:,:))
  print *,minval(meshdom%nbor(:,:)),maxval(meshdom%nbor(:,:))
  print *,minval(meshdom%octlevel(:)),maxval(meshdom%octlevel(:))
  print *,minval(meshdom%son(:)),maxval(meshdom%son(:))
  print *,minval(meshdom%father(:)),maxval(meshdom%father(:))
  !!!print *,'level of leaves =',minval(meshdom%octlevel(:), mask=(meshdom%son(:)<0)),maxval(meshdom%octlevel(:), mask=(meshdom%son(:)<0))
#endif

  call cpu_time(tmptime)
  if (verbose) print '(" --> Time = ",f12.3," seconds.")',tmptime-start


  ! do the Monte Carlo Radiative Transfer
  if (verbose) print *,'--> starting RT...'
  call MCRT(nphot,photgrid,meshdom,compute_dom)

  if (verbose) print *,'--> RT done'

  ! some checks & logs
  if(verbose)then
     ! test status of photons
     do i=1,nphot
        if(photgrid(i)%status == 0)print*,'ohoho problem with photon status...',i,photgrid(i)%status
     enddo
     ! Some stats on photon status
     print *,' '
     print *,'--> photon status...'
     print *,'# of photons             =',size(photgrid(:)%status)
     print *,'# of status=1 (escaped)  =',count(mask=(photgrid(:)%status==1))
     print *,'# of status=2 (absorbed) =',count(mask=(photgrid(:)%status==2))
     print *,'# of status=3 (crap, pb with precision/in_cell_finder) =',count(mask=(photgrid(:)%status==3))
     print *,' '
     print *,'--> Some diagnostics...'
     print *,'min max status      =',minval(photgrid%status),maxval(photgrid%status)
     print *,'min max pos x       =',minval(photgrid%xcurr(1)),maxval(photgrid%xcurr(1))
     print *,'min max pos y       =',minval(photgrid%xcurr(2)),maxval(photgrid%xcurr(2))
     print *,'min max pos z       =',minval(photgrid%xcurr(3)),maxval(photgrid%xcurr(3))
     print *,'min max nb scatt    =',minval(photgrid%nb_abs),maxval(photgrid%nb_abs)
     print *,'min max nu          =',minval(photgrid%nu_ext),maxval(photgrid%nu_ext)
     print *,'min max lambda      =',clight/maxval(photgrid%nu_ext)*cmtoA,clight/minval(photgrid%nu_ext)*cmtoA
     print *,'min max travel time =',minval(photgrid%time),maxval(photgrid%time)
     print *,'Last scattering'
     print *,'min max pos x       =',minval(photgrid%xlast(1)),maxval(photgrid%xlast(1))
     print *,'min max pos y       =',minval(photgrid%xlast(2)),maxval(photgrid%xlast(2))
     print *,'min max pos z       =',minval(photgrid%xlast(3)),maxval(photgrid%xlast(3))
  endif


  ! write results
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
       if (line(1:8) == '[serial]') then
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
       write(unit,'(a)')             '[serial]'
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
       write(*,'(a)')             '[serial]'
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
