program FescFromStars

  ! Use RASCAS to compute column densities along sight-lines. 

  use module_gray_ray 
  use module_mesh
  use module_domain
  use module_constants

  implicit none

  type(ray_type),dimension(:),allocatable  :: rays
  type(mesh)                               :: meshdom
  type(domain)                             :: compute_dom
  integer                                  :: nrays
  real(kind=8)                             :: start, tmptime, finish

  character(2000) :: parameter_file, line, file_compute_dom
  character(2000),dimension(:),allocatable :: mesh_file_list 
  integer(kind=4) :: narg, i, j, ndomain
  
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [ColumnDensity] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- inputs 
  character(2000)           :: DataDir      = './'                      ! where input files below are 
  character(2000)           :: RaysICFile = 'Rays_IC_file.dat'      ! the file containing photons to cast.
  character(2000)           :: DomDumpFile  = 'domain_decomposition_params.dat' ! the file describing the outputs of CreateDomDump.
  ! --- outputs
  character(2000)           :: fileout = 'fescs.dat'   ! output file ... 
  ! --- parameters
  integer(kind=4)           :: ndirections=500 ! Number of directions (rays) from each source
  real(kind=8)              :: maxdist = -1    ! stop rays after this distance [cm] negative => ignored)
  real(kind=8)              :: maxtau  = -1    ! stop rays after this tau (overrides maxdist) (negative => ignored)
  ! --- halos - for escape fractions out of virial radii
  integer(kind=4)           :: nhalos=0 ! Number of halos
  character(2000)           :: HaloFile = 'halos.dat' ! the file containing halo ids, rvir, and pos
  ! --- miscelaneous
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------

  call cpu_time(start)

  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: FescFromStars params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_FescFromStars_params(parameter_file)
  if(verbose) call print_FescFromStars_params
  ! ------------------------------------------------------------

  

  ! -------------------- read ICs photons --------------------
  if (verbose) print *,'--> reading ICs photons in file: ',trim(RaysICFile)
  call init_rays_from_file(RaysICFile,rays)
  nrays = size(rays)
  if (verbose) print *,'--> Nb of ray sources =',nrays
  if (nrays .eq. 0) then
     if (verbose) print *,' '
     if (verbose) print*,'--> writing results in empty file: ',trim(fileout)
     call dump_rays(fileout,rays)
     stop
  endif
  ! ------------------------------------------------------------

  ! -------------------- read halo IDs and domains -------------
  if(use_halos) then
     if (verbose) print *,'--> reading halos in file: ',trim(HaloFile)
     call init_halos_from_file(HaloFile,halos)
     nhalos = size(halos)
     if (verbose) print *,'--> Nb of halos =',nhalos
  endif
  ! -----------------------------------------------------------

  
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
  call ComputeFesc(nrays,rays,meshdom,compute_dom,maxdist,maxtau,ndirections)

  !if (verbose) print *,'--> RT done'
  !if (verbose) print *,' '
  !if (verbose) print*,'--> writing results in file: ',trim(fileout)
  print*,''
  print*,'RT done, writing results in file: ',trim(fileout)
  call dump_rays(fileout,rays)

  call cpu_time(finish)
  if (verbose) then
     print*,'--> work done'
     print '(" --> Time = ",f12.3," seconds.")',finish-start
     print*,' '
  endif

contains

  subroutine read_FescFromStars_params(pfile)

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
       if (line(1:15) == '[FescFromStars]') then
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
          case ('RaysICFile')
             write(RaysICFile,'(a)') trim(value)
          case ('Halos')
             read(value,*) use_halos
          case ('HaloFile')
             write(HaloFile,'(a)') trim(value)
          case ('DomDumpFile')
             write(DomDumpFile,'(a)') trim(value)
          case ('fileout')
             write(fileout,'(a)') trim(value)
          case ('ndirections')
             read(value,*) ndirections
          case ('maxdist')
             read(value,*) maxdist
          case ('maxtau')
             read(value,*) maxtau
          case ('const_directions')
             read(value,*) const_directions
          end select
       end do
    end if
    close(10)

    if (maxtau > 0) then
       print*,'Will stop ray propagation at tau = ',maxtau
       maxdist = -1 
    end if
    
    ! add path (datadir) to input files 
    RaysICFile = trim(DataDir)//trim(RaysICFile)
    HaloFile = trim(DataDir)//trim(HaloFile)
    DomDumpFile  = trim(DataDir)//trim(DomDumpFile)
    
    call read_mesh_params(pfile)
    
    return

  end subroutine read_FescFromStars_params

  
  subroutine print_FescFromStars_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a)')             '[FescFromStars]'
       write(unit,'(a,a)')           '  DataDir     = ',trim(DataDir)
       write(unit,'(a,a)')           '  RaysICFile  = ',trim(RaysICFile)
       write(unit,'(a,a)')           '  DomDumpFile = ',trim(DomDumpFile)
       write(unit,'(a,a)')           '  fileout     = ',trim(fileout)
       write(unit,'(a,I10)')         '  ndirections = ',ndirections
       write(unit,'(a,L1)')          '  const_dir   = ',const_directions
       write(unit,'(a,ES10.3)')      '  maxdist     = ',maxdist
       write(unit,'(a,ES10.3)')      '  maxtau      = ',maxtau
       write(unit,'(a,L1)')          '  verbose     = ',verbose
       write(unit,'(a)')             ' '
       call print_mesh_params(unit)
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ''
       write(*,'(a)')             '[FescFromStars]'
       write(*,'(a,a)')           '  DataDir     = ',trim(DataDir)
       write(*,'(a,a)')           '  RaysICFile  = ',trim(RaysICFile)
       write(*,'(a,a)')           '  DomDumpFile = ',trim(DomDumpFile)
       write(*,'(a,a)')           '  fileout     = ',trim(fileout)
       write(*,'(a,I10)')         '  ndirections = ',ndirections
       write(*,'(a,L1)')           '  const_dir   = ',const_directions
       write(*,'(a,ES10.3)')      '  maxdist     = ',maxdist
       write(*,'(a,ES10.3)')      '  maxtau      = ',maxtau
       write(*,'(a,L1)')          '  verbose     = ',verbose
       write(*,'(a)')             ' '       
       call print_mesh_params
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
    end if

    return

  end subroutine print_FescFromStars_params
  
end program FescFromStars
