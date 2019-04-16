program main

  ! Use RASCAS to compute column densities along sight-lines. 

  use module_CD
  use module_mesh
  use module_domain
  use module_constants

  implicit none

  type(ray_type),dimension(:),allocatable  :: rays
  type(mesh)                               :: meshdom
  type(domain)                             :: compute_dom
  integer                                  :: nrays, n
  real(kind=8)                             :: start, tmptime, finish

  character(2000) :: parameter_file, line, file_compute_dom, DomDumpFile
  character(2000),dimension(:),allocatable :: mesh_file_list 
  integer(kind=4) :: narg, i, j, ndomain
  real(kind=8),allocatable  :: kdir(:,:)

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [ColumnDensity] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- inputs 
  character(2000)           :: DomDumpDir = '/DATA_simu/HI_D_SMC/00164/' ! where outputs of CreateDomDump are 
  character(2000)           :: RaysICFile = '/DATA_simu/HI_D_SMC/00164/star1/raysIC.dat' ! the file containing initial conditions for rays
  ! --- outputs
  character(2000)           :: fileout    = '/DATA_simu/HI_D_SMC/00164/star1/result.dat' ! output file ... 
  ! --- miscelaneous
  logical                   :: verbose    = .true.
  !For the histogram Column_Density VS gas velocity
  integer(kind=4)           :: nDirections    = 1
  character(2000)           :: direction_file = 'params_mock.dat'
  real(kind=8)              :: vmin = -1d3       !velocities between which we want to plot NSiII vs velocity of gas cell
  real(kind=8)              :: vmax = 1d3        
  integer(kind=4)           :: nBins = 100        ! number of bins of velocities,  if <1, don't do the histo 2d
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
  call read_ColumnDensity_params(parameter_file)
  if(verbose) call print_ColumnDensity_params
  ! ------------------------------------------------------------



  ! -------------------- read ICs photons --------------------
  if (verbose) print *,'--> reading ICs photons in file: ',trim(RaysICFile)
  call init_rays_from_file(RaysICFile,rays)
  nrays = size(rays)
  if (verbose) print *,'--> Nb of rays =',nrays
  ! ------------------------------------------------------------

  ! -------------------- read directions -----------------------
  ! n=3
  ! allocate(kdir(n*2*n,3))
  ! do i=1,n
  !    do j=1,2*n
  !       kdir(2*n*(i-1)+j,:) = (/ sin(pi*i/n)*cos(pi*j/n), sin(pi*i/n)*sin(pi*j/n), cos(pi*i/n) /)
  !    end do
  ! end do
  allocate(kdir(nDirections,3))
  open(unit=10,file=direction_file,status='old',action='read',form='formatted')
  do i=1,nDirections
     read(10,*) kdir(i,:)
     kdir(i,:) = kdir(i,:)/norm2(kdir(i,:))
     read(10,*) 
     read(10,*) 
     read(10,*) 
     read(10,*) 
     read(10,*)
  end do
  close(10)
  ! ------------------------------------------------------------


  ! -------------------- Get domain properties --------------------
  ! here, we parse the file written out by CreateDomDump, which contains all file names and nb of domains.
  if (verbose) print *,'--> reading domain and mesh...'
  write(DomDumpFile,'(a,a)') trim(DomDumpDir),'/domain_decomposition_params.dat'
  open(unit=18,file=DomDumpFile,status='old',form='formatted')
  read(18,'(a)') line ; i = scan(line,'=') ; file_compute_dom = trim(DomDumpDir)//trim(adjustl(line(i+1:)))
  read(18,'(a)') line ; i = scan(line,'=') ; read(line(i+1:),*) ndomain
  if(ndomain/=1)then
     print *,'ndomain /= 1 use the MPI version'
     stop
  endif
  allocate(mesh_file_list(ndomain))
  do j = 1, ndomain
     read(18,'(a)') line ! this is .dom files -> we want the .mesh next line 
     read(18,'(a)') line ; i = scan(line,'=') ; mesh_file_list(j) = trim(DomDumpDir)//trim(adjustl(line(i+1:)))
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

  if(nBins>0) call init_histo_v_N(vmin,vmax,nBins,nrays,nDirections)

  call cpu_time(tmptime)
  if (verbose) print '(" --> Time = ",f12.3," seconds.")',tmptime-start

  ! do the RT stuff
  if (verbose) print *,'--> starting RT...'
  call ComputeCD(nrays,rays,meshdom,compute_dom,nDirections,kdir)
  if (verbose) print *,'--> RT done'

  if (verbose) print *,' '
  if (verbose) print*,'--> writing results in file: ',trim(fileout)
  call dump(fileout,rays)

  call cpu_time(finish)
  if (verbose) print '(" --> Time = ",f12.3," seconds.")',finish-start


contains

  subroutine read_ColumnDensity_params(pfile)

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
          case('vmin')
             read(value,*) vmin
          case('vmax')
             read(value,*) vmax
          case('nBins')
             read(value,*) nBins
          case ('DomDumpDir')
             write(DomDumpDir,'(a)') trim(value)
          case ('RaysICFile')
             write(RaysICFile,'(a)') trim(value)
          case('nDirections')
             read(value,*) nDirections
          case ('direction_file')
             write(direction_file,'(a)') trim(value)
          case ('fileout')
             write(fileout,'(a)') trim(value)
          end select
       end do
    end if
    close(10)

    call read_mesh_params(pfile)

    return

  end subroutine read_ColumnDensity_params


  subroutine print_ColumnDensity_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a)')             '[ColumnDensity]'
       write(unit,'(a,a)')           '  DomDumpDir  = ',trim(DomDumpDir)
       write(unit,'(a,a)')           '  RaysICFile  = ',trim(RaysICFile)
       write(unit,'(a,a)')           '  fileout     = ',trim(fileout)
       write(unit,'(ES10.3,1x)')     '  vmin        = ', vmin
       write(unit,'(ES10.3,1x)')     '  vmax        = ', vmax
       write(unit,'(a,i5)')            '  nBins       = ', nBins
       write(unit,'(a,i5)')            '  nDirections = ', nBins
       write(unit,'(a,a)')           'direction_file  = ',trim(direction_file)
       write(unit,'(a,L1)')          '  verbose     = ',verbose
       write(unit,'(a)')             ' '
       call print_mesh_params(unit)
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ''
       write(*,'(a)')             '[ColumnDensity]'
       write(*,'(a,a)')           '  DomDumpDir  = ',trim(DomDumpDir)
       write(*,'(a,a)')           '  RaysICFile  = ',trim(RaysICFile)
       write(*,'(a,a)')           '  fileout     = ',trim(fileout)
       write(*,'(a,ES10.3)')        '  vmin        = ', vmin
       write(*,'(a,ES10.3)')        '  vmax        = ', vmax
       write(*,'(a,i5)')            '  nBins       = ', nBins
       write(*,'(a,i5)')            '  nDirections = ', nDirections
       write(*,'(a,a)')           'direction_file  = ',trim(direction_file)
       write(*,'(a,L1)')          '  verbose     = ',verbose
       write(*,'(a)')             ' '       
       call print_mesh_params
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
    end if

    return

  end subroutine print_ColumnDensity_params

end program main
