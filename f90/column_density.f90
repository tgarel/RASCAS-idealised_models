program optical_bin

  use module_mesh
  use module_domain
  use module_gray_ray
  use module_constants

  character(2000)                          :: parameter_file, line, file_compute_dom
  character(2000),allocatable,dimension(:) :: mesh_file_list, domain_file_list
  integer(kind=4)                          :: narg, i, j, n_angles
  integer(kind=4)                          :: icell, ind, ioct, cell_level

  type(mesh)                               :: meshdom
  type(domain)                             :: compute_dom
  real(kind=8),dimension(3)                :: pos, k

  real(kind=8),allocatable                 :: column_densities(:)

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [OPTICAL_BIN] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- inputs 
  character(2000)           :: DomDumpDir   = 'test/'                           ! where outputs of CreateDomDump are
  character(2000)           :: DomDumpFile  = 'domain_decomposition_params.dat' ! the file describing the outputs of CreateDomDump (in DomDumpDir)
  ! --- outputs
  character(2000)           :: fileout = 'column_density'                             ! output file ... 
  ! --- miscelaneous
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------
  

  ! -------------------- read parameters -----------------------------------------------------
  narg = command_argument_count()
  if(narg .lt. 1 .and. rank==master)then
     write(*,*)'You should type: optical_bin params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_optical_params(parameter_file)
  if(verbose .and. rank==master .and. rank==master) call print_optical_params
  ! ------------------------------------------------------------------------------------------


  ! -------------------- Read domain list from CreateDomDump param file --------------------
  if (verbose .and. rank==master) print *,'--> reading domain and mesh...'
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


  ! -------------------- get the computational domain ----------------------------------------
  call domain_constructor_from_file(file_compute_dom,compute_dom)
  ! ------------------------------------------------------------------------------------------

  ! -------------------- get the mesh, including the gas properties (ndust) ------------------
  call mesh_from_file(mesh_file_list(1),meshdom) ! -> overwrites gas props if parameters are set to.
  ! ------------------------------------------------------------------------------------------

  pos = (/ 5d-1, 5d-1, 5d-1 /)

  n_angles = 20
  allocate(column_densities(n_angles))

  do i=1,n_angles
     k = (/ 0d0, sin((i-1)*2*pi/n_angles), cos((i-1)*2*pi/n_angles) /)
     
     column_densities(i) = column_density(meshdom, compute_dom, pos, k, 3d23)
  end do

  print*, column_densities
  

contains



  subroutine read_optical_params(pfile)

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
       if ((line(1:13) == '[OPTICAL_BIN]').or.(line(1:13) == '[optical_bin]')) then
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
         ! case ('star_infos_file')
            ! write(star_infos_file,'(a)') trim(value)
          case ('DomDumpFile')
             write(DomDumpFile,'(a)') trim(value)
          case ('fileout')
             write(fileout,'(a)') trim(value)
          end select
       end do
    end if
    close(10)

    ! add path (DomDumpDir) to files generated by CreateDomDump
    DomDumpFile  = trim(DomDumpDir)//trim(DomDumpFile)

    return

  end subroutine read_optical_params


  subroutine print_optical_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a)')             '[OPTICAL_BIN]'
       write(unit,'(a,a)')           '  DomDumpDir     = ',trim(DomDumpDir)
       !write(unit,'(a,a)')           '  star_infos_file= ',trim(star_infos_file)
       write(unit,'(a,a)')           '  DomDumpFile    = ',trim(DomDumpFile)
       write(unit,'(a,a)')           '  fileout        = ',trim(fileout)
       write(unit,'(a,L1)')          '  verbose .and. rank==master        = ',verbose .and. rank==master
       write(unit,'(a)')             ' '
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ''
       write(*,'(a)')             '[OPTICAL_BIN]'
       write(*,'(a,a)')           '  DomDumpDir     = ',trim(DomDumpDir)
       !write(*,'(a,a)')           '  star_infos_file= ',trim(star_infos_file)
       write(*,'(a,a)')           '  DomDumpFile    = ',trim(DomDumpFile)
       write(*,'(a,a)')           '  fileout        = ',trim(fileout)
       write(*,'(a,L1)')          '  verbose        = ',verbose
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
    end if

    return

  end subroutine print_optical_params

end program optical_bin





 !print*, meshdom%nCell


  ! do icell=1,30
  !    if(meshdom%son(icell) < 0) then

  !       ray%ID = -meshdom%son(icell)
  !       ray%dist = 0d0
  !       ray%tau = 0d0
  !       ray%x_em = star_pos(:,1)
  !       ray%halo_ID = 1

  !       ind   = (icell - meshdom%nCoarse - 1) / meshdom%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
  !       ioct  = icell - meshdom%nCoarse - (ind - 1) * meshdom%nOct
  !       cell_level   = meshdom%octlevel(ioct)      ! level of current cell
  !       posoct(:)    = meshdom%xoct(ioct,:)

  !       cell_center = get_cell_center(posoct,ind,cell_level)
  !       dist_tot = sqrt((cell_center(1)-star_pos(1,1))**2 + (cell_center(2)-star_pos(2,1))**2 + (cell_center(3)-star_pos(3,1))**2)

  !       ray%k_em = (cell_center - star_pos(:,1))/dist_tot

  !       call ray_advance(ray, meshdom, compute_dom, dist_tot*box_size_cm)
  !       !if(ray%halo_ID/=1)counter = counter + 1
  !       print*, abs(dist_tot*box_size_cm - ray%dist)/1d3
  !       !tau_dust(-meshdom%son(icell)) = ray%tau

  !    end if
! end do


  ! icell = 1
  ! ind   = (icell - meshdom%nCoarse - 1) / meshdom%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
  ! ioct  = icell - meshdom%nCoarse - (ind - 1) * meshdom%nOct
  ! cell_level   = meshdom%octlevel(ioct)      ! level of current cell
  ! ison = meshdom%son(icell)
  ! posoct(:) = meshdom%xoct(ison,:)
  !print*, icell, ind, ioct, cell_level, ison, posoct

  ! print*, maxval(cell_flux(1,:))
