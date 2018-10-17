program optical_bin

  use module_mesh
  use module_domain
  use module_gray_ray
  use mpi

  !Variables for mpi
  integer(kind=4)	 	           :: ierr, rank, npsize, status0(MPI_STATUS_SIZE), chunksize
  integer,parameter			   :: master=0
  real(kind=8)	 			   :: t1=0d0, t2=0d0
  character*10     		       	   :: ranktxt

  character(2000)                          :: parameter_file, line, file_compute_dom
  character(2000),allocatable,dimension(:) :: mesh_file_list, domain_file_list
  integer(kind=4)                          :: narg, i, j, counter, ndomain, counter2
  integer(kind=4)                          :: icell, ind, ioct, cell_level

  type(mesh)                               :: meshdom
  type(domain)                             :: compute_dom

  integer(kind=4)                          :: nstars, nSEDgroups, istar
  real(kind=8),allocatable                 :: star_pos(:,:), star_L(:,:)

  real(kind=8),allocatable                 :: cell_flux_cpu(:,:), cell_flux_all(:,:)

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [OPTICAL_BIN] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- inputs 
  character(2000)           :: DomDumpDir   = 'test/'                           ! where outputs of CreateDomDump are
  character(2000)           :: star_infos_file = 'star_infos'                   ! the file containing photons to cast (incl. full path)
  character(2000)           :: DomDumpFile  = 'domain_decomposition_params.dat' ! the file describing the outputs of CreateDomDump (in DomDumpDir)
  ! --- outputs
  character(2000)           :: fileout = 'flux'                             ! output file ... 
  ! --- miscelaneous
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------


  ! --------------------------------------------------------------------------
  !Commands related to mpi
  ! --------------------------------------------------------------------------
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, npsize, ierr)
  t1=MPI_WTIME()						!Initialize time to display the total time of execution
  open(unit=14,status='scratch')
  rewind(14)
  write (14,*) rank
  rewind(14)							!5 lines to create a character of the rank of the core,  to put it in names of output files   (more straightforward way ?)
  read  (14,*) ranktxt
  close (14)
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

  ! -------------------- Read the star properties and dispatch them between the cores --------------------------------------------
  open(unit=10,file=trim(star_infos_file), form='unformatted', action='read')
  read(10) nSEDgroups
  read(10) nstars

  chunksize = (nstars+npsize-1)/npsize						!Size of the array for each core.  A little bit smaller for the last core.
  if (rank==npsize-1) chunksize = nstars - (npsize-1)*chunksize
  
  if(rank==master) then
     allocate(star_pos(3,nstars),star_L(nSEDgroups,nstars))
     read(10) star_pos
     read(10) star_L

     do i=1,npsize-2			!Then the master sends each piece of data to the other threads
        call MPI_SEND(star_pos(:,i*chunksize+1:(i+1)*chunksize),chunksize*3,MPI_DOUBLE_PRECISION,i,11,MPI_COMM_WORLD,ierr)
        call MPI_SEND(star_L(:,i*chunksize+1:(i+1)*chunksize),chunksize*nSEDgroups,MPI_DOUBLE_PRECISION,i,12,MPI_COMM_WORLD,ierr)
     end do
     call MPI_SEND(star_pos(:,(npsize-1)*chunksize+1:nstars),(nstars-(npsize-1)*chunksize)*3,MPI_DOUBLE_PRECISION,npsize-1,11,MPI_COMM_WORLD,ierr)  !The last thread data has another size
     call MPI_SEND(star_L(:,(npsize-1)*chunksize+1:nstars),(nstars-(npsize-1)*chunksize)*nSEDgroups,MPI_DOUBLE_PRECISION,npsize-1,12,MPI_COMM_WORLD,ierr)
  else
     allocate(star_pos(3,chunksize),star_L(nSEDgroups,chunksize))		!The other cores receive their piece of data
     call MPI_RECV(star_pos,chunksize*3,MPI_DOUBLE_PRECISION,master,11,MPI_COMM_WORLD,status0,ierr)
     call MPI_RECV(star_L,chunksize*nSEDgroups,MPI_DOUBLE_PRECISION,master,12,MPI_COMM_WORLD,status0,ierr)
  end if
  close(10)
  ! ------------------------------------------------------------------------------------------


  allocate(cell_flux_cpu(nSEDgroups,meshdom%nLeaf), cell_flux_all(nSEDgroups,meshdom%nLeaf))
  cell_flux_cpu = 0d0 ; cell_flux_all = 0d0
  if(verbose .and. rank==master)print*, 'nLeaf      = ', meshdom%nLeaf
  if(verbose .and. rank==master)print*, 'nSEDgroups = ', nSEDgroups
  counter = 0 ; counter2 = 0

  noct = meshdom%noct
  ncoarse = meshdom%ncoarse
  
  do istar=1,1
     call add_flux(1,istar,influence_radius(istar))
     open(unit=rank,file=trim(fileout)//trim(ranktxt),status='replace',form='unformatted',action='write')
     write(rank) istar
     write(rank) cell_flux_cpu
     close(rank)
  end do
  

  !print*, minval(sum(cell_flux(:,:),1))
  !print*, minloc(cell_flux(1,:))
  !if(verbose)print*, 'Rank : ', rank, ' Number of rays that exitted the domain : ', counter2
  if(verbose)print*, 'Rank : ', rank, ' Number of rays : ', counter

  call MPI_ALLREDUCE(cell_flux_cpu, cell_flux_all, nSEDgroups*meshdom%nLeaf, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)

  if(rank==master) then
     open(unit=10,file=trim(fileout)//trim('_all'),status='replace',form='unformatted',action='write')
     write(10) cell_flux_all
     close(10)
  end if
  
  
  deallocate(star_pos,star_L)
  

  t2=MPI_WTIME()
  if(rank==master) write(*,*) 'Rank : ', rank, ' Execution time : ', t2-t1

  call MPI_FINALIZE(ierr)



contains


  recursive subroutine add_flux(icell, istar, radius)

    implicit none

    integer(kind=4),intent(in) :: icell, istar
    real(kind=8),intent(in)    :: radius
    integer(kind=4)            :: ind, ileaf
    real(kind=8)               :: cell_center(3)

    if(influence(icell, istar, radius, cell_center)) then

       ileaf = meshdom%son(icell)
       if(ileaf < 0 ) then
          call increment_flux(-ileaf, istar, cell_center)

       else
          do ind=0,7
             call add_flux(ncoarse + ind*noct + meshdom%son(icell), istar, radius)
          end do
       end if
    end if

  end subroutine add_flux
  


  subroutine increment_flux(ileaf, istar, cell_center)

    implicit none

    integer(kind=4),intent(in) :: ileaf, istar
    real(kind=8),intent(in)    :: cell_center(3)
    integer(kind=4)            :: ind, ioct, cell_level
    real(kind=8)               :: posoct(3), dist_tot
    type(ray_type)             :: ray

    ! ileaf = -meshdom%son(icell)
    ! if(ileaf < 0) then
    !    print*, 'Problem, not a leaf cell in subroutine increment_flux. Stopping the program.'
    !    stop
    ! end if

    counter = counter + 1
    ray%ID = ileaf
    ray%dist = 0d0
    ray%tau = 0d0
    ray%x_em = star_pos(:,istar)
    ray%halo_ID = 0

    ! ind   = (icell - nCoarse - 1) / nOct + 1   ! JB: should we make a few simple functions to do all this ? 
    ! ioct  = icell - nCoarse - (ind - 1) * nOct
    ! cell_level   = meshdom%octlevel(ioct)      ! level of current cell
    ! posoct(:)    = meshdom%xoct(ioct,:)

    ! cell_center = get_cell_center(posoct,ind,cell_level)
    !print*, cell_center
    dist_tot = sqrt((cell_center(1)-star_pos(1,istar))**2 + (cell_center(2)-star_pos(2,istar))**2 + (cell_center(3)-star_pos(3,istar))**2)

    ray%k_em = (cell_center - star_pos(:,istar))/dist_tot

    call ray_advance(ray, meshdom, compute_dom, dist_tot*box_size_cm)
    
    if(ray%halo_ID/=1)counter2 = counter2 + 1
    !tau_dust(ileaf) = ray%tau

    !if(ileaf==355311) print*, 'infos : ', ray%tau, dist_tot*box_size_cm, star_L(1,istar)

    cell_flux_cpu(:,ileaf) = cell_flux_cpu(:,ileaf) + exp(-ray%tau)/4d0/pi/dist_tot**2/box_size_cm**2 * star_L(:,istar)

  end subroutine increment_flux



  function influence(icell, istar, radius, cell_center)

    implicit none

    integer(kind=4),intent(in) :: icell, istar
    real(kind=8),intent(in)    :: radius
    real(kind=8),intent(out)   :: cell_center(3)
    integer(kind=4)            :: ileaf, ison
    real(kind=8)               :: dist, size, posoct(3)
    logical                    :: influence

    ileaf = meshdom%son(icell)

    if(ileaf == 0) then
       cell_center = 0d0
       influence = .false.
    else

       if(ileaf < 0) then
          ind   = (icell - meshdom%nCoarse - 1) / meshdom%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
          ioct  = icell - meshdom%nCoarse - (ind - 1) * meshdom%nOct
          cell_level   = meshdom%octlevel(ioct)      ! level of current cell
          posoct(:)    = meshdom%xoct(ioct,:)
          cell_center = get_cell_center(posoct,ind,cell_level)
          size = 0.5d0**cell_level*box_size_cm
       else
          ison = meshdom%son(icell)
          cell_center(:) = meshdom%xoct(ison,:)
          cell_level = meshdom%octlevel(ison)
          size = 0.5d0**(cell_level-1)*box_size_cm
       end if


       dist = sqrt((cell_center(1)-star_pos(1,istar))**2 + (cell_center(2)-star_pos(2,istar))**2 + (cell_center(3)-star_pos(3,istar))**2)*box_size_cm

       if(dist - sqrt(3d0)/2d0*size < radius) influence = .true.

    end if

    return

  end function influence


  function influence_radius(istar)

    implicit none

    integer(kind=4)            :: istar
    real(kind=8)               :: influence_radius

    influence_radius = 1./sqrt(10.)*sqrt(sum(star_L(:,istar))/4d0/pi)

    return

  end function influence_radius



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
          case ('star_infos_file')
             write(star_infos_file,'(a)') trim(value)
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
       write(unit,'(a,a)')           '  star_infos_file= ',trim(star_infos_file)
       write(unit,'(a,a)')           '  DomDumpFile    = ',trim(DomDumpFile)
       write(unit,'(a,a)')           '  fileout        = ',trim(fileout)
       write(unit,'(a,L1)')          '  verbose .and. rank==master        = ',verbose .and. rank==master
       write(unit,'(a)')             ' '
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ''
       write(*,'(a)')             '[OPTICAL_BIN]'
       write(*,'(a,a)')           '  DomDumpDir     = ',trim(DomDumpDir)
       write(*,'(a,a)')           '  star_infos_file= ',trim(star_infos_file)
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
