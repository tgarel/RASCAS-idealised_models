program optical_bin

  use module_mesh
  use module_domain
  use module_ramses
  use module_opt_ray
  use module_spectra
  use mpi

  !Variables for mpi
  integer(kind=4)	 	           :: ierr, rank, npsize, status0(MPI_STATUS_SIZE)
  integer(kind=4),allocatable              :: chunksize(:), disp(:)
  integer,parameter			   :: master=0
  real(kind=8)	 			   :: t1=0d0, t2=0d0, t3=0d0
  character*10     		       	   :: ranktxt

  character(2000)                          :: parameter_file, line, file_compute_dom
  character(2000),allocatable,dimension(:) :: mesh_file_list, domain_file_list
  integer(kind=4)                          :: narg, i, j, counter, ndomain, counter2, n, m
  integer(kind=4)                          :: icell, ind, ioct, cell_level

  type(mesh)                               :: meshdom
  type(domain)                             :: compute_dom

  integer(kind=4)                          :: nstars, istar, nSEDgroups
  real(kind=8),allocatable                 :: star_pos(:,:),star_age(:),star_mass(:),star_vel(:,:),star_met(:), star_pos_cpu(:,:)
  real(kind=8),allocatable                 :: star_L(:,:), star_L_cpu(:,:), star_E(:,:), SiI_csn(:,:), star_E_cpu(:,:), SiI_csn_cpu(:,:), photo_rate(:), photo_rate_sum(:), photo_rate_all(:)
  real(kind=8)                             :: distance, direction(3)

  !real(kind=8),allocatable                 :: cell_flux_cpu(:,:), cell_flux_all(:,:)
  type(ray_type)             :: ray
  !real(kind=8) :: help(3)

  ! ------------------------------------------------------------------------------------------------
  ! user-defined parameters - read from section [OPTICAL_BIN] of the parameter file
  ! ------------------------------------------------------------------------------------------------
  ! --- inputs
  character(2000)           :: repository = './'                                ! ramses run directory (where all output_xxxxx dirs are).
  integer(kind=4)           :: snapnum = 1                                      ! ramses output number to use
  character(2000)           :: DomDumpDir   = 'test/'                           ! where outputs of CreateDomDump are
  character(2000)           :: DomDumpFile  = 'domain_decomposition_params.dat' ! the file describing the outputs of CreateDomDump (in DomDumpDir)
  ! --- outputs
  character(2000)           :: fileout = 'flux'                                 ! output file ... 
  ! --- miscelaneous
  logical                   :: verbose = .true.
  logical                   :: skip    = .true.
  ! ------------------------------------------------------------------------------------------------


  ! ------------------------------------------------------------------------------------------------
  !Commands related to mpi
  ! ------------------------------------------------------------------------------------------------
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, npsize, ierr)
  allocate(chunksize(npsize), disp(npsize))
  t1=MPI_WTIME()						!Initialize time to display the total time of execution
  ! open(unit=14,status='scratch')
  ! rewind(14)
  ! write (14,*) rank
  ! rewind(14)							!5 lines to create a character of the rank of the core,  to put it in names of output files   (more straightforward way ?)
  ! read  (14,*) ranktxt
  ! close (14)
  ! ------------------------------------------------------------------------------------------------


  ! -------------------- read parameters -----------------------------------------------------------
  narg = command_argument_count()
  if(narg .lt. 1 .and. rank==master)then
     write(*,*)'You should type: optical_bin params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_optical_params(parameter_file)
  if(verbose .and. rank==master .and. rank==master) call print_optical_params
  ! ------------------------------------------------------------------------------------------------


  ! -------------------- Read domain list from CreateDomDump param file ----------------------------
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
  ! ------------------------------------------------------------------------------------------------


  ! -------------------- get the computational domain ----------------------------------------------
  call domain_constructor_from_file(file_compute_dom,compute_dom)
  ! ------------------------------------------------------------------------------------------------

  ! -------------------- get the mesh, including the gas properties (ndust) ------------------------
  call mesh_from_file(mesh_file_list(1),meshdom) ! -> overwrites gas props if parameters are set to.
  ! ------------------------------------------------------------------------------------------------


  ! -------------------- the master computes the data ----------------------------------------------
  nSEDgroups = get_nOptBins()   !From module_spectra
  if(rank==master) then

     !get star properties in domain
     call ramses_read_stars_in_domain(repository,snapnum,compute_dom,star_pos,star_age,star_mass,star_vel,star_met)
     deallocate(star_vel)
     nstars = size(star_age)

     !initialize SED properties
     call init_SED_table()

     !compute SED luminosities,  in #photons/s
     allocate(star_L(nSEDgroups,nstars), star_E(nSEDgroups,nstars), SiI_csn(nSEDgroups,nstars))
     do i=1,nstars
        call inp_sed_table(star_age(i), star_met(i), 3, .false., star_E(:,i))
        call inp_sed_table(star_age(i), star_met(i), 1, .false., star_L(:,i))    !Third variable :  1 for L[#photons/s],  3 for mean energy in bin,  2+2*Iion for mean cross-section,  Iion: 1:HI,2:HeI,3:HeII,4:SiI, etc
        call inp_sed_table(star_age(i), star_met(i), 10, .false., SiI_csn(:,i))
        star_L(:,i) = star_L(:,i)*star_mass(i)/msun
     end do
     deallocate(star_age, star_met, star_mass)
  end if
  ! ------------------------------------------------------------------------------------------------


  ! -------------------- dispatch the star properties between the cores ----------------------------
  call MPI_BCAST(nstars, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)          !Every process needs nstars
  chunksize = (nstars+npsize-1)/npsize						!Size of the array for each core.  A little bit smaller for the last core.
  chunksize(npsize) = nstars - (npsize-1)*chunksize(1)
  disp(:) = (/ ((i-1)*chunksize(1), i=1,npsize) /)                              !Necessary for scatterv


  allocate(star_L_cpu(nSEDgroups, chunksize(rank+1)), star_E_cpu(nSEDgroups,chunksize(rank+1)), SiI_csn_cpu(nSEDgroups,chunksize(rank+1)), star_pos_cpu(3,chunksize(rank+1)))

  call MPI_SCATTERV(star_L, chunksize*nSEDgroups, disp*nSEDgroups, MPI_DOUBLE_PRECISION, star_L_cpu, chunksize(rank+1)*nSEDgroups, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(star_E, chunksize*nSEDgroups, disp*nSEDgroups, MPI_DOUBLE_PRECISION, star_E_cpu, chunksize(rank+1)*nSEDgroups, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(SiI_csn, chunksize*nSEDgroups, disp*nSEDgroups, MPI_DOUBLE_PRECISION, SiI_csn_cpu, chunksize(rank+1)*nSEDgroups, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
  call MPI_SCATTERV(star_pos, chunksize*3, disp*3, MPI_DOUBLE_PRECISION, star_pos_cpu, chunksize(rank+1)*3, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

  if(rank==master) deallocate(star_L, star_pos, star_E, SiI_csn)
  ! ------------------------------------------------------------------------------------------------



  !allocate(cell_flux_cpu(nSEDgroups,meshdom%nLeaf), cell_flux_all(nSEDgroups,meshdom%nLeaf))
  !cell_flux_cpu = 0d0 ; cell_flux_all = 0d0
  if(verbose .and. rank==master)print*, 'nLeaf      = ', meshdom%nLeaf
  if(verbose .and. rank==master)print*, 'nstars     = ', nstars
  if(verbose .and. rank==master)print*, 'nSEDgroups = ', nSEDgroups
  counter = 0 ; counter2 = 0

  !noct = meshdom%noct
  !ncoarse = meshdom%ncoarse

  allocate(ray%E(nSEDgroups), ray%tau(nSEDgroups))

  t2=MPI_WTIME()

  distance = 0.003 !code units

  n=4 ; allocate(photo_rate(2*n*n), photo_rate_sum(2*n*n), photo_rate_all(2*n*n)) ; photo_rate_sum = 0d0
  do istar=1,chunksize(rank+1)
     do i=1,n
        do j=1,2*n
           ray%ID = 0
           ray%E(:) = star_E_cpu(:,istar)
           ray%dist = 0d0
           ray%tau = 0d0
           ray%x_em = star_pos_cpu(:,istar)
           direction = compute_dom%sp%center(:) + (/ distance*sin(pi*i/n)*cos(pi*j/n), distance*sin(pi*i/n)*sin(pi*j/n), distance*cos(pi*i/n) /) - star_pos_cpu(:,istar)
           ray%k_em = direction(:)/norm2(direction)

           call ray_advance(ray, meshdom, compute_dom, norm2(direction))
           photo_rate(2*n*(i-1)+j) = 1/4d0/pi/ray%dist**2*sum(star_L_cpu(:,istar)*exp(-ray%tau(:))*SiI_csn_cpu(:,istar))
        end do
     end do
     photo_rate_sum(:) = photo_rate_sum(:) + photo_rate(:)
  end do

  print*, rank, photo_rate_sum
  
  call MPI_ALLREDUCE(photo_rate_sum, photo_rate_all, 2*n*n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  open(unit=10, file=trim(fileout), form='formatted', action='write')
  write(10,*) photo_rate_all
  close(10)


  t3=MPI_WTIME()
  write(*,*) 'Rank : ', rank, ' Computation time : ', t3-t2

  call MPI_FINALIZE(ierr)



contains


  ! recursive subroutine add_flux(icell, istar, radius)

  !   implicit none

  !   integer(kind=4),intent(in) :: icell, istar
  !   real(kind=8),intent(in)    :: radius
  !   integer(kind=4)            :: ind, ileaf
  !   real(kind=8)               :: cell_center(3)

  !   if(influence(icell, istar, radius, cell_center)) then

  !      ileaf = meshdom%son(icell)
  !      if(ileaf < 0 ) then
  !         !print*, cell_center
  !         call increment_flux(-ileaf, istar, cell_center)

  !      else
  !         do ind=0,7
  !            call add_flux(ncoarse + ind*noct + meshdom%son(icell), istar, radius)
  !         end do
  !      end if
  !   end if

  ! end subroutine add_flux




  ! subroutine increment_flux(ileaf, istar, cell_center)

  !   implicit none

  !   integer(kind=4),intent(in) :: ileaf, istar
  !   real(kind=8),intent(in)    :: cell_center(3)
  !   integer(kind=4)            :: ind, ioct, cell_level
  !   real(kind=8)               :: posoct(3), dist_tot
  !   type(ray_type)             :: ray

  !   ! ileaf = -meshdom%son(icell)
  !   ! if(ileaf < 0) then
  !   !    print*, 'Problem, not a leaf cell in subroutine increment_flux. Stopping the program.'
  !   !    stop
  !   ! end if

  !   counter = counter + 1
  !   ray%ID = ileaf
  !   ray%dist = 0d0
  !   ray%tau = 0d0
  !   ray%x_em = star_pos_cpu(:,istar)
  !   !ray%halo_ID = 0

  !   ! ind   = (icell - nCoarse - 1) / nOct + 1   ! JB: should we make a few simple functions to do all this ? 
  !   ! ioct  = icell - nCoarse - (ind - 1) * nOct
  !   ! cell_level   = meshdom%octlevel(ioct)      ! level of current cell
  !   ! posoct(:)    = meshdom%xoct(ioct,:)

  !   ! cell_center = get_cell_center(posoct,ind,cell_level)
  !   !print*, cell_center
  !   dist_tot = sqrt((cell_center(1)-star_pos_cpu(1,istar))**2 + (cell_center(2)-star_pos_cpu(2,istar))**2 + (cell_center(3)-star_pos_cpu(3,istar))**2)

  !   ray%k_em = (cell_center - star_pos_cpu(:,istar))/dist_tot

  !   !call ray_advance(ray, meshdom, compute_dom, dist_tot*box_size_cm)

  !   if(ray%halo_ID/=1)counter2 = counter2 + 1
  !   !tau_dust(ileaf) = ray%tau

  !   !if(ileaf==355311) print*, 'infos : ', ray%tau, dist_tot*box_size_cm, star_L(1,istar)

  !   cell_flux_cpu(:,ileaf) = cell_flux_cpu(:,ileaf) + exp(-ray%tau)/4d0/pi/dist_tot**2/box_size_cm**2 * star_L_cpu(:,istar)

  ! end subroutine increment_flux



  ! function influence(icell, istar, radius, cell_center)

  !   implicit none

  !   integer(kind=4),intent(in) :: icell, istar
  !   real(kind=8),intent(in)    :: radius
  !   real(kind=8),intent(out)   :: cell_center(3)
  !   integer(kind=4)            :: ileaf, ison
  !   real(kind=8)               :: dist_star, dist_dom, size, posoct(3)
  !   logical                    :: influence


  !   ileaf = meshdom%son(icell)

  !   if(ileaf == 0) then
  !      cell_center = 0d0
  !      influence = .false.
  !      print*, 'ileaf=0'
  !   else

  !      if(ileaf < 0) then
  !         ind   = (icell - meshdom%nCoarse - 1) / meshdom%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
  !         ioct  = icell - meshdom%nCoarse - (ind - 1) * meshdom%nOct
  !         cell_level   = meshdom%octlevel(ioct)      ! level of current cell
  !         posoct(:)    = meshdom%xoct(ioct,:)
  !         cell_center = get_cell_center(posoct,ind,cell_level)
  !         size = 0.5d0**cell_level*box_size_cm
  !      else
  !         ison = meshdom%son(icell)
  !         cell_center(:) = meshdom%xoct(ison,:)
  !         cell_level = meshdom%octlevel(ison)
  !         size = 0.5d0**(cell_level-1)*box_size_cm
  !      end if
  !   end if

  !   dist_star = sqrt((cell_center(1)-star_pos_cpu(1,istar))**2 + (cell_center(2)-star_pos_cpu(2,istar))**2 + (cell_center(3)-star_pos_cpu(3,istar))**2)*box_size_cm
  !   dist_dom = sqrt((cell_center(1)-compute_dom%sp%center(1))**2 + (cell_center(2)-compute_dom%sp%center(2))**2 + (cell_center(3)-compute_dom%sp%center(3))**2)*box_size_cm

  !   if(dist_star - sqrt(3d0)/2d0*size < radius .and. dist_dom - sqrt(3d0)/2d0*size < compute_dom%sp%radius) influence = .true.

  !   return

  ! end function influence


  ! function influence_radius(istar)

  !   implicit none

  !   integer(kind=4)            :: istar
  !   real(kind=8)               :: influence_radius

  !   influence_radius = 1./sqrt(10.)*sqrt(sum(star_L_cpu(:,istar))/4d0/pi)

  !   return

  ! end function influence_radius



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
          case ('skip')
             read(value,*) skip
          case ('DomDumpDir')
             write(DomDumpDir,'(a,"/")') trim(value)
          case ('repository')
             write(repository,'(a)') trim(value)
          case('snapnum')
             read(value,*) snapnum
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


    call read_ramses_params(pfile)
    call read_spectra_params(pfile)

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
       write(unit,'(a,a)')           '  repository     = ',trim(repository)
       write(unit,'(a,i5)')          '  snapnum        = ',snapnum
       write(unit,'(a,a)')           '  DomDumpFile    = ',trim(DomDumpFile)
       write(unit,'(a,a)')           '  fileout        = ',trim(fileout)
       write(unit,'(a,L1)')          '  verbose        = ',verbose
       write(unit,'(a,L1)')          '  skip           = ',skip
       write(unit,'(a)')             ' '
       call print_ramses_params(unit)
       call print_spectra_params(unit)
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ''
       write(*,'(a)')             '[OPTICAL_BIN]'
       write(*,'(a,a)')           '  DomDumpDir     = ',trim(DomDumpDir)
       write(*,'(a,a)')           '  repository     = ',trim(repository)
       write(*,'(a,i5)')          '  snapnum        = ',snapnum
       write(*,'(a,a)')           '  DomDumpFile    = ',trim(DomDumpFile)
       write(*,'(a,a)')           '  fileout        = ',trim(fileout)
       write(*,'(a,L1)')          '  verbose        = ',verbose
       write(*,'(a,L1)')          '  skip           = ',skip
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       call print_ramses_params
       call print_spectra_params
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
