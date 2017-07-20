!--PEEL--
program Analyse

  use module_photon
  use module_mesh
  use module_domain
  use module_gas_composition
  
  implicit none

  type(photon_current),allocatable :: pgrid(:)
  logical,allocatable        :: pdone(:)
  real(kind=8),allocatable   :: nhi(:)
  
  integer(kind=4)            :: i,narg,n,ifile,j, ndomain
  integer(kind=4)            :: npeels
  character(2000)            :: parameter_file, line, file_compute_dom
  integer(kind=4)            :: icell, ileaf, iran , cnt
  real(kind=8)               :: nu
  character(2000),dimension(:),allocatable :: mesh_file_list 
  type(mesh)                              :: meshdom
  type(domain)                            :: compute_dom
  logical :: verbose = .true.

  character(2000)           :: DataDir,PhotonFile, fichier
  character(2000)           :: DomDumpFile  = 'domain_decomposition_params.dat' ! the file describing the outputs of CreateDomDump.

  ! ------------------------------------------------------------------------------------------------------------------------

  write(DataDir,'(a)') "/Users/blaizot/Documents/Astro/SIMULATIONS/P13-20h-1Mpc-MUSIC/Zoom-7-10508/hdRun-LSS150pkpc-new-new-prime-jeans-fbk/rascas/HI_D_SMCdust/116/"

  
  ! ------------------------------------------------------------------------------------------------------------------------
  ! Read ICs ...
  ! ------------------------------------------------------------------------------------------------------------------------
  write(PhotonFile,'(a,a)') trim(DataDir),'LyaLineMilesICs_2.dat'
  call init_photons_from_file(PhotonFile,pgrid)
  open(unit=33,file='Snap116-analyseICs.txt',status='unknown',form='formatted')
  do i = 1,size(pgrid)
     write(33,'(3(e18.10,1x))') pgrid(i)%xlast(:)
  end do
  close(33)
  deallocate(pgrid)
  ! ------------------------------------------------------------------------------------------------------------------------


  
  ! ------------------------------------------------------------------------------------------------------------------------
  ! Read photons (result of a rascas run)
  ! ------------------------------------------------------------------------------------------------------------------------
  if (verbose) write(*,*) 'reading photons ...'
  write(PhotonFile,'(a,a)') trim(DataDir),'LyaLineMiles_results_2.dat'
  call read_photon_dump(PhotonFile,pgrid)
  ! ------------------------------------------------------------------------------------------------------------------------

  
  ! ------------------------------------------------------------------------------------------------------------------------
  ! read domains and mesh information 
  ! ------------------------------------------------------------------------------------------------------------------------
  ! here, we parse the file written out by CreateDomDump, which contains all file names and nb of domains.
  if (verbose) print *,'--> reading domain and mesh...'
  write(fichier,'(a,a)') trim(DataDir),trim(DomDumpFile)
  open(unit=18,file=fichier,status='old',form='formatted')
  read(18,'(a)') line ; i = scan(line,'=') ; file_compute_dom = trim(DataDir)//trim(adjustl(line(i+1:)))
  read(18,'(a)') line ; i = scan(line,'=') ; read(line(i+1:),*) ndomain
  allocate(mesh_file_list(ndomain))
  do j = 1, ndomain
     read(18,'(a)') line ! this is .dom files -> we want the .mesh next line 
     read(18,'(a)') line ; i = scan(line,'=') ; mesh_file_list(j) = trim(DataDir)//trim(adjustl(line(i+1:)))
  end do
  close(18)
  call domain_constructor_from_file(file_compute_dom,compute_dom)
  if (verbose) then
     print *,'--> Ndomain =',ndomain
     print *,'    |_ ',trim(file_compute_dom)
     do j=1,ndomain
        print *,'    |_ ',trim(mesh_file_list(j))
     end do
  endif
  ! ------------------------------------------------------------------------------------------------------------------------



  ! ------------------------------------------------------------------------------------------------------------------------
  ! spawn rays from these peels in a chosen direction
  ! -> loop on domains:
  ! ---- read mesh
  ! ---- locate peels in mesh, find their cell and draw their weights.
  ! ---- clear domain
  ! ------------------------------------------------------------------------------------------------------------------------
  npeels = size(pgrid)
  allocate(nhi(npeels),pdone(npeels))
  pdone = .false.
  open(unit=33,file='Snap116-analyse.txt',status='unknown',form='formatted')
  do j = 1, ndomain
     if (verbose) print*,"Reading domain ", j
     call mesh_from_file(mesh_file_list(j),meshdom)
     if (verbose) print*,"-- done reading "
     cnt = 0
     do i = 1,npeels
        if (.not. pdone(i)) then ! if peel was not treated yet
            if (domain_contains_point(pgrid(i)%xlast,meshdom%domain)) then ! if peel is in current domain 
               icell   = in_cell_finder(meshdom,pgrid(i)%xlast)
               ileaf   = - meshdom%son(icell)
               nhi(i)  = meshdom%gas(ileaf)%nHI
               write(33,'(9(e18.10,1x))') nhi(i),pgrid(i)%xlast(:),pgrid(i)%nu_ext,real(pgrid(i)%nb_abs,8),meshdom%gas(ileaf)%v(:)
               pdone(i) = .true.
               cnt = cnt + 1
            end if
        end if
     end do
     print *,'npeels in domain : ',cnt
     call mesh_destructor(meshdom)
  end do
  close(33)
  ! check that all rays are initialised
  do i = 1,npeels
     if (.not. pdone(i)) then
        print*,'Some rays are not defined ... '
        stop
     end if
  end do
  ! ------------------------------------------------------------------------------------------------------------------------



end program Analyse
!--LEEP--
