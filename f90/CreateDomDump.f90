program CreateDomDump

  use module_domain
  use module_ramses
  use module_mesh
  use module_gas_composition
  use module_select

  implicit none

  type(domain)                             :: domaine_de_calcul
  type(domain),dimension(:),allocatable    :: domain_list
  type(mesh)                               :: domain_mesh
  type(gas),dimension(:),allocatable       :: gas_leaves,selected_leaves

  real(kind=8),dimension(:,:),allocatable  :: x_leaf, xleaf_sel
  real(kind=8),dimension(:,:),allocatable  :: ramses_var
  integer,dimension(:),allocatable         :: leaf_level, leaflevel_sel, ind_sel

  integer :: noctsnap,nleaftot,nvar,nleaf_sel,i
  character(2000) :: toto,meshroot
  character(2000),dimension(:),allocatable :: domain_file_list, mesh_file_list

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [CreateDomDump] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: fichier = 'compute_domain.dom'  ! file to which computational domain will be written.
  character(2000)           :: repository = './'               ! ramses run directory (where all output_xxxxx dirs are).
  integer(kind=4)           :: snapnum = 1                     ! ramses output number to use
  ! --- computational domain  
  character(10)             :: typechar = 'sphere'    ! shape type of domain  // only sphere allowed now. 
  real(kind=8),dimension(3) :: pos = (/0.5,0.5,0.5/)  ! center of domain [code units]
  real(kind=8)              :: rvir = 0.3             ! radius of domain [code units]
  ! --- domain decomposition // Only works with a number of shells now, sharing center wiht sphere above.
  character(10)             :: decomp_dom_type = 'shell' ! shape type of domain  // only shell allowed now.
  integer(kind=4)           :: ndomain = 1               ! nb of domains in decomposition
  real(kind=8),allocatable  :: rin(:)                    ! inner radius of shell [code units], default value 0
  real(kind=8),allocatable  :: rout(:)                   ! outer radius of shell [code units], default value = min(rvir+0.05,0.5)
  ! --- miscelaneous
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------

  ! read user-defined parameters 
  call read_CreateDomDump_params('parameters.dat')
  if (verbose) call print_CreateDomDump_params
  
  ! Define a spherical computational domain. This domain describes the volume in which photons fly.
  call domain_constructor_from_scratch(domaine_de_calcul,typechar,xc=pos(1),yc=pos(2),zc=pos(3),r=rvir) 

  
  ! lecture de toutes les feuilles de la simu
  call read_leaf_cells(repository, snapnum, nleaftot, nvar, x_leaf, ramses_var, leaf_level)
  nOctSnap = get_nGridTot(repository,snapnum)

  ! conversion des feuilles en proprietes voulues... -> construction du type gas
  call gas_from_ramses_leaves(repository,snapnum,nleaftot,nvar,ramses_var, gas_leaves)

  ! domain decomposition 
  meshroot = 'domain_'
  allocate(domain_list(ndomain))
  allocate(domain_file_list(ndomain),mesh_file_list(ndomain))
  do i = 1, ndomain
     call domain_constructor_from_scratch(domain_list(i),decomp_dom_type,xc=pos(1),yc=pos(2),zc=pos(3),r_inbound=rin(i),r_outbound=rout(i))
     write(toto,'(i8)') i
     write(toto,'(a)') adjustl(toto)  ! remove leading spaces
     write(domain_file_list(i),'(a,a,a)') trim(meshroot),trim(toto),'.dom'
     write(mesh_file_list(i),'(a,a,a)') trim(meshroot),trim(toto),'.mesh'
     if (verbose) print *,i,trim(domain_file_list(i)),' ',trim(mesh_file_list(i))
  end do

  ! write master info
  call domain_write_file(fichier,domaine_de_calcul)
  open(unit=10, file="MCLya_domain_params.dat")
  write(10,*) 'computational_domain_file = ',trim(fichier)
  write(10,*) 'Ndomain = ',ndomain
  do i=1,ndomain
     write(10,*) 'mesh_domain_file = ',trim(domain_file_list(i))
     call domain_write_file(domain_file_list(i),domain_list(i))
  end do
  close(10)

  ! creation des mesh
  do i = 1,ndomain
     call select_in_domain(domain_list(i), nleaftot, x_leaf, ind_sel)
     call select_from_domain(arr_in=x_leaf,     ind_sel=ind_sel, arr_out=xleaf_sel)
     call select_from_domain(arr_in=leaf_level, ind_sel=ind_sel, arr_out=leaflevel_sel)
     call select_from_domain(arr_in=gas_leaves, ind_sel=ind_sel, arr_out=selected_leaves)
     nleaf_sel = size(ind_sel)
     call mesh_from_leaves(nOctSnap,domain_list(i),nleaf_sel, &
          selected_leaves,xleaf_sel,leaflevel_sel,domain_mesh)
     call dump_mesh(domain_mesh, mesh_file_list(i))
     call mesh_destructor(domain_mesh)
  enddo


contains
  
  subroutine read_CreateDomDump_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! ALSO read parameter form used modules (mesh)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present
    logical         :: ndomain_present 
    
    section_present = .false.
    ndomain_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:15) == '[CreateDomDump]') then
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
          case ('pos')
             read(value,*) pos(1),pos(2),pos(3)
          case ('rvir')
             read(value,*) rvir
          case ('typechar')
             write(typechar,'(a)') trim(value)
          case ('verbose')
             read(value,*) verbose
          case ('fichier')
             write(fichier,'(a)') trim(value)
          case ('repository')
             write(repository,'(a)') trim(value)
          case ('snapnum')
             read(value,*) snapnum
          case('decomp_dom_type')
             write(decomp_dom_type,'(a)') trim(value)
          case ('ndomain')
             ndomain_present = .true.
             read(value,*) ndomain
             allocate(rin(ndomain),rout(ndomain))
          case ('rin')
             read(value,*) rin(:)
          case ('rout')
             read(value,*) rout(:)
          end select
       end do
    end if
    close(10)

    if (.not. ndomain_present) then ! assign default values
       allocate(rin(ndomain),rout(ndomain))
       rin(1)  = 0.0d0
       rout(1) = min(rvir + 0.05d0, 0.5d0)
    end if
    
    call read_mesh_params(pfile)
    
    return

  end subroutine read_CreateDomDump_params

  
  subroutine print_CreateDomDump_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit
    character(100) :: fmt

    if (present(unit)) then 
       write(unit,'(a,a,a)')         '[CreateDomDump]'
       write(unit,'(a)')             '# input / output parameters'
       write(unit,'(a,a)')           '  fichier         = ',trim(fichier)
       write(unit,'(a,a)')           '  repository      = ',trim(repository)
       write(unit,'(a,i5)')          '  snapnum         = ',snapnum
       write(unit,'(a)')             '# computational domain parameters'
       write(unit,'(a,a)')           '  typechar        = ',trim(typechar)
       write(unit,'(a,3(ES9.3,1x))') '  pos             = ',pos(1),pos(2),pos(3)
       write(unit,'(a,ES9.3)')       '  rvir            = ',rvir
       write(unit,'(a)')             '# domain decomposition parameters'
       write(unit,'(a,a)')           '  decomp_dom_type = ',trim(decomp_dom_type)
       write(unit,'(a,i5)')          '  ndomain         = ',ndomain
       write(fmt,'(a,i3,a)') '(a,',ndomain,'(ES9.3,1x))'
       write(unit,fmt)               '  rin             = ',rin(:)
       write(unit,fmt)               '  rout            = ',rout(:)
       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,L1)')          '  verbose         = ',verbose
       write(unit,'(a)')             ' '
       call print_mesh_params(unit)
    else
       write(*,'(a,a,a)')         '[CreateDomDump]'
       write(*,'(a)')             '# input / output parameters'
       write(*,'(a,a)')           '  fichier    = ',trim(fichier)
       write(*,'(a,a)')           '  repository = ',trim(repository)
       write(*,'(a,i5)')          '  snapnum    = ',snapnum
       write(*,'(a)')             '# computational domain parameters'
       write(*,'(a,a)')           '  typechar   = ',trim(typechar)
       write(*,'(a,3(ES9.3,1x))') '  pos        = ',pos(1),pos(2),pos(3)
       write(*,'(a,ES9.3)')       '  rvir       = ',rvir
       write(*,'(a)')             '# domain decomposition parameters'
       write(*,'(a,a)')           '  decomp_dom_type = ',trim(decomp_dom_type)
       write(*,'(a,i5)')          '  ndomain         = ',ndomain
       write(fmt,'(a,i3,a)') '(a,',ndomain,'(ES9.3,1x))'
       write(*,fmt)               '  rin             = ',rin(:)
       write(*,fmt)               '  rout            = ',rout(:)
       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,L1)')          '  verbose    = ',verbose
       write(*,'(a)')             ' '       
       call print_mesh_params
    end if

    return

  end subroutine print_CreateDomDump_params



end program CreateDomDump
