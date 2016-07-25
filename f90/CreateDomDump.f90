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

  integer :: noctsnap,nleaftot,nvar,nleaf_sel,ndomain,i
  integer,parameter :: snapnum = 23 !100
  character(2000),parameter :: repository='/cral2/blaizot/SIMULATIONS/P13-20h-1Mpc-MUSIC/Zoom-7-1290/hdRun-LSS150pkpc-new-new-prime/'!'/Users/leo/ASTROPHYSICS/data/Disks/Run31/'
  character(2000) :: fichier,toto,meshroot
  character(2000),dimension(:),allocatable :: domain_file_list, mesh_file_list
  real(kind=8) :: rin,rout
  
  ! user-defined parameters - read from section [CreateDomDump] of the parameter file
  real(kind=8),dimension(3) :: pos = (/0.5,0.5,0.5/)  ! center of domain [code units]
  real(kind=8)              :: rvir = 0.3             ! radius of domain [code units]
  character(10)             :: typechar = 'sphere'    ! shape type of domain

  fichier = 'compute_domain.dom'

  ! read user-defined parameters 
  call read_CreateDomDump_params('parameters.dat')
  call print_CreateDomDump_params
  stop
  
  ! def du domaine de calcul : sphere
  call domain_constructor_from_scratch(domaine_de_calcul,typechar,xc=pos(1),yc=pos(2),zc=pos(3),r=rvir) 

  ! lecture de toutes les feuilles de la simu
  call read_leaf_cells(repository, snapnum, nleaftot, nvar, x_leaf, ramses_var, leaf_level)
  nOctSnap = get_nGridTot(repository,snapnum)

  ! conversion des feuilles en proprietes voulues... -> construction du type gas
  call gas_from_ramses_leaves(repository,snapnum,nleaftot,nvar,ramses_var, gas_leaves)


  ! decomposition en domaine (geometrie)
  meshroot = 'domain_'
  ndomain = 2
  allocate(domain_list(ndomain))
  allocate(domain_file_list(ndomain),mesh_file_list(ndomain))
  !call domain_constructor_from_scratch(liste_domaines(1),spherechar,xc=pos(1),yc=pos(2),zc=pos(3),r=rvir)
  !pos = (/0.25,0.25,0.25/)
  !call domain_constructor_from_scratch(liste_domaines(2),spherechar,xc=pos(1),yc=pos(2),zc=pos(3),r=rvir)
  ! 2 shells with overlap
  typechar ='shell'
  pos = (/0.5,0.5,0.5/)
  rin  = 0.00000
  rout = 0.03000
  call domain_constructor_from_scratch(domain_list(1),typechar,xc=pos(1),yc=pos(2),zc=pos(3),r_inbound=rin,r_outbound=rout)
  rin  = 0.02000
  rout = 0.30000
  call domain_constructor_from_scratch(domain_list(2),typechar,xc=pos(1),yc=pos(2),zc=pos(3),r_inbound=rin,r_outbound=rout)
  !do idomain = 1,ndomains
  ! def. du domaine courant
  !rin  = ...
  !rout = ... 
  !liste_domaines(idomain) = domain_constructor('shell',xc=0.5,yc=0.5,zc=0.5,rin=rin,rout=rout) 
  !end do

  ! write master info ... a formater... 
  do i=1,ndomain
     write(domain_file_list(i),'(a,i2.2,a)') trim(meshroot),i,'.dom'
     write(mesh_file_list(i),'(a,i2.2,a)') trim(meshroot),i,'.mesh'
!!!toto = 'domain_'//trim(char(i))//'.dat'
     !!print *,toto
     !!domain_file_list(i)=trim(toto)
     print *,i,trim(domain_file_list(i)),' ',trim(mesh_file_list(i))
  end do
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
    ! ALSO read parameter form used modules (gas_composition)
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
             read(value,*) typechar
          end select
       end do
    end if
    close(10)

    call read_gas_composition_params(pfile)
    
    return

  end subroutine read_CreateDomDump_params

  
  subroutine print_CreateDomDump_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)')         '[CreateDomDump]'
       write(unit,'(a,3(ES9.3,1x))') '  pos      : ',pos(1),pos(2),pos(3)
       write(unit,'(a,ES9.3)')       '  rvir     : ',rvir
       write(unit,'(a,a)')           '  typechar : ',trim(typechar)
       call print_gas_composition_params(unit)
    else
       write(*,'(a,a,a)')         '[CreateDomDump]'
       write(*,'(a,3(ES9.3,1x))') '  pos      : ',pos(1),pos(2),pos(3)
       write(*,'(a,ES9.3)')       '  rvir     : ',rvir
       write(*,'(a,a)')           '  typechar : ',trim(typechar)
       call print_gas_composition_params
    end if

    return

  end subroutine print_CreateDomDump_params



end program CreateDomDump
