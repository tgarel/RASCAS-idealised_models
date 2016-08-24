program output2subvol

  use module_domain
  use module_ramses
  use module_mesh
  use module_gas_composition
  use module_select

  implicit none

  type(domain)                            :: subvol
  !type(gas)                :: gas_leaves,selected_leaves
  type(gas),dimension(:),allocatable      :: gas_leaves,selected_leaves
  type(mesh)                              :: meshout
  real(kind=8),dimension(:,:),allocatable :: x_leaf, xleaves
  real(kind=8),dimension(:,:),allocatable :: ramses_var
  integer,dimension(:),allocatable        :: leaf_level, leaveslevel, ind_sel

  integer                   :: noctsnap,nleaftot,nvar,nleaves
  integer,parameter         :: snapnum = 100
  character(2000),parameter :: repository='/Users/leo/ASTROPHYSICS/data/Disks/Run31/'
  character(2000)           :: fichier_subvol
  real(kind=8),dimension(3) :: pos
  real(kind=8)              :: rvir
  character(10)             :: testtype

  testtype='sphere'
  fichier_subvol = 'test.dat'


  ! halomaker connection
  call get_halo(pos,rvir)
  
  ! def du domaine de calcul : sphere
  call domain_constructor_from_scratch(subvol,testtype,xc=pos(1),yc=pos(2),zc=pos(3),r=rvir) 

  ! lecture de toutes les feuilles de la simu
  ! peut-etre get_nleaf d'abord ... 
  call read_leaf_cells(repository, snapnum, nleaftot, nvar, x_leaf, ramses_var, leaf_level)
  nOctSnap = get_nGridTot(repository,snapnum)

  ! conversion des feuilles en proprietes voulues... -> construction du type gas
  call gas_from_ramses_leaves(repository, snapnum, nleaftot, nvar, ramses_var, gas_leaves)

  ! creation de mesh
  !call select_leaves_in_domain(subvol, x_leaf, ramses_var, leaf_level, selected_leaves)
  ! domain_gas = gas(selected_leaves)
  ! leaves = leaves(selected_leaves)  ! eqv. resize function
  ! en fait ici
  call select_in_domain(subvol, nleaftot, x_leaf, ind_sel)
  nleaves = size(ind_sel)
  call select_from_domain(arr_in=x_leaf, ind_sel=ind_sel, arr_out=xleaves)
  call select_from_domain(leaf_level, ind_sel, leaveslevel)
  call select_from_domain(gas_leaves, ind_sel, selected_leaves)

  call mesh_from_leaves(nOctSnap,subvol,nleaves,selected_leaves,xleaves,leaveslevel,meshout)

  ! It would be nice to add a header to the file with some info related to units and halos, to be self-contained...
  ! how to do that?
  call dump_mesh(meshout, fichier_subvol)

  call mesh_destructor(meshout)

contains

  subroutine get_halo(pos,rvir)
    ! fake routine for test...
    real(kind=8),dimension(3),intent(out) :: pos
    real(kind=8),intent(out) :: rvir

    pos = (/0.5,0.5,0.5/)
    rvir = 0.1

  end subroutine get_halo

end program Output2subvol
