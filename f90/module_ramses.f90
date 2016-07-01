module module_ramses

  implicit none

  ! private stuff ... 
  ! stuff read from AMR files
  integer(kind=4),private          :: ncell,ncoarse,ngridmax
  real(kind=8),allocatable,private :: xg(:,:)      ! grids position
  integer,allocatable,private      :: nbor(:,:)    ! neighboring father cells
  integer,allocatable,private      :: next(:)      ! next grid in list
  integer,allocatable,private      :: son(:)       ! sons grids
  integer,allocatable,private      :: cpu_map(:)  ! domain decomposition
  integer,allocatable,private      :: headl(:,:),taill(:,:),numbl(:,:),numbtot(:,:)
  integer,allocatable,private      :: headb(:,:),tailb(:,:),numbb(:,:)
  real(KIND=8),dimension(1:3)::xbound=(/0d0,0d0,0d0/)

  ! Stop pretending this would work in 2D : 
  integer(kind=4),parameter :: ndim = 3
  integer(kind=4),parameter :: twondim = 6
  integer(kind=4),parameter :: twotondim= 8 
  
  ! stuff read from the HYDRO files
  real(kind=8),allocatable,private :: var(:,:)
  real(kind=8),allocatable,private :: cell_x(:),cell_y(:),cell_z(:)
  integer,allocatable,private      :: cell_level(:)

  integer(kind=4),private :: ncpu,nvar


  public  :: read_leaf_cells, get_ngridtot
  private :: read_hydro, read_amr, get_nleaf, get_nvar, clear_amr, get_ncpu, get_param_real

  !==================================================================================
contains

  ! ----------------
  ! public functions 
  ! ----------------
  
  subroutine read_leaf_cells(repository, snapnum, nleaftot, nvar, &
       & xleaf, ramses_var, leaf_level)

    ! read all leaf cell from a simulation snapshot. Return standard 
    ! ramses ramses variables through ramses_var(nvar,nleaftot) and
    ! positions (xleaf(3,nleaftot)) and levels (leaf_level).

    implicit none 
    
    character(2000),intent(in)                :: repository
    integer(kind=4),intent(in)                :: snapnum
    integer(kind=4),intent(inout)             :: nleaftot, nvar
    real(kind=8),allocatable, intent(inout)   :: ramses_var(:,:)
    real(kind=8),allocatable,intent(inout)    :: xleaf(:,:)
    integer(kind=4),allocatable,intent(inout) :: leaf_level(:)

    logical                                   :: do_allocs
    integer(kind=4)                           :: icpu, ileaf, icell, ivar

    nleaftot = get_nleaf(repository,snapnum)  ! sets ncpu too 
    nvar     = get_nvar(repository,snapnum)
    allocate(ramses_var(nvar,nleaftot), xleaf(nleaftot,3), leaf_level(nleaftot))
    ncpu = get_ncpu(repository,snapnum)

    print *,'nleaftot, nvar, ncpu =',nleaftot,nvar,ncpu

    do_allocs = .true.
    ileaf = 0
    do icpu = 1,ncpu
       call read_amr(repository,snapnum,icpu,do_allocs)
       call read_hydro(repository,snapnum,icpu,do_allocs)
       do_allocs = .false.
       ! collect leaf cells
       do icell = 1,ncell
          if (son(icell)==0 .and. cpu_map(icell) == icpu) then
             ileaf = ileaf + 1
             do ivar = 1,nvar
                ramses_var(ivar,ileaf) = var(icell,ivar)
             end do
             xleaf(ileaf,1)    = cell_x(icell)
             xleaf(ileaf,2)    = cell_y(icell)
             xleaf(ileaf,3)    = cell_z(icell)
             leaf_level(ileaf) = cell_level(icell)
          end if
       end do
    end do
    call clear_amr

    return
    
  end subroutine read_leaf_cells

  
  function get_nGridTot(repository,snapnum)
    
    ! get total number of grids in the simulation 

    implicit none 

    integer(kind=4),intent(in)  :: snapnum
    character(1000),intent(in)  :: repository
    integer(kind=4)             :: get_nGridTot
    character(1000)             :: nomfich
    logical                     :: ok
    integer(kind=4)             :: icpu,ngrid_current

    ncpu = get_ncpu(repository,snapnum)
    get_nGridTot = 0
    do icpu = 1,ncpu
       write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
       inquire(file=nomfich, exist=ok)
       if(.not. ok)then
          write(*,*)'File '//TRIM(nomfich)//' not found'    
          stop
       end if
       open(unit=10,file=nomfich,form='unformatted',status='old',action='read')
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)ngrid_current
       close(10)
       get_nGridTot = get_nGridTot + ngrid_current
    end do

    return
    
  end function get_nGridTot


  function get_hi_density(nvar,ramses_var)
    ! return array with HI density (at/cc), computed from ramses raw variables
    integer(kind=4),intent(in) :: nvar
    real(kind=8),intent(in)    :: ramses_var(nvar) ! one cell only
    real(kind=8)               :: get_hi_density

    stop
    return
  end function get_hi_density

  function get_velocity_cgs(nvar,ramses_var)
    integer(kind=4),intent(in) :: nvar
    real(kind=8),intent(in)    :: ramses_var(nvar) ! one cell only
    real(kind=8),dimension(3)  :: get_velocity_cgs

    stop
    return
  end function get_velocity_cgs

    
  !subroutine ramses_read_star_particles(repository, snapnum, nstars, stars)
  !  character(2000),intent(in)                     :: repository
  !  integer(kind=4),intent(in)                     :: snapnum
  !  integer(kind=4),intent(inout)                  :: nstars
  !  type(star_particle),allocatable, intent(inout) :: stars(:)
  !  ! ...
  !  ! allocate(stars(nstars))
  !  ! ... 
  !end subroutine ramses_read_star_particles


  
  !==================================================================================
  ! ----------------
  ! private functions 
  ! ----------------
  
  subroutine read_hydro(repository,snapnum,icpu,do_allocs)

    implicit none

    integer(kind=4),intent(in)  :: snapnum,icpu
    character(1000),intent(in)  :: repository
    logical,intent(in)          :: do_allocs
    real(kind=8)                :: dx
    character(1000)             :: nomfich
    integer(kind=4)             :: i,nlevelmax,nboundary,ix,iy,iz,ind,ilevel,ibound,ncache,istart,ivar,iskip,igrid
    real(kind=8),allocatable    :: xc(:,:),xx(:)
    integer(kind=4),allocatable :: ind_grid(:)

    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/hydro_',snapnum,'.out',icpu
    open(unit=10,file=nomfich,form='unformatted',status='old',action='read')
    read(10)
    read(10)nvar
    read(10)
    read(10)nlevelmax
    read(10)nboundary
    read(10)
    if (do_allocs) then
       allocate(var(1:ncell,1:nvar))
       allocate(cell_x(1:ncell),cell_y(1:ncell),cell_z(1:ncell))
       allocate(cell_level(1:ncell))
    end if
    allocate(xc(1:twotondim,1:ndim))

    cell_level = -1
    
    do ilevel=1,nlevelmax

       dx=0.5d0**ilevel
       do ind=1,twotondim
          iz=(ind-1)/4
          iy=(ind-1-4*iz)/2
          ix=(ind-1-2*iy-4*iz)
          xc(ind,1)=(dble(ix)-0.5D0)*dx
          xc(ind,2)=(dble(iy)-0.5D0)*dx
          xc(ind,3)=(dble(iz)-0.5D0)*dx
       end do

       do ibound=1,nboundary+ncpu
          if(ibound<=ncpu)then  ! in the box 
             ncache=numbl(ibound,ilevel)   ! nb of grids in the simulated box. 
             istart=headl(ibound,ilevel)   ! head of grid list of simulated box
          else                  ! boundaries of simulated volume (aka useless)
             ncache=numbb(ibound-ncpu,ilevel)
             istart=headb(ibound-ncpu,ilevel)
          end if
          read(10)!ilevel2
          read(10)!numbl2
          if(ncache>0)then
             allocate(ind_grid(1:ncache))
             allocate(xx(1:ncache))
             ! Loop over level grids
             igrid=istart
             do i=1,ncache
                ind_grid(i)=igrid
                igrid=next(igrid)
             end do
             ! Loop over cells
             do ind=1,twotondim
                iskip=ncoarse+(ind-1)*ngridmax
                ! Loop over conservative variables
                do ivar=1,nvar
                   read(10) xx
                   if (ibound > ncpu) cycle  ! dont bother with boundaries
                   do i = 1, ncache
                      var(ind_grid(i)+iskip,ivar) = xx(i)
                   end do
                end do
                do i = 1,ncache
                   cell_x(ind_grid(i)+iskip) = xc(ind,1) + xg(ind_grid(i),1) -xbound(1)
                   cell_y(ind_grid(i)+iskip) = xc(ind,2) + xg(ind_grid(i),2) -xbound(2)
                   cell_z(ind_grid(i)+iskip) = xc(ind,3) + xg(ind_grid(i),3) -xbound(3)
                   cell_level(ind_grid(i)+iskip)      = ilevel
                end do
             end do
             deallocate(ind_grid,xx)
          end if
       end do
    end do
    deallocate(xc)
    close(10)

    return

  end subroutine read_hydro

  
  subroutine read_amr(repository,snapnum,icpu,do_allocs)
    
    implicit none 

    integer(kind=4),intent(in)  :: snapnum,icpu
    character(1000),intent(in)  :: repository
    logical,intent(in)          :: do_allocs
    character(1000)             :: nomfich 
    integer,allocatable         :: ind_grid(:),iig(:),grid(:)
    real(kind=8),allocatable    :: xxg(:)
    logical                     :: ok
    integer(kind=4)             :: i,nx,ny,nz,nlevelmax,nboundary
    integer(kind=4)             :: ilevel,ncache,ibound,idim,ind,iskip

    ! Vérification de l'existence des fichiers AMR
    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
    inquire(file=nomfich, exist=ok)
    if(.not. ok)then
       write(*,*)'File '//TRIM(nomfich)//' not found'    
       stop
    end if
    open(unit=10,file=nomfich,form='unformatted',status='old',action='read')
    ! Read grid variables
    read(10)
    read(10)
    read(10)nx,ny,nz
    xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

    ! Critical parameter: define the root level of the tree
    ncoarse=nx*ny*nz
    read(10)nlevelmax
    read(10)ngridmax
    read(10)nboundary
    read(10)!ngrid_current
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    if (do_allocs) allocate( headl(1:ncpu,1:nlevelmax),taill(1:ncpu,1:nlevelmax), &
         & numbl(1:ncpu,1:nlevelmax),numbtot(1:10,1:nlevelmax), &
         & headb(1:nboundary,1:nlevelmax),tailb(1:nboundary,1:nlevelmax), &
         & numbb(1:nboundary,1:nlevelmax) )
    headl=0;taill=0;numbl=0;numbtot=0;headb=0;tailb=0;numbb=0
    ! Allocate tree arrays
    if (do_allocs) then
       allocate(next(1:ngridmax))
       allocate(nbor(1:ngridmax,1:twondim))
    end if
    nbor=0; next=0
    ! Allocate grid center coordinates
    if (do_allocs) allocate(xg(1:ngridmax,1:ndim))
    xg=0.0D0
    ! Read levels variables
    read(10)headl(1:ncpu,1:nlevelmax)
    read(10)taill(1:ncpu,1:nlevelmax)
    read(10)numbl(1:ncpu,1:nlevelmax)
    read(10)numbtot(1:10,1:nlevelmax)
    ! Read boundary linked list
    if(nboundary>0)then
       read(10)headb(1:nboundary,1:nlevelmax)
       read(10)tailb(1:nboundary,1:nlevelmax)
       read(10)numbb(1:nboundary,1:nlevelmax)
    end if
    !  Read free memory
    read(10)
    next(ngridmax) = 0
    ! Read cpu boundaries
    read(10)
    read(10)
    ncell=ncoarse+twotondim*ngridmax
    if (do_allocs) allocate(son(1:ncell),cpu_map(1:ncell))
    son=0; cpu_map=0
    ! Read coarse level
    read(10)son(1:ncoarse)       
    read(10)
    read(10)cpu_map(1:ncoarse)
    do ilevel=1,nlevelmax
       do ibound=1,nboundary+ncpu
          if(ibound<=ncpu)then
             ncache=numbl(ibound,ilevel)
          else
             ncache=numbb(ibound-ncpu,ilevel)
          end if
          if(ncache>0)then
             allocate(ind_grid(1:ncache))
             allocate(xxg(1:ncache))
             allocate(iig(1:ncache))
             allocate(grid(1:ncache))
             ! Read grid index
             read(10)ind_grid
             ! Read next index
             read(10)iig
             do i=1,ncache
                next(ind_grid(i))=iig(i)
             end do
             ! Read prev index (skip)
             read(10)iig
             ! Read grid center
             do idim=1,ndim
                read(10)xxg
                do i=1,ncache
                   xg(ind_grid(i),idim)=xxg(i)
                end do
             end do
             ! Read father index (skip)
             read(10)iig
             ! Read nbor index
             do ind=1,twondim
                read(10)iig
                do i=1,ncache
                   nbor(ind_grid(i),ind)=iig(i)
                end do
             end do
             ! Read son index
             do ind=1,twotondim
                iskip=ncoarse+(ind-1)*ngridmax
                read(10)iig
                do i=1,ncache
                   son(ind_grid(i)+iskip)=iig(i)
                end do
             end do
             ! Read cpu map
             do ind=1,twotondim
                iskip=ncoarse+(ind-1)*ngridmax
                read(10)iig
                do i=1,ncache
                   cpu_map(ind_grid(i)+iskip)=iig(i)
                end do
             end do
             ! Read refinement map (skip)
             do ind=1,twotondim
                read(10)!iig 
             end do
             deallocate(xxg,iig,grid,ind_grid)
          end if
       end do
    end do
    close(10)
    return


  end subroutine read_amr

  

  subroutine clear_amr

    implicit none

    deallocate(son,cpu_map,xg,nbor,next)
    deallocate(headl,taill,numbl,numbtot,headb,tailb,numbb)
    deallocate(var,cell_x,cell_y,cell_z,cell_level)
    
    return
    
  end subroutine clear_amr



  function get_nleaf(repository,snapnum)

    implicit none
    
    integer(kind=4),intent(in)  :: snapnum
    character(1000),intent(in)  :: repository
    integer(kind=4)             :: get_nleaf
    integer(kind=4)             :: icpu,icell
    logical                     :: do_allocs
    character(1000)             :: nomfich 
    integer,allocatable         :: ind_grid(:),iig(:)
    integer,allocatable         :: son(:)
    integer,allocatable         :: cpu_map(:)
    integer,allocatable         :: numbl(:,:)
    integer,allocatable         :: numbb(:,:)
    logical                     :: ok
    integer(kind=4)             :: i,nx,ny,nz,nlevelmax,nboundary,ncell,ncoarse,ngridmax
    integer(kind=4)             :: ilevel,ncache,ibound,idim,ind,iskip
    integer(kind=4)             :: ndim,twondim,twotondim
    
    ncpu = get_ncpu(repository,snapnum)
    get_nleaf = 0
    ndim = 3
    twondim = 2*ndim
    twotondim = 2**ndim

    do icpu = 1,ncpu
       do_allocs = (icpu == 1) 
       write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/amr_',snapnum,'.out',icpu
       inquire(file=nomfich, exist=ok)
       if(.not. ok)then
          write(*,*)'File '//TRIM(nomfich)//' not found'    
          stop
       end if
       open(unit=10,file=nomfich,form='unformatted',status='old',action='read')
       read(10)
       read(10)
       read(10)nx,ny,nz
       ncoarse=nx*ny*nz ! Critical parameter: define the root level of the tree
       read(10)nlevelmax
       read(10)ngridmax
       read(10)nboundary
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       read(10)
       if (do_allocs) allocate(numbl(1:ncpu,1:nlevelmax),numbb(1:nboundary,1:nlevelmax))
       numbl=0;numbb=0
       ! Read levels variables
       read(10)!headl(1:ncpu,1:nlevelmax)
       read(10)!taill(1:ncpu,1:nlevelmax)
       read(10)numbl(1:ncpu,1:nlevelmax)
       read(10)!numbtot(1:10,1:nlevelmax)
       ! Read boundary linked list
       if(nboundary>0)then
          read(10)!headb(1:nboundary,1:nlevelmax)
          read(10)!tailb(1:nboundary,1:nlevelmax)
          read(10)numbb(1:nboundary,1:nlevelmax)
       end if
       read(10) ! Read free memory
       read(10) ! Read cpu boundaries
       read(10)
       ncell=ncoarse+twotondim*ngridmax
       if (do_allocs) allocate(son(1:ncell),cpu_map(1:ncell))
       son=0; cpu_map=0
       ! Read coarse level
       read(10)son(1:ncoarse)       
       read(10)
       read(10)cpu_map(1:ncoarse)
       do ilevel=1,nlevelmax
          do ibound=1,nboundary+ncpu
             if(ibound<=ncpu)then
                ncache=numbl(ibound,ilevel)
             else
                ncache=numbb(ibound-ncpu,ilevel)
             end if
             if(ncache>0)then
                !!!if (ilevel < fg_levelmin) print*,'lev < levmin',ilevel
                allocate(ind_grid(1:ncache))
                allocate(iig(1:ncache))
                read(10)ind_grid ! Read grid index
                read(10) ! Read next index
                read(10) ! Read prev index
                do idim=1,ndim
                   read(10) ! Read grid center
                end do
                read(10) ! Read father index
                do ind=1,twondim
                   read(10) ! Read nbor index
                end do
                do ind=1,twotondim
                   iskip=ncoarse+(ind-1)*ngridmax
                   read(10)iig  ! Read son index 
                   do i=1,ncache
                      son(ind_grid(i)+iskip)=iig(i) 
                   end do
                end do
                do ind=1,twotondim
                   iskip=ncoarse+(ind-1)*ngridmax
                   read(10)iig ! Read cpu map
                   do i=1,ncache
                      cpu_map(ind_grid(i)+iskip)=iig(i)
                   end do
                end do
              
                do ind=1,twotondim
                   read(10) ! Read refinement map (skip)
                end do
                deallocate(iig,ind_grid)
             end if
          end do
       end do
       close(10)
       ! count leaf cells
       do icell = 1,ncell
          if (son(icell)==0 .and. cpu_map(icell) == icpu) then
             get_nleaf = get_nleaf + 1
          end if
       end do
    end do
    deallocate(son,cpu_map,numbl,numbb)
        
    return
  end function get_nleaf


  function get_nvar(repository,snapnum)
    implicit none 
    integer(kind=4),intent(in)  :: snapnum
    character(1000),intent(in)  :: repository
    character(1000)             :: nomfich
    integer(kind=4)             :: get_nvar,icpu
  
    icpu = 1
    write(nomfich,'(a,a,i5.5,a,i5.5,a,i5.5)') trim(repository),'/output_',snapnum,'/hydro_',snapnum,'.out',icpu
    open(unit=10,file=nomfich,form='unformatted',status='old',action='read')
    read(10)
    read(10)get_nvar
    close(10)
    return
  end function get_nvar


  function get_ncpu(repository,snapnum)

    implicit none
    
    integer(kind=4)            :: get_ncpu
    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum
    
    get_ncpu = nint(get_param_real(repository,snapnum,'ncpu'))
    
    return
  end function get_ncpu



  function get_param_real(repository,snapnum,param)

    implicit none 
    
    real(kind=8)               :: get_param_real
    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum
    character(*),intent(in)    :: param
    logical(kind=4)            :: not_ok
    character(512)             :: nomfich
    character(512)             :: line,name,value
    integer(kind=4)            :: i
    integer(kind=4),parameter  :: param_unit = 13

    not_ok = .true.
    write(nomfich,'(a,a,i5.5,a,i5.5,a)') trim(repository),'/output_',snapnum,'/info_',snapnum,'.txt'
    open(unit=param_unit,file=nomfich,status='old',form='formatted')
    do 
       read(param_unit,'(a)',end=2) line
       i = scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))
       ! check for a comment at end of line !
       i = scan(value,'!')
       if (i /= 0) value = trim(adjustl(value(:i-1)))

       if (trim(name) .eq. trim(param)) then 
          read(value,*) get_param_real
          not_ok = .false.
       end if

    end do
2   close (param_unit)  

    if (not_ok) then 
       write(6,*) '> parameter not found in infoxxx.txt :',trim(param)
       stop
    end if
    
    return

  end function get_param_real
  

    
end module module_ramses
!==================================================================================
!==================================================================================
