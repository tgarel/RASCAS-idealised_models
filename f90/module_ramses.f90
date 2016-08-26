module module_ramses

  use module_constants, only : kB, mp, XH

  implicit none

  private 
  
  ! stuff read from AMR files
  integer(kind=4)                  :: ncell,ncoarse,ngridmax
  real(kind=8),allocatable         :: xg(:,:)      ! grids position
  integer,allocatable              :: nbor(:,:)    ! neighboring father cells
  integer,allocatable              :: next(:)      ! next grid in list
  integer,allocatable              :: son(:)       ! sons grids
  integer,allocatable              :: cpu_map(:)  ! domain decomposition
  integer,allocatable              :: headl(:,:),taill(:,:),numbl(:,:),numbtot(:,:)
  integer,allocatable              :: headb(:,:),tailb(:,:),numbb(:,:)
  real(KIND=8),dimension(1:3)      :: xbound=(/0d0,0d0,0d0/)  

  ! Stop pretending this would work in 2D 
  integer(kind=4),parameter :: ndim = 3
  integer(kind=4),parameter :: twondim = 6
  integer(kind=4),parameter :: twotondim= 8 
  
  ! stuff read from the HYDRO files
  real(kind=8),allocatable         :: var(:,:)
  real(kind=8),allocatable         :: cell_x(:),cell_y(:),cell_z(:)
  integer,allocatable              :: cell_level(:)

  integer(kind=4)                  :: ncpu,nvar

  ! conversion factors (units)
  logical                        :: conversion_scales_are_known = .False. 
  real(kind=8)                   :: dp_scale_l,dp_scale_d,dp_scale_t,dp_scale_T2,dp_scale_zsun,dp_scale_nh,dp_scale_v

  ! cooling-related stuff -------------------------------------------------------------
  type cooling_table
     integer(kind=4)          :: n11
     integer(kind=4)          :: n22
     real(kind=8),allocatable :: nH(:)    
     real(kind=8),allocatable :: T2(:)    
     real(kind=8),allocatable :: metal(:,:)  
     real(kind=8),allocatable :: cool(:,:)  
     real(kind=8),allocatable :: heat(:,:)  
     real(kind=8),allocatable :: cool_com(:,:)  
     real(kind=8),allocatable :: heat_com(:,:)  
     real(kind=8),allocatable :: cool_com_prime(:,:)  
     real(kind=8),allocatable :: heat_com_prime(:,:)  
     real(kind=8),allocatable :: metal_prime(:,:)  
     real(kind=8),allocatable :: cool_prime(:,:)  
     real(kind=8),allocatable :: heat_prime(:,:)  
     real(kind=8),allocatable :: mu(:,:)  
     real(kind=8),allocatable :: spec(:,:,:)  ! last dimension (6) is n_e, n_HI, n_HII, n_HeI, n_HeII, n_HeIII
  end type cooling_table
  type(cooling_table) :: cooling
  type cool_interp
     integer(kind=4)  :: n_nH
     real(kind=8)     :: nH_start,nH_step
     integer(kind=4)  :: n_T2
     real(kind=8)     :: T2_start,T2_step
  end type cool_interp
  type(cool_interp)  :: cool_int
  logical            :: cooling_is_read = .False. 
  ! ----------------------------------------------------------------------------------

  public  :: read_leaf_cells, get_ngridtot, ramses_get_velocity_cgs, ramses_get_T_nhi_cgs, ramses_get_metallicity, ramses_get_box_size_cm
  ! default is private now ... !! private :: read_hydro, read_amr, get_nleaf, get_nvar, clear_amr, get_ncpu, get_param_real

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

    print *,' '
    print *,'...reading cells...'
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


  
  subroutine ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,temp,nhi)

    implicit none 

    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: nhi(nleaf), temp(nleaf)
    integer(kind=4) :: ihx,ihy,i
    real(kind=8)    :: xx,yy,dxx1,dxx2,dyy1,dyy2,f
    integer(kind=4) :: if1,if2,jf1,jf2

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if

    if (.not. cooling_is_read) then
       call read_cooling(repository,snapnum)
       cooling_is_read = .True.
    end if
    
    nhi  = ramses_var(1,:) * dp_scale_nh  ! nb of H atoms per cm^3
    temp = ramses_var(5,:) / ramses_var(1,:) * dp_scale_T2  ! T/mu [ K ]
    
    ! compute the ionization state and temperature using the 'cooling' tables
    do i = 1, nleaf 
       xx  = log10(nhi(i))
       ihx = int((xx - cool_int%nh_start)/cool_int%nh_step) + 1
       if (ihx < 1) then 
          ihx = 1 
       else if (ihx > cool_int%n_nh) then
          ihx = cool_int%n_nh
       end if
       yy  = log10(temp(i))
       ihy = int((yy - cool_int%t2_start)/cool_int%t2_step) + 1
       if (ihy < 1) then 
          ihy = 1 
       else if (ihy > cool_int%n_t2) then
          ihy = cool_int%n_t2
       end if
       ! 2D linear interpolation:
       if (ihx < cool_int%n_nh) then 
          dxx1  = max(xx - cooling%nh(ihx),0.0d0) / cool_int%nh_step 
          dxx2  = min(cooling%nh(ihx+1) - xx,cool_int%nh_step) / cool_int%nh_step
          if1  = ihx
          if2  = ihx+1
       else
          dxx1  = 0.0d0
          dxx2  = 1.0d0
          if1  = ihx
          if2  = ihx
       end if
       if (ihy < cool_int%n_t2) then 
          dyy1  = max(yy - cooling%t2(ihy),0.0d0) / cool_int%t2_step
          dyy2  = min(cooling%t2(ihy+1) - yy,cool_int%t2_step) / cool_int%t2_step
          jf1  = ihy
          jf2  = ihy + 1
       else
          dyy1  = 0.0d0
          dyy2  = 1.0d0
          jf1  = ihy
          jf2  = ihy
       end if
       if (abs(dxx1+dxx2-1.0d0) > 1.0d-6 .or. abs(dyy1+dyy2-1.0d0) > 1.0d-6) then 
          write(*,*) 'Fucked up the interpolation ... '
          print*,dxx1+dxx2,dyy1+dyy2
          stop
       end if
       ! neutral H density 
       f = dxx1 * dyy1 * cooling%spec(if2,jf2,2) + dxx2 * dyy1 * cooling%spec(if1,jf2,2) &
            & + dxx1 * dyy2 * cooling%spec(if2,jf1,2) + dxx2 * dyy2 * cooling%spec(if1,jf1,2)
       nhi(i) = 10.0d0**f   ! nHI (cm^-3)
       ! GET MU to convert T/MU into T ... 
       f = dxx1 * dyy1 * cooling%mu(if2,jf2) + dxx2 * dyy1 * cooling%mu(if1,jf2) &
            & + dxx1 * dyy2 * cooling%mu(if2,jf1) + dxx2 * dyy2 * cooling%mu(if1,jf1)
       temp(i) = temp(i) * f   ! This is now T (in K) with no bloody mu ... 
    end do
       
    return

  end subroutine ramses_get_T_nhi_cgs
  
  subroutine ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,velocity_cgs)

    implicit none
    
    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: velocity_cgs(3,nleaf)

    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if
    
    velocity_cgs = ramses_var(2:4,:) * dp_scale_v ! [ cm / s ]
    
    return
    
  end subroutine ramses_get_velocity_cgs

  subroutine ramses_get_metallicity(nleaf,nvar,ramses_var,metallicity)

    implicit none
    
    integer(kind=4),intent(in)  :: nleaf, nvar
    real(kind=8),intent(in)     :: ramses_var(nvar,nleaf) ! one cell only
    real(kind=8),intent(inout)  :: metallicity(nleaf)

    if (nvar < 6) then
       print*,'No metals !!! '
       stop
    end if
    metallicity = ramses_var(6,:) 
    
    return

  end subroutine ramses_get_metallicity

  
  function ramses_get_box_size_cm(repository,snapnum)

    implicit none
    
    character(1000),intent(in)  :: repository
    integer(kind=4),intent(in)  :: snapnum
    real(kind=8)                :: ramses_get_box_size_cm 
    
    ! get conversion factors if necessary
    if (.not. conversion_scales_are_known) then 
       call read_conversion_scales(repository,snapnum)
       conversion_scales_are_known = .True.
    end if
    
    ramses_get_box_size_cm = get_param_real(repository,snapnum,'boxlen') * dp_scale_l  ! [ cm ] 
    
    return
    
  end function ramses_get_box_size_cm

    
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
  
  subroutine read_conversion_scales(repository,snapnum)
    
    implicit none 

    character(512),intent(in)  :: repository
    integer(kind=4),intent(in) :: snapnum

    ! set global variables 
    dp_scale_l    = get_param_real(repository,snapnum,'unit_l')
    dp_scale_d    = get_param_real(repository,snapnum,'unit_d')
    dp_scale_t    = get_param_real(repository,snapnum,'unit_t')
    dp_scale_nH   = XH/mp * dp_scale_d      ! convert mass density (code units) to numerical density of H atoms
    dp_scale_v    = dp_scale_l/dp_scale_t   ! -> converts velocities into cm/s
    dp_scale_T2   = mp/kB * dp_scale_v**2   ! -> converts P/rho to T/mu, in K
    dp_scale_zsun = 1.d0/0.0127

    return

  end subroutine read_conversion_scales
  
!*****************************************************************************************************************

  subroutine read_cooling(repository,snapnum)

    implicit none
    
    character(1000),intent(in) :: repository
    integer(kind=4),intent(in) :: snapnum
    character(1024)            :: filename
    integer(kind=4)            :: n1,n2

    ! initialize cooling variables
    call clear_cooling

    ! read cooling variables from cooling.out file
    write(filename,'(a,a,i5.5,a,i5.5,a)') trim(repository),'/output_',snapnum,"/cooling_",snapnum,".out"
    open(unit=44,file=filename,form='unformatted')
    read(44) n1,n2
    cooling%n11 = n1
    cooling%n22 = n2
    allocate(cooling%nH(n1),cooling%T2(n2))
    allocate(cooling%cool_com(n1,n2),cooling%heat_com(n1,n2),cooling%cool_com_prime(n1,n2),cooling%heat_com_prime(n1,n2))
    allocate(cooling%cool(n1,n2),cooling%heat(n1,n2),cooling%mu(n1,n2))
    allocate(cooling%cool_prime(n1,n2),cooling%heat_prime(n1,n2),cooling%metal_prime(n1,n2))
    allocate(cooling%metal(n1,n2),cooling%spec(n1,n2,6))
    read(44)cooling%nH
    read(44)cooling%T2
    read(44)cooling%cool
    read(44)cooling%heat
    read(44)cooling%cool_com
    read(44)cooling%heat_com
    read(44)cooling%metal
    read(44)cooling%cool_prime
    read(44)cooling%heat_prime
    read(44)cooling%cool_com_prime
    read(44)cooling%heat_com_prime
    read(44)cooling%metal_prime
    read(44)cooling%mu
    read(44)cooling%spec
    close(44)
    
    ! define useful quantities for interpolation 
    cool_int%n_nh     = n1
    cool_int%nh_start = minval(cooling%nh)
    cool_int%nh_step  = cooling%nh(2) - cooling%nh(1)
    cool_int%n_t2     = n2
    cool_int%t2_start = minval(cooling%t2)
    cool_int%t2_step  = cooling%t2(2) - cooling%t2(1)
    
    return

  end subroutine read_cooling

!*****************************************************************************************************************

  subroutine clear_cooling

    implicit none

    if (cooling%n11 > 0 .or. cooling%n22 > 0) then 
       cooling%n11 = 0
       cooling%n22 = 0
       deallocate(cooling%nH,cooling%T2,cooling%metal)
       deallocate(cooling%heat_com,cooling%cool_com,cooling%heat_com_prime,cooling%cool_com_prime)
       deallocate(cooling%cool,cooling%heat,cooling%metal_prime,cooling%cool_prime)
       deallocate(cooling%heat_prime,cooling%mu,cooling%spec)
    end if
    
    cool_int%n_nh = 0
    cool_int%nh_start = 0.0d0
    cool_int%nh_step  = 0.0d0
    cool_int%n_t2 = 0
    cool_int%t2_start = 0.0d0
    cool_int%t2_step  = 0.0d0
    
    return

  end subroutine clear_cooling

    
end module module_ramses
!==================================================================================
!==================================================================================
