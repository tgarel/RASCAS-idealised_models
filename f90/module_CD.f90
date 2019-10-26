module module_CD

  ! This is a trimmed and slightly modified version of module_photon. It is quite redondant with module_opt_ray.f90,  have to put order in all that

  use module_gas_composition
  use module_mesh
  use module_constants
  use module_random
  use module_domain
  use module_utils, only: path, isotropic_direction

  implicit none

  type ray_type
     real(kind=8),dimension(3) :: x_em     ! emission location (box units)
     real(kind=8)              :: weight   ! Weight of the photon,  useful for example to implement stellar luminosities.   Not in jbv's version
  end type ray_type

  type histo_v_N
     real(kind=8)              :: vmin, vmax    ! velocities between which we want to plot NSiII vs velocity of gas cell,   in km/s  (be careful,  velocity of gas cell is in cm/s)
     integer(kind=4)           :: nBins         ! number of bins of velocities
     integer(kind=4)           :: nrays
     integer(kind=4)           :: nGas
     integer(kind=4)           :: nDirections
     real(kind=8), allocatable :: list(:,:,:,:)       ! list of size nBins+2, with the sum of the column densities for a given velocity of gas cell. The nBins+2 intervals are (-infinity,vmin), (vmin, vmin + (vmax-vmin)/nBins), (vmin + (vmax-vmin)/nBins,  vmin + 2*(vmax-vmin)/nBins), .... , (vmin + (nBins-1)*(vmax-vmin)/nBins, vmax), (vmax, infinity)
  end type histo_v_N

  type(histo_v_N) :: histo
  logical         :: do_histo = .false.

  public  :: ComputeCD, ray_advance, init_rays_from_file, dump, init_histo_v_N
  private :: path

contains

  subroutine ComputeCD(nrays,rays,mesh_dom,compute_dom,nDirections,kdir)

    integer(kind=4), intent(in)                         :: nrays,nDirections
    real(kind=8), dimension(nDirections,3), intent(in)  :: kdir
    type(ray_type), dimension(nrays), intent(inout)     :: rays
    type(mesh) ,intent(in)                              :: mesh_dom
    type(domain), intent(in)                            :: compute_dom
    integer(kind=4)                                     :: i,j
  
    !$OMP PARALLEL &
    !$OMP DEFAULT(shared)           
    !$OMP DO SCHEDULE(DYNAMIC, 100)
    do i=1,nrays 
       do j=1,nDirections
          call ray_advance(rays(i),mesh_dom,compute_dom,kdir(j,:),j,i)
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine ComputeCD


  subroutine ray_advance(ray,domesh,domaine_calcul,kray,idir,ID)

    type(ray_type),intent(inout) :: ray
    type(mesh),intent(in)        :: domesh         ! mesh
    type(domain),intent(in)      :: domaine_calcul ! domaine dans lequel on propage les photons...
    integer(kind=4),intent(in)   :: ID, idir
    real(kind=8),intent(in)      :: kray(3)
    type(gas)                    :: cell_gas       ! gas in the current cell 
    integer(kind=4)              :: icell, ioct, ind, ileaf, cell_level  ! current cell indices and level
    real(kind=8)                 :: cell_size, cell_size_cm
    real(kind=8),dimension(3)    :: ppos,ppos_cell ! working coordinates of photon (in box and in cell units)
    real(kind=8)                 :: distance_to_border,distance_to_border_cm
    real(kind=8)                 :: dborder
    integer(kind=4)              :: i, icellnew, npush
    real(kind=8),dimension(3)    :: cell_corner, posoct
    logical                      :: flagoutvol, in_domain, escape_domain_before_cell

    ! initialise ray tracing 
    ppos  = ray%x_em   ! emission position 

    ! check that the ray starts in the domain
    if (.not. domain_contains_point(ppos,domaine_calcul)) then
       print * ,'Ray outside domain at start '
       print*, ppos
       stop
    end if

    ! find the (leaf) cell in which the photon is, and define all its indices
    icell = in_cell_finder(domesh,ppos)
    if(domesh%son(icell)>=0)then
       print*,'ERROR: not a leaf cell ',ppos
       stop
    endif
    ileaf = - domesh%son(icell)
    ind   = (icell - domesh%nCoarse - 1) / domesh%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
    ioct  = icell - domesh%nCoarse - (ind - 1) * domesh%nOct

    ! advance ray until it escapes the computational domain ... 
    ray_propagation : do

       ! gather properties properties of current cell
       cell_level   = domesh%octlevel(ioct)      ! level of current cell
       cell_size    = 0.5d0**cell_level          ! size of current cell in box units
       cell_size_cm = cell_size * box_size_cm    ! size of the current cell in cm
       cell_gas     = domesh%gas(ileaf)  

       ! compute position of photon in current-cell units
       posoct(:)    = domesh%xoct(ioct,:)
       cell_corner  = get_cell_corner(posoct,ind,cell_level)   ! position of cell corner, in box units.
       ppos_cell    = (ppos - cell_corner) / cell_size         ! position of photon in cell units (x,y,z in [0,1] within cell)
       if((ppos_cell(1)>1.0d0).or.(ppos_cell(2)>1.0d0).or.(ppos_cell(3)>1.0d0).or. &
            (ppos_cell(1)<0.0d0).or.(ppos_cell(2)<0.0d0).or.(ppos_cell(3)<0.0d0))then
          print*,"ERROR: problem in computing ppos_cell ",ppos_cell(1),ppos_cell(2),ppos_cell(3)
          stop
       endif

       ! compute distance of photon to border of cell or domain along propagation direction
       distance_to_border    = path(ppos_cell,kray) * cell_size            ! in box units

       ! Check if ray escapes domain before cell 
       escape_domain_before_cell = .false.
       dborder = domain_distance_to_border_along_k(ppos,kray,domaine_calcul)
       if (dborder < distance_to_border) then
          escape_domain_before_cell = .true.
          distance_to_border = dborder
       end if
       distance_to_border_cm = distance_to_border * box_size_cm ! cm

       !Increment histo_v_N
       if(do_histo) then
          i = min( max( floor( ( dot_product( -1d-5*cell_gas%v(:),kray(:) ) - histo%vmin )*histo%nBins/(histo%vmax-histo%vmin) ) + 2, 1 ), histo%nBins+2 )  !computes which index of the list should be incremented
          histo%list(:,idir,ID,i) = histo%list(:,idir,ID,i) + gas_get_CD(cell_gas, distance_to_border_cm)  !increment the corresponding entry
       end if
       !print*,dot_product( -1d-5*cell_gas%v(:),kray(:) ), log10(gas_get_CD(cell_gas, distance_to_border_cm))
       !print*, gas_get_CD(cell_gas, distance_to_border_cm), cell_gas%dopwidth**2 * 15.999 * amu / 2 / kb

       ! update head of ray position
       ppos = ppos + kray * distance_to_border *(1.0d0 + epsilon(1.0d0))
       ! correct for periodicity
       do i=1,3
          if (ppos(i) < 0.0d0) ppos(i)=ppos(i)+1.0d0
          if (ppos(i) > 1.0d0) ppos(i)=ppos(i)-1.0d0
       enddo


       if (escape_domain_before_cell) then ! ray escapes domain -> we are done.
          exit ray_propagation  ! ray is done 
       end if

       ! check if photon still in computational domain after position update
       in_domain = domain_contains_point(ppos,domaine_calcul)
       if (.not.(in_domain)) then 
          exit ray_propagation  ! ray is done 
       end if

       ! Ray moves to next cell : find it
       call whereIsPhotonGoing(domesh,icell,ppos,icellnew,flagoutvol)
       ! It may happen due to numerical precision that the photon is still in the current cell (i.e. icell == icellnew).
       ! -> give it an extra push untill it is out. 
       npush = 0
       do while (icell==icellnew)
          npush = npush + 1
          if (npush>100) then
             print*,ppos
             print*,kray
             print*,'Too many pushes, npush>100 '
             stop
          endif
          ppos(1) = ppos(1) + merge(-1.0d0,1.0d0,kray(1)<0.0d0) * epsilon(ppos(1))
          ppos(2) = ppos(2) + merge(-1.0d0,1.0d0,kray(2)<0.0d0) * epsilon(ppos(2))
          ppos(3) = ppos(3) + merge(-1.0d0,1.0d0,kray(3)<0.0d0) * epsilon(ppos(3))

          ! correct for periodicity
          do i=1,3
             if (ppos(i) < 0.0d0) ppos(i)=ppos(i)+1.0d0
             if (ppos(i) > 1.0d0) ppos(i)=ppos(i)-1.0d0
          enddo

          ! test that we are still in domain before calling WhereIsPhotonGoing... 
          in_domain = domain_contains_point(ppos,domaine_calcul)
          if (.not.(in_domain)) then
             exit ray_propagation  ! ray is done 
          end if
          call whereIsPhotonGoing(domesh,icell,ppos,icellnew,flagoutvol)
       end do
       if (npush > 1) print*,'WARNING : npush > 1 needed in module_gray_ray:propagate.',npush

       ! check if photon outside of cpu domain (flagoutvol)
       if(flagoutvol)then
          ! photon out of cpu domain -> we have a problem... 
          print*,'ERROR: photon out of CPU domain when it should not ... '
          stop
       endif

       ! else, new cell is in the cpu domain and photon goes to this new cell
       icell = icellnew
       if(domesh%son(icell)>=0)then
          print*,'ERROR: not a leaf cell',icell,flagoutvol
          print*,'This should not happen in module_ray (on single domains). '
          stop
       endif
       ileaf = - domesh%son(icell)
       ind   = (icell - domesh%nCoarse - 1) / domesh%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
       ioct  = icell - domesh%nCoarse - (ind - 1) * domesh%nOct       

    end do ray_propagation

  end subroutine ray_advance


  subroutine init_rays_from_file(file,rays)

    character(2000),intent(in)                           :: file    
    type(ray_type),dimension(:),allocatable, intent(out) :: rays
    integer(kind=4)                                      :: i, n_rays

    ! read ICs
    open(unit=14, file=trim(file), status='unknown', form='unformatted', action='read')
    read(14) n_rays
    allocate(rays(n_rays))
    read(14) (rays(i)%x_em(:), i=1,n_rays)
    read(14) (rays(i)%weight, i=1,n_rays)
    
    close(14)

  end subroutine init_rays_from_file


  subroutine dump(file,rays)

    character(2000),intent(in)             :: file
    type(ray_type),dimension(:),intent(in) :: rays
    integer(kind=4)                        :: i,j,k,l,np
    character(2000)                        :: nomfich

    np = size(rays)

    do i=1,max(1,histo%nGas)
       write(nomfich,'(a,a,i1)') trim(file),'_',i
       open(unit=14, file=nomfich, status='unknown', form='unformatted', action='write')
       write(14) np
       write(14) histo%nDirections
       write(14) histo%nGas
       if(do_histo) then
          write(14) histo%vmin, histo%vmax
          write(14) histo%nBins
          do j=1,histo%nDirections
             write(14) ((histo%list(i,j,k,l), k=1,histo%nrays), l=1,histo%nBins+2)
          end do
       end if
       close(14)
    end do

  end subroutine dump


  subroutine init_histo_v_N(vmin,vmax,nBins,nrays,nDirections)

    real(kind=8),intent(in)    :: vmin, vmax
    integer(kind=4),intent(in) :: nBins,nrays,nDirections
    integer(kind=4)            :: nGas

    histo%vmin = vmin
    histo%vmax = vmax
    histo%nBins = nBins
    histo%nrays = nrays
    histo%nDirections = nDirections
    nGas = gas_get_n_CD()
    histo%nGas = nGas
    allocate(histo%list(nGas,nDirections,nrays,nBins+2))
    histo%list = 0d0

    do_histo = .true.

  end subroutine init_histo_v_N

end module module_CD
