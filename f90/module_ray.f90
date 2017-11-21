module module_ray

  ! This module assumes a single domain and will not work with domain decomposition... 
  ! This is a trimmed and slightly modified version of module_photon.
  
  use module_gas_composition
  use module_mesh
  use module_constants
  use module_random
  use module_domain
  use module_utils, only: path 

  implicit none

  ! todonext define accuracy
  real(kind=8),parameter :: accuracy=1.d-15

  type ray_type
     integer                   :: ID       ! a positive unique ID 
     real(kind=8)              :: nu_ext   ! external frame frequency (Hz)
     real(kind=8)              :: dist     ! distance traveled along ray (box units)
     real(kind=8)              :: tau      ! integrated opacity along ray
     real(kind=8),dimension(3) :: x_em     ! emission location (box units)
     real(kind=8),dimension(3) :: k_em     ! emission direction == propagation direction (normalised vector)
  end type ray_type

  public  :: RayTracing, ray_advance, init_rays_from_file, dump_rays
  private :: path

contains

  subroutine RayTracing(nrays,rays,mesh_dom,compute_dom,maxdist,maxtau)

    integer, intent(in)                             :: nrays
    type(ray_type), dimension(nrays), intent(inout) :: rays
    type(mesh), intent(in)                          :: mesh_dom
    type(domain), intent(in)                        :: compute_dom
    real(kind=8),intent(in)                         :: maxdist,maxtau
    integer                                         :: i

    do i=1,nrays
       call ray_advance(rays(i),mesh_dom,compute_dom,maxdist,maxtau)
    enddo

  end subroutine RayTracing



  subroutine ray_advance(ray,domesh,domaine_calcul,maxdist,maxtau)

    type(ray_type),intent(inout)   :: ray            ! a ray 
    type(mesh),intent(in)          :: domesh         ! mesh
    type(domain),intent(in)        :: domaine_calcul ! domaine dans lequel on propage les photons...
    real(kind=8),intent(in)        :: maxdist,maxtau ! stop propagation at either maxdist or maxtau (the one which is positive). 
    type(gas)                      :: cell_gas       ! gas in the current cell 
    integer(kind=4)                :: icell, ioct, ind, ileaf, cell_level  ! current cell indices and level
    real(kind=8)                   :: cell_size, cell_size_cm, scalar, nu_cell, maxdist_cm
    real(kind=8),dimension(3)      :: ppos,ppos_cell ! working coordinates of photon (in box and in cell units)
    real(kind=8)                   :: distance_to_border,distance_to_border_cm
    real(kind=8)                   :: dist, tau_cell, tau, dborder
    integer(kind=4)                :: i, icellnew, npush
    real(kind=8),dimension(3)      :: vgas, kray, cell_corner, posoct
    logical                        :: flagoutvol, in_domain
    real(kind=8)                   :: epsilon_cell
    
    ! initialise ray tracing 
    ppos  = ray%x_em   ! emission position 
    kray  = ray%k_em   ! propagation direction
    dist  = 0.0d0      ! distance covered
    tau   = 0.0d0      ! corresponding optical depth
    maxdist_cm = maxdist * box_size_cm 
    
    ! find the (leaf) cell in which the photon is, and define all its indices
    icell = in_cell_finder(domesh,ppos)
    if(domesh%son(icell)>=0)then
       print*,'ERROR: not a leaf cell'
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
          print*,"ERROR: problem in computing ppos_cell"
          stop
       endif
       
       ! get gas velocity (in cgs units)
       vgas         = get_gas_velocity(cell_gas)
       ! compute photon's frequency in cell's moving frame
       scalar       = kray(1) * vgas(1) + kray(2) * vgas(2) + kray(3) * vgas(3)
       nu_cell      = (1.d0 - scalar/clight) * ray%nu_ext  

       ! compute distance of photon to border of cell or domain along propagation direction
       distance_to_border    = path(ppos_cell,kray) / cell_size            ! in box units
       distance_to_border    = min(distance_to_border, &
            & domain_distance_to_border_along_k(ppos,kray,domaine_calcul)) ! in box units
       distance_to_border_cm = distance_to_border * box_size_cm ! cm
       
       ! compute (total) optical depth along ray in cell 
       tau_cell = gas_get_tau(cell_gas, distance_to_border_cm, nu_cell)
       
       ! update head of ray position
       ppos = ppos + kray * distance_to_border *(1.0d0 + epsilon(1.0d0))
       ! correct for periodicity
       do i=1,3
          if (ppos(i) < 0.0d0) ppos(i)=ppos(i)+1.0d0
          if (ppos(i) > 1.0d0) ppos(i)=ppos(i)-1.0d0
       enddo
       
       ! update traveled distance and optical depth
       dist = dist + distance_to_border_cm
       tau  = tau + tau_cell 
       
       ! check if photon still in computational domain
       in_domain = domain_contains_point(ppos,domaine_calcul)
       if (.not.(in_domain)) exit ray_propagation  ! ray is done 
       ! check if we reached tau or distance limits
       if (dist > maxdist_cm .and. tau > maxtau) then ! dist or tau exceeding boundary -> correct excess and exit. 
          if (maxdist_cm > 0) then
             dborder  = maxdist_cm - dist
             ray%dist = maxdist_cm
             ray%tau  = tau + (dborder / distance_to_border_cm)*tau_cell
          else
             dborder  = maxtau - tau
             ray%dist = dist + dborder / tau_cell * distance_to_border_cm 
             ray%tau  = maxtau
          end if
          exit ray_propagation
       end if
       
       ! Ray moves to next cell : find it
       call whereIsPhotonGoing(domesh,icell,ppos,icellnew,flagoutvol)
       ! It may happen due to numerical precision that the photon is still in the current cell (i.e. icell == icellnew).
       ! -> give it an extra push untill it is out. 
       npush = 0
       do while (icell==icellnew)
          npush = npush + 1
          ppos(1) = ppos(1) + merge(-1.0d0,1.0d0,kray(1)<0.0d0) * epsilon(ppos(1))
          ppos(2) = ppos(2) + merge(-1.0d0,1.0d0,kray(2)<0.0d0) * epsilon(ppos(2))
          ppos(3) = ppos(3) + merge(-1.0d0,1.0d0,kray(3)<0.0d0) * epsilon(ppos(3))
          call whereIsPhotonGoing(domesh,icell,ppos,icellnew,flagoutvol)
       end do
       if (npush > 1) print*,'WARNING : npush > 1 needed in module_photon:propagate.'
       ! test whether photon was pushed out of domain with the extra pushes
       ! (and in that case, call it done). 
       if (npush > 0) then 
          in_domain = domain_contains_point(ppos,domaine_calcul)
          if (.not. in_domain) then
             print*,'WARNING: pushed photon outside domain ... '
             exit ray_propagation
          end if
       end if
       
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
    read(14) (rays(i)%ID,i=1,n_rays)
    read(14) (rays(i)%nu_ext,i=1,n_rays)
    read(14) (rays(i)%x_em(:),i=1,n_rays)
    read(14) (rays(i)%k_em(:),i=1,n_rays)
    close(14)
    ! initialise other properties. 
    do i=1,n_rays
       rays(i)%tau  = 0.d0
       rays(i)%dist = 0.d0
    enddo
    
  end subroutine init_rays_from_file


  subroutine dump_rays(file,rays)

    character(2000),intent(in)             :: file
    type(ray_type),dimension(:),intent(in) :: rays
    integer(kind=4)                        :: i,np

    np = size(rays)
    open(unit=14, file=trim(file), status='unknown', form='unformatted', action='write')
    write(14) np
    write(14) (rays(i)%ID,i=1,np)
    write(14) (rays(i)%dist,i=1,np)
    write(14) (rays(i)%tau,i=1,np)
    close(14)

  end subroutine dump_rays

end module module_ray
