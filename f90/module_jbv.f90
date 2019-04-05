module module_jbv

  ! This is a trimmed and slightly modified version of module_photon.
    
  use module_gas_composition
  use module_mesh
  use module_constants
  use module_random
  use module_domain
  use module_utils, only: path, isotropic_direction

  implicit none

  type ray_type
     real(kind=8),dimension(3) :: x_em     ! emission location (box units)
     real(kind=8),dimension(3) :: k_em     ! emission direction == propagation direction (normalised vector)
     real(kind=8)              :: NHI      ! HI column density [cm-2]
     real(kind=8)              :: dist     ! distance traveled [cm]
  end type ray_type


  public  :: ComputeNHI, ray_advance, init_rays_from_file, dump_rays
  private :: path

contains

  subroutine ComputeNHI(nrays,rays,mesh_dom,compute_dom)

    integer(kind=4), intent(in)                     :: nrays
    type(ray_type), dimension(nrays), intent(inout) :: rays
    type(mesh), intent(in)                          :: mesh_dom
    type(domain), intent(in)                        :: compute_dom
    integer(kind=4)                                 :: i
    real(kind=8)                                    :: x_em(3),k_em(3), NHI, dist, vide

    vide = 20000.0d0 * 3.0857d18  ! in cm
    vide = vide / box_size_cm ! in code units 

!$OMP PARALLEL &
!$OMP DEFAULT(shared) &
!$OMP PRIVATE(NHI, dist, k_em, x_em)               
!$OMP DO SCHEDULE(DYNAMIC, 100) 
    do i=1,nrays 
       ! initialisation 
       x_em(:) = rays(i)%x_em(:) + vide * rays(i)%x_em(:)  ! start vide away from star
       k_em(:) = rays(i)%k_em(:)
       NHI     = 0.0d0
       dist    = 0.0d0 
       call ray_advance(x_em,k_em,dist,NHI,mesh_dom,compute_dom)
       rays(i)%nhi  = NHI
       rays(i)%dist = dist
    end do
!$OMP END DO
!$OMP END PARALLEL

  end subroutine ComputeNHI


  subroutine ray_advance(x_em,k_em,dist,NHI,domesh,domaine_calcul)

    real(kind=8),intent(in)   :: x_em(3),k_em(3)
    real(kind=8),intent(out)  :: dist, NHI
    type(mesh),intent(in)     :: domesh         ! mesh
    type(domain),intent(in)   :: domaine_calcul ! domaine dans lequel on propage les photons...
    type(gas)                 :: cell_gas       ! gas in the current cell 
    integer(kind=4)           :: icell, ioct, ind, ileaf, cell_level  ! current cell indices and level
    real(kind=8)              :: cell_size, cell_size_cm
    real(kind=8),dimension(3) :: ppos,ppos_cell ! working coordinates of photon (in box and in cell units)
    real(kind=8)              :: distance_to_border,distance_to_border_cm
    real(kind=8)              :: dborder, vide
    integer(kind=4)           :: i, icellnew, npush
    real(kind=8),dimension(3) :: kray, cell_corner, posoct
    logical                   :: flagoutvol, in_domain, escape_domain_before_cell
    
    ! initialise ray tracing 
    ppos  = x_em   ! emission position 
    kray  = k_em   ! propagation direction
    dist  = 0.0d0  ! distance covered
    NHI   = 0.0d0  !initialise nhtot density along the ray
    
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
       
       ! gather properties of current cell
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
       dborder = domain_distance_to_border_along_k(ppos,kray,domaine_calcul) ! in module_domain

       if (dborder < distance_to_border) then
          escape_domain_before_cell = .true.
          distance_to_border = dborder
       end if
       distance_to_border_cm = distance_to_border * box_size_cm ! cm
       
       
       ! increment column density and traveled distance
       NHI  = NHI  + distance_to_border_cm * cell_gas%nHI 
       dist = dist + distance_to_border_cm

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
    open(unit=14, file=trim(file), status='unknown', form='formatted', action='read')
    read(14,*) n_rays
    allocate(rays(n_rays))
    do i = 1, n_rays
       read(14,*) rays(i)%x_em(1),rays(i)%x_em(2),rays(i)%x_em(3),rays(i)%k_em(1),rays(i)%k_em(2),rays(i)%k_em(3)
    end do
    close(14)
    
  end subroutine init_rays_from_file 


  subroutine dump_rays(file,rays)

    character(2000),intent(in)             :: file
    type(ray_type),dimension(:),intent(in) :: rays
    integer(kind=4)                        :: i,np

    np = size(rays)
    open(unit=14, file=trim(file), status='unknown', form='formatted', action='write')
    write(14,*) np
    do i=1,np
       write(14,'(8(e14.6,1x))') rays(i)%x_em(1),rays(i)%x_em(2),rays(i)%x_em(3),rays(i)%k_em(1), &
       rays(i)%k_em(2),rays(i)%k_em(3),rays(i)%NHI,rays(i)%dist
    end do
    close(14)

  end subroutine dump_rays

  
end module module_jbv
