module module_domain

  implicit none

  public
  
  !JB: les types ci-dessous pourraient etre prives (et meme dans des modules) ?
  type shell
     real(kind=8),dimension(3) :: center
     real(kind=8)              :: r_inbound,r_outbound
  end type shell

  type cube
     real(kind=8),dimension(3) :: center
     real(kind=8)              :: size            ! convention: size is the full size of the cube, corners are x+-size/2
  end type cube

  type sphere
     real(kind=8),dimension(3) :: center
     real(kind=8)              :: radius
  end type sphere

  type slab                                   ! infinite slab in xy direction (make use of periodic boundaries)
     real(kind=8)              :: zc          ! position of the slab in the z direction
     real(kind=8)              :: thickness   ! thickness of the slab in the z direction 
  end type slab

  type domain
     character(10) :: type  ! one of the possible shapes ('shell', 'cube', 'sphere', 'slab')
     type(shell)   :: sh
     type(cube)    :: cu
     type(sphere)  :: sp
     type(slab)    :: sl
  end type domain
  

contains
  
  ! ================
  ! public (tout...) 
  ! ================

  ! constructeurs
  ! -------------
  subroutine domain_constructor_from_scratch(dom,type,xc,yc,zc,r,r_inbound,r_outbound,size,thickness)
    character(10),intent(in)         :: type
    real(kind=8),intent(in),optional :: xc,yc,zc
    real(kind=8),intent(in),optional :: r                     ! parameters for sphere
    real(kind=8),intent(in),optional :: r_inbound,r_outbound  ! parameters for shell
    real(kind=8),intent(in),optional :: size                  ! parameters for cube
    real(kind=8),intent(in),optional :: thickness             ! parameters for slab
    type(domain),intent(out)         :: dom
    logical :: ok

    select case(type)

    case('sphere')
       ! check if optional argument required for sphere are present
       ok = present(xc).and.present(yc).and.present(zc).and.present(r)
       if (.not.ok) then
          print *,'arguments to construct a sphere domain are missing...'
          stop
       endif
       dom%type=type
       dom%sp%center(1)=xc
       dom%sp%center(2)=yc
       dom%sp%center(3)=zc
       dom%sp%radius=r

    case('shell')
       ! check if optional argument required for sphere are present
       ok = present(xc).and.present(yc).and.present(zc).and.present(r_inbound).and.present(r_outbound)
       if (.not.ok) then
          print *,'arguments to construct a shell domain are missing...'
          stop
       endif
       dom%type=type
       dom%sh%center(1)=xc
       dom%sh%center(2)=yc
       dom%sh%center(3)=zc
       dom%sh%r_inbound=r_inbound
       dom%sh%r_outbound=r_outbound

    case('cube')
       ! check if optional argument required for sphere are present
       ok = present(xc).and.present(yc).and.present(zc).and.present(size)
       if (.not.ok) then
          print *,'arguments to construct a cube domain are missing...'
          stop
       endif
       dom%type=type
       dom%cu%center(1)=xc
       dom%cu%center(2)=yc
       dom%cu%center(3)=zc
       dom%cu%size=size


    case('slab')
       ! check if optional argument required for sphere are present
       ok = present(zc).and.present(thickness)
       if (.not.ok) then
          print *,'arguments to construct a slab domain are missing...'
          stop
       endif
       dom%type=type
       dom%sl%zc=zc
       dom%sl%thickness=thickness

    case default
       print *,'type not defined',type
       stop
    end select

    return
  end subroutine domain_constructor_from_scratch


  
  subroutine domain_constructor_from_file(filename,dom)
    ! read a domain file (filename) and initialise domain dom from it
    character(2000),intent(in) :: filename
    type(domain),intent(inout) :: dom

    open(unit=13, file=trim(filename), status='old', form='formatted', action='read')

    call read_domain(unit=13,dom=dom)

    close(13)

  end subroutine domain_constructor_from_file




  
  ! fonctions utilitaires
  ! ---------------------


  subroutine select_in_domain(dom,n,xp,indsel)
    integer,intent(in)                           :: n
    type(domain),intent(in)                      :: dom 
    real(kind=8),dimension(1:n,1:3),intent(in)   :: xp
    integer,dimension(:),allocatable,intent(out) :: indsel
    integer                                      :: i,ii,nsel
    integer,dimension(:),allocatable             :: tmpi
    real(kind=8)                                 :: dd
    
    select case(dom%type)
    
    case('sphere')
       allocate(indsel(1:n))
       indsel=0
       ii=0
       do i=1,n
          dd = (xp(i,1)-dom%sp%center(1))**2 + (xp(i,2)-dom%sp%center(2))**2 + (xp(i,3)-dom%sp%center(3))**2
          if(sqrt(dd)<=dom%sp%radius)then
             ii=ii+1
             indsel(ii)=i
          endif
       enddo
       nsel=ii
       allocate(tmpi(1:nsel))
       tmpi(1:nsel) = indsel(1:nsel)
       deallocate(indsel)
       allocate(indsel(1:nsel))
       indsel=tmpi
       deallocate(tmpi)

    case('shell')
       allocate(indsel(1:n))
       indsel=0
       ii=0
       do i=1,n
          dd = (xp(i,1)-dom%sh%center(1))**2 + (xp(i,2)-dom%sh%center(2))**2 + (xp(i,3)-dom%sh%center(3))**2
          if((sqrt(dd)>=dom%sh%r_inbound).and.(sqrt(dd)<dom%sh%r_outbound))then
             ii=ii+1
             indsel(ii)=i
          endif
       enddo
       nsel=ii
       allocate(tmpi(1:nsel))
       tmpi(1:nsel) = indsel(1:nsel)
       deallocate(indsel)
       allocate(indsel(1:nsel))
       indsel=tmpi
       deallocate(tmpi)


    case('cube')
       allocate(indsel(1:n))
       indsel=0
       ii=0
       do i=1,n
          if((abs(xp(i,1)-dom%cu%center(1)) <= dom%cu%size/2.).and. & 
             (abs(xp(i,2)-dom%cu%center(2)) <= dom%cu%size/2.).and. &
             (abs(xp(i,3)-dom%cu%center(3)) <= dom%cu%size/2.))then
             ii=ii+1
             indsel(ii)=i
          endif
       enddo
       nsel=ii
       allocate(tmpi(1:nsel))
       tmpi(1:nsel) = indsel(1:nsel)
       deallocate(indsel)
       allocate(indsel(1:nsel))
       indsel=tmpi
       deallocate(tmpi)


    case('slab')
       allocate(indsel(1:n))
       indsel=0
       ii=0
       do i=1,n
          if(abs(xp(i,3)-dom%sl%zc) <= dom%sl%thickness/2.)then
             ii=ii+1
             indsel(ii)=i
          endif
       enddo
       nsel=ii
       allocate(tmpi(1:nsel))
       tmpi(1:nsel) = indsel(1:nsel)
       deallocate(indsel)
       allocate(indsel(1:nsel))
       indsel=tmpi
       deallocate(tmpi)


    case default
       print *,'type not defined',dom%type
       stop
    end select

  end subroutine select_in_domain



  subroutine read_domain(unit,dom)
    integer,intent(in)       :: unit
    type(domain),intent(out) :: dom 

    read(unit,*) dom%type
    select case(dom%type)
    case('sphere')
       read(unit,*) dom%sp%center(:)
       read(unit,*) dom%sp%radius
    case('shell')
       read(unit,*) dom%sh%center(:)
       read(unit,*) dom%sh%r_inbound
       read(unit,*) dom%sh%r_outbound
    case('cube')
       read(unit,*) dom%cu%center(:)
       read(unit,*) dom%cu%size
    case('slab')
       read(unit,*) dom%sl%zc
       read(unit,*) dom%sl%thickness
    case default
       print *,'type not defined',dom%type
       stop
    end select

  end subroutine read_domain


  subroutine read_domain_bin(unit,dom)
    integer,intent(in)       :: unit
    type(domain),intent(out) :: dom 

    read(unit) dom%type
    select case(dom%type)
    case('sphere')
       read(unit) dom%sp%center(:)
       read(unit) dom%sp%radius
    case('shell')
       read(unit) dom%sh%center(:)
       read(unit) dom%sh%r_inbound
       read(unit) dom%sh%r_outbound
    case('cube')
       read(unit) dom%cu%center(:)
       read(unit) dom%cu%size
    case('slab')
       read(unit) dom%sl%zc
       read(unit) dom%sl%thickness
    case default
       print *,'type not defined',dom%type
       stop
    end select

  end subroutine read_domain_bin


  subroutine dump_domain_bin(unit,dom)
    integer,intent(in)      :: unit
    type(domain),intent(in) :: dom 

    write(unit) dom%type
    select case(dom%type)
    case('sphere')
       write(unit) dom%sp%center(:)
       write(unit) dom%sp%radius
    case('shell')
       write(unit) dom%sh%center(:)
       write(unit) dom%sh%r_inbound
       write(unit) dom%sh%r_outbound
    case('cube')
       write(unit) dom%cu%center(:)
       write(unit) dom%cu%size
    case('slab')
       write(unit) dom%sl%zc
       write(unit) dom%sl%thickness
    case default 
       print *,'type not defined',dom%type
       stop
    end select

  end subroutine dump_domain_bin

  
  
  subroutine dump_domain_form(unit,dom)
    integer,intent(in)      :: unit
    type(domain),intent(in) :: dom 

    write(unit,'(a)') trim(dom%type)
    select case(dom%type)
    case('sphere')
       write(unit,'(3(f14.10))') dom%sp%center(:)
       write(unit,'(f14.10)') dom%sp%radius
    case('shell')
       write(unit,'(3(f14.10))') dom%sh%center(:)
       write(unit,'(f14.10)') dom%sh%r_inbound
       write(unit,'(f14.10)') dom%sh%r_outbound
    case('cube')
       write(unit,'(3(f14.10))') dom%cu%center(:)
       write(unit,'(f14.10)') dom%cu%size
    case('slab')
       write(unit,'(f14.10)') dom%sl%zc
       write(unit,'(f14.10)') dom%sl%thickness
    case default 
       print *,'type not defined',dom%type
       stop
    end select

  end subroutine dump_domain_form



  subroutine domain_write_file(file,dom)
    character(*),intent(in) :: file
    type(domain),intent(in)    :: dom 

    !!!write(file ,'(a,a)') trim(fichier),'.dom'
    open(unit=14, file=trim(file), status='unknown', form='formatted', action='write')
    call dump_domain_form(unit=14,dom=dom)
    close(14)

  end subroutine domain_write_file




  function domain_contains_point(x,dom)
    ! -> returns T/F if point xyz is in domain dom.
    type(domain),intent(in)              :: dom
    real(kind=8),dimension(3),intent(in) :: x
    logical                              :: domain_contains_point
    real(kind=8)                         :: rr
    domain_contains_point=.false.
    select case(dom%type)
    case('sphere')
       rr = sqrt((x(1)-dom%sp%center(1))**2 + (x(2)-dom%sp%center(2))**2 + (x(3)-dom%sp%center(3))**2)
       if(rr<=dom%sp%radius)domain_contains_point=.true.
    case('shell')
       rr = sqrt((x(1)-dom%sh%center(1))**2 + (x(2)-dom%sh%center(2))**2 + (x(3)-dom%sh%center(3))**2)
       if((rr>=dom%sh%r_inbound).and.(rr<dom%sh%r_outbound))domain_contains_point=.true.
    case('cube')
       if((abs(x(1)-dom%cu%center(1)) <= dom%cu%size/2.).and. & 
            (abs(x(2)-dom%cu%center(2)) <= dom%cu%size/2.).and. &
            (abs(x(3)-dom%cu%center(3)) <= dom%cu%size/2.))domain_contains_point=.true.
    case('slab')
       if(abs(x(3)-dom%sl%zc) <= dom%sl%thickness/2.)domain_contains_point=.true.
    end select
    return
  end function domain_contains_point

  

  function domain_contains_cell(x,dx,dom)
    ! -> as above with dx tolerance ?
    ! convention: dx is the cell size
    type(domain),intent(in)              :: dom
    real(kind=8),dimension(3),intent(in) :: x
    real(kind=8),intent(in)              :: dx
    logical                              :: domain_contains_cell
    real(kind=8)                         :: rr

    domain_contains_cell=.false.

    select case(dom%type)

    case('sphere')
       rr = sqrt((x(1)-dom%sp%center(1))**2 + (x(2)-dom%sp%center(2))**2 + (x(3)-dom%sp%center(3))**2)
       if((rr+dx/sqrt(2.))<=dom%sp%radius)domain_contains_cell=.true.

    case('shell')
       rr = sqrt((x(1)-dom%sh%center(1))**2 + (x(2)-dom%sh%center(2))**2 + (x(3)-dom%sh%center(3))**2)
       if(((rr-dx/sqrt(2.))>=dom%sh%r_inbound).and.((rr+dx/sqrt(2.))<dom%sh%r_outbound))domain_contains_cell=.true.

    case('cube')
       if((x(1)+dx/2. <= dom%cu%center(1)+dom%cu%size/2.).and. &
          (x(1)-dx/2. <= dom%cu%center(1)-dom%cu%size/2.).and. &
          (x(2)+dx/2. <= dom%cu%center(2)+dom%cu%size/2.).and. &
          (x(2)-dx/2. <= dom%cu%center(2)-dom%cu%size/2.).and. &
          (x(3)+dx/2. <= dom%cu%center(3)+dom%cu%size/2.).and. &
          (x(3)-dx/2. <= dom%cu%center(3)-dom%cu%size/2.)) domain_contains_cell=.true.

    case('slab')
       if((x(3)+dx/2. <= dom%sl%zc+dom%sl%thickness/2.).and. &
          (x(3)-dx/2. <= dom%sl%zc-dom%sl%thickness/2.)) domain_contains_cell=.true.

    end select

    return
  end function domain_contains_cell



  function get_my_new_domain(x,liste_domaines)
    !-> given position xyz of a point, returns the dom where the point is
    ! in case of overlapping domains, should use domain_distance_to_border to choose
    ! how to do that efficiently?
    ! scan each domain and test using domain_contains_point ?
    type(domain),intent(in),dimension(:) :: liste_domaines
    real(kind=8),dimension(3),intent(in) :: x
    integer :: get_my_new_domain, ndom, count_dom, first_dom, i
    logical :: x_in_i
    real(kind=8) :: d1,d2

    ndom = size(liste_domaines)
    count_dom = 0
    first_dom = -1
    do i=1,ndom
       x_in_i = domain_contains_point(x,liste_domaines(i))
       if(x_in_i)then
          count_dom         = count_dom + 1
          if(count_dom>1) first_dom = get_my_new_domain
          get_my_new_domain = i
       endif
    enddo
    if((count_dom > 2).or.(count_dom==0))then
       print *,'--> Problem with get_my_new_domain'
       stop
    endif
    if(count_dom>1)then
       ! point belongs to 2 domains, choose the one for which distance to border is maximum
       d1 = domain_distance_to_border(x,liste_domaines(first_dom))
       d2 = domain_distance_to_border(x,liste_domaines(get_my_new_domain))
       if(d2<d1)then
          get_my_new_domain = first_dom
       endif
       ! else (d2>=d1) then get_my_new_domain is ok
    endif

  end function get_my_new_domain


  
  function domain_distance_to_border(x,dom)
    ! return distance of point xyz to the closest border of domain dom
    real(kind=8),dimension(3),intent(in) :: x
    type(domain),intent(in)              :: dom
    real(kind=8)                         :: domain_distance_to_border, rr, ddx, ddy, ddz

    select case(dom%type)

    case('sphere')
       rr = sqrt((x(1)-dom%sp%center(1))**2 + (x(2)-dom%sp%center(2))**2 + (x(3)-dom%sp%center(3))**2)
       domain_distance_to_border = dom%sp%radius - rr

    case('shell')
       rr = sqrt((x(1)-dom%sh%center(1))**2 + (x(2)-dom%sh%center(2))**2 + (x(3)-dom%sh%center(3))**2)
       domain_distance_to_border = min((rr-dom%sh%r_inbound),(dom%sh%r_outbound-rr))

    case('cube')
       ddx = min((dom%cu%center(1)+dom%cu%size/2.-x(1)), (x(1)-dom%cu%center(1)-dom%cu%size/2.))
       ddy = min((dom%cu%center(2)+dom%cu%size/2.-x(2)), (x(2)-dom%cu%center(2)-dom%cu%size/2.))
       ddz = min((dom%cu%center(3)+dom%cu%size/2.-x(3)), (x(3)-dom%cu%center(3)-dom%cu%size/2.))
       domain_distance_to_border = min(ddx,ddy,ddz)
    
    case('slab')
       domain_distance_to_border = min((dom%sl%zc+dom%sl%thickness/2.-x(3)), (x(3)-dom%sl%zc-dom%sl%thickness)) 

    end select

  end function domain_distance_to_border

  
  
end module module_domain
