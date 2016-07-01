module module_mesh

  use module_domain
  use module_gas_composition

  implicit none

  type, public :: mesh
     type(domain) :: domain
     integer :: nCoarse,nOct,nLeaf,nCell
     integer,allocatable,dimension(:)        :: octlevel,son,father
     integer,allocatable,dimension(:,:)      :: nbor
     real(kind=8),dimension(:,:),allocatable :: xoct
     type(gas),allocatable,dimension(:)      :: gas
     !! Prix a payer de ce "format" pour gas, on doit connaitre sa "forme" dans ce module
  end type mesh

  private

  integer, dimension(1:3,1:8),parameter,private :: nbortest = &
       reshape( (/1,3,5, 2,3,5, 1,4,5, 2,4,5, 1,3,6, 2,3,6, 1,4,6, 2,4,6/), (/3,8/) )  

  real(kind=8), dimension(1:3,1:8), parameter :: offset = &
       reshape((/-.5,-.5,-.5, +.5,-.5,-.5, -.5,+.5,-.5, +.5,+.5,-.5, &
                 -.5,-.5,+.5, +.5,-.5,+.5, -.5,+.5,+.5, +.5,+.5,+.5/), (/3,8/))


  ! TODO: change shape of xoct,xleaf,nbor to optimize access xoct(1:noct,1:3) => xoct(1:3,1:noct)

  ! WORKING VARIABLES (interne au module, pourrait un jour etre modifie)
  
  integer :: nCoarse,nOct,nCell,nleaf
  integer :: countleaf,countempty

  ! oct-tree data
  integer,allocatable,dimension(:)        :: octlevel,son,father
  integer,allocatable,dimension(:,:)      :: nbor
  real(kind=8),dimension(:,:),allocatable :: xoct
  !
  real(kind=8),dimension(:,:),allocatable :: xleaf
  integer,dimension(:),allocatable        :: leaflevel
  

  public :: mesh_from_leaves, mesh_from_file, mesh_destructor, dump_mesh, whereisphotongoing, digincell, &
       in_cell_finder, get_cell_corner, overwrite_mesh
  private :: add_oct_domain, make_nbor_array, xcson, get_nleaflocal, get_ileaflocal, icell2icell


  contains

    !<<< List of routines to include here >>>>
    !
    ! mesh constructor
    !------------------
    ! subroutine mesh_from_leaves(nCoarseSnap,nOctSnap,nCellSnap,domain,nleaves,leaves,xleaves,leaveslevel,m)
    ! subroutine mesh_from_file
    !
    ! public use
    !------------
    ! subroutine mesh_destructor
    ! subroutine dump_mesh
    ! subroutine where_is_photon_going(mesh,icellold,xpnew,icellnew,flagoutvol)
    ! function digincell(mesh,xp,icell)
    ! function in_cell_finder(mesh,xp)
    !
    ! private use
    !-------------
    ! recursive subroutine add_oct_domain(ilastoct,nLeaf,ileaf,xnew,lfather,ifathercell)
    ! subroutine make_nbor_array
    ! function xcson(ind,fx,level)
    ! function get_nleaflocal(n,ileaf,xc,level)
    ! function get_ileaflocal(n,ileaf,xc,level,nleaflocal)
    ! function icell2icell(icell,ngridmax,ngridmin,ncoarse)
    ! some routines for checking & verbosity
    !
    !<><><><><><><><><><><><><><><><><><><><><>


    !===============================================================================================
    ! mesh constructor
    !===============================================================================================
    
    subroutine mesh_from_leaves(nOctMax,dom,nleaves,leaves,xleaves,leaveslevel,m)

      type(domain),intent(in)  :: dom
      type(gas),dimension(:),intent(in)     :: leaves
      integer,intent(in)       :: noctmax,nleaves
      real(kind=8),dimension(1:nleaves,1:3),intent(in) :: xleaves
      integer,dimension(1:nleaves),intent(in) :: leaveslevel
      type(mesh),intent(out)   :: m

      integer                                 :: ilastoct,lfather,ifathercell,i
      real(kind=8),dimension(3)               :: xnew
      !integer                                 :: ind,icount,oct_status,icount2,ioct
      integer, dimension(:),allocatable       :: ileaf
      integer :: nocttrue

      ! INIT WORKING VARIABLES
      ncoarse = 1 ! by construction
      noct = nOctMax
      ncell = ncoarse + 8*noct
      nleaf = nleaves

      ! allocate & initialize working arrays
      allocate(xleaf(1:nleaf,1:3),leaflevel(1:nleaf))
      xleaf = xleaves
      leaflevel = leaveslevel
      allocate(ileaf(1:nLeaf))
      allocate(son(1:nCell))
      allocate(father(1:nOct))
      allocate(octlevel(1:nOct))
      allocate(xoct(1:nOct,1:3))
      ! allocate(nbor(1:nOct,1:6)) ! will be done later
      son=0
      father=-1
      octlevel=-1
      xoct=-99.
      !nbor=0
  
      countleaf=0
      countempty=0

      ! premier passage pour la cellule coarse
      ilastoct=0
      xnew=(/0.5,0.5,0.5/)
      lfather=0
      ifathercell=1
      do i=1,nLeaf
         ileaf(i)=i
      enddo
      son(1)=1

      write(*,*)
      write(*,*)'...entering recursive loop...'

      ! build the mesh
      call add_oct_domain(ilastoct,nLeaf,ileaf,xnew,lfather,ifathercell)

      ! resize arrays
      nOctTrue = ilastoct
      call resize_octtree(nOct,nOctTrue)
      nOct  = nOctTrue
      nCell = nCoarse + 8*nOct

      write(*,*)
      write(*,*) 'nCoarse = ',nCoarse
      write(*,*) 'nOct    = ',nOct
      write(*,*) 'nLeaf   = ',nLeaf
      write(*,*) 'nCell   = ',nCell

      ! make nbor array
      write(*,*)
      write(*,*)'...building nbor array...'
      allocate(nbor(nOct,6))
      nbor = 0
      call make_nbor_array

      ! work should be done
      ! == son, father, xoct, octlevel should be filled
      write(*,*)'work done!'
      write(*,*)'ilastoct = ',ilastoct
      write(*,*)'cellules feuilles comptees =',countleaf
      write(*,*)'cellules ~vides~ =',countempty
      ! do some checks if you want
      call check_octtree(noctmax)


      ! fill the mesh type structure
      allocate(m%son(ncell),m%father(noct),m%octlevel(noct),m%xoct(noct,3),m%nbor(noct,6))

      ! allocate & fill m%domain%stuff

      m%domain = dom

      m%ncoarse = ncoarse
      m%noct    = noct
      m%ncell   = ncell
      m%nleaf   = nleaf

      m%son = son
      m%father = father
      m%octlevel = octlevel
      m%xoct = xoct
      m%nbor = nbor

      ! allocate & fill m%gas
      m%gas = leaves
      
      ! check code
      call check_struct(nleaf,ileaf,m)
      ! clean-up
      ! deallocate arrays
      deallocate(father,son,nbor,octlevel,xoct,ileaf,xleaf,leaflevel)

    end subroutine mesh_from_leaves



    subroutine mesh_from_file(file,m)

      character(2000),intent(in) :: file
      type(mesh),intent(out)   :: m
      !integer                  :: i1,i2,i3,i4
    
      integer,allocatable,dimension(:)        :: iarr1d
      integer,allocatable,dimension(:,:)      :: iarr2d
      !real(kind=8),dimension(:),allocatable   :: farr1d
      real(kind=8),dimension(:,:),allocatable :: farr2d
      type(domain) :: dom
      type(gas),dimension(:),allocatable :: g

      open(unit=13, file=file, status='old', form='unformatted', action='read')

      ! + read type domain
      call read_domain_bin(unit=13,dom=dom)
      m%domain = dom

      read(13) m%ncoarse, m%noct, m%ncell, m%nleaf

      ! allocate m
      allocate(m%son(m%ncell),m%father(m%noct),m%octlevel(m%noct),m%xoct(m%noct,3),m%nbor(m%noct,6))

      ! father
      allocate(iarr1d(m%noct))
      read(13) iarr1d
      m%father = iarr1d
      deallocate(iarr1d)

      ! son
      allocate(iarr1d(m%ncell))
      read(13) iarr1d
      m%son = iarr1d
      deallocate(iarr1d)

      ! nbor
      allocate(iarr2d(m%noct,6))
      read(13) iarr2d
      m%nbor = iarr2d
      deallocate(iarr2d)

      ! octlevel
      allocate(iarr1d(m%noct))
      read(13) iarr1d
      m%octlevel = iarr1d
      deallocate(iarr1d)

      ! xoct
      allocate(farr2d(m%noct,3))
      read(13) farr2d
      m%xoct = farr2d
      deallocate(farr2d)

      ! + read type gas
      call read_gas(unit=13,n=m%nleaf,g=g)
      allocate(m%gas(1:m%nleaf))
      m%gas = g

      !print*,'in load mesh',shape(g),shape(m%gas)

      close(13)
      
    end subroutine mesh_from_file



    subroutine overwrite_mesh(m,nhi,vth)

      type(mesh),intent(inout)  :: m
      real(kind=8),intent(in)   :: nhi,vth

      call overwrite_gas(m%gas,nhi,vth)

#ifdef DEBUG
      print*,'in overwrite_mesh: ',nhi,vth
      print*,'in overwrite_mesh: ',shape(m%gas)
      print*,'in overwrite_mesh: ',m%ncoarse, m%noct, m%ncell, m%nleaf
      print*,'in overwrite_mesh: ',minval(m%gas%nhi),maxval(m%gas%nhi)
#endif

    end subroutine overwrite_mesh


    !===============================================================================================
    ! public use
    !===============================================================================================

    subroutine whereIsPhotonGoing(m,icellold,xpnew,icellnew,flagoutvol)

      type(mesh),intent(in)                :: m
      integer, intent(in)                  :: icellold
      real(kind=8),dimension(3),intent(in) :: xpnew
      integer, intent(out)                 :: icellnew
      logical, intent(out)                 :: flagoutvol
      logical                              :: inside
      integer                              :: ind,ioct,ifather,level,ison,indleaf,ioctleaf,icell,i
      real(kind=8),dimension(3)            :: xc
      real(kind=8)                         :: dx

      ! 1/ test if photon in father(myoct)
      !     if yes -> dig
      !
      ! 2/ if no, test 3 neighbours according to ind in oct
      !
      ! 3/ if no, (photon probably escaped by edges or corners) then find where is photon digging from the top


      !icellnew = in_cell_finder(m,xpnew)
      !! test if icellnew==0 , dans ce cas sortie de domaine...
      !! check if icellnew is inside the domain
      !if(m%son(icellnew)==0) flagoutvol=.true.
      !return

      flagoutvol=.false.

      ! find ind and ioct
      ind  = (icellold-m%ncoarse-1)/m%noct + 1
      ioct = icellold - m%ncoarse - (ind-1)*m%noct

      ! test if photon is in father(ioct)
      ifather = m%father(ioct)    ! is a cell
      xc      = m%xoct(ioct,1:3)
      level   = m%octlevel(ioct)  ! == l
      dx      = 0.5d0**(level)   ! taille des cellules contenues dans l'oct, donc moitie de la cellule father at level l-1
      inside  = ((abs(xc(1)-xpnew(1))<dx).and.(abs(xc(2)-xpnew(2))<dx).and.(abs(xc(3)-xpnew(3))<dx))
      
      if(inside)then
#ifdef DEBUG
         print*,'-WIPG case 1'
#endif
         ! if yes -> dig
         icellnew = digincell(m,xpnew,ifather)
         ! check if icellnew is inside the domain
         if(m%son(icellnew)==0) flagoutvol=.true.
         return

      else ! if no -> test 3 neighbors
         do i=1,3
            icell = m%nbor(ioct,nbortest(i,ind))
            ! les nbor de l'oct sont les cellules voisines de la cellule father de l'oct
            ! 3 cases: 
            ! 1/ if son(icell)>0 -> oct -> donc xoct et on peut tester si dedans et si oui digincell
            ! 2/ if son(icell)<0 -> leaf cell, on recupere ind/ioct, xoct, avec offset on recupere xleaf, si dedans icellnew=ileaf
            ! 3/ if son(icell)=0 -> leaf cell, on recupere ind/ioct, xoct, avec offset on recupere xcell, si dedans (attention au changement de niveau et dx) => renvoyer message comme quoi photon change de domaine

            ison=m%son(icell)
            if(ison>0)then
               xc     = m%xoct(ison,1:3)
               !level  = m%octlevel(ison)    ! = l
               !dx     = 0.5d0**(level)
               inside = ((abs(xc(1)-xpnew(1))<dx).and.(abs(xc(2)-xpnew(2))<dx).and.(abs(xc(3)-xpnew(3))<dx))      
               if(inside)then
#ifdef DEBUG
                  print*,'-WIPG case 2'
#endif
                  ! if yes -> dig
                  icellnew = digincell(m,xpnew,icell)
                  ! check if icellnew is inside the domain
                  if(m%son(icellnew)==0) flagoutvol=.true.
                  return
               endif
            else
               indleaf  = (icell-m%ncoarse-1)/m%noct + 1
               ioctleaf = icell - m%ncoarse - (indleaf-1)*m%noct
               xc       = m%xoct(ioctleaf,1:3)
               !level    = m%octlevel(ioctleaf)  ! = l-1
               !dx       = 0.5d0**(level)
               xc       = xc + offset(1:3,indleaf)*dx*2 ! on garde dx at level l => x2
               inside   = ((abs(xc(1)-xpnew(1))<dx).and.(abs(xc(2)-xpnew(2))<dx).and.(abs(xc(3)-xpnew(3))<dx))      
               if(inside)then
#ifdef DEBUG
                  print*,'-WIPG case 2'
#endif
                  icellnew=icell
                  if(ison==0)flagoutvol=.true.
                  ! else ! ison<0 ok
                  return
               endif
            endif
         enddo
         ! if not in the 3 neighbors, then dig from the top (escaped by an edge or a corner...)
#ifdef DEBUG
         print*,'-WIPG case 3'
#endif
         icellnew = in_cell_finder(m,xpnew)
         ! test if icellnew==0 , dans ce cas sortie de domaine...
         ! check if icellnew is inside the domain
         if(m%son(icellnew)==0) flagoutvol=.true.
         return
      endif

      print *,'oh my god, where is the photon...'
      stop

    end subroutine whereIsPhotonGoing



    function digincell(m,xp,icell)

      type(mesh),               intent(in) :: m
      real(kind=8),dimension(3),intent(in) :: xp
      integer,                  intent(in) :: icell
      integer                              :: digincell
      integer                              :: ison,ix,iy,iz,ind,icellint
      real(kind=8),dimension(3)            :: x

      ! find the leaf cell containing xp, starting from icell
      ! this assumes that xp is contained inside icell... 

      ison = m%son(icell)               ! its (oct) son
      do while (ison > 0)               ! there is an oct in current cell 
         x  = m%xoct(ison,1:3)          ! oct position
         ix = merge(0,1,xp(1) < x(1)) 
         iy = merge(0,1,xp(2) < x(2))
         iz = merge(0,1,xp(3) < x(3))
         ind = ix + 2*iy + 4*iz 
         icellint = m%ncoarse + ind * m%nOct + ison
         ison  = m%son(icellint)
      end do
      !digincell=-ison
      digincell=icellint

      return

    end function digincell
    



    function in_cell_finder(m,xp)

      type(mesh),intent(in)                :: m
      real(kind=8),dimension(3),intent(in) :: xp
      integer                              :: in_cell_finder
      integer                              :: ison,ix,iy,iz,ind,icell
      real(kind=8),dimension(3)            :: x

      ! find leaf cell containing point xp
      icell = 1 ! the coarse cell containing the whole box is always nb 1
      ison  = m%son(icell) ! its (oct) son
      do while (ison > 0)  ! there is an oct in current cell 
         x  = m%xoct(ison,1:3)  ! oct position
         ix = merge(0,1,xp(1) < x(1)) 
         iy = merge(0,1,xp(2) < x(2))
         iz = merge(0,1,xp(3) < x(3))
         ind = ix + 2*iy + 4*iz
         icell = m%ncoarse + ind * m%nOct + ison
         ison  = m%son(icell)
      end do
      !in_cell_finder=-ison
      in_cell_finder=icell

      return
    end function in_cell_finder



    function get_cell_corner(xoct, ind, cell_level) 
      ! return position of cell corner, in box units.
      
      real(kind=8),dimension(3),intent(in) :: xoct
      integer, intent(in)                  :: ind, cell_level
      real(kind=8),dimension(3)            :: get_cell_corner, xcell
      real(kind=8)                         :: dx

      dx              = 0.5d0**cell_level
      xcell           = xoct + offset(1:3,ind)*dx
      get_cell_corner = xcell - dx/2.d0

!#ifdef DEBUG
!      print *,'--> in function get_cell_corner'
!      print *,cell_level, dx, ind
!      print *,xoct
!      print *,xcell
!#endif

    end function get_cell_corner



    subroutine dump_mesh(m,file)

      type(mesh),intent(in)    :: m
      character(2000),intent(in) :: file

      ! dump data in binary format

      write(*,*)
      write(*,*) '...dump file...'
      write(*,*)


#ifdef DEBUG
      print *,'--> check mesh dom'
      print *,m%domain
      print *,m%nCoarse,m%nOct,m%nLeaf,m%nCell
      print *,minval(m%xoct(:,:)),maxval(m%xoct(:,:))
      print *,minval(m%nbor(:,:)),maxval(m%nbor(:,:))
      print *,minval(m%octlevel(:)),maxval(m%octlevel(:))
      print *,minval(m%son(:)),maxval(m%son(:))
      print *,minval(m%father(:)),maxval(m%father(:))
#endif


      open(unit=13, file=file, status='unknown', form='unformatted', action='write')

      ! dump domain
      call dump_domain_bin(unit=13,dom=m%domain)

      write(13) m%nCoarse, m%nOct, m%nCell, m%nLeaf
      write(13) m%father(1:nOct)
      write(13) m%son(1:nCell)
      write(13) m%nbor(1:nOct,1:6)
      write(13) m%octlevel(1:nOct)
      write(13) m%xoct(1:nOct,1:3)

      ! dump gas
      call dump_gas(unit=13,g=m%gas)

      close(13)

    end subroutine dump_mesh



    subroutine mesh_destructor(m)
      type(mesh),intent(inout) :: m
      call gas_destructor(m%gas)
      ! maybe not needed...
      !call domain_destructor(m%domain)
      deallocate(m%father,m%son,m%nbor,m%xoct,m%octlevel)
    end subroutine mesh_destructor


    !===============================================================================================
    ! PRIVATE routines
    !===============================================================================================


    recursive subroutine add_oct_domain(ilastoct,n,id,fcx,fcl,ifc)
    
      integer,intent(inout)                :: ilastoct
      integer,intent(in)                   :: n
      integer,dimension(n),intent(in)      :: id
      integer,intent(in)                   :: fcl,ifc
      real(kind=8),dimension(3),intent(in) :: fcx
      integer                              :: ind,level,nleaflocal,icell,ilastoctLevel
      integer,dimension(:),allocatable     :: ileaflocal
      real(kind=8),dimension(3)            :: xcentre
      logical                              :: same_level

      ilastoct = ilastoct + 1 
      octlevel(ilastoct)=fcl+1
      xoct(ilastoct,1:3) = fcx(1:3)
      father(ilastoct) = ifc
      level=octlevel(ilastoct)
      ilastoctLevel = ilastoct

      do ind=1,8

         if(level==1)then
            write(*,'(a,i1,a)')' level=1 cell # ',ind,' done'
         endif

         ! define cell centre and dx 
         xcentre = xcson(ind,fcx,level)

         nleaflocal = get_nleaflocal(n,id,xcentre,level)

         icell = nCoarse+ilastoctLevel+(ind-1)*nOct

         !add_a_new_oct = .false.
         !add_a_new_oct=((neaflocal>1).or.((nleaflocal==1).and.(
         ! define new condition:
         ! nleaflocal>1 or nleaflocal==1 && level(oct)==level(leaf)

         ! 3 cases:
         if(nleaflocal>1)then
            ! this is the old case, you add a new oct

            allocate(ileaflocal(nleaflocal))
            ileaflocal = get_ileaflocal(n,id,xcentre,level,nleaflocal)                   ! indices dans xl(:,1:3)
            son(icell) = ilastoct+1
            call add_oct_domain(ilastoct,nleaflocal,ileaflocal,xcentre,level,icell)
            deallocate(ileaflocal)

         else

            if(nleaflocal==1)then

               ! this is the true leaf OR one of the leaves (the other are not present in the selection)
               allocate(ileaflocal(nleaflocal))
               ileaflocal = get_ileaflocal(n,id,xcentre,level,nleaflocal)

               if(level/=octlevel(ilastoctLevel))then
                  print*,'test level false...'
               end if

               !same_level=all(octlevel(ilastoctLevel)==leaflevel(ileaflocal))
               same_level=all(level==leaflevel(ileaflocal))

               if(same_level)then
                  son(icell)=-ileaflocal(1)
                  countleaf = countleaf+1
                  deallocate(ileaflocal)
               else
                  ! add a new oct
                  son(icell) = ilastoct+1
                  call add_oct_domain(ilastoct,nleaflocal,ileaflocal,xcentre,level,icell)
                  deallocate(ileaflocal)
               endif

            else
               ! nleaflocal=0
               ! there is no leaf, flag son to 0

               if(nleaflocal/=0)then
                  print*,'Oh no, what happens',nleaflocal
                  stop
               endif

               son(icell)=0
               countempty=countempty+1

            endif

         endif

      enddo

    end subroutine add_oct_domain

    
    
    function xcson(ind,fx,level)

      integer,intent(in)                    :: ind,level
      real(kind=8),dimension(3),intent(in)  :: fx
      real(kind=8)                          :: dx
      real(kind=8),dimension(3)             :: xcson

      dx = 0.5d0**level / 2.
      select case(ind) 
      case(1) ! 1st cell in oct
         xcson(1) = fx(1) - dx
         xcson(2) = fx(2) - dx
         xcson(3) = fx(3) - dx
      case(2) ! 2nd cell in oct
         xcson(1) = fx(1) + dx
         xcson(2) = fx(2) - dx
         xcson(3) = fx(3) - dx
      case(3) ! 3rd cell in oct
         xcson(1) = fx(1) - dx
         xcson(2) = fx(2) + dx
         xcson(3) = fx(3) - dx
      case(4) ! 4th cell in oct
         xcson(1) = fx(1) + dx
         xcson(2) = fx(2) + dx
         xcson(3) = fx(3) - dx
      case(5) ! 5th cell in oct
         xcson(1) = fx(1) - dx
         xcson(2) = fx(2) - dx
         xcson(3) = fx(3) + dx
      case(6) ! 6th cell in oct
         xcson(1) = fx(1) + dx
         xcson(2) = fx(2) - dx
         xcson(3) = fx(3) + dx
      case(7) ! 7th cell in oct
         xcson(1) = fx(1) - dx
         xcson(2) = fx(2) + dx
         xcson(3) = fx(3) + dx
      case(8) ! 8th cell in oct
         xcson(1) = fx(1) + dx
         xcson(2) = fx(2) + dx
         xcson(3) = fx(3) + dx
      end select

      return

    end function xcson


    
    function get_nleaflocal(n,ileaf,xc,level)

      integer,intent(in)                   :: n,level
      integer,dimension(n),intent(in)      :: ileaf
      real(kind=8),dimension(3),intent(in) :: xc
      integer                              :: icell,jcell,kcell,ii,jj,kk,i,count,get_nleaflocal

      count=0
      icell=int(xc(1)*2**level)+1 
      jcell=int(xc(2)*2**level)+1 
      kcell=int(xc(3)*2**level)+1 

      do i=1,n
         ii=int(xleaf(ileaf(i),1)*2**level)+1
         jj=int(xleaf(ileaf(i),2)*2**level)+1
         kk=int(xleaf(ileaf(i),3)*2**level)+1
         if(ii==icell.and.jj==jcell.and.kk==kcell)then
            ! ileaf(i) in cell
            count=count+1
         endif
      enddo
      get_nleaflocal=count

      return

    end function get_nleaflocal




    function get_ileaflocal(n,ileaf,xc,level,nleaflocal)

      integer,intent(in)                        :: n,level,nleaflocal
      integer,dimension(n),intent(in)           :: ileaf
      real(kind=8),dimension(3),intent(in)      :: xc
      integer                                   :: icell,jcell,kcell,ii,jj,kk,i,count
      integer,dimension(nleaflocal)             :: get_ileaflocal

      count=0
      icell=int(xc(1)*2**level)+1 
      jcell=int(xc(2)*2**level)+1 
      kcell=int(xc(3)*2**level)+1 

      do i=1,n
         ii=int(xleaf(ileaf(i),1)*2**level)+1
         jj=int(xleaf(ileaf(i),2)*2**level)+1
         kk=int(xleaf(ileaf(i),3)*2**level)+1
         if(ii==icell.and.jj==jcell.and.kk==kcell)then
            ! ileaf(i) in cell
            count=count+1
            get_ileaflocal(count)=ileaf(i)
         endif
      enddo

      return
    end function get_ileaflocal

    
    !JB-
    function icell2icell(icell,ngridmax,ngridmin,ncoarse)

      ! returns icell indexed with ngridmin instead of ngridmax

      integer(kind=4),intent(in) :: icell,ngridmax,ngridmin,ncoarse
      integer(kind=4)            :: icell2icell
      integer(kind=4)            :: ind,jcell

      ind   = (icell - ncoarse - 1) / ngridmax + 1
      jcell = icell - ncoarse - (ind - 1) * ngridmax

      icell2icell = ncoarse + jcell + (ind-1) * ngridmin

      return

    end function icell2icell
    !-JB


    subroutine make_nbor_array
      
      ! make nbor array
      !JB-
      
      ! need octlevel, xoct, ncoarse, noct, son 
      
      integer                                 :: inbor,level,ix,iy,iz,icell,l,ison
      integer                                 :: i,ind,ioct
      real(kind=8)                            :: xc(3),xnbor(3),x(3),dx
      integer                                 :: noops1,noops2


      write(*,*)'into make nbor array...'
      write(*,*)ncoarse,noct
      write(*,*)minval(son),maxval(son)
      write(*,*)minval(octlevel),maxval(octlevel)

      if (ncoarse /= 1) then ! on a vraiment besoin de ca en fait... 
         Print*,'Oh no, ncoarse /= 1... '
         stop
      end if
      noops1=0
      noops2=0
      do ioct=1,nOct  ! loop over all octs. 
         level = octlevel(ioct)  ! level of oct -> 6 neighbors are cells at level-1
         dx    = 0.5d0**(level-1) ! size of neighbor cells (or cell containing the oct)
         xc    = xoct(ioct,1:3)  ! position of current oct.
         if (level == 1) then   ! this is the oct within the unique coarse cell (at l=0)/
            nbor(ioct,1:6) = 1  ! -> define its neighbors as the coarse cell itself.
            cycle
         end if
         if (level <= 0) then  ! this should not happen, but it does : what are these octs ? 
            if (level == 0) noops1 = noops1 + 1
            !if (level == -1) noops2 = noops2 + 1
            if (level < -1) print*,'oops'
            cycle
         end if
         do inbor=1,6   ! loop over 6 neighbor cells and compute the coords. of their center
            select case(inbor) 
            case (1) ! 1st neighbour (-dx)
               xnbor(1) = xc(1) - dx
               xnbor(2) = xc(2)
               xnbor(3) = xc(3)
            case (2) ! 2nd neighbour (+dx)
               xnbor(1) = xc(1) + dx
               xnbor(2) = xc(2)
               xnbor(3) = xc(3)
            case (3) ! 3rd neighbour (-dy)
               xnbor(1) = xc(1)
               xnbor(2) = xc(2) - dx
               xnbor(3) = xc(3)
            case (4) ! 4th neighbour (+dy)
               xnbor(1) = xc(1)
               xnbor(2) = xc(2) + dx
               xnbor(3) = xc(3)
            case (5) ! 5th neighbour (-dz)
               xnbor(1) = xc(1)
               xnbor(2) = xc(2)
               xnbor(3) = xc(3) - dx
            case (6) ! 6th neighbour (+dz)
               xnbor(1) = xc(1)
               xnbor(2) = xc(2)
               xnbor(3) = xc(3) + dx
            end select
            ! correct for periodic boundary conditions :
            do i = 1,3
               if (xnbor(i) > 1.0d0) xnbor(i) = xnbor(i) - 1.0d0
               if (xnbor(i) < 0.0d0) xnbor(i) = xnbor(i) + 1.0d0
            end do
            ! find cell of level "level-1" containing xnbor
            icell = 1 ! the coarse cell containing the whole box is always nb 1
            l     = 0 ! its level 
            ison  = son(icell) ! its (oct) son
            do while (ison > 0)  ! there is an oct in current cell 
               x  = xoct(ison,1:3)  ! oct position
               ix = merge(0,1,xnbor(1) < x(1)) 
               iy = merge(0,1,xnbor(2) < x(2))
               iz = merge(0,1,xnbor(3) < x(3))
               ind = ix + 2*iy + 4*iz 
               icell = ncoarse + ind * nOct + ison
               ison  = son(icell)
               l   = l + 1
               if (l == level-1) exit
            end do
            nbor(ioct,inbor) = icell 
            if((ison/=0).and.(l/=level-1))then
               print*,'Ouh la la',ioct,ison,l,level,inbor
               stop
               !noops2=noops2+1
            end if
         end do
      end do
      print*,noops1,noops2
      !-JB

    end subroutine make_nbor_array
  


    subroutine resize_octtree(nold,nnew)

      integer, intent(in) :: nold,nnew

      real(kind=8),dimension(:,:),allocatable :: tmpr
      integer,dimension(:),allocatable        :: tmpi
      integer                                 :: ncellnew,i

      ncellnew = nCoarse+8*nnew

      write(*,*)
      write(*,*) '...resize arrays...'

      ! transform & resize father
      allocate(tmpi(nnew))
      do i=1,nnew
         tmpi(i) = icell2icell(father(i),nold,nnew,nCoarse) ! change father cell num from indexation with nOctTot octs to indexation with ilastoct octs.
      end do
      deallocate(father)
      allocate(father(nnew))
      father = tmpi
      deallocate(tmpi)

      ! transform & resize son
      allocate(tmpi(ncellnew))
      tmpi(1)=son(1)
      do i=1,8  ! gather all cells in begining of array (don't change values which are oct indexes).
         tmpi(nCoarse+1+(i-1)*nnew:nCoarse+i*nnew) = son(nCoarse+1+(i-1)*nold:nCoarse+(i-1)*nold+nnew)
      end do
      deallocate(son)
      allocate(son(ncellnew))
      son = tmpi
      deallocate(tmpi)

      ! resize octlevel
      allocate(tmpi(1:nnew))
      tmpi = octlevel(1:nnew)
      deallocate(octlevel)
      allocate(octlevel(nnew))
      octlevel = tmpi
      deallocate(tmpi)
      
      ! resize nbor
      !deallocate(nbor)
      !allocate(nbor(nOctTot,6))
      !nbor = 0

      ! resize xoct
      allocate(tmpr(nnew,3))
      tmpr = xoct(1:nnew,:)
      deallocate(xoct)
      allocate(xoct(nnew,3))
      xoct = tmpr
      deallocate(tmpr)
      
    end subroutine resize_octtree
    


    subroutine check_octtree(noctmax)

      integer,intent(in) :: noctmax
      integer :: oct_status,i,ind,ioct,icount,icount2

      ! check results

      write(*,*)
      write(*,*)'...checking arrays...'
      write(*,*)'min max son      ',minval(son(:)),maxval(son(:))
      write(*,*)'min max father   ',minval(father(:)),maxval(father(:))
      write(*,*)'min max octlevel ',minval(octlevel(:)),maxval(octlevel(:))
      write(*,*)'min max xoct',minval(xoct(:,1)),maxval(xoct(:,1))
      write(*,*)'min max xoct',minval(xoct(:,2)),maxval(xoct(:,3))
      write(*,*)'min max xoct',minval(xoct(:,3)),maxval(xoct(:,2))

      write(*,*)'min max octlevel (>0)',minval(octlevel, mask=(octlevel >= 0)), maxval(octlevel, mask=(octlevel >= 0))
      write(*,*)'min max xoct (>0)',minval(xoct(:,1), mask=(xoct(:,1)>=0)),maxval(xoct(:,1), mask=(xoct(:,1)>=0))
      write(*,*)'min max xoct (>0)',minval(xoct(:,2), mask=(xoct(:,2)>=0)),maxval(xoct(:,2), mask=(xoct(:,2)>=0))
      write(*,*)'min max xoct (>0)',minval(xoct(:,3), mask=(xoct(:,3)>=0)),maxval(xoct(:,3), mask=(xoct(:,3)>=0))

      !write(*,*)'min max ileaf      ',minval(ileaf(:)),maxval(ileaf(:))

      countleaf=0
      icount=0
      do i=1,nCoarse
         icount=icount+1
         ioct=i
         if(son(ioct)<0)then
            countleaf=countleaf+1
         endif
      enddo
      write(*,*) 'ncoarse, son(ncoarse) =',ncoarse,son(1:ncoarse)
      do i=1,nOct
         do ind=1,8
            icount=icount+1
            ioct = nCoarse + i + (ind-1)*nOct
            if(son(ioct)<0)then
               countleaf=countleaf+1
            endif
         enddo
      enddo
      write(*,*)'scan son array # of elements = ',icount
      write(*,*)'nleaf in son                 = ',countleaf
      write(*,*)'diff                         = ',icount-countleaf


      write(*,*)
      ! count incomplete octs...
      icount=0
      icount2=0
      do i=1,nOct
         oct_status=0
         do ind=1,8
            ioct = nCoarse + i + (ind-1)*nOct
            if(son(ioct)/=0)then
               oct_status=oct_status+1
            endif
         enddo
         if (oct_status/=8)then
            icount=icount+1
         endif
         if (oct_status==0)then
            icount2=icount2+1
         endif
         !!if(oct_status>0 .and. oct_status<8)then
         !!   print*,octlevel(i)
         !!endif
      enddo
  
      write(*,*)'incomplete & empty octs =',icount
      write(*,*)'empty octs              =',icount2

      write(*,*)'nOctTot - Noct          =',nOctmax-noct

      write(*,*)

      write(*,*)'Ncoarse + Nleaf + (Noct-1) + 8*NemptyOct =',nCoarse + nLeaf + (noct-1) + icount2*8
      ! Leo: Noct-1 car la cellule coarse principale est aussi un oct maintenant !
      !JB-
      !ouh la la, vraiment ?
      !-JB
      write(*,*)'Ncell                   =',nCell
      write(*,*)'nCell - (nOct + Nleaf)  =',nCell - nLeaf - nOct
      
      write(*,*)
      write(*,*)'...check father array...'
      !icount=0
      !do i=1,nOcttot
      !   if(father(i)==-1)then
      !      icount=icount+1
      !   endif
      !enddo
      !write(*,*)'# of -1 in father array =',icount
      write(*,*) '# of -1 in father array =', count(mask=(father==-1))
      write(*,*) '# of  0 in father array =', count(mask=(father==0))
      write(*,*) '# of >0 in father array =', count(mask=(father>0))
      
    end subroutine check_octtree



    subroutine check_struct(n,ileaf,m)
      ! pour chaque cellule, parcourir l'arbre, retrouver la cellule feuille correspondante et checker valeur de rho
      ! use do while loop as in make_nbor...

      integer, intent(in) :: n
      integer,dimension(1:n), intent(in) :: ileaf
      type(mesh),intent(in)     :: m

      integer                   :: i,ic,ind,icell
      real(kind=8),dimension(3) :: xc

      do i=1,n
         xc(1:3) = xleaf(i,1:3)
         ic = in_cell_finder(m,xc)
         ! some check
         if(m%son(ic)/=-ileaf(i))then
            print*,'Problem with in_cell_finder...',i,ic,ileaf(i),m%son(ic)
            stop
         endif
         !!!if(m%gas%density(ic)/=m%gas%density(i))then
         !if(m%gas(ic)%density/=m%gas(i)%density)then
         !   print*,'Problem with in_cell_finder...',i,ic,ileaf(i)
         !   stop
         !end if
      end do


      ! for each leaf cell check level
      countleaf=0
      do i=1,m%nOct
         do ind=1,8
            icell = m%nCoarse + i + (ind-1)*m%nOct
            if(m%son(icell)<0)then
               countleaf=countleaf+1
               if(m%octlevel(i)/=8)then
                  print*,'pb with mesh',icell,i,m%son(icell),m%octlevel(i)
                  stop
               endif
            endif
         enddo
      enddo
      print*,'nleaves =',countleaf

      print*,'end of check: result = ok'

    end subroutine check_struct


  end module module_mesh
