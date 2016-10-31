module module_photon
 
  use module_gas_composition
  use module_mesh
  use module_constants
  use module_random
  use module_domain

  implicit none

  !use ramses_info_utils, only box_size_cm   ! => not defined in this module....
  ! todonext
!  real(kind=8) :: box_size_cm  -> moved to module_gas_composition
  ! also todonext define accuracy
  real(kind=8),parameter :: accuracy=1.d-15

  ! 2 types for photons, one for the initial properties and one for the properties evolving during RT
  type photon_current
     integer                   :: ID
     integer                   :: status       ! =0 if flying, =1 if escape, =2 if absorption (by dust)
     real(kind=8),dimension(3) :: xlast        ! coordinates of last interaction in box units
     real(kind=8),dimension(3) :: xcurr        ! current position of the photon in box units
     real(kind=8)              :: nu_ext       ! external frame frequency (Hz)
     real(kind=8),dimension(3) :: k            ! normalised propagation vector 
     integer                   :: nb_abs       ! number of interactions before escape
     real(kind=8)              :: time         ! time in [s] from emission to escape/absorption        
     real(kind=8)              :: tau_abs_curr ! current optical depth (useful when photon change mesh domain)
     integer                   :: iran         ! state of the random generator
  end type photon_current

  !real(kind=8)              ::  theta       ! direction of escape
  !real(kind=8)              ::  phi
  !#ifdef POLARIZATION
  !     matrix(26,nInBuffer) = x_pol
  !     matrix(27,nInBuffer) = y_pol
  !     matrix(28,nInBuffer) = z_pol
  !#endif
  !     ! lies a MPI ???
  !     matrix(29,nInBuffer) = iran_em
  !     matrix(30,nInBuffer) = ix_ran_em
  !     matrix(31,nInBuffer) = iy_ran_em
  !     matrix(32,nInBuffer) = am_ran_em


  type photon_init
     integer                   :: ID
     real(kind=8)              :: nu_em    ! emitted frequency in external frame
     real(kind=8),dimension(3) :: x_em
     real(kind=8),dimension(3) :: k_em
     integer                   :: iran         ! state of the random generator
  end type photon_init

  !    real(kind=8) :: nu_star   ! source-frame emission frequency (Hz), encore necessaire?
  !    real(kind=8) :: nu_em    ! emitted frequency in external frame
  !    real(kind=8) :: log_lum_lya_em   ! poid des photons, #ifdef GAS_EMITS case
  !    real(kind=8) :: xx_em    ! original location of the photon        
  !    real(kind=8) :: yy_em    
  !    real(kind=8) :: zz_em
  !    real(kind=8) :: a0_em    ! direction along which the photon was emitted.
  !    real(kind=8) :: b0_em
  !    real(kind=8) :: c0_em
  !    real(kind=8) :: theta_em
  !    real(kind=8) :: phi_em
  !    ! veut-on garder les props initiales de la polarisation??? Oui
  !    real(kind=8) :: x_pol
  !    real(kind=8) :: y_pol
  !    real(kind=8) :: z_pol
  !    ! random generator stuff
  !    real(kind=8) :: iran_em
  !    real(kind=8) :: ix_ran_em
  !    real(kind=8) :: iy_ran_em
  !    real(kind=8) :: am_ran_em
  ! end type photon_init
  !------------------------

  public  :: MCRT, propagate, init_photons_from_file, dump_photons
  private :: path

contains

  ! subroutine initialisation
  ! 1 / calcul des proprietes initiales
  !     plusieurs cas selon les sources, les types de simu (RAMSES-RT ou pas)
  !     positions, mais aussi la distribution spectrale, et distribution angulaire (si pas isotrope) 
  ! 2 / initialiser les proprietes qui evoluent a ces valeurs initiales

  ! subroutine dump_photons => fin de la sim
  ! subroutine backup_photons
  ! subroutine restart
  ! subroutine propagate_photon
  ! subroutine transform_theta_phi ?
  ! subroutine transform_a_b_c

  subroutine MCRT(nbuffer,photpacket,mesh_dom,compute_dom)

    ! this is THE MCLya routine doing the radiative transfer... 
    ! but also a single loop over photons...

    !$ use OMP_LIB

    integer, intent(in)                                     :: nbuffer
    type(photon_current), dimension(nbuffer), intent(inout) :: photpacket
    type(mesh), intent(in)                                  :: mesh_dom
    type(domain), intent(in)                                :: compute_dom
    integer                                                 :: i,counter
    !integer                                                 :: rank,nthreads

!$OMP PARALLEL
    !$ rank = OMP_GET_THREAD_NUM()
    !$ if (rank == 0) then
    !$    nthreads = OMP_GET_NUM_THREADS()
    !$   print *, "Number of threads = ", nthreads
    !$ end if    
!$OMP END PARALLEL


    counter = 0
!$OMP PARALLEL DO &
    !$OMP DEFAULT(SHARED) &
    !$OMP PRIVATE(i,rank)
    do i=1,nbuffer

#ifdef DEBUG
       print *,'--> MCRT launching new photon number =',i,photpacket(i)
#endif

       !if(mod(i,nbuffer/10)==0)print *,'--> MCRT launching new photon number =',i

       !$ rank = OMP_GET_THREAD_NUM()
       !$ print *,'Rank:',rank,' --> MCRT launching new photon number =',i!!!,photpacket(i)

       !counter = counter+1
       !if(mod(counter,nbuffer/10)==0)print *,'--> MCRT processing photon number =',i,counter,' /',nbuffer

       !call propagate(photpacket(i),mesh_dom,compute_dom)
       ! case of photpacket not fully filled...
       if (photpacket(i)%ID>0) then
          call propagate(photpacket(i),mesh_dom,compute_dom)
       endif
       !stop
    enddo
!$OMP END PARALLEL DO

  end subroutine MCRT



  subroutine propagate(p,domesh,domaine_calcul)

    type(photon_current),intent(inout)   :: p              ! photon 
    type(mesh),intent(in)                :: domesh         ! mesh
    type(domain),intent(in)              :: domaine_calcul ! domaine dans lequel on propage les photons...
    
    type(gas)                            :: cell_gas  ! gas in the current cell 
    integer(kind=4)                      :: icell, ioct, ind, ileaf, cell_level  ! current cell indices and level
    real(kind=8)                         :: cell_size, cell_size_cm, tau_abs, scalar, nu_cell, nu_ext, rtau
    real(kind=8),dimension(3)            :: ppos,ppos_cell ! working coordinates of photon (in box and in cell units)
    real(kind=8)                         :: distance_to_border,distance_to_border_cm, d
    real(kind=8)                         :: time
    integer(kind=4)                      :: scatter_flag, i, icellnew, iran
    real(kind=8),dimension(3)            :: vgas, k, cell_corner, posoct, ppos_old
    logical                              :: cell_fully_in_domain, flagoutvol, in_domain
    real(kind=8)                         :: epsilon_cell
    real(kind=8)                         :: dtime,dborder
    
    ! initialise working props of photon
    ppos    = p%xcurr        ! position within full simulation box, in box units.
    time    = p%time
    tau_abs = p%tau_abs_curr ! a initialiser a -1 au debut avec time=0 et xcurr=xinit
    iran    = p%iran

#ifdef DEBUG
    print *,'--> propagating photon from',ppos, iran
    print *,'--> in mesh domain',domesh%nCoarse,domesh%nOct,domesh%nLeaf,domesh%nCell
    print *,minval(domesh%xoct(:,1)),maxval(domesh%xoct(:,1))
    print *,minval(domesh%xoct(:,2)),maxval(domesh%xoct(:,2))
    print *,minval(domesh%xoct(:,3)),maxval(domesh%xoct(:,3))
    print *,minval(domesh%gas%nHI),maxval(domesh%gas%nHI)
#endif

    ! find cell in which the photon is, and define all its indices
    icell = in_cell_finder(domesh,ppos)
    if(domesh%son(icell)>=0)then
       print*,'not a leaf cell'
       stop
    endif
    ileaf = - domesh%son(icell)
    ind   = (icell - domesh%nCoarse - 1) / domesh%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
    ioct  = icell - domesh%nCoarse - (ind - 1) * domesh%nOct

!JB--
!!$    if((domesh%son(icell)<0).and.(domesh%octlevel(ioct)/=8))then
!!$       print*,'>>>>> leaf cell at level lower than 8!!!!!!'
!!$       stop
!!$    endif
!--JB

#ifdef DEBUG
    print *,'--> cell where photon starts',icell,ileaf,ind,ioct
    print *,'--> starting photon_propagation loop'
#endif

    ! propagate photon until escape or death ... 
    photon_propagation : do 

       ! gather properties properties of current cell
       cell_level   = domesh%octlevel(ioct)      ! level of current cell
       cell_size    = 0.5d0**cell_level          ! size of current cell in box units
       cell_size_cm = cell_size * box_size_cm    ! size of the current cell in cm
       cell_gas     = domesh%gas(ileaf)
       ! compute position of photon in current-cell units
       posoct(:)    = domesh%xoct(ioct,:)
       cell_corner  = get_cell_corner(posoct,ind,cell_level)   ! position of cell corner, in box units.
       ppos_cell    = (ppos - cell_corner) / cell_size       ! position of photon in cell units (x,y,z in [0,1] within cell)
       if((ppos_cell(1)>1.).or.(ppos_cell(2)>1.).or.(ppos_cell(3)>1.))then
          print*,"--> Problem in computing ppos_cell"
          stop
       endif

       ! get gas velocity (in cgs units)
       vgas         = get_gas_velocity(cell_gas)
       ! compute photon's frequency in cell's moving frame
       scalar       = p%k(1) * vgas(1) + p%k(2) * vgas(2) + p%k(3) * vgas(3)
       nu_cell      = (1.d0 - scalar/clight) * p%nu_ext  

       ! define epsilon according to cell size & numerical accuracy, a creuser plus tard et tester
       epsilon_cell = 2.d0*accuracy/cell_size * 1.d5
       ! /cell_size car 
       ! on veut epsilon_box > tiny
       ! => epsilon_cell > tiny/cell_size

       ! define/update flag_cell_fully_in_comp_dom pour eviter tous les tests par la suite ET les updates de ppos (??)
       cell_fully_in_domain =  domain_contains_cell(ppos,cell_size,domaine_calcul)


#ifdef DEBUG
       print *,'--> cell properties',cell_level,cell_size,cell_size_cm
       print *,'        cell_gas  =',cell_gas
       print *,'        pos corner=',cell_corner
       print *,'                   ',(ppos-cell_corner)
       print *,'       pos in cell=',ppos_cell
       print *,'                   ',cell_fully_in_domain
       print *,'                   ',epsilon_cell
       print *,'                k= ',p%k
       print *,'          pos oct =',posoct
       print *,'--> current photon position',ppos
#endif

       propag_in_cell : do
          
          ! generate the opt depth where the photon is scattered/absorbed
          if (tau_abs <= 0.0d0) then
             rtau    = ran3(iran)
             tau_abs = -log(1.0d0-rtau+1.d-30)
          end if

          ! compute distance of photon to border of cell along propagation direction
          distance_to_border    = path(ppos_cell,p%k)               ! in cell units
          distance_to_border_cm = distance_to_border * cell_size_cm ! cm

          ! check whether scattering occurs within cell (scatter_flag > 0) or not (scatter_flag==0)
          scatter_flag = gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs, iran)

#ifdef DEBUG
          print*,'distance_to_border_cm =',distance_to_border_cm
          print*,'distance_to_border    =',distance_to_border
          print*,'scatter_flag          =',scatter_flag
          print*,'box_size_cm           =',box_size_cm
          print*,'cell crossing time    =',cell_size_cm/clight
          print*,'epsilon crossing time =',epsilon_cell*cell_size_cm/clight
#endif

          if (scatter_flag == 0) then
             ! next scattering event will not occur in the cell, then move photon to next cell

             ! update ppos_cell with distance_to_border (distance_to_border_cm has not been modified if flag==0)
             ! also add epsilon to ensure finding next cell
             ppos_cell = ppos_cell + p%k * (distance_to_border + epsilon_cell)

             ! update ppos (in box units)
             ppos_old = ppos
             ppos = ppos_cell * cell_size + cell_corner
#ifdef DEBUG
             print*,'diff ppos =',ppos-ppos_old
#endif
             ! correct for periodicity
             do i=1,3
                if (ppos(i) < 0.d0) ppos(i)=ppos(i)+1.d0
                if (ppos(i) > 1.d0) ppos(i)=ppos(i)-1.d0
             enddo

             ! update travel time
             time = time + distance_to_border_cm/clight

             ! check if photon still in computational domain
             in_domain = domain_contains_point(ppos,domaine_calcul)
             if(.not.(in_domain))then     ! => photon done
                p%status       = 1
                p%xcurr        = ppos
                ! correct time
                dborder = domain_distance_to_border(ppos,domaine_calcul)
                dtime   = dborder*box_size_cm/clight ! should be negative
                time    = time+dtime
                p%time         = time
                p%tau_abs_curr = tau_abs
                p%iran         = iran
#ifdef DEBUG
                print*,'-exit propagation, photon escaped comp. domain'
#endif
                exit photon_propagation
             endif

             call whereIsPhotonGoing(domesh,icell,ppos,icellnew,flagoutvol)
             ! check routine
             if(icell==icellnew)then
                print*,'Problem with routine WhereIsPhotonGoing'
                print*,'epsilon    =',accuracy, epsilon_cell, distance_to_border
                print*,'ppos_cell  =',ppos_cell
                print*,'delta_ppos =', p%k * (distance_to_border + epsilon_cell)
                print*,'ppos       =',ppos
                print*,'cell_size  =',cell_size
#ifdef DEBUG
                stop
#endif
                ! does not stop/kill the job anymore in this case, just flag this photon as crap (=3)
                p%status       = 3
                p%xcurr        = ppos
                p%time         = time
                p%tau_abs_curr = tau_abs
                p%iran         = iran
                exit photon_propagation
             endif
             ! check if photon outside of cpu domain (flagoutvol)
             if(flagoutvol)then
                ! photon out of cpu domain, to be sent back to master
!!#ifdef DEBUG
!!                p%status       = 1
!!#endif
                p%xcurr        = ppos
                p%time         = time
                p%tau_abs_curr = tau_abs
                p%iran         = iran
#ifdef DEBUG
                print*,'-exit propagation, photon escaped mesh domain'
#endif
                exit photon_propagation
             endif
             ! else, new cell is in the cpu domain and photon goes to this new cell
             icell = icellnew
             if(domesh%son(icell)>=0)then
                print*,'not a leaf cell',icell,flagoutvol
                !stop
             endif
             ileaf = - domesh%son(icell)
             ind   = (icell - domesh%nCoarse - 1) / domesh%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
             ioct  = icell - domesh%nCoarse - (ind - 1) * domesh%nOct
!JB-- 
!!$             if((domesh%son(icell)<0).and.(domesh%octlevel(ioct)/=8))then
!!$                print*,'>>>>> leaf cell at level lower than 8!!!!!!'
!!$                stop
!!$             endif
!--JB

             ! there has been no interaction in the cell, tau_abs has been updated in gas_get_scatter_flag
             ! => move to next cell
#ifdef DEBUG
             print*,'--> exit cell, photon moves to cell',icellnew,ileaf,ind,ioct
#endif

             exit propag_in_cell

          else
             ! Next event happens inside this cell

             ! length and time travelled by the photon before event
             d    = distance_to_border_cm   ! NB: at this point, distance_to_border became "distance_to_interaction" in gas_get_scatter_flag
             time = time + d/clight
             d    = d / cell_size_cm        ! in box units

             ! update ppos_cell
             do i=1,3
                ppos_cell(i) = ppos_cell(i) + p%k(i) * d
             enddo
             ! update ppos according to ppos_cell
             ppos = ppos_cell * cell_size + cell_corner

             ! Check if photon is still in the computational domain
             if(.not.(cell_fully_in_domain))then   ! this flag allows to not check at each scattering the new position 
                                                   ! as long as the cell is fully contained in the computational domain
                in_domain = domain_contains_point(ppos,domaine_calcul)
                if(.not.(in_domain))then    ! photon done, nothing else to do
                   p%status       = 1
                   p%xcurr        = ppos
                   ! correct time
                   dborder = domain_distance_to_border(ppos,domaine_calcul)
                   dtime   = dborder*box_size_cm/clight ! should be negative
                   time    = time+dtime
                   p%time         = time
                   p%tau_abs_curr = tau_abs
                   p%iran         = iran
#ifdef DEBUG
                   print*,'-exit propagation, photon escaped comp. domain'
#endif
                   exit photon_propagation
                endif
             endif
             !------------
             ! scattering
             !------------

             p%nb_abs = p%nb_abs + 1     ! increment nb of scatterings
             p%xlast = ppos              ! memorize the location of "potential" last interaction
             ! scattering => modify nu_cell, k, nu_ext
             ! so it needs to transport nu_cell, k, nu_ext, not p (because type p not known in gas)
             nu_ext = p%nu_ext
             k = p%k
             call gas_scatter(scatter_flag, cell_gas, nu_cell, k, nu_ext, iran)    ! NB: nu_cell, k, nu_ext, and iran en inout

#ifdef DEBUG
!             print*,'scattering impossible...'
             print*,'stop after one scattering...'
             stop
#endif

             ! if(p%nb_abs==3)then
             !    print*,ppos
             !    print*,ppos_cell
             !    print*,k
             !    print*,nu_cell
             !    print*,nu_ext
             !    print*,p%nb_abs
             !    print*,tau_abs
             !    !stop
             ! endif

             p%nu_ext = nu_ext
             ! for TEST case, to have photons propagating straight on, comment the following line
             p%k = k
             ! there has been an interaction -> reset tau_abs
             tau_abs = -1.0d0

             ! scatter_flag allows to know the status (aborbed or not) of the photon in case of dust
             if(scatter_flag==4)then
                ! photon has been absorbed by dust, photon done, nothing else to do
                p%status       = 2
                p%xcurr        = ppos
                p%time         = time
                p%tau_abs_curr = tau_abs
                p%iran         = iran
#ifdef DEBUG
                print*,'-exit propagation, photon is dead'
#endif
                exit photon_propagation
             endif

             !if(mod(p%nb_abs,100)==0)then
             !   print*,p%nb_abs
             !   print*,ppos
             !   print*,ppos_cell
             !endif
          
          end if

       end do propag_in_cell
#ifdef DEBUG
       print*,'nb scattering in cell',p%nb_abs
#endif
    end do photon_propagation

    ! fin de la propagation du photon p => soit le photon sort du domaine (du cpu ou de calcul), soit le photon est mort
    ! un simple test 
    ! soit out of computational domain -> in_domain=.false. == p%status=1
    ! soit out of mesh-cpu domain      -> flagoutvol==.true.
    ! soit death                       -> p%status=2
    ! sanity check
    !if(.not.(flagoutvol).or..not.(p%status==1).or..not.(p%status==2))then
    if(.not.(flagoutvol).and.(p%status==0))then
       print *,'--> problem 1 with photon propagation in module_photon.f90!',flagoutvol,p%status
       stop
    endif

    if(.not.(flagoutvol.or.(p%status==1).or.(p%status==2).or.(p%status==3)))then
       ! not that status=3 is still problematic...
       print *,'--> problem 2 with photon propagation in module_photon.f90!',flagoutvol,p%status
       stop
    endif

  end subroutine propagate
  




  function path(pos,dir)

    ! compute distance to border of a cell (in cell units), from position
    ! pos (in cell units) and in direction dir. 
    
    implicit none

    real(kind=8),intent(in) :: pos(3)   ! position of photon in cell units
    real(kind=8),intent(in) :: dir(3)   ! propagation direction of photon
    integer(kind=4)         :: i
    real(kind=8)            :: dx(3)
    real(kind=8)            :: path     ! distance from pos to exit point

    do i = 1,3
       if(dir(i) < 0.) then
          dx(i) = -pos(i) / dir(i)
       else if (dir(i) > 0.) then
          dx(i) = (1.0d0 - pos(i)) / dir(i)
       else ! dir(i) == 0
          dx(i) = 10.  ! larger than maximum extent of cell (sqrt(3)) in cell units
       end if
    end do
    path = minval(dx)

    return
    
  end function path



  subroutine init_photons_from_file(file,pgrid)

    character(2000),intent(in)                                 :: file
    type(photon_current),dimension(:),allocatable, intent(out) :: pgrid
    type(photon_init),dimension(:),allocatable                 :: pgridinit
    integer(kind=4)                                            :: i, n_photon, iseed

    ! read ICs
    open(unit=14, file=trim(file), status='unknown', form='unformatted', action='read')
    read(14) n_photon
    allocate(pgridinit(n_photon))
    read(14) iseed
    read(14) (pgridinit(i)%ID,i=1,n_photon)
    read(14) (pgridinit(i)%nu_em,i=1,n_photon)
    read(14) (pgridinit(i)%x_em(:),i=1,n_photon)
    read(14) (pgridinit(i)%k_em(:),i=1,n_photon)
    read(14) (pgridinit(i)%iran,i=1,n_photon)
    close(14)

    ! build photgrid current
    allocate(pgrid(n_photon))
    do i=1,n_photon
       pgrid(i)%ID           = pgridinit(i)%ID
       pgrid(i)%status       = 0
       !pgrid(i)%xlast        = (/-99.,-99.,-99./)
       pgrid(i)%xlast        = pgridinit(i)%x_em
       pgrid(i)%xcurr        = pgridinit(i)%x_em
       pgrid(i)%nu_ext       = pgridinit(i)%nu_em
       pgrid(i)%k            = pgridinit(i)%k_em
       pgrid(i)%nb_abs       = 0
       pgrid(i)%time         = 0.d0
       pgrid(i)%tau_abs_curr = -1.
       pgrid(i)%iran         = pgridinit(i)%iran
    enddo
    deallocate(pgridinit)

  end subroutine init_photons_from_file



  subroutine dump_photons(file,pgrid)

    character(2000),intent(in)                   :: file
    type(photon_current),dimension(:),intent(in) :: pgrid
    integer(kind=4)                              :: i,np

    np = size(pgrid)
    open(unit=14, file=trim(file), status='unknown', form='unformatted', action='write')
    write(14) np
    write(14) (pgrid(i)%ID,      i=1,np)
    write(14) (pgrid(i)%status,  i=1,np)
    write(14) (pgrid(i)%xlast(:),i=1,np)
    write(14) (pgrid(i)%nu_ext,  i=1,np)
    write(14) (pgrid(i)%k(:),    i=1,np)
    write(14) (pgrid(i)%nb_abs,  i=1,np)
    write(14) (pgrid(i)%time,    i=1,np)
    close(14)

  end subroutine dump_photons

end module module_photon
