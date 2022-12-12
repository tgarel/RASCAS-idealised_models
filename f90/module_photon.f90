module module_photon
 
  use module_gas_composition
  use module_mesh
  use module_constants
  use module_random
  use module_domain
  use module_utils, only: path 
  !--PEEL--
  use module_mock
  !--LEEP--
  
  implicit none

  ! 2 types for photons, one for the initial properties is called photon_init
  ! and one for the properties that evolve during the RT called photon_current

  type photon_current
     integer(kind=4)           :: ID
     integer(kind=4)           :: status       ! =0 if flying, =1 if escape, =2 if absorption (by dust)
     real(kind=8),dimension(3) :: xlast        ! coordinates of last interaction in box units
     real(kind=8),dimension(3) :: xcurr        ! current position of the photon in box units
     real(kind=8)              :: nu_ext       ! external frame frequency (Hz)
     real(kind=8),dimension(3) :: k            ! normalised propagation vector 
     integer(kind=4)           :: nb_abs       ! number of interactions before escape
     real(kind=8)              :: time         ! time in [s] from emission to escape/absorption        
     real(kind=8)              :: tau_abs_curr ! current optical depth (useful when photon change mesh domain)
     integer(kind=4)           :: iran         ! state of the random generator
     real(kind=8),dimension(3) :: v_src        ! velocity of the source -- for peeling off
  end type photon_current
  ! Note: if you change something in photon_current, don't forget to update the mpi_photon_type in module_parallel_mpi.f90


  type photon_init
     integer(kind=4)           :: ID
     real(kind=8)              :: nu_em    ! emitted frequency in external frame
     real(kind=8),dimension(3) :: x_em
     real(kind=8),dimension(3) :: k_em
     integer(kind=4)           :: iran     ! state of the random generator
     real(kind=8),dimension(3) :: v_em     ! velocity of the source -- for peeling off
  end type photon_init

  !--PEEL--
  type peel
     real(kind=8)              :: nu     ! frequency before scattering, converted to freq. after virtual scat. towards dir. of observation. 
     real(kind=8),dimension(3) :: x      ! position at scattering
     real(kind=8),dimension(3) :: kin    ! incoming photon direcrtion
     integer(kind=4)           :: icell  ! cell in which scattering occurs, saves one search... 
     real(kind=8)              :: weight ! probability of being re-emitted in obs direction
     integer(kind=4)           :: scatter_flag ! (if negative, this is a emission peel-> specific processing)
     real(kind=8),dimension(3) :: v_src  ! velocity of the source -- for peeling off
  end type peel
  integer(kind=4),parameter            :: PeelBufferSize = 1000000
  type(peel),dimension(PeelBufferSize) :: PeelBuffer
  integer(kind=4)                      :: nPeeled
  !--LEEP--

  public  :: MCRT, propagate, init_photons_from_file, dump_photons

contains


  subroutine MCRT(npp,photpacket,mesh_dom,compute_dom)

    ! this is the Monte Carlo Radiative transfer routine... but also a single loop over photons...

    integer(kind=4),intent(in)                        :: npp
    type(photon_current),dimension(npp),intent(inout) :: photpacket
    type(mesh),intent(in)                             :: mesh_dom
    type(domain),intent(in)                           :: compute_dom
    integer(kind=4)                                   :: i

    do i=1,npp
       ! case of photpacket not fully filled...
       if (photpacket(i)%ID>0) then
          call propagate(photpacket(i),mesh_dom,compute_dom)
       endif
    enddo
 
  end subroutine MCRT



  subroutine propagate(p,domesh,domaine_calcul)

    type(photon_current),intent(inout)   :: p              ! photon 
    type(mesh),intent(in)                :: domesh         ! mesh
    type(domain),intent(in)              :: domaine_calcul ! computational domain in which photons are propagating
    
    type(gas)                            :: cell_gas                   ! gas in the current cell 
    integer(kind=4)                      :: icell, ioct, ind, ileaf, cell_level  ! current cell indices and level
    real(kind=8)                         :: cell_size, cell_size_cm, tau_abs, scalar, nu_cell, nu_ext, rtau
    real(kind=8),dimension(3)            :: ppos,ppos_cell             ! working coordinates of photon (in box and in cell units)
    real(kind=8)                         :: distance_to_border,distance_to_border_cm, d,distance_to_border_box_units
    real(kind=8)                         :: time
    integer(kind=4)                      :: scatter_flag, i, icellnew, iran, npush
    real(kind=8),dimension(3)            :: vgas, k, cell_corner, posoct, pcell
    logical                              :: cell_fully_in_domain, flagoutvol, in_domain, OutOfDomainBeforeCell
    real(kind=8)                         :: dborder, dborder_cm, error
    !--CORESKIP--
    real(kind=8)                         :: xcrit,dist_cm
    !--PIKSEROC--
    
    
    ! initialise working props of photon
    ppos    = p%xcurr        ! position within full simulation box, in box units.
    time    = p%time
    tau_abs = p%tau_abs_curr
    iran    = p%iran

    ! find cell in which the photon is, and define all its indices
    icell = in_cell_finder(domesh,ppos)
    if(domesh%son(icell)>=0)then
       print*,'ERROR: not a leaf cell'
       stop
    endif
    ileaf = - domesh%son(icell)
    ind   = (icell - domesh%nCoarse - 1) / domesh%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
    ioct  = icell - domesh%nCoarse - (ind - 1) * domesh%nOct
    flagoutvol = .false.

    !--PEEL-- initialise counter of peels ...
    if (peeling_off) then
       ! Start with a peel off initial photon 
       nPeeled = 1
       PeelBuffer(nPeeled)%x      = ppos  ! position (at scattering)
       PeelBuffer(nPeeled)%icell  = icell
       ! compute first term of ray's weights : the peeling-off strategy
       PeelBuffer(nPeeled)%weight = 0.5      ! assume isotropy for this particular (non-)event
       PeelBuffer(nPeeled)%nu     = p%nu_ext ! frequency in the direction of observation
       PeelBuffer(nPeeled)%scatter_flag = -1 ! flag peel as an initial peel
       ! JB-
       PeelBuffer(nPeeled)%kin = p%k ! emission direction (nu is in this direction, not in the direction to mock observer)
       ! -JB 
       PeelBuffer(nPeeled)%v_src = p%v_src
       if (nPeeled == PeelBufferSize) then ! buffer is full -> process.
          call process_peels(domesh,domaine_calcul,iran)
          nPeeled=0
       endif
    end if
    !--LEEP-- 

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
       ppos_cell    = (ppos - cell_corner) / cell_size         ! position of photon in cell units (x,y,z in [0,1] within cell)
       if((ppos_cell(1)>1.0d0).or.(ppos_cell(2)>1.0d0).or.(ppos_cell(3)>1.0d0).or. &
            (ppos_cell(1)<0.0d0).or.(ppos_cell(2)<0.0d0).or.(ppos_cell(3)<0.0d0))then
          print*,"ERROR: problem in computing ppos_cell"
          stop
       endif

       ! get gas velocity (in cgs units)
       vgas         = get_gas_velocity(cell_gas)
       ! compute photon's frequency in cell's moving frame
       scalar       = p%k(1) * vgas(1) + p%k(2) * vgas(2) + p%k(3) * vgas(3)
       nu_cell      = (1.0d0 - scalar/clight) * p%nu_ext  

       ! define/update flag_cell_fully_in_comp_dom to avoid various tests in the following
       pcell = cell_corner + 0.5d0*cell_size
       cell_fully_in_domain = domain_fully_contains_cell(pcell,cell_size,domaine_calcul)

       propag_in_cell : do
          
          ! generate the opt depth where the photon is scattered/absorbed
          if (tau_abs <= 0.0d0) then
             rtau    = ran3(iran)
             tau_abs = -log(1.0d0-rtau+1.0d-30)
          end if

          ! compute distance of photon to border of cell along propagation direction
          distance_to_border           = path(ppos_cell,p%k)                   ! in cell units
          distance_to_border_cm        = distance_to_border * cell_size_cm     ! cm
          distance_to_border_box_units = distance_to_border * cell_size        ! in box units
          ! if cell not fully in domain, modify distance_to_border to "distance_to_domain_border" if relevant
          OutOfDomainBeforeCell = .False.
          if(.not.(cell_fully_in_domain))then
             dborder    = domain_distance_to_border_along_k(ppos,p%k,domaine_calcul)  ! in box units
             dborder_cm = dborder * box_size_cm                                       ! from box units to cm
             ! compare distance to cell border and distance to domain border and take the min
             if (dborder_cm < distance_to_border_cm) then
                OutOfDomainBeforeCell        = .True.
                distance_to_border_cm        = dborder_cm
                distance_to_border_box_units = dborder
                distance_to_border           = distance_to_border_cm / cell_size_cm
             end if
          endif
          
          ! check whether scattering occurs within cell or domain (scatter_flag > 0) or not (scatter_flag==0)
          !--CORESKIP--
          ! also compute xcrit for core-skipping if needed
          ! for core-skipping we are interested in the true distance to cell border, NOT the distance along k 
          dist_cm = min(ppos_cell(1),1.0d0-ppos_cell(1),ppos_cell(2),1.0d0-ppos_cell(2),ppos_cell(3),1.0d0-ppos_cell(3)) * cell_size_cm
          !scatter_flag = gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs, iran)
          scatter_flag = gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs, iran, dist_cm, xcrit)
          !--PIKSEROC--

          if (scatter_flag == 0) then   ! next scattering event will not occur in the cell or in the domain

             ! move photon out of cell or domain
             ppos = ppos + p%k * distance_to_border_box_units *(1.0d0 + epsilon(1.0d0))

             ! correct for periodicity
             do i=1,3
                if (ppos(i) < 0.0d0) ppos(i)=ppos(i)+1.0d0
                if (ppos(i) > 1.0d0) ppos(i)=ppos(i)-1.0d0
             enddo
             ! update travel time
             time = time + distance_to_border_cm/clight

             if (OutOfDomainBeforeCell) then ! photon exits computational domain and is done 
                ! it may happen due to numerical precision that the photon is still in the domain despite epsilon above.
                ! -> check and issue warning if it is the case. The error should not be larger than a few times epsilon. 
                in_domain = domain_contains_point(ppos,domaine_calcul)
                if (in_domain) then
                   if (domain_distance_to_border_along_k(ppos,p%k,domaine_calcul)>3.d0*epsilon(distance_to_border)) then  
                      print*,'WARNING : photon still in domain when it should not ... '
                      error = nint(domain_distance_to_border_along_k(ppos,p%k,domaine_calcul)/epsilon(distance_to_border))
                      print*,'          (error ~ ',error,' times num. prec.) '
                   end if
                end if
                p%status       = 1
                p%xcurr        = ppos
                p%time         = time
                p%tau_abs_curr = tau_abs
                p%iran         = iran
                exit photon_propagation

             else
                
                ! photon exits current cell -> find into which new cell it goes
                call whereIsPhotonGoing(domesh,icell,ppos,icellnew,flagoutvol)
                ! It may happen due to numerical precision that the photon is still in the current cell (i.e. icell == icellnew).
                ! -> give it an extra push untill it is out. 
                npush = 0
                do while (icell==icellnew)
                   npush = npush + 1
                   ppos(1) = ppos(1) + merge(-1.0d0,1.0d0,p%k(1)<0.0d0) * epsilon(ppos(1))
                   ppos(2) = ppos(2) + merge(-1.0d0,1.0d0,p%k(2)<0.0d0) * epsilon(ppos(2))
                   ppos(3) = ppos(3) + merge(-1.0d0,1.0d0,p%k(3)<0.0d0) * epsilon(ppos(3))
                   ! correct for periodicity
                   do i=1,3
                      if (ppos(i) < 0.0d0) ppos(i)=ppos(i)+1.0d0
                      if (ppos(i) > 1.0d0) ppos(i)=ppos(i)-1.0d0
                   enddo
                   call whereIsPhotonGoing(domesh,icell,ppos,icellnew,flagoutvol)
                end do
                if (npush > 1) print*,'WARNING : npush > 1 needed in module_photon:propagate.'
                ! test whether photon was pushed out of domain with the extra pushes
                ! (and in that case, call it done). 
                if (npush > 0) then 
                   in_domain = domain_contains_point(ppos,domaine_calcul)
                   if (.not. in_domain) then
                      print*,'WARNING: pushed photon outside domain ... '
                      p%status       = 1
                      p%xcurr        = ppos
                      p%time         = time
                      p%tau_abs_curr = tau_abs
                      p%iran         = iran
                      exit photon_propagation
                   end if
                end if
                ! check if the new cell is outside of the current cpu domain (flagoutvol)
                ! And if so send it back to master
                if(flagoutvol)then
                   p%xcurr        = ppos
                   p%time         = time
                   p%tau_abs_curr = tau_abs
                   p%iran         = iran
                   exit photon_propagation
                endif
                ! Finally, if we're here, the photon entered a cell within the current cpu domain so we go on. 
                icell = icellnew
                ileaf = - domesh%son(icell)
                ind   = (icell - domesh%nCoarse - 1) / domesh%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
                ioct  = icell - domesh%nCoarse - (ind - 1) * domesh%nOct

                ! there has been no interaction in the cell, tau_abs has been updated in gas_get_scatter_flag
                ! -> move to next cell
                exit propag_in_cell

             end if


          else
             ! Next event happens inside this cell and in the domain.

             ! length and time travelled by the photon before event
             d    = distance_to_border_cm   ! NB: at this point, distance_to_border became "distance_to_interaction" in gas_get_scatter_flag
             time = time + d/clight
             d    = d / cell_size_cm        ! in cell units

             ! update ppos_cell
             do i=1,3
                ppos_cell(i) = ppos_cell(i) + p%k(i) * d
             enddo
             ! update ppos according to ppos_cell
             ppos = ppos_cell * cell_size + cell_corner

             !------------
             ! scattering
             !------------
             !--PEEL--
             if (peeling_off) then 
                ! save ray for batch processing later.
                nPeeled = nPeeled + 1
                PeelBuffer(nPeeled)%x     = ppos  ! position (at scattering)
                PeelBuffer(nPeeled)%icell = icell    
                PeelBuffer(nPeeled)%kin   = p%k      ! incoming photon direction 
                PeelBuffer(nPeeled)%nu    = p%nu_ext ! incoming photon frequency
                PeelBuffer(nPeeled)%scatter_flag = scatter_flag
                if (nPeeled == PeelBufferSize) then ! buffer is full -> process.
                   call process_peels(domesh,domaine_calcul,iran)
                   nPeeled=0
                end if
             end if
             !--LEEP-- 

             p%nb_abs = p%nb_abs + 1     ! increment nb of scatterings
             p%xlast = ppos              ! memorize the location of "potential" last interaction
             ! a scattering event modifies nu_cell, k, and nu_ext
             ! so it needs to transport nu_cell, k, nu_ext, but not type p (because type p not known in gas)
             nu_ext = p%nu_ext
             k = p%k
             !--CORESKIP--
             call gas_scatter(scatter_flag, cell_gas, nu_cell, k, nu_ext, iran, xcrit)    ! NB: nu_cell, k, nu_ext, and iran en inout
             !call gas_scatter(scatter_flag, cell_gas, nu_cell, k, nu_ext, iran)    ! NB: nu_cell, k, nu_ext, and iran en inout             
             !--PIKSEROC--
             p%nu_ext = nu_ext
             ! NB: for TEST case, to have photons propagating straight on, comment the following line
             p%k = k
             ! there has been an interaction -> reset tau_abs
             tau_abs = -1.0d0
             ! scatter_flag allows to know the status (aborbed or not) of the photon in case of dust
             ! new convention: negative if absorbed
             if(scatter_flag<0)then
                ! photon has been absorbed by dust, photon done, nothing else to do
                p%status       = 2
                p%xcurr        = ppos
                p%time         = time
                p%tau_abs_curr = tau_abs
                p%iran         = iran
                exit photon_propagation
             endif

          end if

       end do propag_in_cell

    end do photon_propagation
    !--PEEL--
    ! finish processing peel buffer before moving to next photon packet. 
    if (peeling_off) then 
       if (nPeeled > 0) then ! buffer is not empty -> process.
          call process_peels(domesh,domaine_calcul,iran)
          nPeeled=0
       endif
    end if
    !--LEEP-- 

    
    ! End of the photon propagation. There are 3 possible cases:
    !   1/ photon is out of the computational domain == escaped           -> in_domain=.false. && p%status=1
    !   2/ photon is out of mesh-cpu domain -> sent back to master, etc.  -> flagoutvol==.true.
    !   3/ photon is dead                                                 -> p%status=2
    
    ! some simple sanity checks
    if(.not.(flagoutvol).and.(p%status==0))then
       print*,'ERROR: problem 1 with photon propagation in module_photon.f90!',flagoutvol,p%status
       stop
    endif

    if(.not.(flagoutvol.or.(p%status==1).or.(p%status==2)))then
       print *,'ERROR: problem 2 with photon propagation in module_photon.f90!',flagoutvol,p%status
       stop
    endif

  end subroutine propagate
  
  !--PEEL--
  function tau_to_border(p,domesh,domaine_calcul,tau_max,kobs)

    type(peel),intent(inout)             :: p              ! peel 
    type(mesh),intent(in)                :: domesh         ! mesh
    type(domain),intent(in)              :: domaine_calcul ! computational domain in which photons are propagating
    real(kind=8),intent(in)              :: tau_max    ! stop computation when tau reaches tau_max
    real(kind=8),intent(in)              :: kobs(3) ! direction to observer 
    real(kind=8)                         :: tau_to_border
    real(kind=8)                         :: tau_cell
    type(gas)                            :: cell_gas                   ! gas in the current cell 
    integer(kind=4)                      :: icell, ioct, ind, ileaf, cell_level  ! current cell indices and level
    real(kind=8)                         :: cell_size, cell_size_cm, scalar, nu_cell
    real(kind=8),dimension(3)            :: ppos,ppos_cell             ! working coordinates of photon (in box and in cell units)
    real(kind=8)                         :: distance_to_border,distance_to_border_cm, distance_to_border_box_units
    integer(kind=4)                      :: i, icellnew, npush
    real(kind=8),dimension(3)            :: vgas, cell_corner, posoct, pcell
    logical                              :: cell_fully_in_domain, flagoutvol, in_domain, OutOfDomainBeforeCell
    real(kind=8)                         :: dborder, dborder_cm

    tau_to_border = 0.0d0 
    
    ! initialise working props of peel (aka ray)
    ppos  = p%x        ! position within full simulation box, in box units.
    icell = p%icell 

    ! define cell indices
    ileaf = - domesh%son(icell)
    ind   = (icell - domesh%nCoarse - 1) / domesh%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
    ioct  = icell - domesh%nCoarse - (ind - 1) * domesh%nOct
    flagoutvol = .false.

    ! propagate peel to border of computational domain
    photon_propagation : do 

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
          print*,"ERROR: problem in computing ppos_cell",ppos_cell
          stop
       endif
       
       ! get gas velocity (in cgs units)
       vgas         = get_gas_velocity(cell_gas)
       ! compute photon's frequency in cell's moving frame
       scalar       = kobs(1) * vgas(1) + kobs(2) * vgas(2) + kobs(3) * vgas(3)
       nu_cell      = (1.0d0 - scalar/clight) * p%nu 

       ! define/update flag_cell_fully_in_comp_dom to avoid various tests in the following
       pcell = cell_corner + 0.5d0*cell_size
       cell_fully_in_domain = domain_fully_contains_cell(pcell,cell_size,domaine_calcul)
       
       ! compute distance of photon to border of cell along propagation direction
       distance_to_border           = path(ppos_cell,kobs)                   ! in cell units
       distance_to_border_cm        = distance_to_border * cell_size_cm     ! cm
       distance_to_border_box_units = distance_to_border * cell_size        ! in box units
       ! if cell not fully in domain, modify distance_to_border to "distance_to_domain_border" if relevant
       OutOfDomainBeforeCell = .False.
       if(.not.(cell_fully_in_domain))then
          dborder    = domain_distance_to_border_along_k(ppos,kobs,domaine_calcul)  ! in box units
          dborder_cm = dborder * box_size_cm                                        ! from box units to cm
          ! compare distance to cell border and distance to domain border and take the min
          if (dborder_cm < distance_to_border_cm) then
             OutOfDomainBeforeCell        = .True.
             distance_to_border_cm        = dborder_cm
             distance_to_border_box_units = dborder
             distance_to_border           = distance_to_border_cm / cell_size_cm
          end if
       endif

       ! compute (total) optical depth along ray in cell 
       tau_cell = gas_get_tau(cell_gas, distance_to_border_cm, nu_cell)

       tau_to_border = tau_to_border + tau_cell
       if (tau_to_border > tau_max) then
          exit photon_propagation
       end if

       ! advance photon
       ppos = ppos + kobs * distance_to_border_box_units *(1.0d0 + epsilon(1.0d0))
       ! correct for periodicity
       do i=1,3
          if (ppos(i) < 0.0d0) ppos(i)=ppos(i)+1.0d0
          if (ppos(i) > 1.0d0) ppos(i)=ppos(i)-1.0d0
       enddo


       if (OutOfDomainBeforeCell) then ! photon exits computational domain and is done 
          exit photon_propagation
       end if

       ! photon moves to next cell 
       call whereIsPhotonGoing(domesh,icell,ppos,icellnew,flagoutvol)

       ! It may happen due to numerical precision that the photon is still in the current cell (i.e. icell == icellnew).
       ! -> give it an extra push untill it is out. 
       npush = 0
       do while (icell==icellnew)
          npush = npush + 1
          ppos(1) = ppos(1) + merge(-1.0d0,1.0d0,kobs(1)<0.0d0) * epsilon(ppos(1))
          ppos(2) = ppos(2) + merge(-1.0d0,1.0d0,kobs(2)<0.0d0) * epsilon(ppos(2))
          ppos(3) = ppos(3) + merge(-1.0d0,1.0d0,kobs(3)<0.0d0) * epsilon(ppos(3))
          ! correct for periodicity
          do i=1,3
             if (ppos(i) < 0.0d0) ppos(i)=ppos(i)+1.0d0
             if (ppos(i) > 1.0d0) ppos(i)=ppos(i)-1.0d0
          enddo
          call whereIsPhotonGoing(domesh,icell,ppos,icellnew,flagoutvol)
          if (npush == 10) then
             print*,'npush == 10 ... using les grands moyens ... '
             ppos(1) = ppos(1) + kobs(1)*(1000.*epsilon(ppos(1)))
             ppos(2) = ppos(2) + kobs(2)*(1000.*epsilon(ppos(2)))
             ppos(3) = ppos(3) + kobs(3)*(1000.*epsilon(ppos(3)))
          end if
          if (npush == 11) then
             print*,'npush == 11 ... using les very grands moyens ... '
             ppos(1) = ppos(1) + kobs(1)*(1e5*epsilon(ppos(1)))
             ppos(2) = ppos(2) + kobs(2)*(1e5*epsilon(ppos(2)))
             ppos(3) = ppos(3) + kobs(3)*(1e5*epsilon(ppos(3)))
          end if
          if (npush > 11) then
             print*,'oh well ...'
             stop
          end if
       end do
       if (npush > 1) print*,'WARNING : npush > 1 needed in module_photon:propagate.'
       ! test whether photon was pushed out of domain with the extra pushes
       ! (and in that case, call it done). 
       if (npush > 0) then 
          in_domain = domain_contains_point(ppos,domaine_calcul)
          if (.not. in_domain) then
             print*,'WARNING: pushed photon outside domain ... '
             exit photon_propagation
          end if
       end if
       ! check if the new cell is outside of the current cpu domain (flagoutvol)
       ! And if so send it back to master
       if(flagoutvol)then
          print*,'ERROR: peeling off works with a single domain for now... '
          stop
       endif
       ! Finally, if we're here, the photon entered a cell within the current cpu domain so we go on. 
       icell = icellnew
       ileaf = - domesh%son(icell)
       ind   = (icell - domesh%nCoarse - 1) / domesh%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
       ioct  = icell - domesh%nCoarse - (ind - 1) * domesh%nOct
       
    end do photon_propagation

  end function tau_to_border
  !--LEEP--
  
  !--PEEL--
  subroutine process_peels(domesh,domaine_calcul,iran)
    implicit none
    type(mesh),intent(in)         :: domesh         ! mesh
    type(domain),intent(in)       :: domaine_calcul ! computational domain in which photons are propagating
    integer(kind=4),intent(inout) :: iran 
    real(kind=8),parameter        :: tau_max = 60   ! stop computation when tau reaches tau_max ... 
    integer(kind=4)               :: ipeel,idir,ileaf
    real(kind=8)                  :: tau,peel_contrib
    real(kind=8)                  :: projpos(2),kobs(3)
    logical                       :: increment_flux, increment_spec, increment_image, increment_cube
    type(gas)                     :: cell_gas       ! gas in the current cell 
    real(kind=8)                  :: nupeel, nu_src, nu_ext, scalar

    peels_count = peels_count+nPeeled
    do idir = 1,nDirections
       kobs = mock_line_of_sight(idir)
       do ipeel = 1,nPeeled
          ! if projected peel falls into one of the detector compute tau to detector.
          call mock_projected_pos(PeelBuffer(ipeel)%x,projpos,idir)
          increment_flux  = mock_point_in_flux_aperture(projpos,idir)
          increment_spec  = mock(idir)%compute_spectrum .and. mock_point_in_spectral_aperture(projpos,idir)
          increment_image = mock(idir)%compute_image .and. mock_point_in_image(projpos,idir)
          increment_cube  = mock(idir)%compute_cube .and. mock_point_in_cube(projpos,idir)
          tau = tau_max
          nupeel = PeelBuffer(ipeel)%nu ! save to restore for next directions
          if (increment_flux .or. increment_spec .or. increment_image .or. increment_cube) then
             rays_count = rays_count+1
             if (PeelBuffer(ipeel)%scatter_flag > 0) then 
                ileaf    = - domesh%son(PeelBuffer(ipeel)%icell)
                cell_gas = domesh%gas(ileaf)
                ! NB: the following line updates peel%nu to the frequency in the direction of observation. 
                PeelBuffer(ipeel)%weight = gas_peeloff_weight(PeelBuffer(ipeel)%scatter_flag, cell_gas, PeelBuffer(ipeel)%nu, PeelBuffer(ipeel)%kin, kobs, iran)
             end if
             ! JB-
             if (PeelBuffer(ipeel)%scatter_flag < 0) then ! this is initialisation (i.e. a peel from emission site)
                ! -> re-compute frequency in the direction of mock observation...
                ! ---> nu is initialised as frequency in external frame, in the direction of emission.
                ! ---> Go back to emitting source frame (cell, stars, or model) and then to external frame using direction of mock instead. 
                ! get nu_source from nu_ext
                scalar = PeelBuffer(ipeel)%kin(1) * PeelBuffer(ipeel)%v_src(1) + &
                     PeelBuffer(ipeel)%kin(2) * PeelBuffer(ipeel)%v_src(2) + PeelBuffer(ipeel)%kin(3) * PeelBuffer(ipeel)%v_src(3)
                nu_src = PeelBuffer(ipeel)%nu * (1.0d0 - scalar/clight) 
                ! get nu_ext from nu_source now in the direction of observation
                scalar = PeelBuffer(ipeel)%v_src(1)*kobs(1) + PeelBuffer(ipeel)%v_src(2)*kobs(2) + PeelBuffer(ipeel)%v_src(3)*kobs(3)
                nu_ext = (1.0d0 + scalar/clight) * nu_src
                PeelBuffer(ipeel)%nu = nu_ext
             end if
             ! -JB
             tau = tau_to_border(PeelBuffer(ipeel),domesh,domaine_calcul,tau_max,kobs)
          end if
          ! if tau is not absurdly large, increment detectors 
          if (tau < tau_max) then
             detectors_count(idir) = detectors_count(idir) + 1
             peel_contrib = PeelBuffer(ipeel)%weight * exp(-tau) * 2.0d0 
             if (increment_flux)  call peel_to_flux(PeelBuffer(ipeel)%nu,peel_contrib,idir) 
             if (increment_spec)  call peel_to_spec(PeelBuffer(ipeel)%nu,peel_contrib,idir)
             if (increment_image) call peel_to_map(projpos,PeelBuffer(ipeel)%nu,peel_contrib,idir)
             if (increment_cube)  call peel_to_cube(projpos,PeelBuffer(ipeel)%nu,peel_contrib,idir)
          end if
          PeelBuffer(ipeel)%nu = nupeel ! restore for next directions
       end do
    end do
    
  end subroutine process_peels
  !--LEEP--


  subroutine init_photons_from_file(file,pgrid)

    character(2000),intent(in)                                 :: file
    type(photon_current),dimension(:),allocatable, intent(out) :: pgrid
    type(photon_init),dimension(:),allocatable                 :: pgridinit
    integer(kind=4)                                            :: i, n_photon, iseed
    real(kind=8)                                               :: total_flux,knorm

    ! read ICs
    open(unit=14, file=trim(file), status='unknown', form='unformatted', action='read')
    read(14) n_photon
    read(14) total_flux    ! nb of real photons [# / s]
    allocate(pgridinit(n_photon))
    read(14) iseed
    read(14) (pgridinit(i)%ID,i=1,n_photon)
    read(14) (pgridinit(i)%nu_em,i=1,n_photon)
    read(14) (pgridinit(i)%x_em(:),i=1,n_photon)
    read(14) (pgridinit(i)%k_em(:),i=1,n_photon)
    read(14) (pgridinit(i)%iran,i=1,n_photon)
    read(14) (pgridinit(i)%v_em(:),i=1,n_photon)
    close(14)

    ! build photgrid current
    allocate(pgrid(n_photon))
    do i=1,n_photon
       pgrid(i)%ID           = pgridinit(i)%ID
       pgrid(i)%status       = 0
       pgrid(i)%xlast        = pgridinit(i)%x_em
       pgrid(i)%xcurr        = pgridinit(i)%x_em
       pgrid(i)%nu_ext       = pgridinit(i)%nu_em
       ! make sure k is normalised
       knorm                 = sqrt(pgridinit(i)%k_em(1)**2 + pgridinit(i)%k_em(2)**2 + pgridinit(i)%k_em(3)**2)
       pgrid(i)%k            = pgridinit(i)%k_em / knorm
       pgrid(i)%nb_abs       = 0
       pgrid(i)%time         = 0.0d0
       pgrid(i)%tau_abs_curr = -1.0d0
       pgrid(i)%iran         = pgridinit(i)%iran
       pgrid(i)%v_src        = pgridinit(i)%v_em
    enddo
    deallocate(pgridinit)

  end subroutine init_photons_from_file



  subroutine restore_photons(file,pgrid)

    character(2000),intent(in)                                 :: file
    type(photon_current),dimension(:),allocatable, intent(out) :: pgrid
    integer(kind=4)                                            :: i, np

    ! restore photons from saved file
    open(unit=16, file=trim(file), status='old', form='unformatted', action='read')
    read(16) np

    allocate(pgrid(np))

    read(16) (pgrid(i)%ID,           i=1,np)
    read(16) (pgrid(i)%status,       i=1,np)
    read(16) (pgrid(i)%xlast(:),     i=1,np)
    read(16) (pgrid(i)%xcurr(:),     i=1,np)
    read(16) (pgrid(i)%nu_ext,       i=1,np)
    read(16) (pgrid(i)%k(:),         i=1,np)
    read(16) (pgrid(i)%nb_abs,       i=1,np)
    read(16) (pgrid(i)%time,         i=1,np)
    read(16) (pgrid(i)%tau_abs_curr, i=1,np)
    read(16) (pgrid(i)%iran,         i=1,np)
    read(16) (pgrid(i)%v_src(:),     i=1,np)

    close(16)

  end subroutine restore_photons



  subroutine save_photons(file,pgrid)

    character(2000),intent(in)                   :: file
    type(photon_current),dimension(:),intent(in) :: pgrid
    integer(kind=4)                              :: i,np

    np = size(pgrid)
    open(unit=16, file=trim(file), status='unknown', form='unformatted', action='write')
    write(16) np
    write(16) (pgrid(i)%ID,           i=1,np)
    write(16) (pgrid(i)%status,       i=1,np)
    write(16) (pgrid(i)%xlast(:),     i=1,np)
    write(16) (pgrid(i)%xcurr(:),     i=1,np)
    write(16) (pgrid(i)%nu_ext,       i=1,np)
    write(16) (pgrid(i)%k(:),         i=1,np)
    write(16) (pgrid(i)%nb_abs,       i=1,np)
    write(16) (pgrid(i)%time,         i=1,np)
    write(16) (pgrid(i)%tau_abs_curr, i=1,np)
    write(16) (pgrid(i)%iran,         i=1,np)
    write(16) (pgrid(i)%v_src(:),     i=1,np)
    close(16)

  end subroutine save_photons



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
