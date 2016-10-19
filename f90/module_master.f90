module module_master

  use module_parallel_mpi
  use module_domain
  use module_photon

  implicit none

  private
  
  integer,dimension(:),allocatable              :: first,last,nqueue,ncpuperdom
  integer,dimension(:,:),allocatable            :: next
  integer,dimension(:), allocatable             :: cpu
  real(kind=8), dimension(:), allocatable       :: delta
  type(photon_current),dimension(:),allocatable :: photgrid, photpacket
  integer                                       :: nphot
  type(domain)                                  :: compute_dom
  type(mesh)                                    :: meshdom
  type(domain),dimension(:),allocatable         :: domain_list
  
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [worker] of the parameter file
  ! --------------------------------------------------------------------------
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------

  public :: master, read_master_params, print_master_params
  
contains

  subroutine master(file_compute_dom, ndomain, domain_file_list, file_ICs, nbuffer, fileout)

    implicit none
    
    character(2000),intent(in)                    :: file_compute_dom
    integer(kind=4),intent(in)                    :: ndomain
    character(2000),dimension(ndomain),intent(in) :: domain_file_list
    character(2000),intent(in)                    :: file_ICs
    integer(kind=4),intent(in)                    :: nbuffer
    character(2000),intent(in)                    :: fileout
    
    integer :: i,j,icpu,idcpu,jnewdom,nphottodo,count,ncpuended,ntest
    logical :: survivor
    logical :: init,everything_not_done
    real(kind=8) :: start_initphot, end_initphot

    call cpu_time(start_initphot)

    ! read ICs photons
    if (verbose) print *,'[master] --> reading ICs photons in file: ',trim(file_ICs)
    call init_photons_from_file(file_ICs,photgrid)
    nphot = size(photgrid)
    if (verbose) print *,'[master] --> Nphoton =',nphot

    allocate(photpacket(nbuffer))
    allocate(cpu(1:nslave))
    allocate(first(ndomain),last(ndomain),nqueue(ndomain),ncpuperdom(ndomain))
    allocate(next(nphot,ndomain))
    allocate(delta(ndomain))
    allocate(domain_list(ndomain))


    if (verbose) print *,'[master] --> reading domain and mesh...'
    ! Get domain properties
    call domain_constructor_from_file(file_compute_dom,compute_dom)
    if (verbose) print *,'[master] --> Ndomain =',ndomain
    if (verbose) print *,'[master] --> |_ ',trim(file_compute_dom)

    ! get list of domains and their properties
    do i=1,ndomain
       call domain_constructor_from_file(domain_file_list(i),domain_list(i))
       if (verbose) print *,'[master] --> |_ ',trim(domain_file_list(i))
    end do

    ! initialize the queue for each domain
    next(:,:)=-1
    first(:)=-1
    last(:)=-1
    nqueue(:)=0
    do i=1,nphot
       j = get_my_new_domain(photgrid(i)%xcurr,domain_list)
       call add_photon_to_domain(i,j)
    enddo

    ! according to the distribution of photons in domains, ditribute cpus to each domain
    call init_loadb(nbuffer,ndomain)

    call cpu_time(end_initphot)
    if (verbose) print '("[master] --> time to initialize photons in master = ",f12.3," seconds.")',end_initphot-start_initphot
    if (verbose) print*,'[master] send a first chunk of photons to each worker'

    ! send a first chunk of photons to each worker
    do icpu=1,nslave

       j=cpu(icpu)

       if(verbose) print*,'icpu / idom =',icpu,j
       
       ! build an array/list/buffer of photons to send
       call fill_buffer(j,photpacket,nbuffer)

       ! Send the photon buffer to the workers
       ! TODO : definir le MPI_TYPE correspondant
       ! and identify the target cpu

       !if(verbose)print*,'send domain number',j
       call MPI_SEND(j, 1, MPI_INTEGER, icpu, tag , MPI_COMM_WORLD, code)

       !if(verbose)print*,'send photpacket to cpu',icpu
       call MPI_SEND(photpacket(1)%id, nbuffer, MPI_TYPE_PHOTON, icpu, tag , MPI_COMM_WORLD, code)
       !!call MPI_SEND (nbuff_local, 1, MPI_INTEGER, proc, WORK_TAG , MPI_COMM_WORLD, IERROR)
       !!call MPI_SEND (photarray, nbuff_local, MPI_PHOTON, icpu, WORK_TAG, MPI_COMM_WORLD, IERROR)
       
    end do

    everything_not_done=.true.

    count=0
    ncpuended=0
    
    if(verbose) print*,'[master] starting loop...'

    ! Receive and append to pertinent domain list
    do while(everything_not_done)

       ! First, receive information from a given CPU and identify the CPU
       if(verbose) then
          print*,'[master] waiting for worker...'
          print*,'[master] Nqueue =',nqueue(:)
          print*,'[master] NCpuPerDom = ',ncpuperdom(:)
          print*,' '
       endif

       call MPI_RECV(ntest, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, IERROR)

       idcpu = status(MPI_SOURCE)

       if (verbose) print *,'[master] received from worker...',nbuffer
       !print *,MPI_ANY_SOURCE,status
       !print *,idcpu

       call MPI_RECV(photpacket(1)%id, nbuffer, MPI_TYPE_PHOTON, idcpu, DONE_TAG, MPI_COMM_WORLD, status, IERROR)

       if (verbose) print*,'[master] receive a packet of',nbuffer,' photons from worker',idcpu
       if (verbose) print*,'[master] ', photpacket(1)%id
       count = count+1
       !if(verbose)print*,'count =',count

       ! By construction, photons that arrive here are not going to go back to the same domain, 
       !   i.e they are either dead (status=1or2) either in transit (status=0)
       do i=1,nbuffer

          if((photpacket(i)%status==0).and.(ndomain==1))then
             print*,'pb with photon buffer in master...'
             stop
          endif
          if(photpacket(i)%status==-1)exit

          !!!!write(*,*)'photon update',i,photpacket(i)%ID,photpacket(i)%status   !!!,photpacket(i)%pos(1)

          survivor=.false.

          ! define code or change status into logical...
          if (photpacket(i)%status == 0) then
             survivor=.true.
          endif

          if(survivor)then
             ! find my new domain, photon by photon
             jnewdom = get_my_new_domain(photpacket(i)%xcurr,domain_list)
             !print*,'jnewdom=',jnewdom

             ! append the jnewdom queue
             call add_photon_to_domain(photpacket(i)%ID,jnewdom)
             ! j'ai l'impression que la meme routine peut servir ici, check with JB

          else
             ! it means that photon escaped (status=1) or has been absorbed by dust (status=2)
             !call write_result(photpacket(i)%ID)
             !...
          endif

          ! whatever the status, update grid values for the photon
          call update_grid(i)
          !...

       enddo

       ! empty photpacket...
       photpacket(:)%status=-1

       ! check end of work
       nphottodo=sum(nqueue)
       if(nphottodo <= 0)then
          !!call MPI_SEND(idcpu, 1, MPI_INTEGER, idcpu, exi_tag , MPI_COMM_WORLD, code)
          ! first count ended cpu, to not skip last working cpu...
          ncpuended=ncpuended+1
          if(ncpuended==nslave)then
             everything_not_done=.false.
          endif
       else
          ! keep sending photons

          ! top instead of bottom (last -> first)

          ! here, it should be checked wether the load-balancing is good...
          ! if not re-affect CPU to another domain..
          ! how to do that ?
          call update_domain(idcpu,nbuffer,ndomain)
          !...

          j=cpu(idcpu)
          call MPI_SEND(j, 1, MPI_INTEGER, idcpu, tag, MPI_COMM_WORLD, code)

          ! Construct a new list of photons
          call fill_buffer(j,photpacket,nbuffer)

          ! send it
          call MPI_SEND(photpacket(1)%id, nbuffer, MPI_TYPE_PHOTON, idcpu, tag , MPI_COMM_WORLD, code)
          !...

       endif

    enddo

    ! need to cleanly finish...???
    do icpu=1,nslave
       if(verbose) print*,'[master] sending exit code to',icpu
       call MPI_SEND(icpu, 1, MPI_INTEGER, icpu, exi_tag , MPI_COMM_WORLD, code)
       !call MPI_SEND(photpacket(1)%id, nbuffer, MPI_TYPE_PHOTON, icpu, exi_tag , MPI_COMM_WORLD, code)
    end do

    ! TODO: 
    ! - need to send end to workers => OK
    ! - check results...

    call check_results


    print *,' '
    print*,'[master] --> writing results in file: ',trim(fileout)
    call dump_photons(fileout,photgrid)


    deallocate(photgrid)
    deallocate(photpacket)
    deallocate(cpu)
    deallocate(first,last,nqueue,ncpuperdom)
    deallocate(next)
    deallocate(delta)

    if(verbose) print *,'[master] end'

  end subroutine master
  !===================================================================================================




  subroutine check_results

    implicit none
    integer :: i

    write(*,*)'--> check results routine...'

    ! test status of photons
    do i=1,nphot
       if(photgrid(i)%status == 0)print*,'ohoho problem with photon status...',i,photgrid(i)%status
    enddo
    ! Some stats on photon status
    print *,' '
    print *,'--> photon status...'
    print *,'# of photons             =',size(photgrid(:)%status)
    print *,'# of status=1 (escaped)  =',count(mask=(photgrid(:)%status==1))
    print *,'# of status=2 (absorbed) =',count(mask=(photgrid(:)%status==2))
    print *,'# of status=3 (crap, pb with precision/in_cell_finder) =',count(mask=(photgrid(:)%status==3))
    print *,' '
    ! write results                                                                                               
    ! some check
    if(verbose)then
       print *,'--> Some diagnostics...'
       print *,'min max status      =',minval(photgrid%status),maxval(photgrid%status)
       print *,'min max pos x       =',minval(photgrid%xcurr(1)),maxval(photgrid%xcurr(1))
       print *,'min max pos y       =',minval(photgrid%xcurr(2)),maxval(photgrid%xcurr(2))
       print *,'min max pos z       =',minval(photgrid%xcurr(3)),maxval(photgrid%xcurr(3))
       print *,'min max nb scatt    =',minval(photgrid%nb_abs),maxval(photgrid%nb_abs)
       print *,'min max nu          =',minval(photgrid%nu_ext),maxval(photgrid%nu_ext)
       print *,'min max lambda      =',clight/minval(photgrid%nu_ext)*cmtoA,clight/maxval(photgrid%nu_ext)*cmtoA
       print *,'min max travel time =',minval(photgrid%time),maxval(photgrid%time)
       print *,'Last scattering'
       print *,'min max pos x       =',minval(photgrid%xlast(1)),maxval(photgrid%xlast(1))
       print *,'min max pos y       =',minval(photgrid%xlast(2)),maxval(photgrid%xlast(2))
       print *,'min max pos z       =',minval(photgrid%xlast(3)),maxval(photgrid%xlast(3))
       print *,' '
    endif

  end subroutine check_results




  subroutine init_loadb(nbuffer,ndomain)
    ! give cpus to domains, according to the number of photons to deal with in each domain (=nqueue)

    implicit none
    integer, intent(in) :: nbuffer
    integer, intent(in) :: ndomain
    integer::nphottot,icpu,ndom,j
    integer,dimension(1)::jtoo
    
    nphottot=sum(nqueue)
    ncpuperdom(:)=nint(real(nslave*nqueue(:))/nphottot)
    ! check wether rounding is ok
    do while(sum(ncpuperdom)/=nslave)
       if(sum(ncpuperdom)>nslave)then
          ! quit one cpu
          jtoo=maxloc(ncpuperdom)
          ncpuperdom(jtoo)=ncpuperdom(jtoo)-1
       else
          if(sum(ncpuperdom)<nslave)then
             ! add one cpu (to the max or to the min?)
             jtoo=minloc(ncpuperdom)
             ncpuperdom(jtoo)=ncpuperdom(jtoo)+1
          endif
       endif
    end do

    if(verbose) write(*,*)'[master] init load-balancing',nslave,sum(ncpuperdom),ncpuperdom(:)

    icpu=1
    do j=1,ndomain
       ndom = ncpuperdom(j)
       cpu(icpu:icpu+ndom-1)=j
       icpu=icpu+ndom
       !print*,j,ndom,icpu

       ! we also have to initialize delta(j)
       delta(j) = real(nqueue(j))/nbuffer/(ncpuperdom(j)+1.)

    end do

    if(verbose)then
       write(*,*)'[master] init cpu mapping',cpu(:)
       write(*,*)'[master] init delta(j)',delta(:)
       write(*,*)'[master] init nqueue(j)',nqueue(:)
    endif

  end subroutine init_loadb




  subroutine update_domain(icpu,nbuffer,ndomain)

    implicit none
    integer, intent(in)::icpu
    integer, intent(in) :: nbuffer
    integer, intent(in) :: ndomain
    integer :: j,jold,nphot
    integer, dimension(1) :: jtarget, jnew

    ! methode 1: on garde les cpus associes aux domaines jusqu'a ce qu'il n'y ait plus assez de photons 
    ! notion de Nrun == nombre de photons pour un domaine / taille du buffer 

    ! need to update delta for all domains
    do j=1,ndomain
       delta(j) = real(nqueue(j))/nbuffer/(ncpuperdom(j)+1.)
    enddo

    j=cpu(icpu)

    nphot=sum(nqueue)

    !write(*,*)'updating domain',j,delta(j),real(nqueue(j)),ncpuperdom(j)

    if(delta(j)<0.5)then

       jtarget = maxloc(delta)

       if((delta(jtarget(1))>1.).or.(delta(j)<=0.))then

          if(verbose) write(*,*)'[master] updating domain suite',j,jtarget,delta(:)

          jold=j
          jnew=jtarget
          cpu(icpu)=jnew(1)
          ncpuperdom(jnew)=ncpuperdom(jnew)+1
          ncpuperdom(jold)=ncpuperdom(jold)-1

          if(verbose)then
             write(*,*)'[master] updated load-balancing',ncpuperdom(:)
             write(*,*)'[master] updated cpu mapping',cpu(:)
          endif
       endif

    endif


  end subroutine update_domain




  subroutine add_photon_to_domain(i,j)
    ! add one photon to the domain queue

    implicit none
    integer, intent(in) :: i,j
    
    if(first(j)==-1)then
       first(j)=i
       last(j)=i
    else
       next(last(j),j)=i
       last(j)=i
    endif
    nqueue(j)=nqueue(j)+1

  end subroutine add_photon_to_domain




  subroutine fill_buffer(j,photpacket,nbuffer)
    ! fill photpacket(nbuffer)

    implicit none
    integer, intent(in) :: j
    integer, intent(in) :: nbuffer
    type(photon_current), dimension(nbuffer), intent(out) :: photpacket
    integer :: i,fsave

    i=1
    if(first(j)==-1)then
       print*,'no no no, cannot fill buffer...'
       stop
    endif
    photpacket(:)%ID=0
    do
       photpacket(i) = photgrid(first(j))
       fsave = first(j)
       first(j)=next(first(j),j)
       next(fsave,j)=-1
       i=i+1
       nqueue(j)=nqueue(j)-1
       if (i>nbuffer.or.first(j)==-1)exit
    end do

  end subroutine fill_buffer



  subroutine update_grid(i)

    integer :: myID
    integer, intent(in) :: i

    myID = photpacket(i)%ID
    if (photgrid(myID)%ID /= myID) then
       print*,'oh oh id mismatch'
       stop
       return
    endif
    photgrid(myID) = photpacket(i)

  end subroutine update_grid


  subroutine read_master_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! ALSO call read_params of depdencies (mesh)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present
    
    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:8) == '[master]') then
          section_present = .true.
          exit
       end if
    end do
    ! read section if present
    if (section_present) then 
       do
          read (10,'(a)',iostat=err) line
          if(err/=0) exit
          if (line(1:1) == '[') exit ! next section starting... -> leave
          i = scan(line,'=')
          if (i==0 .or. line(1:1)=='#' .or. line(1:1)=='!') cycle  ! skip blank or commented lines
          name=trim(adjustl(line(:i-1)))
          value=trim(adjustl(line(i+1:)))
          i = scan(value,'!')
          if (i /= 0) value = trim(adjustl(value(:i-1)))
          select case (trim(name))
          case ('verbose')
             read(value,*) verbose
          end select
       end do
    end if
    close(10)

    call read_mesh_params(pfile)

    return

  end subroutine read_master_params


  
  subroutine print_master_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a)')             '[master]'
       write(unit,'(a,L1)')          '  verbose        = ',verbose
       write(unit,'(a)')             ' '
       call print_mesh_params(unit)
    else
       write(*,'(a)')             '[master]'
       write(*,'(a,L1)')          '  verbose        = ',verbose
       write(*,'(a)')             ' '       
       call print_mesh_params
    end if

    return

  end subroutine print_master_params

  

end module module_master
