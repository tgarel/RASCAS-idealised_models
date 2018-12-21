module module_master

  use module_parallel_mpi
  use module_domain
  use module_photon

  implicit none

  private
  
  integer(kind=4),dimension(:),allocatable         :: first,last,nqueue,ncpuperdom
  integer(kind=4),dimension(:,:),allocatable       :: next
  integer(kind=4),dimension(:), allocatable        :: cpu
  real(kind=8), dimension(:), allocatable          :: delta
  type(photon_current),dimension(:),allocatable    :: photgrid, photpacket
  integer(kind=4)                                  :: nphot
  type(domain)                                     :: compute_dom
  type(domain),dimension(:),allocatable            :: domain_list
  
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [master] of the parameter file
  ! --------------------------------------------------------------------------
  logical                   :: verbose = .false.
  ! checkpoint/restart
  logical                   :: restart = .false.    ! if true, start the run from backup file PhotonBakFile
  character(2000)           :: PhotonBakFile = 'backup_photons.dat'
  real(kind=8)              :: dt_backup = 7200.    ! time in seconds between 2 backups, default is 7200
  ! --------------------------------------------------------------------------

  public :: master, read_master_params, print_master_params
  
contains

  subroutine master(file_compute_dom, ndomain, domain_file_list, file_ICs, nbundle, fileout)

    implicit none
    
    character(2000),intent(in)                    :: file_compute_dom
    integer(kind=4),intent(in)                    :: ndomain
    character(2000),dimension(ndomain),intent(in) :: domain_file_list
    character(2000),intent(in)                    :: file_ICs
    integer(kind=4),intent(in)                    :: nbundle
    character(2000),intent(in)                    :: fileout
    integer(kind=4)                               :: i,j,icpu,idcpu,jnewdom,nphottodo,ncpuended,ntest,nphotdone
    logical                                       :: everything_not_done
    real(kind=8)                                  :: start_initphot, end_initphot, time_now, dt_since_last_backup, time_last_backup
    real(kind=8)                                  :: percentDone, percentBefore
    
    call cpu_time(start_initphot)

    ! read ICs photons or restore from backup
    if(restart)then
       if (verbose) print *,'[master] restoring photons from file: ',trim(PhotonBakFile)
       call restore_photons(PhotonBakFile,photgrid)
       nphot = size(photgrid)
       nphottodo = count(mask=(photgrid(:)%status==0))
       if (verbose)then
          print *,'[master] Nphoton =',nphot
          print *,'[master] Nphoton to do =',nphottodo
       endif
    else
       if (verbose) print *,'[master] reading ICs photons in file: ',trim(file_ICs)
       call init_photons_from_file(file_ICs,photgrid)
       nphot = size(photgrid)
       nphottodo = nphot
       if (verbose) print *,'[master] Nphoton =',nphot
    endif
    percentBefore = real(nphot-nphottodo)/nphot*100.
    
    ! some sanity checks
    if(nbundle*nworker>nphottodo)then
       print *,'ERROR: decrease nbundle and/or ncpu'
       call stop_mpi
    endif
    ! guidance for a good load-balancing
    if(4*nbundle*nworker>nphottodo)then
       print *,'ERROR: decrease nbundle for a good load-balancing of the code'
       print *,'--> suggested nbundle =', nphottodo/nworker/10
       call stop_mpi
    endif

    allocate(photpacket(nbundle))
    allocate(cpu(1:nworker))
    allocate(first(ndomain),last(ndomain),nqueue(ndomain),ncpuperdom(ndomain))
    allocate(next(nphot,ndomain))
    allocate(delta(ndomain))
    allocate(domain_list(ndomain))

    if (verbose) print *,'[master] reading domains'
    ! Get domain properties
    call domain_constructor_from_file(file_compute_dom,compute_dom)
    if (verbose) then
       print *,'[master] Ndomain =',ndomain
       print *,'   |_ ',trim(file_compute_dom)
    end if
    
    ! get list of domains and their properties
    do i=1,ndomain
       call domain_constructor_from_file(domain_file_list(i),domain_list(i))
       if (verbose) print *,'   |_ ',trim(domain_file_list(i))
    end do

    ! initialize the queue for each domain
    next(:,:)=-1
    first(:)=-1
    last(:)=-1
    nqueue(:)=0
    do i=1,nphot
       if(photgrid(i)%status==0)then  ! for restart
          j = get_my_new_domain(photgrid(i)%xcurr,domain_list)
          call add_photon_to_domain(i,j)
       endif
    enddo
    
    ! according to the distribution of photons in domains, ditribute cpus to each domain
    call init_loadb(nbundle,ndomain)

    call cpu_time(end_initphot)
    if (verbose) print '(" [master] time to initialize photons in master = ",f12.3," seconds.")',end_initphot-start_initphot
    if (verbose) print*,'[master] send a first bundle of photon packets to each worker'
    time_last_backup = end_initphot

    ! send a first bundle of photon packets to each worker
    do icpu=1,nworker

       ! identify the mesh domain of the targeted worker
       j=cpu(icpu)
       if(verbose) print '(" [master] allocates domain ",i5," to cpu ",i5)',j,icpu
       
       ! construct a list of photon packets to send
       call fill_bundle(j,photpacket,nbundle)

       ! Send the mesh domain number to the worker
       call MPI_SEND(j, 1, MPI_INTEGER, icpu, tag , MPI_COMM_WORLD, code)
 
       ! Send the bundle of photon packets to the worker
       !debug
       !call MPI_SEND(photpacket(1)%id, nbundle, MPI_TYPE_PHOTON, icpu, tag , MPI_COMM_WORLD, code)
       do i = 1,nbundle 
          call MPI_SEND(photpacket(i)%id, 1, MPI_INTEGER, icpu, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(photpacket(i)%status, 1, MPI_INTEGER, icpu, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(photpacket(i)%xlast, 3, MPI_DOUBLE_PRECISION, icpu, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(photpacket(i)%xcurr, 3, MPI_DOUBLE_PRECISION, icpu, tag , MPI_COMM_WORLD, code) 
          call MPI_SEND(photpacket(i)%nu_ext, 1, MPI_DOUBLE_PRECISION, icpu, tag , MPI_COMM_WORLD, code) 
          call MPI_SEND(photpacket(i)%k, 3, MPI_DOUBLE_PRECISION, icpu, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(photpacket(i)%nb_abs, 1, MPI_INTEGER, icpu, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(photpacket(i)%time, 1, MPI_DOUBLE_PRECISION, icpu, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(photpacket(i)%tau_abs_curr, 1, MPI_DOUBLE_PRECISION, icpu, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(photpacket(i)%iran, 1, MPI_INTEGER, icpu, tag , MPI_COMM_WORLD, code)
       end do
       !gubed
    end do

    everything_not_done=.true.

    ncpuended=0
    
    if(verbose) print*,'[master] starting waiting/receiving/sending cycle'

    ! Receive and append to pertinent domain list
    do while(everything_not_done)

       !if(verbose) then
       !   print*,'[master] NQueue     =',nqueue(:)
       !   print*,'[master] NCpuPerDom =',ncpuperdom(:)
       !   print*,'[master] waiting for worker...'
       !   print*,' '
       !endif

       ! First, receive information from a given CPU and identify the CPU
       call MPI_RECV(ntest, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, IERROR)

       idcpu = status(MPI_SOURCE)
       
       !debug
       !call MPI_RECV(photpacket(1)%id, nbundle, MPI_TYPE_PHOTON, idcpu, DONE_TAG, MPI_COMM_WORLD, status, IERROR)
       do i = 1,nbundle
          call MPI_RECV(photpacket(i)%id, 1, MPI_INTEGER, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(photpacket(i)%status, 1, MPI_INTEGER, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(photpacket(i)%xlast, 3, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(photpacket(i)%xcurr, 3, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(photpacket(i)%nu_ext, 1, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(photpacket(i)%k, 3, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(photpacket(i)%nb_abs, 1, MPI_INTEGER, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(photpacket(i)%time, 1, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(photpacket(i)%tau_abs_curr, 1, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(photpacket(i)%iran, 1, MPI_INTEGER, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
       end do
       !gubed                                                                                                                                                                                                                                                                                                              \

       !if (verbose) print*,'[master] receive a bundle of',nbundle,' photon packets from worker',idcpu

       ! By construction, photons that arrive here are not going to go back to the same domain, 
       !   i.e they are either dead (status=1or2) either in transit (status=0)
       do i=1,nbundle

          if((photpacket(i)%status==0).and.(ndomain==1))then
             print*,'ERROR: pb with photon packets bundle in master...'
             call stop_mpi
          endif
          ! empty bundle, nothing to do
          if(photpacket(i)%status==-1)exit

          if (photpacket(i)%status==0) then             ! photon in transit
             ! find my new domain
             jnewdom = get_my_new_domain(photpacket(i)%xcurr,domain_list)
             ! append the jnewdom queue
             call add_photon_to_domain(photpacket(i)%ID,jnewdom)
          endif

          ! whatever the status, update grid values for the photon
          call update_grid(i)

       enddo

       ! clear photpacket
       photpacket(:)%status=-1

       ! check if it is time to back up
       call cpu_time(time_now)
       dt_since_last_backup = time_now - time_last_backup
       if(dt_since_last_backup > dt_backup)then
          call backup_run
          time_last_backup = time_now
       endif

       ! check end of work
       nphottodo=sum(nqueue)
       if(nphottodo <= 0)then
          ! no more photon to send

          ! first count ended cpu, to not skip last working cpu...
          ncpuended=ncpuended+1
          if(ncpuended==nworker)then
             everything_not_done=.false.
          endif
          if(verbose) print '(" [master] no more photon to send to worker ",i5," then send exit code.")',idcpu
          !if (verbose) print*,'[master] no more photon to send to worker ',idcpu
          !if (verbose) print*,'[master] sending exit code to',idcpu
          call MPI_SEND(idcpu, 1, MPI_INTEGER, idcpu, exi_tag , MPI_COMM_WORLD, code)

       else
          ! keep sending photons

          ! check load-balancing and re-allocate CPU to a new domain if needed
          call update_domain(idcpu,nbundle,ndomain)

          j=cpu(idcpu)

          call MPI_SEND(j, 1, MPI_INTEGER, idcpu, tag, MPI_COMM_WORLD, code)

          ! Construct a new bundle of photon packets
          call fill_bundle(j,photpacket,nbundle)

          !if (verbose) print*,'[master] sending a new bundle of photpackets to worker/domain ',idcpu,j
          !if (verbose) print*,'[master] status RT ',nphottodo, nphot
          
          ! send it

          !debug
          !call MPI_SEND(photpacket(1)%id, nbundle, MPI_TYPE_PHOTON, idcpu, tag , MPI_COMM_WORLD, code)
          do i = 1,nbundle
             call MPI_SEND(photpacket(i)%id, 1, MPI_INTEGER, idcpu, tag , MPI_COMM_WORLD, code)
             call MPI_SEND(photpacket(i)%status, 1, MPI_INTEGER, idcpu, tag , MPI_COMM_WORLD, code)
             call MPI_SEND(photpacket(i)%xlast, 3, MPI_DOUBLE_PRECISION, idcpu, tag , MPI_COMM_WORLD, code)
             call MPI_SEND(photpacket(i)%xcurr, 3, MPI_DOUBLE_PRECISION, idcpu, tag , MPI_COMM_WORLD, code)
             call MPI_SEND(photpacket(i)%nu_ext, 1, MPI_DOUBLE_PRECISION, idcpu, tag , MPI_COMM_WORLD, code)
             call MPI_SEND(photpacket(i)%k, 3, MPI_DOUBLE_PRECISION, idcpu, tag , MPI_COMM_WORLD, code)
             call MPI_SEND(photpacket(i)%nb_abs, 1, MPI_INTEGER, idcpu, tag , MPI_COMM_WORLD, code)
             call MPI_SEND(photpacket(i)%time, 1, MPI_DOUBLE_PRECISION, idcpu, tag , MPI_COMM_WORLD, code)
             call MPI_SEND(photpacket(i)%tau_abs_curr, 1, MPI_DOUBLE_PRECISION, idcpu, tag , MPI_COMM_WORLD, code)
             call MPI_SEND(photpacket(i)%iran, 1, MPI_INTEGER, idcpu, tag , MPI_COMM_WORLD, code)
          end do
          !gubed                                                                                                                                                                                                                                                                                                                                                             

       endif

       ! print progress
       nphotdone = count(mask=(photgrid(:)%status/=0))
       percentdone = real(nphotdone)/nphot*100.
       !if(nphotdone/=nphotdonebefore)then
       !   !print*,'[master] progress report RT, number of photon packets done = ',nphotdone,real(nphotdone)/nphot*100.
       !   print '(" [master] --> number of photon packets done = ",i8," (",f5.1," %)")',nphotdone,real(nphotdone)/nphot*100.
       !   nphotdonebefore=nphotdone
       !endif
       if(percentdone>=percentBefore+1.)then
          print '(" [master] number of photon packets done = ",i8," (",f4.1," %)")',nphotdone,percentDone
          percentBefore = percentDone + 1.
       endif

    enddo

    ! finale synchronization, for profiling purposes
    call MPI_BARRIER(MPI_COMM_WORLD,code)    

    if(verbose)then
       call print_diagnostics
       !print *,' '
       print *,'[master] writing results in file: ',trim(fileout)
    endif
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


  subroutine backup_run

    character(1000) :: filebak, copyfile
    logical         :: file_exists

    ! first, copy last backup file into filebak
    filebak = trim(PhotonBakFile)//'.bak'
    ! check if file exists
    INQUIRE(FILE=PhotonBakFile, EXIST=file_exists)
    !if(access(PhotonBakFile,' ') == 0) then
    if(file_exists)then
       copyfile = 'mv '//trim(PhotonBakFile)//' '//trim(filebak) 
       !call execute_command_file(copyfile)
       call system(copyfile)
    endif

    ! then, write a new backup file
    ! in this 1st attempt, we save only photon_current grid
    ! => restart from the grid, need to reconstruct all the queues, but can restart with any number of CPU
    call save_photons(PhotonBakFile,photgrid)

    if (verbose) print *,'[master] backup done'

  end subroutine backup_run



  subroutine print_diagnostics

    implicit none
    integer(kind=4) :: i

    write(*,*)'[master] some quick checks on the results...'

    ! test status of photons
    do i=1,nphot
       if(photgrid(i)%status == 0)print*,'ERROR: problem with photon status...',i,photgrid(i)%status
    enddo
    ! Some stats on photon status
    !print *,' '
    print *,'   photon status:'
    print *,'   # of photons             =',size(photgrid(:)%status)
    print *,'   # of status=1 (escaped)  =',count(mask=(photgrid(:)%status==1))
    print *,'   # of status=2 (absorbed) =',count(mask=(photgrid(:)%status==2))
    !print *,' '
    ! write results                                                                                               
    ! some check
    if(verbose)then
       print *,'   some diagnostics:'
       print *,'   min max status      =',minval(photgrid%status),maxval(photgrid%status)
       print *,'   min max pos x       =',minval(photgrid%xcurr(1)),maxval(photgrid%xcurr(1))
       print *,'   min max pos y       =',minval(photgrid%xcurr(2)),maxval(photgrid%xcurr(2))
       print *,'   min max pos z       =',minval(photgrid%xcurr(3)),maxval(photgrid%xcurr(3))
       print *,'   min max nb scatt    =',minval(photgrid%nb_abs),maxval(photgrid%nb_abs)
       print *,'   min max nu          =',minval(photgrid%nu_ext),maxval(photgrid%nu_ext)
       print *,'   min max lambda      =',clight/maxval(photgrid%nu_ext)*cmtoA,clight/minval(photgrid%nu_ext)*cmtoA
       print *,'   min max travel time =',minval(photgrid%time),maxval(photgrid%time)
       print *,'   Last scattering:'
       print *,'   min max pos x       =',minval(photgrid%xlast(1)),maxval(photgrid%xlast(1))
       print *,'   min max pos y       =',minval(photgrid%xlast(2)),maxval(photgrid%xlast(2))
       print *,'   min max pos z       =',minval(photgrid%xlast(3)),maxval(photgrid%xlast(3))
    endif

  end subroutine print_diagnostics




  subroutine init_loadb(nbundle,ndomain)
    ! give cpus to domains, according to the number of photons to deal with in each domain (=nqueue)

    implicit none
    integer(kind=4), intent(in)  :: nbundle
    integer(kind=4), intent(in)  :: ndomain
    integer(kind=4)              :: nphottot,icpu,ndom,j
    integer(kind=4),dimension(1) :: jtoo
    
    nphottot=sum(nqueue)
    ! due to limitation in integer precision (default is kind=4), nworker*nqueue could easily give an overflow...
    ncpuperdom(:)=int(real(nworker)*real(nqueue(:))/nphottot)

    ! assign cpus left over by rounding error to the domain which has the max nb of photons.
    ! (This is best guess and will be balanced dynamically later).
    if (sum(ncpuperdom) < nworker) then
       jtoo=maxloc(ncpuperdom)
       ncpuperdom(jtoo) = ncpuperdom(jtoo) + nworker - sum(ncpuperdom)
    else if (sum(ncpuperdom) > nworker) then 
       write(*,*) '[master] ERROR in init_loadb ... aborting'
       call stop_mpi
    end if

    !if(verbose) write(*,*)'[master] init load-balancing',nworker,sum(ncpuperdom),ncpuperdom(:)

    icpu=1
    do j=1,ndomain
       ndom = ncpuperdom(j)
       cpu(icpu:icpu+ndom-1)=j
       icpu=icpu+ndom

       ! we also have to initialize delta(j)
       delta(j) = real(nqueue(j))/nbundle/(ncpuperdom(j)+1.)

    end do

    !if(verbose)then
    !   write(*,*)'[master] init cpu mapping',cpu(:)
    !   write(*,*)'[master] init delta(j)',delta(:)
    !   write(*,*)'[master] init nqueue(j)',nqueue(:)
    !endif
    if(verbose)then
       write(*,*)'[master] initial mapping workers <--> domains'
       write(*,*)'   ------------------------------'
       write(*,*)'   domain   Nqueue   Ncpu   delta'
       write(*,*)'   ------------------------------'
       do j=1,ndomain
          write(*,'(i8,i10,i6,f10.3)') j, nqueue(j), ncpuperdom(j), delta(j)
       enddo
       write(*,*)'   ------------------------------'
    endif
    
  end subroutine init_loadb




  subroutine update_domain(icpu,nbundle,ndomain)

    implicit none
    integer(kind=4), intent(in)   :: icpu
    integer(kind=4), intent(in)   :: nbundle
    integer(kind=4), intent(in)   :: ndomain
    integer(kind=4)               :: j,jold,nphot
    integer(kind=4), dimension(1) :: jtarget, jnew

    ! Note: the method used is that we keep workers allocated to their domains until there is no more photon in the domain queue.

    ! need to update delta for all domains
    do j=1,ndomain
       delta(j) = real(nqueue(j))/nbundle/(ncpuperdom(j)+1.)
    enddo

    j=cpu(icpu)

    nphot=sum(nqueue)

    if(delta(j)<0.5)then

       jtarget = maxloc(delta)

       if((delta(jtarget(1))>1.).or.(delta(j)<=0.))then

          !if(verbose) write(*,*)'[master] updating domain suite',j,jtarget,delta(:)

          jold=j
          jnew=jtarget
          cpu(icpu)=jnew(1)
          ncpuperdom(jnew)=ncpuperdom(jnew)+1
          ncpuperdom(jold)=ncpuperdom(jold)-1

          !if(verbose)then
          !   write(*,*)'[master] updated load-balancing',ncpuperdom(:)
          !   write(*,*)'[master] updated cpu mapping',cpu(:)
          !endif
          if(verbose)then
             write(*,*)'[master] updated mapping workers <--> domains'
             write(*,*)'   ------------------------------'
             write(*,*)'   domain   Nqueue   Ncpu   delta'
             write(*,*)'   ------------------------------'
             do j=1,ndomain
                write(*,'(i8,i10,i6,f10.3)') j, nqueue(j), ncpuperdom(j), delta(j)
             enddo
             write(*,*)'   ------------------------------'
           endif
       endif

    endif

  end subroutine update_domain




  subroutine add_photon_to_domain(i,j)
    ! add one photon to the domain queue

    implicit none
    integer(kind=4), intent(in) :: i,j
    
    if(first(j)==-1)then
       first(j)=i
       last(j)=i
    else
       next(last(j),j)=i
       last(j)=i
    endif
    nqueue(j)=nqueue(j)+1

  end subroutine add_photon_to_domain




  subroutine fill_bundle(j,photpacket,nbundle)
    ! fill bundle of photon packets photpacket(nbundle)

    implicit none
    integer(kind=4), intent(in)                           :: j
    integer(kind=4), intent(in)                           :: nbundle
    type(photon_current), dimension(nbundle), intent(out) :: photpacket
    integer(kind=4)                                       :: i,fsave

    i=1
    if(first(j)==-1)then
       print*,'ERROR: cannot fill bundle...'
       call stop_mpi
    endif
    photpacket(:)%ID=0
    do
       photpacket(i) = photgrid(first(j))
       fsave = first(j)
       first(j)=next(first(j),j)
       next(fsave,j)=-1
       i=i+1
       nqueue(j)=nqueue(j)-1
       if (i>nbundle.or.first(j)==-1)exit
    end do

  end subroutine fill_bundle



  subroutine update_grid(i)
    ! update photon grid with photpacket
    ! NB: it would be clearer to pass photpakcet in argument

    implicit none
    integer(kind=4)             :: myID
    integer(kind=4), intent(in) :: i

    myID = photpacket(i)%ID
    if (photgrid(myID)%ID /= myID) then
       print*,'ERROR: id mismatch updating the grid'
       call stop_mpi
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
    character(1000)         :: line,name,value
    integer(kind=4)         :: err,i
    logical                 :: section_present
    
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
          case ('restart')
             read(value,*) restart
          case ('PhotonBakFile')
             write(PhotonBakFile,'(a)') trim(value)
          case ('dt_backup')
             read(value,*) dt_backup
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
       write(unit,'(a,L1)')          '  restart        = ',restart
       write(unit,'(a,a)')           '  PhotonBakFile  = ',trim(PhotonBakFile)
       write(unit,'(a,f12.3)')       '  dt_backup      = ',dt_backup
       write(unit,'(a)')             ' '
       call print_mesh_params(unit)
    else
       write(*,'(a)')             '[master]'
       write(*,'(a,L1)')          '  verbose        = ',verbose
       write(*,'(a,L1)')          '  restart        = ',restart
       write(*,'(a,a)')           '  PhotonBakFile  = ',trim(PhotonBakFile)
       write(*,'(a,f12.3)')       '  dt_backup      = ',dt_backup
       write(*,'(a)')             ' '       
       call print_mesh_params
    end if

    return

  end subroutine print_master_params

  

end module module_master
