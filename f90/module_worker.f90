module module_worker

  use module_parallel_mpi
  use module_domain
  use module_photon
  use module_mesh
  !--PEEL--
  use module_mock
  !--LEEP--

  private 

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [worker] of the parameter file
  ! --------------------------------------------------------------------------
  logical                   :: verbose = .false.
  ! --------------------------------------------------------------------------

  public :: worker, read_worker_params, print_worker_params
  
contains

  subroutine worker(file_compute_dom, ndomain, mesh_file_list, nbuffer)

    implicit none
    
    character(2000),intent(in)                    :: file_compute_dom
    integer(kind=4),intent(in)                    :: ndomain
    character(2000),dimension(ndomain),intent(in) :: mesh_file_list
    integer(kind=4),intent(in)                    :: nbuffer
    
    integer(kind=4)                               :: jdomain, mydom
    type(photon_current),dimension(nbuffer)       :: photpacket
    type(mesh)                                    :: meshdom
    type(domain)                                  :: compute_dom
    real(kind=8)                                  :: start_photpacket,end_photpacket

    ! get the computational domain
    call domain_constructor_from_file(file_compute_dom,compute_dom)

    mydom=-1

    do while (status(MPI_TAG) .ne. EXI_TAG)

       ! receive my domain number
       call MPI_RECV(jdomain, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, IERROR)

       if(status(MPI_TAG) == EXI_TAG) exit

       ! if necessary read domain data
       if(jdomain/=mydom)then

          if(verbose) write(*,'(a,i4.4,a,i4)') ' [w',rank,'] : loading new domain',jdomain

          call mesh_from_file(mesh_file_list(jdomain),meshdom) ! -> overwrites gas props if parameters are set to. 

          mydom=jdomain

          if(verbose) write(*,'(a,i4.4,a)') ' [w',rank,'] : done loading new domain'
          
       endif

       ! receive my list of photons to propagate
       call MPI_RECV(photpacket(1)%id, nbuffer, MPI_TYPE_PHOTON, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, IERROR)

       if(verbose) write(*,'(a,i4.4,a)') ' [w',rank,'] : just receive a photpacket and start processing it...' 

       call cpu_time(start_photpacket)

       ! do the RT stuff. This is a single loop over photons in photpacket.
       call MCRT(nbuffer,photpacket,meshdom,compute_dom)
       
       call cpu_time(end_photpacket)

       if(verbose)then
          write(*,'(a,i4.4,a,f12.3,a)') ' [w',rank,'] : time to propagate photpacket = ',end_photpacket-start_photpacket,' seconds.'
          write(*,'(a,i4.4,a)') ' [w',rank,'] : finish packet of photons, sending back to master ' 
       endif

       ! send my results
       call MPI_SEND(nbuffer, 1, MPI_INTEGER, 0, tag , MPI_COMM_WORLD, code)

       call MPI_SEND(photpacket(1)%id, nbuffer, MPI_TYPE_PHOTON, 0, tag , MPI_COMM_WORLD, code)

    enddo

    !--PEEL--
    if (peeling_off) then 
       call dump_mocks(rank)
    end if
    !--LEEP--
    
    if(verbose) write(*,'(a,i4.4,a)') ' [w',rank,'] : exit of loop...'

    ! synchronization, for profiling purposes
    call MPI_BARRIER(MPI_COMM_WORLD,code)


  end subroutine worker


  
  subroutine read_worker_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! NB: does not call read_params of depdencies (master module does that). 
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
       if (line(1:8) == '[worker]') then
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

    return

  end subroutine read_worker_params


  
  subroutine print_worker_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a)')             '[worker]'
       write(unit,'(a,L1)')          '  verbose        = ',verbose
       write(unit,'(a)')             ' '
    else
       write(*,'(a)')             '[worker]'
       write(*,'(a,L1)')          '  verbose        = ',verbose
       write(*,'(a)')             ' '       
    end if

    return

  end subroutine print_worker_params
  
end module module_worker
