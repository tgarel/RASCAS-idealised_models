module module_worker

  use module_parallel_mpi
  use module_domain
  use module_photon
  use module_mesh
  use module_params
  !!!use module_uparallel

contains

  subroutine worker

    implicit none
    
    integer                                 :: jdomain, mydom
    logical                                 :: mpi_ok, newdomain
    type(photon_current),dimension(nbuffer) :: photpacket
    type(mesh)                              :: meshdom
    type(domain)                            :: compute_dom

    ! Question: est-ce qu'il faut d'abord initialiser quelque chose, un premier paquet de photos?

    ! get the compute domain
    ! how to get filename?
    call domain_constructor_from_file(file_compute_dom,compute_dom)

    ! load uparallel table
    !!print*,'--> loading uparallel tables...'
    !!call init_uparallel_tables

    mydom=-1
    MPI_ok=.true.

!!    do while (MPI_ok)
    do while (status(MPI_TAG) .ne. EXI_TAG)

       ! receive my domain
       ! Q: how to retreive the domain? through MPI or an integer jdomain is a key? 
       call MPI_RECV(jdomain, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, IERROR)

       !print*,DONE_TAG, MPI_TAG
       !print*,'[worker] status(MPI_TAG) =',status(MPI_TAG),EXI_TAG
       if(status(MPI_TAG) == EXI_TAG) exit

       ! if necessary read domain data
       if(jdomain/=mydom)then

          call mesh_from_file(mesh_file_list(jdomain),meshdom) ! -> overwrites gas props if parameters are set to. 
          mydom=jdomain
!!$          if(overwritegas)then
!!$             call overwrite_mesh(meshdom,nhi_new,vth_new)
!!$          endif

          if (rank==1 .and. verbose)then
             print*,'[worker 1] read mesh domain in file: ',trim(mesh_file_list(jdomain))
             if(overwritegas)then
                print*,'  & overwrite gas with nHI =',nhi_new
                print*,'  & overwrite gas with vth =',vth_new
             endif
          endif

       endif

       ! receive my list of photons to propagate
       call MPI_RECV(photpacket(1)%id, nbuffer, MPI_TYPE_PHOTON, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, IERROR)

       !if(verbose)then
       !   print*,'[worker] my rank',rank
       !   print*,'[worker] just receive a photpacket and start processing it...'
       !endif

       ! do the RT stuff
       call MCRT(nbuffer,photpacket,meshdom,compute_dom)
       ! this is a single loop over photons
       ! but with all the RT stuff <= plugin with McLya
       !...

       !print*,'[worker] finish packet of photons, sending back to master...',nbuffer,code

       ! send my results
       call MPI_SEND(nbuffer, 1, MPI_INTEGER, 0, tag , MPI_COMM_WORLD, code)

       !print*,'[worker] finish packet of photons, sending back to master...',photpacket(1)%id,code

       call MPI_SEND(photpacket(1)%id, nbuffer, MPI_TYPE_PHOTON, 0, tag , MPI_COMM_WORLD, code)

       !print *,'[worker] ',status(MPI_TAG),EXI_TAG

    enddo

    if(verbose) print *,'[worker] number',rank,' exit of loop...'

  end subroutine worker

end module module_worker
