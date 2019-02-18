program main

  use module_parallel_mpi
  use module_photon

  implicit none

  integer(kind=4) :: i, j
  type(photon_current),dimension(10) :: p, temp_p

  call start_mpi
  call define_mpi_type

  ! initialize p_old
  do i=1,10
     p(i)%id = i
     p(i)%status = 0
     p(i)%xlast = (/0.,0.,0./)
     p(i)%xcurr = (/0.,0.,0./)
     p(i)%nu_ext = 1216.
     p(i)%k = (/0.,0.,1./)
     p(i)%nb_abs = 0
     p(i)%time = 0.
     p(i)%tau_abs_curr = 0.
     p(i)%iran = -99
  enddo
  
  ! Envoi des particules de 0 vers 1 
  if (rank == 0) then
     call MPI_SEND(p(1)%id, 10, MPI_TYPE_PHOTON, 1,  tag , MPI_COMM_WORLD, code)
  else
     call MPI_RECV(temp_p(1)%id, 10,  MPI_TYPE_PHOTON, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, IERROR) 

     do i=1,10
        if(p(i)%id /= temp_p(i)%id) print*,'MPI problem'
        if(p(i)%status /= temp_p(i)%status) print*,'MPI problem'
        do j=1,3
           if(p(i)%xlast(j) /= temp_p(i)%xlast(j)) print*,'MPI problem'
           if(p(i)%xcurr(j) /= temp_p(i)%xcurr(j)) print*,'MPI problem'
        enddo
        if(p(i)%nu_ext /= temp_p(i)%nu_ext) print*,'MPI problem'
        do j=1,3
           if(p(i)%k(j) /= temp_p(i)%k(j)) print*,'MPI problem'
        enddo
        if(p(i)%nb_abs /= temp_p(i)%nb_abs) print*,'MPI problem'
        if(p(i)%time /= temp_p(i)%time) print*,'MPI problem'
        if(p(i)%tau_abs_curr /= temp_p(i)%tau_abs_curr) print*,'MPI problem'
        if(p(i)%iran /= temp_p(i)%iran) print*,'MPI problem'
     enddo
     
  endif

  !!!call MPI_BARRIER(MPI_COMM_WORLD,code)

  ! Libération du type
  call MPI_TYPE_FREE(MPI_TYPE_PHOTON,code)
  call MPI_FINALIZE(code)
  
  
  
end program main
