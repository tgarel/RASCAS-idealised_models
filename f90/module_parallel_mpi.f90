module module_parallel_mpi

  use mpi
  use module_photon
  
  implicit none

  ! MPI constants
  integer(kind=4)                                 :: rank          ! cpu number
  integer(kind=4)                                 :: nb_cpus       ! total number of cpus
  integer(kind=4)                                 :: code, error
  integer(kind=4)                                 :: tag, done_tag, ierror, exi_tag=2
  integer(kind=4)                                 :: nworker
  integer(kind=4),dimension(MPI_STATUS_SIZE)      :: status

  ! define a MPI derived type for photons
  integer(kind=4),parameter                       :: nbloc=11 ! nbloc=10 ! TIBO
  integer(kind=4),dimension(nbloc)                :: types, longueurs_blocs
  integer(kind=MPI_ADDRESS_KIND),dimension(nbloc) :: deplacements, adresses
  integer(kind=MPI_ADDRESS_KIND)                  :: lb,extent,nextelement,lbound,asize
  integer(kind=4)                                 :: mpi_type_photon,temp


contains

  !//=========================================
  !||
  !|| subroutine start_mpi
  !|| subroutine finish_mpi
  !|| subroutine stop_mpi
  !|| subroutine define_mpi_type
  !||
  !\\==========================================



  subroutine start_mpi

    ! Initialize the MPI execution environment
    call mpi_init(code)

    ! Determines the rank of the cpu
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

    ! Determines the number of cpus
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_cpus, code)
    
  end subroutine start_mpi



  subroutine finish_mpi

    ! Terminates MPI execution environment
    call mpi_finalize(code)
    
  end subroutine finish_mpi



  subroutine stop_mpi

    call MPI_ABORT(MPI_COMM_WORLD,error,code)

  end subroutine stop_mpi



  subroutine define_mpi_type
    ! Create an MPI datatype for sending photons

    integer(kind=4)      :: i
    ! TIBO
    !type(photon_current),dimension(10) :: p
    type(photon_current),dimension(11) :: p
    ! OBIT
    
    ! type photon_current
    !    integer(kind=4)           :: ID
    !    integer(kind=4)           :: status       ! =0 if flying, =1 if escape, =2 if absorption (by dust)
    !    real(kind=8),dimension(3) :: xlast        ! coordinates of last interaction in box units
    !    real(kind=8),dimension(3) :: xcurr        ! current position of the photon in box units
    !    real(kind=8)              :: nu_ext       ! external frame frequency (Hz)
    !    real(kind=8),dimension(3) :: k            ! normalised propagation vector 
    !    integer(kind=4)           :: nb_abs       ! number of interactions before escape
    !    real(kind=8)              :: time         ! time in [s] from emission to escape/absorption        
    !    real(kind=8)              :: tau_abs_curr ! current optical depth (useful when photon change mesh domain)
    !    integer(kind=4)           :: iran         ! state of the random generator
    ! end type photon_current


    ! Create array of types
    types = (/ MPI_INTEGER, &
         MPI_INTEGER, &
         MPI_DOUBLE_PRECISION, &
         MPI_DOUBLE_PRECISION, &
         MPI_DOUBLE_PRECISION, &
         MPI_DOUBLE_PRECISION, &
         MPI_INTEGER, &
         MPI_DOUBLE_PRECISION, &
         MPI_DOUBLE_PRECISION, &
         MPI_INTEGER, &
         MPI_INTEGER /)

    ! block lengths
    ! TIBO
    longueurs_blocs = (/1,1,3,3,1,3,1,1,1,1,1/)
    ! OBIT
    
    call MPI_GET_ADDRESS(p(1)%id,           adresses(1),code)
    call MPI_GET_ADDRESS(p(1)%status,       adresses(2),code)
    call MPI_GET_ADDRESS(p(1)%xlast,        adresses(3),code)
    call MPI_GET_ADDRESS(p(1)%xcurr,        adresses(4),code)
    call MPI_GET_ADDRESS(p(1)%nu_ext,       adresses(5),code)
    call MPI_GET_ADDRESS(p(1)%k,            adresses(6),code)
    call MPI_GET_ADDRESS(p(1)%nb_abs,       adresses(7),code)
    call MPI_GET_ADDRESS(p(1)%time,         adresses(8),code)
    call MPI_GET_ADDRESS(p(1)%tau_abs_curr, adresses(9),code)
    call MPI_GET_ADDRESS(p(1)%iran,         adresses(10),code)
    ! TIBO
    call MPI_GET_ADDRESS(p(1)%n_backscatt,  adresses(11),code)
    ! OBIT
    
    ! Compute array of displacements
    do i=1,nbloc
       deplacements(i)=adresses(i) - adresses(1)
       !if (i>1)then
       !   print*,'displacement',i,deplacements(i),adresses(i)-adresses(i-1)
       !else
       !   print*,'displacement',i,deplacements(i)
       !end if
    end do
    ! Create mpi_type_photon
    call MPI_TYPE_CREATE_STRUCT(nbloc,longueurs_blocs,deplacements,types,temp,code)
    ! Extent correct
    !call MPI_TYPE_GET_EXTENT(temp, lbound, asize, code)
    !print*,'MPI type get extent -> ',lbound, asize, code

    call MPI_GET_ADDRESS(p(2)%id,nextelement,code)
    lb = 0
    extent = nextelement-adresses(1)
    !print*,'extent =',extent, nextelement-adresses(10)
    call MPI_TYPE_CREATE_RESIZED(temp,lb,extent,mpi_type_photon,code)

    !call MPI_TYPE_GET_EXTENT(mpi_type_photon, lbound, asize, code)
    !print*,'MPI type get extent -> ',lbound, asize, code

    ! Commits the new type
    call MPI_TYPE_COMMIT(mpi_type_photon,code)

  end subroutine define_mpi_type



  subroutine test_mpi_type
    ! With some hardware and/or MPI libs we had problem with MPI_TYPE_CREATE_STRUCT.
    ! This routine tests the correct behaviour of MPI with the type mpi_type_photon used in rascas.
    
    use module_photon
    
    implicit none
    
    integer(kind=4) :: i, j
    type(photon_current),dimension(10) :: p_test, temp_p
    
    ! initialize p_test
    do i=1,10
       p_test(i)%id = i
       p_test(i)%status = 0
       p_test(i)%xlast = (/1.*i,0.5*i,2.*i/)
       p_test(i)%xcurr = (/0.1*i,0.2*i,0.3*i/)
       p_test(i)%nu_ext = 1216.
       p_test(i)%k = (/0.,0.,1./)
       p_test(i)%nb_abs = 100*i
       p_test(i)%time = 1.23456789 * i
       p_test(i)%tau_abs_curr = 0.
       p_test(i)%iran = -99
    enddo
    
    ! Send p_test from 0 to 1
    if (rank == 0) then
       call MPI_SEND(p_test(1)%id, 10, MPI_TYPE_PHOTON, 1, tag , MPI_COMM_WORLD, code)
    endif
    if (rank == 1) then
       call MPI_RECV(temp_p(1)%id, 10, MPI_TYPE_PHOTON, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, IERROR)
       do i=1,10
          if(p_test(i)%id /= temp_p(i)%id) stop 'problem with MPI type'
          if(p_test(i)%status /= temp_p(i)%status) stop 'problem with MPI type'
          do j=1,3
             if(p_test(i)%xlast(j) /= temp_p(i)%xlast(j)) stop 'problem with MPI type'
             if(p_test(i)%xcurr(j) /= temp_p(i)%xcurr(j)) stop 'problem with MPI type'
          enddo
          if(p_test(i)%nu_ext /= temp_p(i)%nu_ext) stop 'problem with MPI type'
          do j=1,3
             if(p_test(i)%k(j) /= temp_p(i)%k(j)) stop 'problem with MPI type'
          enddo
          if(p_test(i)%nb_abs /= temp_p(i)%nb_abs) stop 'problem with MPI type'
          if(p_test(i)%time /= temp_p(i)%time) stop 'problem with MPI type'
          if(p_test(i)%tau_abs_curr /= temp_p(i)%tau_abs_curr) stop 'problem with MPI type'
          if(p_test(i)%iran /= temp_p(i)%iran) stop 'problem with MPI type'
       enddo
    endif
  end subroutine test_mpi_type
  
  
end module module_parallel_mpi

