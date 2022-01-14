module module_parallel_mpi

  use mpi
  use module_photon
  !--PEEL-- 
  use module_mock
  !--LEEP--

  implicit none

  ! MPI constants
  integer(kind=4)                                 :: rank          ! cpu number
  integer(kind=4)                                 :: nb_cpus       ! total number of cpus
  integer(kind=4)                                 :: code, error
  integer(kind=4)                                 :: tag, done_tag, ierror, exi_tag=2
  integer(kind=4)                                 :: nworker
  integer(kind=4),dimension(MPI_STATUS_SIZE)      :: status

  ! define a MPI derived type for photons
  integer(kind=4),parameter                       :: nbloc=11
  integer(kind=4),dimension(nbloc)                :: types, longueurs_blocs
  integer(kind=MPI_ADDRESS_KIND),dimension(nbloc) :: deplacements, adresses
  integer(kind=MPI_ADDRESS_KIND)                  :: lb,extent,nextelement,lbound,asize
  integer(kind=4)                                 :: mpi_type_photon,temp

  !--PEEL--
  public :: master_receives_mock, send_mock_to_master
  !--LEEP-- 

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
    type(photon_current),dimension(10) :: p

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
    !    real(kind=8),dimension(3) :: v_src        ! velocity of the source -- for peeling off
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
         MPI_DOUBLE_PRECISION /)

    ! block lengths
    longueurs_blocs = (/1,1,3,3,1,3,1,1,1,1,3/)
  
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
    call MPI_GET_ADDRESS(p(1)%v_src,        adresses(11),code)

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
       p_test(i)%v_src = (/20.e5*i,-10.e5*i,5.e5*i/)
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
          do j=1,3
             if(p_test(i)%v_src(j) /= temp_p(i)%v_src(j)) stop 'problem with MPI type'
          enddo
       enddo
    endif
  end subroutine test_mpi_type

  !--PEEL--
  subroutine send_mock_to_master(idcpu)
    ! send mocks to master

    integer(kind=4),intent(in) :: idcpu
    integer(kind=4) :: idir, n

    n = idcpu
    print*,'cpu ',idcpu,' is sending info ... '
    call MPI_SEND(n, 1, MPI_INTEGER, 0, tag , MPI_COMM_WORLD, code)

    do idir = 1,nDirections
       ! flux
       call MPI_SEND(mock(idir)%flux_aperture, 1, MPI_DOUBLE_PRECISION, 0, tag , MPI_COMM_WORLD, code)
       call MPI_SEND(mock(idir)%flux, 1, MPI_DOUBLE_PRECISION, 0, tag , MPI_COMM_WORLD, code)
       ! spectrum
       call MPI_SEND(mock(idir)%compute_spectrum, 1, MPI_LOGICAL, 0, tag , MPI_COMM_WORLD, code)
       if (mock(idir)%compute_spectrum) then 
          call MPI_SEND(mock(idir)%spec_npix, 1, MPI_INTEGER, 0, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(mock(idir)%spec_aperture, 1, MPI_DOUBLE_PRECISION, 0, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(mock(idir)%spec_lmin, 1, MPI_DOUBLE_PRECISION, 0, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(mock(idir)%spec_lmax, 1, MPI_DOUBLE_PRECISION, 0, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(mock(idir)%spectrum, mock(idir)%spec_npix, MPI_DOUBLE_PRECISION, 0, tag , MPI_COMM_WORLD, code)
       end if
       ! image
       call MPI_SEND(mock(idir)%compute_image, 1, MPI_LOGICAL, 0, tag , MPI_COMM_WORLD, code)
       if (mock(idir)%compute_image) then 
          call MPI_SEND(mock(idir)%image_npix, 1, MPI_INTEGER, 0, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(mock(idir)%image_side, 1, MPI_DOUBLE_PRECISION, 0, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(mock(idir)%center, 3, MPI_DOUBLE_PRECISION, 0, tag , MPI_COMM_WORLD, code)
          n = mock(idir)%image_npix*mock(idir)%image_npix
          call MPI_SEND(mock(idir)%image, n, MPI_DOUBLE_PRECISION, 0, tag , MPI_COMM_WORLD, code)
       end if
       ! cube
       call MPI_SEND(mock(idir)%compute_cube, 1, MPI_LOGICAL, 0, tag , MPI_COMM_WORLD, code)
       if (mock(idir)%compute_cube) then 
          call MPI_SEND(mock(idir)%cube_lbda_npix, 1, MPI_INTEGER, 0, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(mock(idir)%cube_image_npix, 1, MPI_INTEGER, 0, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(mock(idir)%cube_lmin, 1, MPI_DOUBLE_PRECISION, 0, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(mock(idir)%cube_lmax, 1, MPI_DOUBLE_PRECISION, 0, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(mock(idir)%cube_side, 1, MPI_DOUBLE_PRECISION, 0, tag , MPI_COMM_WORLD, code)
          call MPI_SEND(mock(idir)%center, 3, MPI_DOUBLE_PRECISION, 0, tag , MPI_COMM_WORLD, code)
          n = mock(idir)%cube_lbda_npix*mock(idir)%cube_image_npix*mock(idir)%cube_image_npix
          call MPI_SEND(mock(idir)%cube, n, MPI_DOUBLE_PRECISION, 0, tag , MPI_COMM_WORLD, code)
       end if
    end do
    
    return
    
  end subroutine send_mock_to_master

  subroutine master_receives_mock

    ! Receive mocks from a worker and add it to master's

    implicit none 
    integer(kind=4) :: idir, n, idcpu
    real(kind=8)             :: flux
    real(kind=8),allocatable :: spec(:), image(:,:), cube(:,:,:)
    
    
    ! First, receive information from a given CPU and identify the CPU
    print*,'master about to receive ... '
    call MPI_RECV(n, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, IERROR)
    print*,'master received idcpu = ',n
    idcpu = status(MPI_SOURCE)
    if (idcpu /= n) then
       print*,'Jeje has to learn MPI ... '
    end if
    
    do idir = 1,nDirections
       ! flux
       call MPI_RECV(mock(idir)%flux_aperture, 1, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
       call MPI_RECV(flux, 1, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
       mock(idir)%flux = mock(idir)%flux + flux
       ! spectrum
       call MPI_RECV(mock(idir)%compute_spectrum, 1, MPI_LOGICAL, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
       if (mock(idir)%compute_spectrum) then
          call MPI_RECV(mock(idir)%spec_npix, 1, MPI_INTEGER, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(mock(idir)%spec_aperture, 1, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(mock(idir)%spec_lmin, 1, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(mock(idir)%spec_lmax, 1, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          allocate(spec(mock(idir)%spec_npix))
          call MPI_RECV(spec, mock(idir)%spec_npix, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          mock(idir)%spectrum = mock(idir)%spectrum + spec
          deallocate(spec)
       end if
       ! image
       call MPI_RECV(mock(idir)%compute_image, 1, MPI_LOGICAL, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
       if (mock(idir)%compute_image) then 
          call MPI_RECV(mock(idir)%image_npix, 1, MPI_INTEGER, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(mock(idir)%image_side, 1, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(mock(idir)%center, 3, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          n = mock(idir)%image_npix*mock(idir)%image_npix
          allocate(image(mock(idir)%image_npix,mock(idir)%image_npix))
          call MPI_RECV(image, n, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          mock(idir)%image = mock(idir)%image + image
          deallocate(image)
       end if
       ! cube
       call MPI_RECV(mock(idir)%compute_cube, 1, MPI_LOGICAL, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
       if (mock(idir)%compute_cube) then 
          call MPI_RECV(mock(idir)%cube_lbda_npix, 1, MPI_INTEGER, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(mock(idir)%cube_image_npix, 1, MPI_INTEGER, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(mock(idir)%cube_lmin, 1, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(mock(idir)%cube_lmax, 1, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(mock(idir)%cube_side, 1, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          call MPI_RECV(mock(idir)%center, 3, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          n = mock(idir)%cube_lbda_npix*mock(idir)%cube_image_npix*mock(idir)%cube_image_npix
          allocate(cube(mock(idir)%cube_lbda_npix,mock(idir)%cube_image_npix,mock(idir)%cube_image_npix))
          call MPI_RECV(cube, n, MPI_DOUBLE_PRECISION, idcpu, DONE_TAG , MPI_COMM_WORLD, status,IERROR)
          mock(idir)%cube = mock(idir)%cube + cube
          deallocate(cube)
       end if
    end do
    
    return
    
  end subroutine master_receives_mock
  !--LEEP--
  
end module module_parallel_mpi

