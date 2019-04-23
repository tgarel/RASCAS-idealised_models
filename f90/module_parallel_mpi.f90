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
  integer(kind=4),parameter                       :: nbloc=10
  integer(kind=4),dimension(nbloc)                :: types, longueurs_blocs
  integer(kind=MPI_ADDRESS_KIND),dimension(nbloc) :: deplacements, adresses
  integer(kind=4)                                 :: mpi_type_photon


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
    type(photon_current) :: p

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
         MPI_INTEGER /)

    ! block lengths
    longueurs_blocs = (/1,1,3,3,1,3,1,1,1,1/)
  
    call MPI_GET_ADDRESS(p%id,           adresses(1),code)
    call MPI_GET_ADDRESS(p%status,       adresses(2),code)
    call MPI_GET_ADDRESS(p%xlast,        adresses(3),code)
    call MPI_GET_ADDRESS(p%xcurr,        adresses(4),code)
    call MPI_GET_ADDRESS(p%nu_ext,       adresses(5),code)
    call MPI_GET_ADDRESS(p%k,            adresses(6),code)
    call MPI_GET_ADDRESS(p%nb_abs,       adresses(7),code)
    call MPI_GET_ADDRESS(p%time,         adresses(8),code)
    call MPI_GET_ADDRESS(p%tau_abs_curr, adresses(9),code)
    call MPI_GET_ADDRESS(p%iran,         adresses(10),code)

    ! Compute array of displacements
    do i=1,nbloc
       deplacements(i)=adresses(i) - adresses(1)
    end do
    ! Create mpi_type_photon
    call MPI_TYPE_CREATE_STRUCT (nbloc,longueurs_blocs,deplacements,types,mpi_type_photon,&
         code)
    ! Commits the new type
    call MPI_TYPE_COMMIT(mpi_type_photon,code)

  end subroutine define_mpi_type



end module module_parallel_mpi

