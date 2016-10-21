!//=========================================
!||module_parallel_mpi.f90
!||
!|| subroutine initialize_mpi
!|| subroutine finalisation_mpi
!|| subroutine type_derive
!||
!\\==========================================
module module_parallel_mpi

  use mpi
  use module_photon
  
  implicit none

  !Numero du sous domaine
  INTEGER                                   :: rank
  !Nombre de processus
  INTEGER                                   :: nb_cpus
  !Constantes MPI
  INTEGER                                   :: code  

  integer :: tag, done_tag, ierror, exi_tag=2
  integer :: nslave

  integer, dimension(MPI_STATUS_SIZE) :: status

  integer, parameter :: nbloc=10
  integer, dimension(nbloc) :: types, longueurs_blocs
  integer(kind= MPI_ADDRESS_KIND), dimension(nbloc) :: deplacements, adresses
  integer :: mpi_type_photon


contains


  subroutine initialisation_mpi
    !************
    !Initialisation pour chaque processus de son rang et du 
    !nombre total de processus nb_procs
    !************

    !Initialisation de MPI
    call mpi_init(code)

    !savoir quel processus je suis
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

    !Connaitre le nombre total de processus
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_cpus, code)
    
  end subroutine initialisation_mpi


  subroutine finalisation_mpi
    !************
    !Desactivation de l'environnement MPI
    !************

    ! Desactivation de MPI
    call mpi_finalize(code)
    
  end subroutine finalisation_mpi


  subroutine clean_stop

    call MPI_FINALIZE(code)
    stop

  end subroutine clean_stop


  subroutine type_derive

    integer :: i
    type(photon_current) :: p

    ! type photon_current
    !    integer                   :: ID
    !    integer                   :: status       ! =0 if flying, =1 if escape, =2 if absorption (by dust)
    !    real(kind=8),dimension(3) :: xlast        ! coordinates of last interaction in box units
    !    real(kind=8),dimension(3) :: xcurr        ! current position of the photon in box units
    !    real(kind=8)              :: nu_ext       ! external frame frequency (Hz)
    !    real(kind=8),dimension(3) :: k            ! normalised propagation vector 
    !    integer                   :: nb_abs       ! number of interactions before escape
    !    real(kind=8)              :: time         ! time in [s] from emission to escape/absorption        
    !    real(kind=8)              :: tau_abs_curr ! current optical depth (useful when photon change mesh domain)
    !    integer                   :: iran         ! state of the random generator


    ! create photon type for MPI send or list of photons ?


    ! Construction du type de donnees
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

    ! Calcul des deplacements relatifs a l'adresse de depart
    do i=1,nbloc
       deplacements(i)=adresses(i) - adresses(1)
    end do
    call MPI_TYPE_CREATE_STRUCT (nbloc,longueurs_blocs,deplacements,types,mpi_type_photon,&
         code)
    ! Validation du type structure
    call MPI_TYPE_COMMIT(mpi_type_photon,code)


  end subroutine type_derive



end module module_parallel_mpi

