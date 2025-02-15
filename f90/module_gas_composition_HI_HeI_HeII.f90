module module_gas_composition

  ! mix of HI, HeI, HeII
  ! NB: no lines here, this is for escape fraction computations ... -> this module will not work for regular rascas runs
  !
  ! WARNING : no velocity here ... no frequencies ... 

  use module_random
  use module_ramses
  use module_constants

  implicit none

  private

  character(100),parameter :: moduleName = 'module_gas_composition_HI_HeI_HeII.f90'
  
  type, public :: gas
     real(kind=8) :: nH    ! H  numerical density [H/cm3]
     real(kind=8) :: nHI   ! HI numerical density [HI/cm3]
     real(kind=8) :: nHeI  ! HeI numerical density [HeI/cm3]
     real(kind=8) :: nHeII ! HeII numerical density [HeII/cm3]
  end type gas
  
  real(kind=8),public :: box_size_cm   ! size of simulation box in cm. 

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [gas_composition] of the parameter file
  ! --------------------------------------------------------------------------
  ! grey cross sections
  real(kind=8)             :: sigma_HI
  real(kind=8)             :: sigma_HeI
  real(kind=8)             :: sigma_HeII
  ! miscelaneous
  logical                  :: verbose             = .false. ! display some run-time info on this module
  ! --------------------------------------------------------------------------

  ! public functions:
  public :: gas_from_ramses_leaves,dump_gas,gas_get_tau
  public :: read_gas,gas_destructor,read_gas_composition_params,print_gas_composition_params

contains
  

  subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var, g)

    ! define gas contents from ramses raw data

    character(2000),intent(in)                     :: repository 
    integer(kind=4),intent(in)                     :: snapnum
    integer(kind=4),intent(in)                     :: nleaf,nvar
    real(kind=8),intent(in),dimension(nvar,nleaf)  :: ramses_var
    type(gas),dimension(:),allocatable,intent(out) :: g
    real(kind=8),dimension(:),allocatable          :: nh, nhi, nhei, nheii

    ! allocate gas-element array
    allocate(g(nleaf))

    box_size_cm = ramses_get_box_size_cm(repository,snapnum)
    ! get nH, nHI, nHeI, nHeII
    if (verbose) write(*,*) '-- module_gas_composition_HI_HeI_HeII : extracting nHI, nHeI, and nHeII from ramses'
    allocate(nh(nleaf),nhi(nleaf),nhei(nleaf),nheii(nleaf))
    call ramses_get_nh_nhi_nhei_nehii_cgs(repository,snapnum,nleaf,nvar,ramses_var,nh,nhi,nhei,nheii)
    g(:)%nH    = nh(:)
    g(:)%nHI   = nhi(:)
    g(:)%nHeI  = nhei(:)
    g(:)%nHeII = nheii(:)
    deallocate(nh,nhi,nhei,nheii)

    return

  end subroutine gas_from_ramses_leaves

  
  function  gas_get_tau(cell_gas, distance_cm)

    ! --------------------------------------------------------------------------
    ! compute total opacity of gas accross distance_cm 
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - cell_gas : HI, HeI, HeII
    ! - distance_cm : the distance along which to compute tau [cm]
    ! OUTPUTS:
    ! - gas_get_tau : the total optical depth
    ! --------------------------------------------------------------------------

    type(gas),intent(in)    :: cell_gas
    real(kind=8),intent(in) :: distance_cm
    real(kind=8)            :: gas_get_tau
    real(kind=8)            :: tau_HI, tau_HeI, tau_HeII

    ! compute optical depths for different components of the gas.
    tau_HI = cell_gas%nHI * distance_cm * sigma_HI
    tau_HeI = cell_gas%nHeI * distance_cm * sigma_HeI
    tau_HeII = cell_gas%nHeII * distance_cm * sigma_HeII    
    gas_get_tau = tau_HI + tau_HeI + tau_HeII

    return
    
  end function gas_get_tau


  subroutine dump_gas(unit,g)
    type(gas),dimension(:),intent(in) :: g
    integer(kind=4),intent(in)        :: unit
    integer(kind=4)                   :: i,nleaf
    nleaf = size(g)
    write(unit) (g(i)%nH,  i=1,nleaf)
    write(unit) (g(i)%nHI, i=1,nleaf)
    write(unit) (g(i)%nHeI, i=1,nleaf)
    write(unit) (g(i)%nHeII, i=1,nleaf)
    write(unit) box_size_cm 
  end subroutine dump_gas
  


  subroutine read_gas(unit,n,g)
    integer(kind=4),intent(in)                     :: unit,n
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                                :: i
    allocate(g(1:n))
    read(unit) (g(i)%nH,i=1,n)
    read(unit) (g(i)%nHI,i=1,n)
    read(unit) (g(i)%nHeI,i=1,n)
    read(unit) (g(i)%nHeII,i=1,n)
    read(unit) box_size_cm 

    !print*,'minmax densities : '
    !print*,minval(g(:)%nH), maxval(g(:)%nH)
    !print*,minval(g(:)%nHI), maxval(g(:)%nHI)
    !print*,minval(g(:)%nHeI), maxval(g(:)%nHeI)
    !print*,minval(g(:)%nHeII), maxval(g(:)%nHeII)
    

  end subroutine read_gas

  

  subroutine gas_destructor(g)
    type(gas),dimension(:),allocatable,intent(inout) :: g
    deallocate(g)
  end subroutine gas_destructor


  subroutine read_gas_composition_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! ALSO read parameter form used modules (HI_model)
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
       if (line(1:17) == '[gas_composition]') then
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
          case ('sigma_HI')
             read(value,*) sigma_HI
          case ('sigma_HeI')
             read(value,*) sigma_HeI
          case ('sigma_HeII')
             read(value,*) sigma_HeII
          case ('verbose')
             read(value,*) verbose
          end select
       end do
    end if
    close(10)

    call read_ramses_params(pfile)

    return

  end subroutine read_gas_composition_params



  subroutine print_gas_composition_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)')    '[gas_composition]'
       write(unit,'(a,a)')     '# code compiled with: ',trim(moduleName)
       write(unit,'(a)')        '# cross sections '
       write(unit,'(a,ES10.3)') '  sigma_HI            = ',sigma_HI
       write(unit,'(a,ES10.3)') '  sigma_HeI           = ',sigma_HeI
       write(unit,'(a,ES10.3)') '  sigma_HeII          = ',sigma_HeII
       write(unit,'(a)')        '# miscelaneous parameters'
       write(unit,'(a,L1)')     '  verbose             = ',verbose
       write(unit,'(a)')             ' '
       call print_ramses_params(unit)
       write(unit,'(a)')             ' '
    else
       write(*,'(a,a,a)') '[gas_composition]'
       write(*,'(a,a,a)')    '[gas_composition]'
       write(*,'(a)')        '# cross sections '
       write(*,'(a,ES10.3)') '  sigma_HI            = ',sigma_HI
       write(*,'(a,ES10.3)') '  sigma_HeI           = ',sigma_HeI
       write(*,'(a,ES10.3)') '  sigma_HeII          = ',sigma_HeII
       write(*,'(a)')       '# miscelaneous parameters'
       write(*,'(a,L1)')    '  verbose             = ',verbose
       write(*,'(a)')             ' '
       call print_ramses_params
       write(*,'(a)')             ' '
    end if

    return

  end subroutine print_gas_composition_params


end module module_gas_composition
