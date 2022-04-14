module module_gas_composition

  use module_dust_model
  use module_random
  use module_constants
  use module_lyc_model

  implicit none

  private

  character(100),parameter :: moduleName = 'module_gas_composition_HI_HeI_HeII_dust.f90'

  type, public :: gas
     real(kind=8)                :: nHI       ! HI numerical density [HI/cm3]
     real(kind=8)                :: nHeI      ! HeI numerical density [HeI/cm3]
     real(kind=8)                :: nHeII     ! HeII numerical density [HeII/cm3]
     real(kind=8)                :: ndust     ! pseudo-numerical density of dust particles [#/cm3]
     real(kind=8)                :: v(3)      ! gas velocity [cm/s]
  end type gas
  real(kind=8),public :: box_size_cm   ! size of simulation box in cm. 

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [gas_composition] of the parameter file
  ! --------------------------------------------------------------------------
  ! mixture parameters
  real(kind=8)             :: f_ion           = 0.01   ! ndust = (n_HI + f_ion*n_HII) * Z/Zref [Laursen+09]
  real(kind=8)             :: Zref            = 0.005  ! reference metallicity. Should be ~ 0.005 for SMC and ~ 0.01 for LMC.
  logical                  :: use_dust        = .true.
  logical                  :: use_helium      = .true.
  logical                  :: use_v           = .true.
  ! --------------------------------------------------------------------------

  ! public functions:
  public :: gas_from_ramses_leaves,gas_get_scatter_flag,gas_scatter,dump_gas,get_gas_velocity
  public :: read_gas,gas_destructor,read_gas_composition_params,print_gas_composition_params
  !--PEEL--
  public :: gas_peeloff_weight,gas_get_tau
  !--LEEP--

  !--CORESKIP-- push variable from module_lines_model up so that module_photon knows about it... 
  public :: HI_core_skip 
  !--PIKSEROC-- 

contains


  subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var,g)

    ! define gas contents from ramses raw data

    use module_ramses

    character(2000),intent(in)                     :: repository 
    integer(kind=4),intent(in)                     :: snapnum
    integer(kind=4),intent(in)                     :: nleaf,nvar
    real(kind=8),intent(in),dimension(nvar,nleaf)  :: ramses_var
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                                :: ileaf
    real(kind=8),dimension(:),allocatable          :: T, nhi, metallicity, nh, nhei, nheii, nO
    real(kind=8),dimension(:,:),allocatable        :: v
    
    ! allocate gas-element array
    allocate(g(nleaf))

    box_size_cm = ramses_get_box_size_cm(repository,snapnum)

    ! compute velocities in cm / s
    write(*,*) '-- module_gas_composition_HI_HeI_HeII_dust : extracting velocities from ramses '
    allocate(v(3,nleaf))
    call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
    do ileaf = 1,nleaf
       g(ileaf)%v = v(:,ileaf)
    end do
    deallocate(v)

    ! get HI, HeI and HeII from ramses
    write(*,*) '-- module_gas_composition_HI_HeI_HeII_dust : extracting nHI, nHeI, and nHeII from ramses '
    allocate(nhi(nleaf),nh(nleaf),nhei(nleaf),nheii(nleaf))
    call ramses_get_nh_nhi_nhei_nehii_cgs(repository,snapnum,nleaf,nvar,ramses_var,nh,nhi,nhei,nheii)
    g(:)%nHI   = nhi(:)
    g(:)%nHeI  = nhei(:)
    g(:)%nHeII = nheii(:)
    deallocate(nhei,nheii)

    ! get ndust (pseudo dust density)
    write(*,*) '-- module_gas_composition_HI_HeI_HeII_dust, hilbert mode : extracting ndust from ramses '
    allocate(metallicity(nleaf),T(nleaf),nO(nleaf))
    call ramses_get_metallicity(nleaf,nvar,ramses_var,metallicity)
    call ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nhi)
    nO = 4.9d-4*nh*metallicity/0.0139
    do ileaf = 1,nleaf
       g(ileaf)%ndust = get_n_dust(nh(ileaf),nhi(ileaf),nh(ileaf)-nhi(ileaf),metallicity(ileaf),Zref,f_ion,nO(ileaf),T(ileaf)) ! [ /cm3 ]
    end do
    deallocate(metallicity,nh,nhi,nO,T)

    return

  end subroutine gas_from_ramses_leaves


  function get_gas_velocity(cell_gas)
    type(gas),intent(in)      :: cell_gas
    real(kind=8),dimension(3) :: get_gas_velocity
    get_gas_velocity(:) = cell_gas%v(:)
    return
  end function get_gas_velocity


  function  gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs, iran, CS_dist_cm, CS_xcrit)

    ! --------------------------------------------------------------------------
    ! Decide whether a scattering event occurs, and if so, on which element
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - cell_gas : nHI, nHeI, nHeII and dust
    ! - distance_to_border_cm : the maximum distance the photon may travel (before leaving the cell)
    ! - nu_cell : photon frequency in cell's frame [ Hz ]
    ! - tau_abs : optical depth at which the next scattering event will occur
    ! - iran    : random generator state of the photon
    ! - CS_dist_cm: not useful here, for HI core-skipping only.
    ! OUTPUTS:
    ! - distance_to_border_cm : comes out as the distance to scattering event (if there is an event)
    ! - tau_abs : 0 if a scatter occurs, decremented by tau_cell if photon escapes cell. 
    ! - gas_get_scatter_flag : 0 [no scatter], 1 [MgII-2796 scatter], 2 [MgII-2804 scatter], 3 [dust] 
    ! - CS_xcrit: not useful here, for HI core-skipping only.
    ! --------------------------------------------------------------------------

    type(gas),intent(in)                  :: cell_gas
    real(kind=8),intent(inout)            :: distance_to_border_cm
    real(kind=8),intent(in)               :: nu_cell
    real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur.
    integer,intent(inout)                 :: iran 
    integer(kind=4)                       :: gas_get_scatter_flag 
    real(kind=8)                          :: tau_gas, tau_cell, tau_dust, x, proba_gas
    real(kind=8),intent(in)               :: CS_dist_cm
    real(kind=8),intent(inout)            :: CS_xcrit

    ! compute optical depths for different components of the gas.
    tau_gas       = get_tau_gas(cell_gas%nHI, cell_gas%nHeI, cell_gas%nHeII, distance_to_border_cm, nu_cell)
    tau_dust      = get_tau_dust(cell_gas%ndust, distance_to_border_cm, nu_cell)
    tau_cell      = tau_gas + tau_dust

    if (tau_abs > tau_cell) then  ! photon is due for absorption outside the cell 
       gas_get_scatter_flag = 0
       tau_abs = tau_abs - tau_cell
       if (tau_abs.lt.0.0d0) then
          print*, 'tau_abs est negatif'
          stop
       endif
    else  ! the scattering happens inside the cell
       ! decide if it is gas or dust
       proba_gas = tau_gas /  tau_cell
       x = ran3(iran)
       if (x <= proba_gas) then
          gas_get_scatter_flag = 1 ! absorption by gas (HI, HeI or HeII, doesn't matter)
       else
          gas_get_scatter_flag = 2 ! absorption by dust
       end if
       ! and transform "distance_to_border_cm" in "distance_to_absorption_cm"
       distance_to_border_cm = distance_to_border_cm * (tau_abs / tau_cell)
    end if

    return

  end function gas_get_scatter_flag


  subroutine gas_scatter(flag,cell_gas,nu_cell,k,nu_ext,iran,xcrit)

    integer(kind=4),intent(inout)            :: flag
    type(gas),intent(in)                     :: cell_gas
    real(kind=8),intent(inout)               :: nu_cell, nu_ext
    real(kind=8),dimension(3), intent(inout) :: k
    integer(kind=4),intent(inout)            :: iran
    !--CORESKIP--
    real(kind=8),intent(in)                  :: xcrit
    !--PIKSEROC--
    integer(kind=4)                          :: ilost

    select case(flag)
    case(1)            ! Absorption by gas, photon always destroyed
       flag = -2       ! New definition of flag. If -2, photon was killed by gas. if -1, was killed by dust
    case(2)
       call scatter_dust(cell_gas%v, nu_cell, k, nu_ext, iran, ilost)
       if(ilost==1) flag=-1
    end select

  end subroutine gas_scatter


  !--PEEL--
  function gas_get_tau(cell_gas, distance_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute total opacity of gas accross distance_cm at freq. nu_cell
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - cell_gas : a mix of stuff
    ! - distance_cm : the distance along which to compute tau [cm]
    ! - nu_cell : photon frequency in cell's frame [ Hz ]
    ! OUTPUTS:
    ! - gas_get_tau : the total optical depth
    ! --------------------------------------------------------------------------

    ! check whether scattering occurs within cell (scatter_flag > 0) or not (scatter_flag==0)
    type(gas),intent(in)    :: cell_gas
    real(kind=8),intent(in) :: distance_cm
    real(kind=8),intent(in) :: nu_cell
    real(kind=8)            :: gas_get_tau
    real(kind=8)            :: tau_gas, tau_dust

    ! compute optical depths for different components of the gas.
    tau_gas = get_tau_gas(cell_gas%nHI ,cell_gas%nHeI, cell_gas%nHeII, distance_cm, nu_cell)
    tau_dust = get_tau_dust(cell_gas%ndust, distance_cm, nu_cell)
    gas_get_tau = tau_gas + tau_dust

    return

  end function gas_get_tau
  ! --------------------------------------------------------------------------
  !--LEEP--


  !--PEEL--
  function gas_peeloff_weight(flag,cell_gas,nu_ext,kin,kout,iran)

    integer(kind=4),intent(in)            :: flag
    type(gas),intent(in)                  :: cell_gas
    real(kind=8),intent(inout)            :: nu_ext
    real(kind=8),dimension(3), intent(in) :: kin, kout
    integer(kind=4),intent(inout)         :: iran
    real(kind=8)                          :: gas_peeloff_weight

    select case(flag)
    case(1)                   ! Absorption by gas, no contribution to the peeling-off
       gas_peeloff_weight = 0d0
    case(2)
       gas_peeloff_weight = dust_peeloff_weight(cell_gas%v, nu_ext, kin, kout)
    end select

  end function gas_peeloff_weight
  !--LEEP--



  subroutine dump_gas(unit,g)
    type(gas),dimension(:),intent(in) :: g
    integer(kind=4),intent(in)        :: unit
    integer(kind=4)                   :: i,j,nleaf
    real(kind=8)                      :: mean_T, n_tot
    nleaf = size(g)
    write(unit) (g(i)%v(:), i=1,nleaf)
    write(unit) (g(i)%nHI, i=1,nleaf)
    write(unit) (g(i)%nHeI, i=1,nleaf)
    write(unit) (g(i)%nHeII, i=1,nleaf)
    write(unit) (g(i)%ndust, i=1,nleaf)
    write(unit) box_size_cm
  end subroutine dump_gas
  

  subroutine read_gas(unit,n,g)
    integer(kind=4),intent(in)                     :: unit,n
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                                :: i
    allocate(g(1:n))

    if(use_v) then
       read(unit) (g(i)%v(:),i=1,n)
    else
       do i=1,3
          g(:)%v(i) = 0d0
       end do
       read(unit)
    end if

    read(unit) (g(i)%nHI,i=1,n)
    if(use_helium) then
       read(unit) (g(i)%nHeI,i=1,n)
       read(unit) (g(i)%nHeII,i=1,n)
    else
       g(:)%nHeI = 0d0
       g(:)%nHeII = 0d0
       read(unit)
       read(unit)
    end if

    if(use_dust) then
       read(unit) (g(i)%ndust,i=1,n)
    else
       g(:)%ndust = 0d0
       read(unit)
    end if
    
    read(unit) box_size_cm

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
    ! ALSO read parameter form used modules (HI, D, dust models)
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
          case ('f_ion')
             read(value,*) f_ion
          case ('Zref')
             read(value,*) Zref
          case ('use_v')
             read(value,*) use_v
          case ('use_dust')
             read(value,*) use_dust
          case ('use_helium')
             read(value,*) use_helium
          end select
       end do
    end if
    close(10)

    call read_line_params(pfile)
    call read_dust_params(pfile)

    return

  end subroutine read_gas_composition_params



  subroutine print_gas_composition_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)') '[gas_composition]'
       write(unit,'(a,a)')      '# code compiled with: ',trim(moduleName)
       write(unit,'(a)')        '# mixture parameters'
       write(unit,'(a,ES10.3)') '  f_ion           = ',f_ion
       write(unit,'(a,ES10.3)') '  Zref            = ',Zref
       write(unit,'(a,L1)')     '  use_v           = ',use_v
       write(unit,'(a,L1)')     '  use_dust        = ',use_dust
       write(unit,'(a,L1)')     '  use_helium      = ',use_helium
       write(unit,'(a)')             ' '
       call print_line_params(unit)
       write(unit,'(a)')             ' '
       call print_dust_params(unit)
    else
       write(*,'(a,a,a)') '[gas_composition]'
       write(*,'(a,a)')      '# code compiled with: ',trim(moduleName)
       write(*,'(a)')        '# mixture parameters'
       write(*,'(a,ES10.3)') '  f_ion           = ',f_ion
       write(*,'(a,ES10.3)') '  Zref            = ',Zref
       write(*,'(a,L1)')     '  use_v           = ',use_v
       write(*,'(a,L1)')     '  use_dust        = ',use_dust
       write(*,'(a,L1)')     '  use_helium      = ',use_helium
       write(*,'(a)')             ' '
       call print_line_params
       write(*,'(a)')             ' '
       call print_dust_params
    end if

    return

  end subroutine print_gas_composition_params


end module module_gas_composition
