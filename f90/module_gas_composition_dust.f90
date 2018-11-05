module module_gas_composition

  ! Dust only. 
  ! Useful for computing the intensity of non-ionizing radiation in a domain, since it is only affected by dust.
  ! - The dust content is a function of metallicity, HI and HII, following Laursen+09.

  use module_dust_model
  use module_random
  use module_ramses
  use module_constants

  implicit none

  private

  type, public :: gas
     ! fluid
     real(kind=8) :: v(3)      ! gas velocity [cm/s]
     ! Hydrogen 
     !real(kind=8) :: nHI       ! HI numerical density [HI/cm3]
     real(kind=8) :: dopwidth  ! Doppler width [cm/s]
     ! DUST -> model of Laursen, Sommer-Larsen and Andersen 2009.
     ! ->  ndust = (nHI + f_ion nHII)*Z/Zref
     ! f_ion and Zref are two free parameters . 
     real(kind=8) :: ndust     ! pseudo-numerical density of dust particles [#/cm3]
  end type gas
  real(kind=8),public :: box_size_cm   ! size of simulation box in cm. 

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [gas_composition] of the parameter file
  ! --------------------------------------------------------------------------
  ! mixture parameters 
  real(kind=8)             :: f_ion           = 0.01   ! ndust = (n_HI + f_ion*n_HII) * Z/Zsun [Laursen+09]
  real(kind=8)             :: Zref            = 0.005  ! reference metallicity. Should be ~ 0.005 for SMC and ~ 0.01 for LMC. 
  ! possibility to overwrite ramses values with an ad-hoc model 
  logical                  :: gas_overwrite       = .false. ! if true, define cell values from following parameters 
  real(kind=8)             :: fix_nhi             = 0.0d0   ! ad-hoc HI density (H/cm3)
  real(kind=8)             :: fix_vth             = 1.0d5   ! ad-hoc thermal velocity (cm/s)
  real(kind=8)             :: fix_ndust           = 0.0d0   ! ad-hoc dust number density (/cm3)
  real(kind=8)             :: fix_vel             = 0.0d0   ! ad-hoc cell velocity (cm/s) -> NEED BETTER PARAMETERIZATION for more than static... 
  real(kind=8)             :: fix_box_size_cm     = 1.0d8   ! ad-hoc box size in cm. 
  ! miscelaneous
  logical                  :: verbose             = .false. ! display some run-time info on this module
  ! --------------------------------------------------------------------------

  ! public functions:
  public :: gas_from_ramses_leaves,get_gas_velocity,gas_get_scatter_flag,gas_scatter,dump_gas
  public :: read_gas,gas_destructor,read_gas_composition_params,print_gas_composition_params

  !--PEEL--
  public :: gas_peeloff_weight,gas_get_tau,gas_get_tau_gray
  !--LEEP--

contains


  subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var, g)

    ! define gas contents from ramses raw data

    character(2000),intent(in)                     :: repository 
    integer(kind=4),intent(in)                     :: snapnum
    integer(kind=4),intent(in)                     :: nleaf,nvar
    real(kind=8),intent(in),dimension(nvar,nleaf)  :: ramses_var
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                                :: ileaf
    real(kind=8),dimension(:),allocatable          :: T, nhi, metallicity, nhii
    real(kind=8),dimension(:,:),allocatable        :: v

    ! allocate gas-element array
    allocate(g(nleaf))

    if (gas_overwrite) then
       call overwrite_gas(g)
    else

       box_size_cm = ramses_get_box_size_cm(repository,snapnum)

       ! compute velocities in cm / s
       if (verbose) write(*,*) '-- module_gas_composition_dust : extracting velocities from ramses '
       allocate(v(3,nleaf))
       call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
       do ileaf = 1,nleaf
          g(ileaf)%v = v(:,ileaf)
       end do
       deallocate(v)

       ! get nHI and temperature from ramses
       if (verbose) write(*,*) '-- module_gas_composition_dust : extracting nHI and T from ramses '
       allocate(T(nleaf),nhi(nleaf))
       call ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nhi)
       !g(:)%nHI = nhi(:)
       ! compute thermal velocity 
       ! ++++++ TURBULENT VELOCITY >>>>> parameter to add and use here
       !g(:)%dopwidth = sqrt((2.0d0*kb/mp)*T) ! [ cm/s ]

       ! get ndust (pseudo dust density from Laursen, Sommer-Larsen, Andersen 2009)
       if (verbose) write(*,*) '-- module_gas_composition_dust : extracting ndust from ramses '
       allocate(metallicity(nleaf),nhii(nleaf))
       call ramses_get_metallicity(nleaf,nvar,ramses_var,metallicity)
       call ramses_get_nh_cgs(repository,snapnum,nleaf,nvar,ramses_var,nhii)
       nhii = nhii - nhi
       do ileaf = 1,nleaf
          g(ileaf)%ndust = metallicity(ileaf) / Zref * ( nhi(ileaf) + f_ion*nhii(ileaf) )   ! [ /cm3 ]
       end do
       deallocate(metallicity,T,nhi,nhii)
    end if

    if (verbose) print*,'min/max of ndust : ',minval(g(:)%ndust),maxval(g(:)%ndust)

    return

  end subroutine gas_from_ramses_leaves



  subroutine overwrite_gas(g)
    ! overwrite ramses values with an ad-hoc model

    type(gas),dimension(:),intent(inout) :: g

    box_size_cm   = fix_box_size_cm

    g(:)%v(1)     = fix_vel
    g(:)%v(2)     = fix_vel
    g(:)%v(3)     = fix_vel
    !g(:)%nHI      = fix_nhi
    !g(:)%dopwidth = fix_vth
    g(:)%ndust    = fix_ndust

  end subroutine overwrite_gas



  function get_gas_velocity(cell_gas)
    type(gas),intent(in)      :: cell_gas
    real(kind=8),dimension(3) :: get_gas_velocity
    get_gas_velocity(:) = cell_gas%v(:)
    return
  end function get_gas_velocity


  function  gas_get_tau(cell_gas, distance_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute total opacity of gas accross distance_cm at freq. nu_cell
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - cell_gas : a mix of H, D, and dust
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

    ! compute optical depths for different components of the gas.
    gas_get_tau = get_tau_dust(cell_gas%ndust, distance_cm, nu_cell)

    return

  end function gas_get_tau
  ! --------------------------------------------------------------------------


  function  gas_get_tau_gray(cell_gas, distance_cm)

    ! --------------------------------------------------------------------------
    ! compute opacity due to dust accross distance_cm
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - cell_gas : only dust
    ! - distance_cm : the distance along which to compute tau [cm]
    ! OUTPUTS:
    ! - gas_get_tau_gray : the total optical depth
    ! --------------------------------------------------------------------------

    ! check whether scattering occurs within cell (scatter_flag > 0) or not (scatter_flag==0)
    type(gas),intent(in)    :: cell_gas
    real(kind=8),intent(in) :: distance_cm
    real(kind=8)            :: gas_get_tau_gray

    ! compute optical depths for different components of the gas.
    gas_get_tau_gray = get_tau_dust_gray(cell_gas%ndust, distance_cm)

    return

  end function gas_get_tau_gray
  ! --------------------------------------------------------------------------


  !--PEEL--
  function gas_peeloff_weight(flag,cell_gas,nu_ext,kin,kout,iran)

    integer(kind=4),intent(in)            :: flag
    type(gas),intent(in)                  :: cell_gas
    real(kind=8),intent(inout)            :: nu_ext
    real(kind=8),dimension(3), intent(in) :: kin, kout
    integer(kind=4),intent(inout)         :: iran
    real(kind=8)                          :: gas_peeloff_weight

    select case(flag)
    case(1)
       gas_peeloff_weight = dust_peeloff_weight(cell_gas%v, nu_ext, kin, kout)
    case default
       print*,'ERROR in module_gas_composition_dust.f90:gas_peeloff_weight - unknown case : ',flag 
       stop
    end select

  end function gas_peeloff_weight
  !--LEEP--



  function  gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs,iran)

    ! --------------------------------------------------------------------------
    ! Decide whether a scattering event occurs, and if so, on which element
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - cell_gas : a mix of H and dust
    ! - distance_to_border_cm : the maximum distance the photon may travel (before leaving the cell)
    ! - nu_cell : photon frequency in cell's frame [ Hz ]
    ! - tau_abs : optical depth at which the next scattering event will occur
    ! - iran    : random generator state of the photon 
    ! OUTPUTS:
    ! - distance_to_border_cm : comes out as the distance to scattering event (if there is an event)
    ! - tau_abs : 0 if a scatter occurs, decremented by tau_cell if photon escapes cell. 
    ! - gas_get_scatter_flag : 0 [no scatter], 1 [H scatter], 2 [dust]
    ! --------------------------------------------------------------------------

    ! check whether scattering occurs within cell (scatter_flag > 0) or not (scatter_flag==0)
    type(gas),intent(in)                  :: cell_gas
    real(kind=8),intent(inout)            :: distance_to_border_cm
    real(kind=8),intent(in)               :: nu_cell
    real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur.
    integer(kind=4),intent(inout)         :: iran
    integer(kind=4)                       :: gas_get_scatter_flag 
    real(kind=8)                          :: tau_cell !, tirage, proba1 !tau_HI, tau_dust, 

    ! compute optical depths for different components of the gas.
    !tau_HI   = get_tau_HI(cell_gas%nHI, cell_gas%dopwidth, distance_to_border_cm, nu_cell)
    !tau_dust = get_tau_dust(cell_gas%ndust, distance_to_border_cm, nu_cell)
    tau_cell = get_tau_dust(cell_gas%ndust, distance_to_border_cm, nu_cell)

    if (tau_abs > tau_cell) then  ! photon is due for absorption outside the cell 
       gas_get_scatter_flag = 0   ! no scatter
       tau_abs = tau_abs - tau_cell
       if (tau_abs.lt.0.0d0) then
          print*, 'tau_abs est negatif'
          stop
       endif
    else  ! the scattering happens inside the cell. 
       !  decide whether scattering is due to HI or dust
       !proba1 = tau_HI / tau_cell         
       !tirage = ran3(iran)
       !if(tirage <= proba1)then
       !gas_get_scatter_flag = 1 ! HI scatter
       ! else 
       gas_get_scatter_flag = 1 ! interaction with dust
       ! endif
       ! and transform "distance_to_border_cm" in "distance_to_absorption_cm"
       distance_to_border_cm = distance_to_border_cm * (tau_abs / tau_cell)
    end if

    return
  end function gas_get_scatter_flag



  subroutine gas_scatter(flag,cell_gas,nu_cell,k,nu_ext,iran)

    integer(kind=4),intent(inout)            :: flag
    type(gas),intent(in)                     :: cell_gas
    real(kind=8),intent(inout)               :: nu_cell, nu_ext
    real(kind=8),dimension(3), intent(inout) :: k
    integer(kind=4),intent(inout)            :: iran
    integer(kind=4)                          :: ilost

    select case(flag)
    case(1)
       ! call scatter_HI(cell_gas%v, cell_gas%dopwidth, nu_cell, k, nu_ext, iran)
       ! case(2)
       call scatter_dust(cell_gas%v, nu_cell, k, nu_ext, iran, ilost)
       if(ilost==1)flag=-1
    end select

  end subroutine gas_scatter



  subroutine dump_gas(unit,g)

    type(gas),dimension(:),intent(in) :: g
    integer(kind=4),intent(in)        :: unit
    integer(kind=4)                   :: i,nleaf
    nleaf = size(g)
    write(unit) (g(i)%v(:), i=1,nleaf)
    !write(unit) (g(i)%nHI, i=1,nleaf)
    !write(unit) (g(i)%dopwidth, i=1,nleaf)
    write(unit) (g(i)%ndust, i=1,nleaf)
    write(unit) box_size_cm 
  end subroutine dump_gas



  subroutine read_gas(unit,n,g)
    integer(kind=4),intent(in)                     :: unit,n
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                                :: i
    allocate(g(1:n))
    if (gas_overwrite) then
       call overwrite_gas(g)
    else
       read(unit) (g(i)%v(:),i=1,n)
       ! read(unit) (g(i)%nHI,i=1,n)
       !read(unit) (g(i)%dopwidth,i=1,n)
       read(unit) (g(i)%ndust,i=1,n)
       read(unit) box_size_cm 
    end if
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
          case ('gas_overwrite')
             read(value,*) gas_overwrite
          case ('fix_nhi')
             read(value,*) fix_nhi
          case ('fix_vth')
             read(value,*) fix_vth
          case ('fix_ndust')
             read(value,*) fix_ndust
          case ('fix_vel')
             read(value,*) fix_vel
          case ('verbose')
             read(value,*) verbose
          case ('fix_box_size_cm')
             read(value,*) fix_box_size_cm
          end select
       end do
    end if
    close(10)

    call read_ramses_params(pfile)
    !call read_HI_params(pfile)
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
       write(unit,'(a)')        '# mixture parameters'
       write(unit,'(a,ES10.3)') '  f_ion               = ',f_ion
       write(unit,'(a,ES10.3)') '  Zref                = ',Zref
       write(unit,'(a)')        '# overwrite parameters'
       write(unit,'(a,L1)')     '  gas_overwrite       = ',gas_overwrite
       write(unit,'(a,ES10.3)') '  fix_nhi             = ',fix_nhi
       write(unit,'(a,ES10.3)') '  fix_vth             = ',fix_vth
       write(unit,'(a,ES10.3)') '  fix_ndust           = ',fix_ndust
       write(unit,'(a,ES10.3)') '  fix_vel             = ',fix_vel
       write(unit,'(a,ES10.3)') '  fix_box_size_cm     = ',fix_box_size_cm
       write(unit,'(a)')        '# miscelaneous parameters'
       write(unit,'(a,L1)')     '  verbose             = ',verbose
       write(unit,'(a)')             ' '
       call print_ramses_params(unit)
       !write(unit,'(a)')             ' '
       !call print_HI_params(unit)
       write(unit,'(a)')             ' '
       call print_dust_params(unit)
    else
       write(*,'(a,a,a)') '[gas_composition]'
       write(*,'(a)')        '# mixture parameters'
       write(*,'(a,ES10.3)') '  f_ion               = ',f_ion
       write(*,'(a,ES10.3)') '  Zref                = ',Zref
       write(*,'(a)')        '# overwrite parameters'
       write(*,'(a,L1)')     '  gas_overwrite       = ',gas_overwrite
       write(*,'(a,ES10.3)') '  fix_nhi             = ',fix_nhi
       write(*,'(a,ES10.3)') '  fix_vth             = ',fix_vth
       write(*,'(a,ES10.3)') '  fix_ndust           = ',fix_ndust
       write(*,'(a,ES10.3)') '  fix_vel             = ',fix_vel
       write(*,'(a,ES10.3)') '  fix_box_size_cm     = ',fix_box_size_cm
       write(*,'(a)')        '# miscelaneous parameters'
       write(*,'(a,L1)')     '  verbose             = ',verbose
       write(*,'(a)')             ' '
       call print_ramses_params
       !write(*,'(a)')             ' '
       !call print_HI_params
       write(*,'(a)')             ' '
       call print_dust_params
    end if

    return

  end subroutine print_gas_composition_params


end module module_gas_composition