module module_gas_composition

  ! mix of HI, Deuterium, and dust.
  ! - The Deuterium density is set as a fraction of H density through the parameter deut2H_nb_ratio.
  ! - The thermal velocity of D is derived from the one of H using the mass ratio (sqrt_H2Deut_mass_ratio in module_constants.f90).
  ! - The HI content is from RAMSES
  ! - The dust content is a function of metallicity, HI and HII, following Laursen+09.

  use module_HI_1216_model
  use module_D_1215_model
  use module_dust_model
  use module_random
  use module_constants

  implicit none

  private

  character(100),parameter :: moduleName = 'module_gas_composition_HI_D_dust.f90'
  
  type, public :: gas
     ! fluid
     real(kind=8) :: v(3)      ! gas velocity [cm/s]
     ! Hydrogen 
     real(kind=8) :: nHI       ! HI numerical density [HI/cm3]
     real(kind=8) :: dopwidth  ! Doppler width [cm/s]
     ! Deuterium
     ! -> density is computed as nHI * deut2H_nb_ratio
     ! -> dopwidth is computed as dopwidth * sqrt_H2Deut_mass_ratio.
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
  real(kind=8)             :: deut2H_nb_ratio = 3.d-5  ! D to H number ratio. Default is primordial value 3e-5
  real(kind=8)             :: f_ion           = 0.01   ! ndust = (n_HI + f_ion*n_HII) * Z/Zsun [Laursen+09]
  real(kind=8)             :: Zref            = 0.005  ! reference metallicity. Should be ~ 0.005 for SMC and ~ 0.01 for LMC. 
  ! possibility to overwrite ramses values with an ad-hoc model 
  logical                  :: gas_overwrite       = .false. ! if true, define cell values from following parameters 
  real(kind=8)             :: fix_nhi             = 0.0d0   ! ad-hoc HI density (H/cm3)
  real(kind=8)             :: fix_vth             = 1.0d5   ! ad-hoc thermal velocity (cm/s)
  real(kind=8)             :: fix_ndust           = 0.0d0   ! ad-hoc dust number density (/cm3)
  real(kind=8)             :: fix_vel             = 0.0d0   ! ad-hoc cell velocity (cm/s) -> NEED BETTER PARAMETERIZATION for more than static... 
  real(kind=8)             :: fix_box_size_cm     = 1.0d8   ! ad-hoc box size in cm. 
  ! --------------------------------------------------------------------------

  ! public functions:
  public :: gas_from_ramses_leaves,get_gas_velocity,gas_get_scatter_flag,gas_scatter,dump_gas
  public :: read_gas,gas_destructor,read_gas_composition_params,print_gas_composition_params

contains

  
  subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var, g)

    ! define gas contents from ramses raw data

    use module_ramses
    
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
       write(*,*) '-- module_gas_composition_HI_D_dust : extracting velocities from ramses '
       allocate(v(3,nleaf))
       call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
       do ileaf = 1,nleaf
          g(ileaf)%v = v(:,ileaf)
       end do
       deallocate(v)

       ! get nHI and temperature from ramses
       write(*,*) '-- module_gas_composition_HI_D_dust : extracting nHI and T from ramses '
       allocate(T(nleaf),nhi(nleaf))
       call ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nhi)
       g(:)%nHI = nhi(:)
       ! compute thermal velocity 
       ! ++++++ TURBULENT VELOCITY >>>>> parameter to add and use here
       g(:)%dopwidth = sqrt((2.0d0*kb/mp)*T) ! [ cm/s ]

       ! get ndust (pseudo dust density from Laursen, Sommer-Larsen, Andersen 2009)
       write(*,*) '-- module_gas_composition_HI_D_dust : extracting ndust from ramses '
       allocate(metallicity(nleaf),nhii(nleaf))
       call ramses_get_metallicity(nleaf,nvar,ramses_var,metallicity)
       call ramses_get_nh_cgs(repository,snapnum,nleaf,nvar,ramses_var,nhii)
       nhii = nhii - nhi
       do ileaf = 1,nleaf
          g(ileaf)%ndust = metallicity(ileaf) / Zref * ( nhi(ileaf) + f_ion*nhii(ileaf) )   ! [ /cm3 ]
       end do
       deallocate(metallicity,T,nhi,nhii)
    end if

    return

  end subroutine gas_from_ramses_leaves


  
  subroutine overwrite_gas(g)
    ! overwrite ramses values with an ad-hoc model

    type(gas),dimension(:),intent(inout) :: g

    box_size_cm   = fix_box_size_cm
    
    g(:)%v(1)     = fix_vel
    g(:)%v(2)     = fix_vel
    g(:)%v(3)     = fix_vel
    g(:)%nHI      = fix_nhi
    g(:)%dopwidth = fix_vth
    g(:)%ndust    = fix_ndust
    
  end subroutine overwrite_gas



  function get_gas_velocity(cell_gas)
    type(gas),intent(in)      :: cell_gas
    real(kind=8),dimension(3) :: get_gas_velocity
    get_gas_velocity(:) = cell_gas%v(:)
    return
  end function get_gas_velocity



  !--CORESKIP--
  function  gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs, iran, CS_dist_cm, CS_xcrit)

    ! --------------------------------------------------------------------------
    ! Decide whether a scattering event occurs, and if so, on which element
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - cell_gas : a mix of H, D, and dust
    ! - distance_to_border_cm : the maximum distance the photon may travel (before leaving the cell)
    ! - nu_cell : photon frequency in cell's frame [ Hz ]
    ! - tau_abs : optical depth at which the next scattering event will occur
    ! - iran    : random generator state of the photon
    ! - CS_dist_cm: true distance to cell border, not along k. For Core-Skipping only.
    ! OUTPUTS:
    ! - distance_to_border_cm : comes out as the distance to scattering event (if there is an event)
    ! - tau_abs : 0 if a scatter occurs, decremented by tau_cell if photon escapes cell. 
    ! - gas_get_scatter_flag : 0 [no scatter], 1 [H scatter], 2 [D scatter], 3 [dust]
    ! - CS_xcrit: critical frequency above which core-skipping applies. For Core-Skipping only.
    ! --------------------------------------------------------------------------

    ! check whether scattering occurs within cell (scatter_flag > 0) or not (scatter_flag==0)
    type(gas),intent(in)                  :: cell_gas
    real(kind=8),intent(inout)            :: distance_to_border_cm
    real(kind=8),intent(in)               :: nu_cell
    real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur.
    integer(kind=4),intent(inout)         :: iran
    integer(kind=4)                       :: gas_get_scatter_flag 
    real(kind=8)                          :: tau_HI, tau_dust, tau_cell, tau_D, tirage, proba1, proba2
    !--CORESKIP--
    real(kind=8),intent(in)               :: CS_dist_cm
    real(kind=8),intent(inout)            :: CS_xcrit
    real(kind=8)                          :: CS_tau_cell,delta_nu_doppler,a,xcw,x
    !--PIKSEROC--
    
    ! compute optical depths for different components of the gas.
    tau_HI   = get_tau_HI_1216(cell_gas%nHI, cell_gas%dopwidth, distance_to_border_cm, nu_cell)
    tau_dust = get_tau_dust(cell_gas%ndust, distance_to_border_cm, nu_cell)
    tau_D    = get_tau_D_1215(cell_gas%nHI * deut2H_nb_ratio, cell_gas%dopwidth * sqrt_H2Deut_mass_ratio,distance_to_border_cm, nu_cell)
    tau_cell = tau_HI + tau_D + tau_dust

    if (tau_abs > tau_cell) then  ! photon is due for absorption outside the cell 
       gas_get_scatter_flag = 0   ! no scatter
       tau_abs = tau_abs - tau_cell
       if (tau_abs.lt.0.0d0) then
          print*, 'tau_abs est negatif'
          stop
       endif
    else  ! the scattering happens inside the cell. 
       ! decide whether scattering is due to HI, D or dust
       proba1 = tau_HI / tau_cell         
       proba2 = proba1 + tau_D / tau_cell
       tirage = ran3(iran)
       if(tirage <= proba1)then
          gas_get_scatter_flag = 1 ! HI scatter
       else if (tirage <= proba2) then 
          gas_get_scatter_flag = 2 ! interaction with a Deuterium atom 
       else 
          gas_get_scatter_flag = 3 ! interaction with dust
       endif
       ! and transform "distance_to_border_cm" in "distance_to_absorption_cm"
       distance_to_border_cm = distance_to_border_cm * (tau_abs / tau_cell)
    end if

    if (HI_core_skip) then 
       delta_nu_doppler = cell_gas%dopwidth/lambda_0_cm
       a = 6.265d8/fourpi/delta_nu_doppler
       xcw = 6.9184721d0 + 81.766279d0 / (log10(a)-14.651253d0)  ! Smith+15, Eq. 21
       x = (nu_cell - nu_0)/delta_nu_doppler
       if (abs(x) < xcw) then ! apply core-skipping
          CS_tau_cell = get_tau_HI_1216(cell_gas%nHI, cell_gas%dopwidth, CS_dist_cm, nu_cell)
          CS_tau_cell = CS_tau_cell + get_tau_D_1215(cell_gas%nHI * deut2H_nb_ratio, cell_gas%dopwidth * sqrt_H2Deut_mass_ratio, CS_dist_cm, nu_cell)
          CS_tau_cell = CS_tau_cell +get_tau_dust(cell_gas%ndust, CS_dist_cm, nu_cell)
          CS_tau_cell = CS_tau_cell * a
          if (CS_tau_cell > 1.0d0) CS_xcrit = CS_tau_cell**(1./3.)/5.  ! Smith+15, Eq. 35
       end if
    end if
    
    
    return
  end function gas_get_scatter_flag



  !--CORESKIP-- 
  subroutine gas_scatter(flag,cell_gas,nu_cell,k,nu_ext,iran,xcrit)
  !subroutine gas_scatter(flag,cell_gas,nu_cell,k,nu_ext,iran)
  !--PIKSEROC--
    
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
    case(1)
       !--CORESKIP--
       call scatter_HI_1216(cell_gas%v, cell_gas%dopwidth, nu_cell, k, nu_ext, iran, xcrit)
       !call scatter_HI_1216(cell_gas%v, cell_gas%dopwidth, nu_cell, k, nu_ext, iran)
       !--PIKSEROC--
    case(2)
       call scatter_D_1215(cell_gas%v,cell_gas%dopwidth*sqrt_H2Deut_mass_ratio, nu_cell, k, nu_ext, iran)
    case(3)
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
    write(unit) (g(i)%nHI, i=1,nleaf)
    write(unit) (g(i)%dopwidth, i=1,nleaf)
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
       read(unit) (g(i)%nHI,i=1,n)
       read(unit) (g(i)%dopwidth,i=1,n)
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
          case ('deut2H_nb_ratio')
             read(value,*) deut2H_nb_ratio
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
          case ('fix_box_size_cm')
             read(value,*) fix_box_size_cm
          end select
       end do
    end if
    close(10)

    call read_HI_1216_params(pfile)
    call read_D_1215_params(pfile)
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
       write(unit,'(a,ES10.3)') '  deut2H_nb_ratio = ',deut2H_nb_ratio
       write(unit,'(a,ES10.3)') '  f_ion           = ',f_ion
       write(unit,'(a,ES10.3)') '  Zref            = ',Zref
       write(unit,'(a)')        '# overwrite parameters'
       write(unit,'(a,L1)')     '  gas_overwrite   = ',gas_overwrite
       if(gas_overwrite)then
          write(unit,'(a,ES10.3)') '  fix_nhi         = ',fix_nhi
          write(unit,'(a,ES10.3)') '  fix_vth         = ',fix_vth
          write(unit,'(a,ES10.3)') '  fix_ndust       = ',fix_ndust
          write(unit,'(a,ES10.3)') '  fix_vel         = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_box_size_cm = ',fix_box_size_cm
       endif
       write(unit,'(a)')             ' '
       call print_HI_1216_params(unit)
       write(unit,'(a)')             ' '
       call print_D_1215_params(unit)
       write(unit,'(a)')             ' '
       call print_dust_params(unit)
    else
       write(*,'(a,a,a)') '[gas_composition]'
       write(*,'(a,a)')      '# code compiled with: ',trim(moduleName)
       write(*,'(a)')        '# mixture parameters'
       write(*,'(a,ES10.3)') '  deut2H_nb_ratio = ',deut2H_nb_ratio
       write(*,'(a,ES10.3)') '  f_ion           = ',f_ion
       write(*,'(a,ES10.3)') '  Zref            = ',Zref
       write(*,'(a)')        '# overwrite parameters'
       write(*,'(a,L1)')     '  gas_overwrite   = ',gas_overwrite
       if(gas_overwrite)then
          write(*,'(a,ES10.3)') '  fix_nhi         = ',fix_nhi
          write(*,'(a,ES10.3)') '  fix_vth         = ',fix_vth
          write(*,'(a,ES10.3)') '  fix_ndust       = ',fix_ndust
          write(*,'(a,ES10.3)') '  fix_vel         = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_box_size_cm = ',fix_box_size_cm
       endif
       write(*,'(a)')             ' '
       call print_HI_1216_params
       write(*,'(a)')             ' '
       call print_D_1215_params
       write(*,'(a)')             ' '
       call print_dust_params
    end if

    return

  end subroutine print_gas_composition_params


end module module_gas_composition
