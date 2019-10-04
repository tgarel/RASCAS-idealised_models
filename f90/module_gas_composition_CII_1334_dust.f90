module module_gas_composition

  ! Mix of CII and dust 
  ! This modules handles the transition in absorption (1334.53) and dust

  !!! Only works for Idealised model:
  ! ramses_get_T_nCII_cgs not implemented if processing a RAMSES simulation
  ! Still use ramses_get_T_nSiII_cgs instead
  
  use module_CII_1334_model
  use module_dust_model
  use module_random
  use module_constants
  use module_idealised_models

  implicit none

  private

  character(100),parameter :: moduleName = 'module_gas_composition_CII_1334_dust.f90'
  type, public :: gas
     ! fluid
     real(kind=8) :: v(3)      ! gas velocity [cm/s]
     ! CII
     ! -> density is computed as 2.8d-5 * nHI * metallicity / solar_metallicity
     real(kind=8) :: nCII     ! numerical density of CII  [#/cm3]
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
  real(kind=8)             :: fix_nCII           = 0.0d0   ! ad-hoc HI density (H/cm3)
  real(kind=8)             :: fix_vth             = 1.0d5   ! ad-hoc thermal velocity (cm/s)
  real(kind=8)             :: fix_vel             = 0.0d0   ! ad-hoc cell velocity (cm/s) -> NEED BETTER PARAMETERIZATION for more than static... 
  real(kind=8)             :: fix_ndust           = 0.0d0
  real(kind=8)             :: fix_box_size_cm     = 1.0d8   ! ad-hoc box size in cm.
  ! --------------------------------------------------------------------------

  ! public functions:
  public :: gas_from_ramses_leaves,get_gas_velocity,gas_get_scatter_flag,gas_scatter,dump_gas
  ! TIBO
  public :: gas_from_idealised_models
  ! OBIT

  public :: read_gas,gas_destructor,read_gas_composition_params,print_gas_composition_params

contains
  
  ! TIBO
  subroutine gas_from_idealised_models(outputdir, nleaf, g, x_leaf, leaf_level)
    
    integer(kind=4),intent(in)                       :: nleaf
    type(gas),dimension(:),allocatable,intent(inout) :: g
    integer(kind=4)                                  :: ileaf
    real(kind=8),intent(in),dimension(nleaf,3)       :: x_leaf
    integer(kind=4),intent(in),dimension(nleaf)      :: leaf_level
    real(kind=8),dimension(:),allocatable            :: ndust_temp,ngas_temp,dopwidth_temp
    real(kind=8),dimension(:,:),allocatable          :: vel_temp
    character(2000),intent(in)                       :: outputdir
    character(300)                                   :: modelprops_file,file
    
    ! allocate gas-element array
    allocate(g(nleaf))

    box_size_cm = box_size_IM_cm

    allocate(ndust_temp(nleaf),ngas_temp(nleaf),dopwidth_temp(nleaf))
    allocate(vel_temp(3,nleaf))
    
    ndust_temp(:)    = 0.0d0
    ngas_temp(:)     = 0.0d0
    dopwidth_temp(:) = 0.0d0
    vel_temp(:,:)    = 0.0d0
    
    call compute_idealised_gas(ndust_temp,ngas_temp,dopwidth_temp,vel_temp,x_leaf,leaf_level,nleaf)

    do ileaf = 1,nleaf
       g(ileaf)%v = vel_temp(:,ileaf)
    end do
    
    ! if no-dust module, comment the line below 
    g(:)%ndust    = ndust_temp(:)
    ! adapt n_species name to each module_gas_composition
    g(:)%nCII    = ngas_temp(:)

    ! Deal with m_atom here: adapt m_species name to each module_gas_composition

    ! HERE, I read true dopwidth (in cm/s) from param file (i.e. vth_norm, read in module_idealised_models)
    ! Attention, 'temp' in dopwidth_temp stands for temporary, not T..... ;-)
    ! g(:)%dopwidth = dopwidth_temp(:) / sqrt(mSi)
    g(:)%dopwidth = dopwidth_temp(:)
    
!!$    modelprops_file = 'modelprops_file'
!!$    file = trim(outputdir)//trim(modelprops_file)
!!$    open(unit=15, file=trim(file), status='unknown', form='unformatted', action='write')
!!$    write(15) nleaf
!!$    do ileaf = 1,nleaf
!!$       write(15) x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),g(ileaf)%v(1),g(ileaf)%v(2),g(ileaf)%v(3),g(ileaf)%dopwidth,g(ileaf)%nCII,g(ileaf)%ndust
!!$    end do
!!$    close(15)
    
    deallocate(ndust_temp,ngas_temp,dopwidth_temp,vel_temp)

    return
    
  end subroutine gas_from_idealised_models
  ! OBIT
  
  subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var, g)

    ! define gas contents from ramses raw data

    use module_ramses

    character(2000),intent(in)        :: repository 
    integer(kind=4),intent(in)        :: snapnum
    integer(kind=4),intent(in)        :: nleaf,nvar
    real(kind=8),intent(in)           :: ramses_var(nvar,nleaf)
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                   :: ileaf
    real(kind=8),allocatable          :: v(:,:), T(:), nCII(:), nHI(:), nHII(:), metallicity(:)

    ! allocate gas-element array
    allocate(g(nleaf))

    if (gas_overwrite) then
       call overwrite_gas(g)
    else
       
       box_size_cm = ramses_get_box_size_cm(repository,snapnum)

       ! compute velocities in cm / s
       write(*,*) '-- module_gas_composition_CII_dust : extracting velocities from ramses '
       allocate(v(3,nleaf))
       call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
       do ileaf = 1,nleaf
          g(ileaf)%v = v(:,ileaf)
       end do
       deallocate(v)

       ! get nCII and temperature from ramses
       write(*,*) '-- module_gas_composition_CII_dust : extracting nCII from ramses '
       allocate(T(nleaf),nCII(nleaf))
       call ramses_get_T_nSiII_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nCII)
       g(:)%nCII = nCII(:)
       ! compute thermal velocity 
       ! ++++++ TURBULENT VELOCITY >>>>> parameter to add and use here
       g(:)%dopwidth = sqrt(2.0d0*kb/mSi*T) ! [ cm/s ]
       deallocate(T,nCII)

       print*,'min/max of nCII : ',minval(g(:)%nCII),maxval(g(:)%nCII)
       
       ! get ndust (pseudo dust density from Laursen, Sommer-Larsen, Andersen 2009)
       write(*,*) '-- module_gas_composition_CII_dust : extracting ndust from ramses '
       allocate(T(nleaf),nhi(nleaf),metallicity(nleaf),nhii(nleaf))
       call ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nhi)
       call ramses_get_metallicity(nleaf,nvar,ramses_var,metallicity)
       call ramses_get_nh_cgs(repository,snapnum,nleaf,nvar,ramses_var,nhii)
       nhii = nhii - nhi
       do ileaf = 1,nleaf
          g(ileaf)%ndust = metallicity(ileaf) / Zref * ( nhi(ileaf) + f_ion*nhii(ileaf) )   ! [ /cm3 ]
       end do
       deallocate(metallicity,T,nhi,nhii)

       print*,'min/max of ndust : ',minval(g(:)%ndust),maxval(g(:)%ndust)
       
    end if

    return

  end subroutine gas_from_ramses_leaves
  

  subroutine overwrite_gas(g)

    type(gas),dimension(:),intent(inout) :: g

    box_size_cm   = fix_box_size_cm
    
    g(:)%v(1)     = fix_vel
    g(:)%v(2)     = fix_vel
    g(:)%v(3)     = fix_vel
    g(:)%nCII    = fix_nCII
    g(:)%dopwidth = fix_vth
    g(:)%ndust    = fix_ndust
    
  end subroutine overwrite_gas



  function get_gas_velocity(cell_gas)
    type(gas),intent(in)      :: cell_gas
    real(kind=8),dimension(3) :: get_gas_velocity
    get_gas_velocity(:) = cell_gas%v(:)
    return
  end function get_gas_velocity


  
  function  gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs, iran)

    ! --------------------------------------------------------------------------
    ! Decide whether a scattering event occurs, and if so, on which element
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - cell_gas : CII (with two absorption channels) and dust
    ! - distance_to_border_cm : the maximum distance the photon may travel (before leaving the cell)
    ! - nu_cell : photon frequency in cell's frame [ Hz ]
    ! - tau_abs : optical depth at which the next scattering event will occur
    ! - iran    : random generator state of the photon 
    ! OUTPUTS:
    ! - distance_to_border_cm : comes out as the distance to scattering event (if there is an event)
    ! - tau_abs : 0 if a scatter occurs, decremented by tau_cell if photon escapes cell. 
    ! - gas_get_scatter_flag : 0 [no scatter], 1 [CII-1334 scatter], 2 [dust] 
    ! --------------------------------------------------------------------------

    type(gas),intent(in)                  :: cell_gas
    real(kind=8),intent(inout)            :: distance_to_border_cm
    real(kind=8),intent(in)               :: nu_cell
    real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur.
    integer,intent(inout)                 :: iran 
    integer(kind=4)                       :: gas_get_scatter_flag 
    real(kind=8)                          :: tau_CII_1334, tau_dust, tau_cell, proba1, x
    
    ! compute optical depths for different components of the gas.
    tau_CII_1334 = get_tau_CII_1334(cell_gas%nCII, cell_gas%dopwidth, distance_to_border_cm, nu_cell)
    tau_dust      = get_tau_dust(cell_gas%ndust, distance_to_border_cm, nu_cell)
    tau_cell      = tau_CII_1334 + tau_dust

    if (tau_abs > tau_cell) then  ! photon is due for absorption outside the cell 
       gas_get_scatter_flag = 0
       tau_abs = tau_abs - tau_cell
       if (tau_abs.lt.0.d0) then
          print*, 'tau_abs est negatif'
          stop
       endif
    else  ! the scattering happens inside the cell
       ! decide if it is 1334 or dust 
       proba1 = tau_CII_1334 /  tau_cell
       x = ran3(iran)
       if (x <= proba1) then
          gas_get_scatter_flag = 1 ! absorption by CII-1334
       else
          gas_get_scatter_flag = 2 ! absorption by dust
       end if
       ! and transform "distance_to_border_cm" in "distance_to_absorption_cm"
       distance_to_border_cm = distance_to_border_cm * (tau_abs / tau_cell)
    end if

    return

  end function gas_get_scatter_flag



  subroutine gas_scatter(flag,cell_gas,nu_cell,k,nu_ext,iran)

    integer, intent(inout)                    :: flag
    type(gas), intent(in)                     :: cell_gas
    real(kind=8), intent(inout)               :: nu_cell, nu_ext
    real(kind=8), dimension(3), intent(inout) :: k
    integer, intent(inout)                    :: iran
    integer(kind=4)                           :: ilost 

    select case(flag)
    case(1)
       call scatter_CII_1334(cell_gas%v, cell_gas%dopwidth, nu_cell, k, nu_ext, iran)
    case(2)
       call scatter_dust(cell_gas%v, nu_cell, k, nu_ext, iran, ilost)
       if(ilost==1)flag=-1
    end select

  end subroutine gas_scatter



  subroutine dump_gas(unit,g)
    type(gas),dimension(:),intent(in) :: g
    integer,intent(in)                :: unit
    integer                           :: i,nleaf
    nleaf = size(g)
    write(unit) (g(i)%v(:), i=1,nleaf)
    write(unit) (g(i)%nCII, i=1,nleaf)
    write(unit) (g(i)%dopwidth, i=1,nleaf)
    write(unit) (g(i)%ndust, i=1,nleaf)
    write(unit) box_size_cm 
  end subroutine dump_gas



  subroutine read_gas(unit,n,g)
    integer,intent(in)                             :: unit,n
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer                                        :: i
    allocate(g(1:n))
    if (gas_overwrite) then
       call overwrite_gas(g)
    else
       read(unit) (g(i)%v(:),i=1,n)
       read(unit) (g(i)%nCII,i=1,n)
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
    ! ALSO read parameter form used modules (HI_model)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present

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
          case ('fix_nCII')
             read(value,*) fix_nCII
          case ('fix_vth')
             read(value,*) fix_vth
          case ('fix_vel')
             read(value,*) fix_vel
          case ('fix_ndust')
             read(value,*) fix_ndust
          case ('fix_box_size_cm')
             read(value,*) fix_box_size_cm
          end select
       end do
    end if
    close(10)

    call read_dust_params(pfile)
    call read_CII_1334_params(pfile)

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
       write(unit,'(a,ES10.3)') '  f_ion                = ',f_ion
       write(unit,'(a,ES10.3)') '  Zref                 = ',Zref
       write(unit,'(a)')       '# overwrite parameters'
       write(unit,'(a,L1)')    '  gas_overwrite         = ',gas_overwrite
       if(gas_overwrite)then
          write(unit,'(a,ES10.3)') '  fix_nCII            = ',fix_nCII
          write(unit,'(a,ES10.3)') '  fix_vth              = ',fix_vth
          write(unit,'(a,ES10.3)') '  fix_vel              = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_box_size_cm      = ',fix_box_size_cm
       endif
       write(unit,'(a)')             ' '
       call print_dust_params(unit)
       write(unit,'(a)')             ' '
       call print_CII_1334_params(unit)
    else
       write(*,'(a,a,a)') '[gas_composition]'
       write(*,'(a,a)')      '# code compiled with: ',trim(moduleName)
       write(*,'(a)')        '# mixture parameters'
       write(*,'(a,ES10.3)') '  f_ion                = ',f_ion
       write(*,'(a,ES10.3)') '  Zref                 = ',Zref
       write(*,'(a)')       '# overwrite parameters'
       write(*,'(a,L1)')    '  gas_overwrite         = ',gas_overwrite
       if(gas_overwrite)then
          write(*,'(a,ES10.3)') '  fix_nCII            = ',fix_nCII
          write(*,'(a,ES10.3)') '  fix_vth              = ',fix_vth
          write(*,'(a,ES10.3)') '  fix_vel              = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_box_size_cm      = ',fix_box_size_cm
       endif
       write(*,'(a)')             ' '
       call print_dust_params
       write(*,'(a)')             ' '
       call print_CII_1334_params
    end if

    return

  end subroutine print_gas_composition_params



end module module_gas_composition
