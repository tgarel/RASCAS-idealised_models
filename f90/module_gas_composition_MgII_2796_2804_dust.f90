module module_gas_composition

  ! Mix of MgII and dust. 
  ! This modules handles two transitions in absorption by MgII (2796 and 2804) and dust 

  use module_MgII_2796_model
  use module_MgII_2804_model
  use module_dust_model
  use module_random
  use module_constants

  implicit none

  private

  character(100),parameter :: moduleName = 'module_gas_composition_MgII_2796_2804_dust.f90'
  type, public :: gas
     ! fluid
     real(kind=8) :: v(3)      ! gas velocity [cm/s]
     ! MgII -> density is computed as abundance_number_MgII * nHI * metallicity / solar_metallicity
     real(kind=8) :: nMgII     ! numerical density of MgII  [#/cm3]
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
  character(2000)          :: input_ramses_file = './ramses'  !Path to the file 'input_ramses', with all the info on the simulation (ncells, box_size, positions, velocities, nHI, nHII, Z, dopwidth
  character(2000)          :: Ion_file          = './ions'    !Path to the file with the MgII densities
  real(kind=8)             :: f_ion           = 0.01   ! ndust = (n_HI + f_ion*n_HII) * Z/Zsun [Laursen+09]
  real(kind=8)             :: Zref            = 0.005  ! reference metallicity. Should be ~ 0.005 for SMC and ~ 0.01 for LMC. 
  ! possibility to overwrite ramses values with an ad-hoc model 
  logical                  :: gas_overwrite       = .false. ! if true, define cell values from following parameters 
  real(kind=8)             :: fix_nMgII           = 0.0d0   ! ad-hoc HI density (H/cm3)
  real(kind=8)             :: fix_vth             = 1.0d5   ! ad-hoc thermal velocity (cm/s)
  real(kind=8)             :: fix_vel             = 0.0d0   ! ad-hoc cell velocity (cm/s) -> NEED BETTER PARAMETERIZATION for more than static... 
  real(kind=8)             :: fix_ndust           = 0.0d0
  real(kind=8)             :: fix_box_size_cm     = 1.0d8   ! ad-hoc box size in cm. 
  logical                  :: HI_core_skip        = .false.
  logical                  :: no_v                = .false.
  logical                  :: no_dust             = .false.
  ! --------------------------------------------------------------------------

  ! public functions:
  public :: gas_from_ramses_leaves,gas_from_ramses_leaves_ions,gas_from_list,get_gas_velocity,gas_get_scatter_flag,gas_scatter,dump_gas
  public :: read_gas,gas_destructor,read_gas_composition_params,print_gas_composition_params
  !Val
  public :: gas_get_n_CD, gas_get_CD
  !laV
  !--PEEL--
  public :: gas_peeloff_weight,gas_get_tau
  !--LEEP--

  !--CORESKIP-- 
  public :: HI_core_skip
  !--PIKSEROC--

contains

  !Val
  function gas_get_n_CD()

    integer(kind=4)           :: gas_get_n_CD

    gas_get_n_CD = 2

  end function gas_get_n_CD
  !laV

  !Val
  function gas_get_CD(cell_gas, distance_cm)

    real(kind=8), intent(in) :: distance_cm
    type(gas),intent(in)     :: cell_gas
    real(kind=8)             :: gas_get_CD(2)

    gas_get_CD = (/ distance_cm*cell_gas%nMgII, distance_cm*cell_gas%ndust /)

  end function gas_get_CD
  !laV

   subroutine gas_from_list(ncells, cell_pos, cell_l, gas_leaves)

    integer(kind=4), intent(out)		:: ncells
    integer(kind=4), allocatable, intent(out)	:: cell_l(:)
    real(kind=8), allocatable, intent(out)	:: cell_pos(:,:)
    type(gas), allocatable, intent(out)		:: gas_leaves(:)

    integer(kind=4)				:: i
    real(kind=8), allocatable			:: nhi(:), nhii(:), metallicity(:)


   open(unit=20, file=input_ramses_file, status='old', form='unformatted')
    open(unit=21, file=Ion_file, status='old', form='unformatted')  !remove for 2*2*2 test

    read(20) ncells

    allocate(cell_l(ncells), cell_pos(ncells,3), gas_leaves(ncells), nhi(ncells), nhii(ncells), metallicity(ncells))

    read(20) box_size_cm
    read(20) cell_l
    read(20) cell_pos(:,1)
    read(20) cell_pos(:,2)
    read(20) cell_pos(:,3)
    read(20) gas_leaves%v(1)
    read(20) gas_leaves%v(2)
    read(20) gas_leaves%v(3)
    read(20) nhi          !remove for 2*2*2 test
    read(20) nhii         !remove for 2*2*2 test
    read(20) metallicity  !remove for 2*2*2 test
    read(20) gas_leaves%dopwidth
     gas_leaves%dopwidth = gas_leaves%dopwidth / sqrt(mMg/amu)  ! remove for 2*2*2 test
    !read(20) gas_leaves%ndust !add for 2*2*2 test
    !read(20) gas_leaves%nMgII !add for 2*2*2 test
    read(21) gas_leaves%nMgII !remove for 2*2*2 test

    close(20)
    close(21)  !remove for 2*2*2 test
    
    gas_leaves%ndust = metallicity / Zref * ( nhi + f_ion*nhii )   ! [ /cm3 ]   !remove for 2*2*2 test


    print*,'boxsize in cm : ', box_size_cm
    print*,'min/max of vth   : ',minval(gas_leaves(:)%dopwidth),maxval(gas_leaves(:)%dopwidth)
    print*,'min/max of nMgII : ',minval(gas_leaves(:)%nMgII),maxval(gas_leaves(:)%nMgII)
    print*,'min/max of ndust : ',minval(gas_leaves(:)%ndust),maxval(gas_leaves(:)%ndust)

    deallocate(nhi, nhii, metallicity)


  end subroutine gas_from_list
  ! --------------------------------------------------------------------------
  !--laV

 !Val--
  ! --------------------------------------------------------------------------
  subroutine gas_from_ramses_leaves_ions(repository,snapnum,nleaf,nvar,ramses_var,ion_number,MgII_density,g)

    use module_ramses

    character(2000),intent(in)                          :: repository 
    integer(kind=4),intent(in)                          :: snapnum, ion_number
    integer(kind=4),intent(in)                          :: nleaf,nvar
    real(kind=8),intent(in),dimension(nvar,nleaf)       :: ramses_var
    real(kind=8),intent(in),dimension(ion_number,nleaf) :: MgII_density
    type(gas),dimension(:),allocatable,intent(out)      :: g
    integer(kind=4)                                     :: ileaf
    real(kind=8),dimension(:),allocatable               :: T, nhi, metallicity, nhii
    real(kind=8),dimension(:,:),allocatable             :: v

    if(ion_number /= 1) then
       print*, 'Error in module_gas_composition_MgII_2796_2804_dust.f90,  the number of ions in ion_parameter_file should be 1.'
       stop
    end if
    
    ! allocate gas-element array
    allocate(g(nleaf))

    if (gas_overwrite) then
       call overwrite_gas(g)
    else
       
       box_size_cm = ramses_get_box_size_cm(repository,snapnum)

       ! compute velocities in cm / s
       write(*,*) '-- module_gas_composition_MgII_2796_2804_dust : extracting velocities from ramses '
       allocate(v(3,nleaf))
       call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
       do ileaf = 1,nleaf
          g(ileaf)%v = v(:,ileaf)
       end do
       deallocate(v)

       ! get nHI and temperature from ramses
       write(*,*) '-- module_gas_composition_MgII_2796_2804_dust : extracting nHI and T from ramses '
       allocate(T(nleaf),nhi(nleaf))
       call ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nhi)
       
       ! compute thermal velocity 
       ! ++++++ TURBULENT VELOCITY >>>>> parameter to add and use here
       g(:)%dopwidth = sqrt((2.0d0*kb/mMg)*T) ! [ cm/s ]  

       ! get ndust (pseudo dust density from Laursen, Sommer-Larsen, Andersen 2009)
       write(*,*) '-- module_gas_composition_MgII_2796_2804_dust : extracting ndust from ramses '
       allocate(metallicity(nleaf),nhii(nleaf))
       call ramses_get_metallicity(nleaf,nvar,ramses_var,metallicity)
       call ramses_get_nh_cgs(repository,snapnum,nleaf,nvar,ramses_var,nhii)
       nhii = nhii - nhi
       do ileaf = 1,nleaf
          g(ileaf)%ndust = metallicity(ileaf) / Zref * ( nhi(ileaf) + f_ion*nhii(ileaf) )   ! [ /cm3 ]
       end do
       deallocate(metallicity,T,nhi,nhii)

       g(:)%nMgII = MgII_density(1,:)
    end if

    return


  end subroutine gas_from_ramses_leaves_ions
  ! --------------------------------------------------------------------------
  

  subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var, g)

    ! define gas contents from ramses raw data

    use module_ramses
    
    character(2000),intent(in)        :: repository 
    integer(kind=4),intent(in)        :: snapnum
    integer(kind=4),intent(in)        :: nleaf,nvar
    real(kind=8),intent(in)           :: ramses_var(nvar,nleaf)
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                   :: ileaf
    real(kind=8),allocatable          :: v(:,:), T(:), nMgII(:), nHI(:), nHII(:), metallicity(:)

    ! allocate gas-element array
    allocate(g(nleaf))

    if (gas_overwrite) then
       call overwrite_gas(g)
    else
       
       box_size_cm = ramses_get_box_size_cm(repository,snapnum)

       ! compute velocities in cm / s
       write(*,*) '-- module_gas_composition_MgII_dust : extracting velocities from ramses '
       allocate(v(3,nleaf))
       call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
       do ileaf = 1,nleaf
          g(ileaf)%v = v(:,ileaf)
       end do
       deallocate(v)

       ! get nMgII and temperature from ramses
       write(*,*) '-- module_gas_composition_MgII_dust : extracting nMgII from ramses '
       allocate(T(nleaf),nMgII(nleaf))
       call ramses_get_T_nMgII_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nMgII)
       g(:)%nMgII = nMgII(:)
       ! compute thermal velocity 
       ! ++++++ TURBULENT VELOCITY >>>>> parameter to add and use here
       g(:)%dopwidth = sqrt(2.0d0*kb/mMg*T) ! [ cm/s ]
       deallocate(T,nMgII)

       print*,'min/max of nMgII : ',minval(g(:)%nMgII),maxval(g(:)%nMgII)
       
       ! get ndust (pseudo dust density from Laursen, Sommer-Larsen, Andersen 2009)
       write(*,*) '-- module_gas_composition_MgII_dust : extracting ndust from ramses '
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
    g(:)%nMgII    = fix_nMgII
    g(:)%dopwidth = fix_vth
    g(:)%ndust    = fix_ndust
       
  end subroutine overwrite_gas



  function get_gas_velocity(cell_gas)
    type(gas),intent(in)      :: cell_gas
    real(kind=8),dimension(3) :: get_gas_velocity
    get_gas_velocity(:) = cell_gas%v(:)
    return
  end function get_gas_velocity


  
  function  gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs,iran)

    ! --------------------------------------------------------------------------
    ! Decide whether a scattering event occurs, and if so, on which element
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - cell_gas : MgII (with two absorption channels) and dust
    ! - distance_to_border_cm : the maximum distance the photon may travel (before leaving the cell)
    ! - nu_cell : photon frequency in cell's frame [ Hz ]
    ! - tau_abs : optical depth at which the next scattering event will occur
    ! - iran    : random generator state of the photon 
    ! OUTPUTS:
    ! - distance_to_border_cm : comes out as the distance to scattering event (if there is an event)
    ! - tau_abs : 0 if a scatter occurs, decremented by tau_cell if photon escapes cell. 
    ! - gas_get_scatter_flag : 0 [no scatter], 1 [MgII-2796 scatter], 2 [MgII-2804 scatter], 3 [dust] 
    ! --------------------------------------------------------------------------

    type(gas),intent(in)                  :: cell_gas
    real(kind=8),intent(inout)            :: distance_to_border_cm
    real(kind=8),intent(in)               :: nu_cell
    real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur.
    integer,intent(inout)                 :: iran 
    integer(kind=4)                       :: gas_get_scatter_flag 
    real(kind=8)                          :: tau_MgII_2796, tau_MgII_2804, tau_cell, tau_dust, proba27, proba28, x 

    ! compute optical depths for different components of the gas.
    tau_MgII_2796 = get_tau_MgII_2796(cell_gas%nMgII, cell_gas%dopwidth, distance_to_border_cm, nu_cell)
    tau_MgII_2804 = get_tau_MgII_2804(cell_gas%nMgII, cell_gas%dopwidth, distance_to_border_cm, nu_cell)
    tau_dust      = get_tau_dust(cell_gas%ndust, distance_to_border_cm, nu_cell)
    tau_cell      = tau_MgII_2796 + tau_MgII_2804 + tau_dust

    if (tau_abs > tau_cell) then  ! photon is due for absorption outside the cell 
       gas_get_scatter_flag = 0
       tau_abs = tau_abs - tau_cell
       if (tau_abs.lt.0.0d0) then
          print*, 'tau_abs est negatif'
          stop
       endif
    else  ! the scattering happens inside the cell
       ! decide if it is 2796 or 2804 or dust
       proba27 = tau_MgII_2796 /  tau_cell
       proba28 = proba27 + tau_MgII_2804 /  tau_cell
       x = ran3(iran)
       if (x <= proba27) then
          gas_get_scatter_flag = 1 ! absorption by MgII-2796
       else if (x <= proba28) then
          gas_get_scatter_flag = 2 ! absorption by MgII-2804
       else
          gas_get_scatter_flag = 3 ! absorption by dust
       end if
       ! and transform "distance_to_border_cm" in "distance_to_absorption_cm"
       distance_to_border_cm = distance_to_border_cm * (tau_abs / tau_cell)
    end if

    return

  end function gas_get_scatter_flag


     !--PEEL--
  function  gas_get_tau(cell_gas, distance_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute total opacity of gas accross distance_cm at freq. nu_cell
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - cell_gas : a mix of MgIII and dust
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
    real(kind=8)            :: tau_MgII_2796, tau_MgII_2804, tau_dust

    ! compute optical depths for different components of the gas.
    tau_MgII_2796    = get_tau_MgII_2796(cell_gas%nMgII, cell_gas%dopwidth, distance_cm, nu_cell)
    tau_MgII_2804   = get_tau_MgII_2804(cell_gas%nMgII, cell_gas%dopwidth, distance_cm, nu_cell)
    tau_dust        = get_tau_dust(cell_gas%ndust, distance_cm, nu_cell)
    gas_get_tau = tau_MgII_2796 + tau_MgII_2804 + tau_dust

    return

  end function gas_get_tau
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
       gas_peeloff_weight = MgII_2796_peeloff_weight(cell_gas%v, cell_gas%dopwidth, nu_ext, kin, kout, iran)
    case(2)
       gas_peeloff_weight = MgII_2804_peeloff_weight(cell_gas%v, cell_gas%dopwidth, nu_ext, kin, kout, iran)
    case(3)
       gas_peeloff_weight = dust_peeloff_weight(cell_gas%v, nu_ext, kin, kout)
    case default
       print*,'ERROR in module_gas_composition_MgII_2796_2804_dust.f90:gas_peeloff_weight - unknown case : ',flag 
       stop
    end select

  end function gas_peeloff_weight
  !--LEEP--



  !--CORESKIP-- 
  !subroutine gas_scatter(flag,cell_gas,nu_cell,k,nu_ext,iran)
  subroutine gas_scatter(flag,cell_gas,nu_cell,k,nu_ext,iran,xcrit)
  !--PIKSEROC--

    integer, intent(inout)                    :: flag
    type(gas), intent(in)                     :: cell_gas
    real(kind=8), intent(inout)               :: nu_cell, nu_ext
    real(kind=8), dimension(3), intent(inout) :: k
    integer, intent(inout)                    :: iran
    integer(kind=4)                           :: ilost
    !--CORESKIP--
    real(kind=8),intent(in)                  :: xcrit
    !--PIKSEROC--

    select case(flag)
    case(1)
       call scatter_MgII_2796(cell_gas%v, cell_gas%dopwidth, nu_cell, k, nu_ext, iran)
    case(2)
       call scatter_MgII_2804(cell_gas%v, cell_gas%dopwidth, nu_cell, k, nu_ext, iran)
    case(3)
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
    write(unit) (g(i)%nMgII, i=1,nleaf)
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
       if(no_v) then
          do i=1,n
             g(i)%v(:) = 0d0
          end do
          read(unit)
       else
          read(unit) (g(i)%v(:),i=1,n)
       end if

       read(unit) (g(i)%nMgII,i=1,n)
       read(unit) (g(i)%dopwidth,i=1,n)

       if(no_dust) then
          do i=1,n
             g(i)%ndust = 0d0
          end do
          read(unit)
       else
          read(unit) (g(i)%ndust,i=1,n)
       end if
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
           !Val---
          case ('input_ramses_file')
             write(input_ramses_file,'(a)') trim(value)
          case ('Ion_file')
             write(Ion_file,'(a)') trim(value)
          !--Val
          case ('f_ion')
             read(value,*) f_ion
          case ('Zref')
             read(value,*) Zref
          case ('gas_overwrite')
             read(value,*) gas_overwrite
          case ('fix_nMgII')
             read(value,*) fix_nMgII
          case ('fix_vth')
             read(value,*) fix_vth
          case ('fix_vel')
             read(value,*) fix_vel
          case ('fix_ndust')
             read(value,*) fix_ndust
          case ('no_v')
             read(value,*) no_v
          case ('no_dust')
             read(value,*) no_dust
          case ('fix_box_size_cm')
             read(value,*) fix_box_size_cm
          end select
       end do
    end if
    close(10)

    call read_dust_params(pfile)
    call read_MgII_2796_params(pfile)
    call read_MgII_2804_params(pfile)

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
       write(unit,'(a)')        '# overwrite parameters'
       write(unit,'(a,L1)')     '  gas_overwrite        = ',gas_overwrite
       if(gas_overwrite)then
          write(unit,'(a,ES10.3)') '  fix_nMgII            = ',fix_nMgII
          write(unit,'(a,ES10.3)') '  fix_vth              = ',fix_vth
          write(unit,'(a,ES10.3)') '  fix_vel              = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_ndust            = ',fix_ndust
          write(unit,'(a,ES10.3)') '  fix_box_size_cm      = ',fix_box_size_cm
       endif
       write(unit,'(a)')             ' '
       call print_dust_params(unit)
       write(unit,'(a)')             ' '
       call print_MgII_2796_params(unit)
       write(unit,'(a)')             ' '
       call print_MgII_2804_params(unit)
    else
       write(*,'(a,a,a)') '[gas_composition]'
       write(*,'(a,a)')      '# code compiled with: ',trim(moduleName)
       write(*,'(a)')        '# mixture parameters'
       write(*,'(a,ES10.3)') '  f_ion                = ',f_ion
       write(*,'(a,ES10.3)') '  Zref                 = ',Zref
       write(*,'(a)')        '# overwrite parameters'
       write(*,'(a,L1)')     '  gas_overwrite        = ',gas_overwrite
       if(gas_overwrite)then
          write(*,'(a,ES10.3)') '  fix_nMgII            = ',fix_nMgII
          write(*,'(a,ES10.3)') '  fix_vth              = ',fix_vth
          write(*,'(a,ES10.3)') '  fix_vel              = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_ndust            = ',fix_ndust
          write(*,'(a,ES10.3)') '  fix_box_size_cm      = ',fix_box_size_cm
       endif
       write(*,'(a)')             ' '
       call print_dust_params
       write(*,'(a)')             ' '
       call print_MgII_2796_params
       write(*,'(a)')             ' '
       call print_MgII_2804_params
    end if

    return

  end subroutine print_gas_composition_params



end module module_gas_composition