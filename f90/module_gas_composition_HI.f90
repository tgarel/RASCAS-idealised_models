module module_gas_composition

  ! Pure HI gas. 
  ! - The HI content is derived from RAMSES directly

  use module_HI_1216_model
  use module_random
  use module_constants

  implicit none

  private

  character(100),parameter :: moduleName = 'module_gas_composition_HI.f90'
  
  type, public :: gas
     ! fluid
     real(kind=8) :: v(3)      ! gas velocity [cm/s]
     ! Hydrogen 
     real(kind=8) :: nHI       ! HI numerical density [HI/cm3]
     real(kind=8) :: dopwidth  ! Doppler width [cm/s]
  end type gas
  real(kind=8),public :: box_size_cm   ! size of simulation box in cm. 

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [gas_composition] of the parameter file
  ! --------------------------------------------------------------------------
  ! possibility to overwrite ramses values with an ad-hoc model
  character(2000)          :: input_ramses_file   = './ramses'  !Path to the file 'input_ramses', with all the info on the simulation (ncells, box_size, positions, velocities, nHI, nHII, Z, dopwidth
  logical                  :: gas_overwrite       = .false. ! if true, define cell values from following parameters 
  real(kind=8)             :: fix_nhi             = 0.0d0   ! ad-hoc HI density (H/cm3)
  real(kind=8)             :: fix_vth             = 1.0d5   ! ad-hoc thermal velocity (cm/s)
  real(kind=8)             :: fix_vel             = 0.0d0   ! ad-hoc cell velocity (cm/s) -> NEED BETTER PARAMETERIZATION for more than static... 
  real(kind=8)             :: fix_box_size_cm     = 1.0d8   ! ad-hoc box size in cm. 
  ! --------------------------------------------------------------------------

  ! public functions:
  public :: gas_from_ramses_leaves,gas_from_list,get_gas_velocity,gas_get_scatter_flag,gas_scatter,dump_gas
  public :: read_gas,gas_destructor,read_gas_composition_params,print_gas_composition_params
  !Val
  public :: gas_get_CD, gas_get_n_CD
  !laV
  !--PEEL--
  public :: gas_peeloff_weight,gas_get_tau
  !--LEEP--

  !--CORESKIP-- push variable from module_HI_model up so that module_photon knows about it...
  public :: HI_core_skip
  !--PIKSEROC--    

contains


  !Val --
  subroutine gas_from_list(ncells, cell_pos, cell_l, gas_leaves)

    integer(kind=4), intent(out)		:: ncells
    integer(kind=4), allocatable, intent(out)	:: cell_l(:)
    real(kind=8), allocatable, intent(out)	:: cell_pos(:,:)
    type(gas), allocatable, intent(out)		:: gas_leaves(:)

    integer(kind=4)				:: i


    open(unit=20, file=input_ramses_file, status='old', form='unformatted')
    read(20) ncells

    allocate(cell_l(ncells), cell_pos(ncells,3), gas_leaves(ncells))

    read(20) box_size_cm
    read(20) cell_l
    read(20) cell_pos(:,1)
    read(20) cell_pos(:,2)
    read(20) cell_pos(:,3)
    read(20) gas_leaves%v(1)
    read(20) gas_leaves%v(2)
    read(20) gas_leaves%v(3)
    read(20) gas_leaves%nHI
    read(20)
    read(20)
    read(20) gas_leaves%dopwidth
    
    close(20)

    print*,'boxsize in cm : ', box_size_cm
    print*,'min/max of vth   : ',minval(gas_leaves(:)%dopwidth),maxval(gas_leaves(:)%dopwidth)
    print*,'min/max of nHI : ',minval(gas_leaves(:)%nHI),maxval(gas_leaves(:)%nHI)


  end subroutine gas_from_list
  ! --------------------------------------------------------------------------
  !--laV
  

  subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var, g)

    ! define gas contents from ramses raw data

    use module_ramses

    character(2000),intent(in)                     :: repository 
    integer(kind=4),intent(in)                     :: snapnum
    integer(kind=4),intent(in)                     :: nleaf,nvar
    real(kind=8),intent(in),dimension(nvar,nleaf)  :: ramses_var
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                                :: ileaf
    real(kind=8),dimension(:),allocatable          :: T, nhi
    real(kind=8),dimension(:,:),allocatable        :: v

    ! allocate gas-element array
    allocate(g(nleaf))

    if (gas_overwrite) then
       call overwrite_gas(g)
    else
       box_size_cm = ramses_get_box_size_cm(repository,snapnum)
       ! compute velocities in cm / s
       write(*,*) '-- module_gas_composition_HI : extracting velocities from ramses '
       allocate(v(3,nleaf))
       call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
       do ileaf = 1,nleaf
          g(ileaf)%v = v(:,ileaf)
       end do
       deallocate(v)
       ! get nHI and temperature from ramses
       write(*,*) '-- module_gas_composition_HI : extracting nHI and T from ramses '
       allocate(T(nleaf),nhi(nleaf))
       call ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nhi)
       g(:)%nHI = nhi(:)
       ! compute thermal velocity 
       ! ++++++ TURBULENT VELOCITY >>>>> parameter to add and use here
       g(:)%dopwidth = sqrt((2.0d0*kb/mp)*T) ! [ cm/s ]
       deallocate(T,nhi)
    end if

    return

  end subroutine gas_from_ramses_leaves



  !--PEEL--
  function  gas_get_tau(cell_gas, distance_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute total opacity of gas accross distance_cm at freq. nu_cell
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - cell_gas : H alone
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
    real(kind=8)            :: tau_HI

    ! compute optical depths for different components of the gas.
    tau_HI   = get_tau_HI_1216(cell_gas%nHI, cell_gas%dopwidth, distance_cm, nu_cell)
    gas_get_tau = tau_HI

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
    case(1)
       gas_peeloff_weight = HI_1216_peeloff_weight(cell_gas%v, cell_gas%dopwidth, nu_ext, kin, kout, iran)
    case default
       print*,'ERROR in module_gas_composition_HI.f90:gas_peeloff_weight - unknown case : ',flag 
       stop
    end select

  end function gas_peeloff_weight
  !--LEEP--
  
  

  subroutine overwrite_gas(g)
    ! overwrite ramses values with an ad-hoc model

    type(gas),dimension(:),intent(inout) :: g

    box_size_cm   = fix_box_size_cm
    
    g(:)%v(1)     = fix_vel
    g(:)%v(2)     = fix_vel
    g(:)%v(3)     = fix_vel
    g(:)%nHI      = fix_nhi
    g(:)%dopwidth = fix_vth
    
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
    ! - cell_gas : pure HI gas
    ! - distance_to_border_cm : the maximum distance the photon may travel (before leaving the cell)
    ! - nu_cell : photon frequency in cell's frame [ Hz ]
    ! - tau_abs : optical depth at which the next scattering event will occur
    ! - iran    : random generator state of the photon 
    ! OUTPUTS:
    ! - distance_to_border_cm : comes out as the distance to scattering event (if there is an event)
    ! - tau_abs : 0 if a scatter occurs, decremented by tau_cell if photon escapes cell. 
    ! - gas_get_scatter_flag : 0 [no scatter], 1 [H scatter]
    ! --------------------------------------------------------------------------

    type(gas),intent(in)                  :: cell_gas
    real(kind=8),intent(inout)            :: distance_to_border_cm
    real(kind=8),intent(in)               :: nu_cell
    real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur.
    integer(kind=4),intent(inout)         :: iran 
    integer(kind=4)                       :: gas_get_scatter_flag 
    real(kind=8)                          :: tau_HI, tau_cell

    ! compute optical depths for different components of the gas.
    tau_HI   = get_tau_HI_1216(cell_gas%nHI, cell_gas%dopwidth, distance_to_border_cm, nu_cell)
    tau_cell = tau_HI

    if (tau_abs > tau_cell) then  ! photon is due for absorption outside the cell 
       gas_get_scatter_flag = 0
       tau_abs = tau_abs - tau_cell
       if (tau_abs.lt.0.0d0) then
          print*, 'tau_abs est negatif'
          stop
       endif
    else  ! the scattering happens inside the cell 
       gas_get_scatter_flag = 1
       ! and transform "distance_to_border_cm" in "distance_to_absorption_cm"
       distance_to_border_cm = distance_to_border_cm * (tau_abs / tau_cell)
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


    select case(flag)
    case(1)
       !--CORESKIP--
       call scatter_HI_1216(cell_gas%v, cell_gas%dopwidth, nu_cell, k, nu_ext, iran,xcrit)
       !call scatter_HI(cell_gas%v, cell_gas%dopwidth, nu_cell, k, nu_ext, iran)
       !--PIKSEROC--

    end select
    
  end subroutine gas_scatter

    !Val
  function gas_get_n_CD()

    integer(kind=4)           :: gas_get_n_CD

    gas_get_n_CD = 1

  end function gas_get_n_CD
  !laV

  !Val
  function gas_get_CD(cell_gas, distance_cm)

    real(kind=8), intent(in) :: distance_cm
    type(gas),intent(in)     :: cell_gas
    real(kind=8)             :: gas_get_CD(1)

    gas_get_CD = (/ distance_cm*cell_gas%nHI /)

  end function gas_get_CD
  !laV



  subroutine dump_gas(unit,g)
    type(gas),dimension(:),intent(in) :: g
    integer(kind=4),intent(in)        :: unit
    integer(kind=4)                   :: i,nleaf
    nleaf = size(g)
    write(unit) (g(i)%v(:), i=1,nleaf)
    write(unit) (g(i)%nHI, i=1,nleaf)
    write(unit) (g(i)%dopwidth, i=1,nleaf)
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
          case ('input_ramses_file')
             write(input_ramses_file,'(a)') trim(value)
          case ('gas_overwrite')
             read(value,*) gas_overwrite
          case ('fix_nhi')
             read(value,*) fix_nhi
          case ('fix_vth')
             read(value,*) fix_vth
          case ('fix_vel')
             read(value,*) fix_vel
          case ('fix_box_size_cm')
             read(value,*) fix_box_size_cm
          end select
       end do
    end if
    close(10)

    call read_HI_1216_params(pfile)

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
       write(unit,'(a,a)')     '# code compiled with: ',trim(moduleName)
       write(unit,'(a)')       '# overwrite parameters'
       write(unit,'(a,L1)')    '  gas_overwrite       = ',gas_overwrite
       if(gas_overwrite)then
          write(unit,'(a,ES10.3)') '  fix_nhi            = ',fix_nhi
          write(unit,'(a,ES10.3)') '  fix_vth            = ',fix_vth
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_box_size_cm    = ',fix_box_size_cm
       endif
       write(unit,'(a)')             ' '
       call print_HI_1216_params(unit)
    else
       write(*,'(a,a,a)') '[gas_composition]'
       write(*,'(a,a)')     '# code compiled with: ',trim(moduleName)
       write(*,'(a)')       '# overwrite parameters'
       write(*,'(a,L1)')    '  gas_overwrite       = ',gas_overwrite
       if(gas_overwrite)then
          write(*,'(a,ES10.3)') '  fix_nhi            = ',fix_nhi
          write(*,'(a,ES10.3)') '  fix_vth            = ',fix_vth
          write(*,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_box_size_cm    = ',fix_box_size_cm
       endif
       write(*,'(a)')             ' '
       call print_HI_1216_params
    end if

    return

  end subroutine print_gas_composition_params



end module module_gas_composition
