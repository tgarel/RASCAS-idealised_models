module module_gas_composition

  ! Pure FeII gas. 
  ! This modules handles two transitions in absorption by FeII (2587, 2600)  
  ! and five decay channels (resonant and fluorescent) at 2587, 2600, 2612, 2626, and 2637. 
  ! This corresponds to the multiplet UV1

  use module_FeII_2587_model
  use module_FeII_2600_model
  use module_random
  use module_constants

  implicit none

  private

  character(100),parameter :: moduleName = 'module_gas_composition_FeII_UV1.f90'
  type, public :: gas
     ! fluid
     real(kind=8) :: v(3)      ! gas velocity [cm/s]
     real(kind=8) :: nFeII     ! numerical density of FeII  [#/cm3]
                               ! FeII -> density is computed as abundance_number_FeII * nHI * metallicity / solar_metallicity
     real(kind=8) :: dopwidth  ! Doppler width [cm/s]
  end type gas
  real(kind=8),public :: box_size_cm   ! size of simulation box in cm. 

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [gas_composition] of the parameter file
  ! --------------------------------------------------------------------------
  ! possibility to overwrite ramses values with an ad-hoc model 
  logical                  :: gas_overwrite       = .false. ! if true, define cell values from following parameters 
  real(kind=8)             :: fix_nFeII           = 0.0d0   ! ad-hoc HI density (H/cm3)
  real(kind=8)             :: fix_vth             = 1.0d5   ! ad-hoc thermal velocity (cm/s)
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

    character(2000),intent(in)        :: repository 
    integer(kind=4),intent(in)        :: snapnum
    integer(kind=4),intent(in)        :: nleaf,nvar
    real(kind=8),intent(in)           :: ramses_var(nvar,nleaf)
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                   :: ileaf
    real(kind=8),allocatable          :: v(:,:), T(:), nFeII(:)

    ! allocate gas-element array
    allocate(g(nleaf))

    if (gas_overwrite) then
       call overwrite_gas(g)
    else
       
       box_size_cm = ramses_get_box_size_cm(repository,snapnum)

       ! compute velocities in cm / s
       write(*,*) '-- module_gas_composition_FeII_UV1 : extracting velocities from ramses '
       allocate(v(3,nleaf))
       call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
       do ileaf = 1,nleaf
          g(ileaf)%v = v(:,ileaf)
       end do
       deallocate(v)

       ! get nFeII and temperature from ramses
       write(*,*) '-- module_gas_composition_FeII_UV1 : extracting nFeII from ramses '
       allocate(T(nleaf),nFeII(nleaf))
       call ramses_get_T_nFeII_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nFeII)
       g(:)%nFeII = nFeII(:)
       ! compute thermal velocity 
       ! ++++++ TURBULENT VELOCITY >>>>> parameter to add and use here
       g(:)%dopwidth = sqrt(2.0d0*kb/mFe*T) ! [ cm/s ]
       deallocate(T,nFeII)

       print*,'min/max of nFeII : ',minval(g(:)%nFeII),maxval(g(:)%nFeII)
       
    end if

    return

  end subroutine gas_from_ramses_leaves
  

  subroutine overwrite_gas(g)

    type(gas),dimension(:),intent(inout) :: g

    box_size_cm   = fix_box_size_cm
    
    g(:)%v(1)     = fix_vel
    g(:)%v(2)     = fix_vel
    g(:)%v(3)     = fix_vel
    g(:)%nFeII    = fix_nFeII
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
    ! - cell_gas : FeII (with two absorption channels)
    ! - distance_to_border_cm : the maximum distance the photon may travel (before leaving the cell)
    ! - nu_cell : photon frequency in cell's frame [ Hz ]
    ! - tau_abs : optical depth at which the next scattering event will occur
    ! - iran    : random generator state of the photon 
    ! OUTPUTS:
    ! - distance_to_border_cm : comes out as the distance to scattering event (if there is an event)
    ! - tau_abs : 0 if a scatter occurs, decremented by tau_cell if photon escapes cell. 
    ! - gas_get_scatter_flag : 0 [no scatter], 1 [FeII-2587 scatter], 2 [FeII-2600 scatter]
    ! --------------------------------------------------------------------------

    type(gas),intent(in)                  :: cell_gas
    real(kind=8),intent(inout)            :: distance_to_border_cm
    real(kind=8),intent(in)               :: nu_cell
    real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur.
    integer,intent(inout)                 :: iran 
    integer(kind=4)                       :: gas_get_scatter_flag 
    real(kind=8)                          :: tau_FeII_2587, tau_FeII_2600 , tau_cell, proba87, x
    
    ! compute optical depths for different components of the gas.
    tau_FeII_2587 = get_tau_FeII_2587(cell_gas%nFeII, cell_gas%dopwidth, distance_to_border_cm, nu_cell)
    tau_FeII_2600 = get_tau_FeII_2600(cell_gas%nFeII, cell_gas%dopwidth, distance_to_border_cm, nu_cell)
    tau_cell      = tau_FeII_2587 + tau_FeII_2600

    if (tau_abs > tau_cell) then  ! photon is due for absorption outside the cell 
       gas_get_scatter_flag = 0
       tau_abs = tau_abs - tau_cell
       if (tau_abs.lt.0.d0) then
          print*, 'tau_abs est negatif'
          stop
       endif
    else  ! the scattering happens inside the cell
       ! decide if it is 2587 or 2600
       proba87 = tau_FeII_2587 / tau_cell
       x = ran3(iran)
       if (x <= proba87) then
          gas_get_scatter_flag = 1 ! absorption by FeII-2587
       else
          gas_get_scatter_flag = 2 ! absorption by FeII-2600
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

    select case(flag)
    case(1)
       call scatter_FeII_2587(cell_gas%v, cell_gas%dopwidth, nu_cell, k, nu_ext, iran)
    case(2)
       call scatter_FeII_2600(cell_gas%v, cell_gas%dopwidth, nu_cell, k, nu_ext, iran)
    end select

  end subroutine gas_scatter



  subroutine dump_gas(unit,g)
    type(gas),dimension(:),intent(in) :: g
    integer,intent(in)                :: unit
    integer                           :: i,nleaf
    nleaf = size(g)
    write(unit) (g(i)%v(:), i=1,nleaf)
    write(unit) (g(i)%nFeII, i=1,nleaf)
    write(unit) (g(i)%dopwidth, i=1,nleaf)
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
       read(unit) (g(i)%nFeII,i=1,n)
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
          case ('gas_overwrite')
             read(value,*) gas_overwrite
          case ('fix_nFeII')
             read(value,*) fix_nFeII
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

    call read_FeII_2587_params(pfile)
    call read_FeII_2600_params(pfile)
    
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
       write(unit,'(a)')        '# overwrite parameters'
       write(unit,'(a,L1)')     '  gas_overwrite        = ',gas_overwrite
       if(gas_overwrite)then
          write(unit,'(a,ES10.3)') '  fix_nFeII            = ',fix_nFeII
          write(unit,'(a,ES10.3)') '  fix_vth              = ',fix_vth
          write(unit,'(a,ES10.3)') '  fix_vel              = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_box_size_cm      = ',fix_box_size_cm
       endif
       write(unit,'(a)')             ' '
       call print_FeII_2587_params(unit)
       write(unit,'(a)')             ' '
       call print_FeII_2600_params(unit)
    else
       write(*,'(a,a,a)') '[gas_composition]'
       write(*,'(a,a)')      '# code compiled with: ',trim(moduleName)
       write(*,'(a)')        '# overwrite parameters'
       write(*,'(a,L1)')     '  gas_overwrite        = ',gas_overwrite
       if(gas_overwrite)then
          write(*,'(a,ES10.3)') '  fix_nFeII            = ',fix_nFeII
          write(*,'(a,ES10.3)') '  fix_vth              = ',fix_vth
          write(*,'(a,ES10.3)') '  fix_vel              = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_box_size_cm      = ',fix_box_size_cm
       endif
       write(*,'(a)')             ' '
       call print_FeII_2587_params
       write(*,'(a)')             ' '
       call print_FeII_2600_params
    end if

    return

  end subroutine print_gas_composition_params



end module module_gas_composition
