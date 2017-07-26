module module_gas_composition

  ! Pure FeII gas. 
  ! This modules handles two transitions in absorption by FeII (2587, 2600)  
  ! and five decay channels (resonant and fluorescent) at 2587, 2600, 2612, 2626, and 2637. 
  ! This corresponds to the multiplet UV1

  use module_FeII_2587_model
  use module_FeII_2600_model
  use module_random
  use module_ramses
  use module_constants
  use module_idealised_models

  implicit none

  private

  type, public :: gas
     ! fluid
     real(kind=8) :: v(3)      ! gas velocity [cm/s]
     ! FeII
     ! -> density is computed as 2.8d-5 * nHI * metallicity / solar_metallicity
     real(kind=8) :: nFeII     ! numerical density of FeII  [#/cm3]
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
  ! possibility to overwrite ramses values with an ad-hoc model 
  logical                  :: gas_overwrite       = .false. ! if true, define cell values from following parameters 
  real(kind=8)             :: fix_nFeII           = 0.0d0   ! ad-hoc HI density (H/cm3)
  real(kind=8)             :: fix_vth             = 1.0d5   ! ad-hoc thermal velocity (cm/s)
  real(kind=8)             :: fix_vel             = 0.0d0   ! ad-hoc cell velocity (cm/s) -> NEED BETTER PARAMETERIZATION for more than static... 
  real(kind=8)             :: fix_box_size_cm     = 1.0d8   ! ad-hoc box size in cm.
  ! miscelaneous
  logical                  :: verbose             = .true. ! display some run-time info on this module
  ! --------------------------------------------------------------------------

  ! public functions:
  public :: gas_from_ramses_leaves,get_gas_velocity,gas_get_scatter_flag,gas_scatter,dump_gas
  public :: read_gas,gas_destructor,read_gas_composition_params,print_gas_composition_params

contains
  

  subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var, g, x_leaf, leaf_level)

    ! define gas contents from ramses raw data

    character(2000),intent(in)        :: repository 
    integer(kind=4),intent(in)        :: snapnum
    integer(kind=4),intent(in)        :: nleaf,nvar
    real(kind=8),intent(in)           :: ramses_var(nvar,nleaf)
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                   :: ileaf
    real(kind=8),allocatable          :: v(:,:), T(:), nFeII(:), nHI(:), nHII(:), metallicity(:)

    character(2000)                                :: file,datadir_path
    real(kind=8),intent(in),dimension(nleaf,3)     :: x_leaf
    integer(kind=4),intent(in),dimension(nleaf)    :: leaf_level
    character(200)                                 :: modelprops_file
    
    ! allocate gas-element array
    allocate(g(nleaf))

    if (gas_overwrite) then
       call overwrite_gas(g,x_leaf,leaf_level,nleaf)
       ! Dump idealised model (nh, temp, pos...)
       call read_datadir(datadir_path)
       modelprops_file = 'modelprops_file'
       file = trim(datadir_path)//trim(modelprops_file)
       open(unit=15, file=trim(file), status='unknown', form='unformatted', action='write')
       write(15) nleaf
       do ileaf = 1,nleaf
          write(15) x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),g(ileaf)%v(1),g(ileaf)%v(2),g(ileaf)%v(3),g(ileaf)%dopwidth,g(ileaf)%nFeII
       end do
       close(15)
    else
       
       box_size_cm = ramses_get_box_size_cm(repository,snapnum)

       ! compute velocities in cm / s
       if (verbose) write(*,*) '-- module_gas_composition_FeII_UV1 : extracting velocities from ramses '
       allocate(v(3,nleaf))
       call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
       do ileaf = 1,nleaf
          g(ileaf)%v = v(:,ileaf)
       end do
       deallocate(v)

       ! get nFeII and temperature from ramses
       if (verbose) write(*,*) '-- module_gas_composition_FeII_UV1 : extracting nFeII from ramses '
       allocate(T(nleaf),nFeII(nleaf))
       call ramses_get_T_nFeII_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nFeII)
       g(:)%nFeII = nFeII(:)
       ! compute thermal velocity 
       ! ++++++ TURBULENT VELOCITY >>>>> parameter to add and use here
       g(:)%dopwidth = sqrt(2.0d0*kb/mFe*T) ! [ cm/s ]
       deallocate(T,nFeII)

       if (verbose) print*,'min/max of nFeII : ',minval(g(:)%nFeII),maxval(g(:)%nFeII)
       
    end if

    return

  end subroutine gas_from_ramses_leaves
  
  
  ! move to idealised_models ?
  subroutine overwrite_gas(g,x_leaf,leaf_level,nleaf)
    
    type(gas),dimension(:),intent(inout) :: g
    real(kind=8),intent(in)              :: x_leaf(nleaf,3)
    integer(kind=4),intent(in)           :: leaf_level(nleaf)
    integer(kind=4)                      :: ileaf
    integer(kind=4),intent(in)           :: nleaf
    real(kind=8)                         :: dx_cell
    character(50)                        :: overwrite_model
    
    call read_overwrite_params(overwrite_model)
    
    box_size_cm   = fix_box_size_cm
    
    select case (overwrite_model)
    case('uniform')
       g(:)%v(1)     = fix_vel
       g(:)%v(2)     = fix_vel
       g(:)%v(3)     = fix_vel
       g(:)%nFeII    = fix_nFeII
       g(:)%dopwidth = fix_vth     
       
#ifdef DEBUG
       print*,'in overwrite_gas: allocated g?',shape(g)
       print*,'in overwrite_gas: ',minval(g%nFeII),maxval(g%nFeII)
       print*,'in overwrite_gas: ',minval(g%dopwidth),maxval(g%dopwidth)
       print*,'in overwrite_gas: ',minval(g%v),maxval(g%v)
       print*,'in overwrite_gas: ',box_size_cm
#endif
    case('sphere_homogen_velfix')      
       do ileaf=1,nleaf
          dx_cell = 0.5d0**leaf_level(ileaf)
          call sphere_homogen_velfix(g(ileaf)%v(1),g(ileaf)%v(2),g(ileaf)%v(3),g(ileaf)%nFeII,g(ileaf)%dopwidth,x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do
       
    case('shell_homogen_velfix')      
       do ileaf=1,nleaf
          dx_cell = 0.5d0**leaf_level(ileaf)
          call shell_homogen_velfix(g(ileaf)%v(1),g(ileaf)%v(2),g(ileaf)%v(3),g(ileaf)%nFeII,g(ileaf)%dopwidth,x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do

    case('sphere_homog_velgrad')      
       do ileaf=1,nleaf
          dx_cell = 0.5d0**leaf_level(ileaf)
          call sphere_homogen_velgrad(g(ileaf)%v(1),g(ileaf)%v(2),g(ileaf)%v(3),g(ileaf)%nFeII,g(ileaf)%dopwidth,x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do

    case('sphere_homogen_steidel')      
       do ileaf=1,nleaf
          dx_cell = 0.5d0**leaf_level(ileaf)
          call sphere_homogen_steidel(g(ileaf)%v(1),g(ileaf)%v(2),g(ileaf)%v(3),g(ileaf)%nFeII,g(ileaf)%dopwidth,x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do
       
    case('sphere_homogen_velgrad_ct_outflow_rate')      
       do ileaf=1,nleaf
          dx_cell = 0.5d0**leaf_level(ileaf)
          call sphere_homogen_velgrad_ct_outflow_rate(g(ileaf)%v(1),g(ileaf)%v(2),g(ileaf)%v(3),g(ileaf)%nFeII,g(ileaf)%dopwidth,x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do

    case('sphere_homogen_velgrad_rad_pressure')      
       do ileaf=1,nleaf
          dx_cell = 0.5d0**leaf_level(ileaf)
          call sphere_homogen_velgrad_rad_pressure(g(ileaf)%v(1),g(ileaf)%v(2),g(ileaf)%v(3),g(ileaf)%nFeII,g(ileaf)%dopwidth,x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do
       
    case('sphere_densgrad_velgrad')      
       do ileaf=1,nleaf
          dx_cell = 0.5d0**leaf_level(ileaf)
          call sphere_densgrad_velgrad(g(ileaf)%v(1),g(ileaf)%v(2),g(ileaf)%v(3),g(ileaf)%nFeII,g(ileaf)%dopwidth,x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do

    case('sphere_densgrad_velfix')      
       do ileaf=1,nleaf
          dx_cell = 0.5d0**leaf_level(ileaf)
          call sphere_densgrad_velfix(g(ileaf)%v(1),g(ileaf)%v(2),g(ileaf)%v(3),g(ileaf)%nFeII,g(ileaf)%dopwidth,x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do

    case('sphere_scarlata15')      
       do ileaf=1,nleaf
          dx_cell = 0.5d0**leaf_level(ileaf)
          call sphere_scarlata15(g(ileaf)%v(1),g(ileaf)%v(2),g(ileaf)%v(3),g(ileaf)%nFeII,g(ileaf)%dopwidth,x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do

    case('sphere_prochaska11')      
       do ileaf=1,nleaf
          dx_cell = 0.5d0**leaf_level(ileaf)
          call sphere_prochaska11(g(ileaf)%v(1),g(ileaf)%v(2),g(ileaf)%v(3),g(ileaf)%nFeII,g(ileaf)%dopwidth,x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell,g(ileaf)%ndust)
       end do
       
    case('disc_thin')      
       do ileaf=1,nleaf
          dx_cell = 0.5d0**leaf_level(ileaf)
          call disc_thin(g(ileaf)%v(1),g(ileaf)%v(2),g(ileaf)%v(3),g(ileaf)%nFeII,g(ileaf)%dopwidth,x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do
    case('disc_thick')      
       do ileaf=1,nleaf
          dx_cell = 0.5d0**leaf_level(ileaf)
          call disc_thick(g(ileaf)%v(1),g(ileaf)%v(2),g(ileaf)%v(3),g(ileaf)%nFeII,g(ileaf)%dopwidth,x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do
       
    end select
      
  end subroutine overwrite_gas
  
!!$  subroutine overwrite_gas(g)
!!$
!!$    type(gas),dimension(:),intent(inout) :: g
!!$
!!$    box_size_cm   = fix_box_size_cm
!!$    
!!$    g(:)%v(1)     = fix_vel
!!$    g(:)%v(2)     = fix_vel
!!$    g(:)%v(3)     = fix_vel
!!$    g(:)%nFeII    = fix_nFeII
!!$    g(:)%dopwidth = fix_vth
!!$
!!$  end subroutine overwrite_gas



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
    ! - cell_gas : FeII (with two absorption channels) and dust
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
    integer(kind=4)                           :: ilost 

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
    integer(kind=4),intent(in)                     :: unit,n
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                                :: i
    
    allocate(g(1:n))
    
    read(unit) (g(i)%v(:),i=1,n)
    read(unit) (g(i)%nFeII,i=1,n)
    read(unit) (g(i)%dopwidth,i=1,n)
    read(unit) box_size_cm 
    
    if (verbose) print*,'min/max of nFeII : ',minval(g(:)%nFeII),maxval(g(:)%nFeII)
    
  end subroutine read_gas
  
  
!!$  subroutine read_gas(unit,n,g)
!!$    integer,intent(in)                             :: unit,n
!!$    type(gas),dimension(:),allocatable,intent(out) :: g
!!$    integer                                        :: i
!!$    allocate(g(1:n))
!!$    if (gas_overwrite) then
!!$       call overwrite_gas(g)
!!$    else
!!$       read(unit) (g(i)%v(:),i=1,n)
!!$       read(unit) (g(i)%nFeII,i=1,n)
!!$       read(unit) (g(i)%dopwidth,i=1,n)
!!$       read(unit) box_size_cm 
!!$    end if
!!$
!!$    if (verbose) print*,'min/max of nFeII : ',minval(g(:)%nFeII),maxval(g(:)%nFeII)
!!$
!!$  end subroutine read_gas

  

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
          case ('verbose')
             read(value,*) verbose
          case ('fix_box_size_cm')
             read(value,*) fix_box_size_cm
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
       write(unit,'(a,a,a)') '[gas_composition]'
       write(unit,'(a)')       '# overwrite parameters'
       write(unit,'(a,L1)')    '  gas_overwrite         = ',gas_overwrite
       write(unit,'(a,ES10.3)') '  fix_nFeII            = ',fix_nFeII
       write(unit,'(a,ES10.3)') '  fix_vth              = ',fix_vth
       write(unit,'(a,ES10.3)') '  fix_vel              = ',fix_vel
       write(unit,'(a,ES10.3)') '  fix_box_size_cm      = ',fix_box_size_cm
       write(unit,'(a)')       '# miscelaneous parameters'
       write(unit,'(a,L1)')    '  verbose               = ',verbose
       write(unit,'(a)')             ' '
       call print_ramses_params(unit)
    else
       write(*,'(a,a,a)') '[gas_composition]'
       write(*,'(a)')       '# overwrite parameters'
       write(*,'(a,L1)')    '  gas_overwrite         = ',gas_overwrite
       write(*,'(a,ES10.3)') '  fix_nFeII            = ',fix_nFeII
       write(*,'(a,ES10.3)') '  fix_vth              = ',fix_vth
       write(*,'(a,ES10.3)') '  fix_vel              = ',fix_vel
       write(*,'(a,ES10.3)') '  fix_box_size_cm      = ',fix_box_size_cm
       write(*,'(a)')       '# miscelaneous parameters'
       write(*,'(a,L1)')    '  verbose               = ',verbose
       write(*,'(a)')             ' '
       call print_ramses_params
    end if

    return

  end subroutine print_gas_composition_params


  subroutine read_datadir(DataDir)
    
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present
    character(2000) :: pfile,DataDir
    
    call get_command_argument(1, pfile)
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    section_present = .false.
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       !if (line(1:70) == '[RASCAS]') then
       if (line(1:70) == '[RASCAS-serial]') then
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
          case ('DataDir')
             write(DataDir,'(a)') trim(value)
          end select
       end do
    end if
    close(10)
    
    return
    
  end subroutine read_datadir
  

end module module_gas_composition
