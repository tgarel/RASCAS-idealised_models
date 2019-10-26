module module_gas_composition

  use module_dust_model
  use module_random
  use module_constants
  use module_lines_model

  implicit none

  private

  character(100),parameter :: moduleName = 'module_gas_composition.f90'
  
  type, public :: gas
     real(kind=8)                :: v(3)      ! gas velocity [cm/s]
     real(kind=8),allocatable    :: n(:)
     real(kind=8)                :: dopwidth  ! Doppler width [cm/s]  sqrt(2*kb*T/amu),  have to divide by sqrt(atomic mass) for each element, for example sqrt(28.085) for Silicone
     real(kind=8)                :: ndust     ! pseudo-numerical density of dust particles [#/cm3]
  end type gas
  real(kind=8),public :: box_size_cm   ! size of simulation box in cm. 
  
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [gas_composition] of the parameter file
  ! --------------------------------------------------------------------------
  ! mixture parameters
  real(kind=8)             :: f_ion           = 0.01   ! ndust = (n_HI + f_ion*n_HII) * Z/Zsun [Laursen+09]
  real(kind=8)             :: Zref            = 0.005  ! reference metallicity. Should be ~ 0.005 for SMC and ~ 0.01 for LMC.
  
  real(kind=8)             :: v_turb          = 0d0    ! Constant turbulent velocity accross the simulation,  in km/s

  logical                  :: no_v                = .false.
  logical                  :: no_dust             = .false.
  logical                  :: mix_HI_metals       = .false.
  logical                  :: file_excited        = .false.
  logical                  :: test_excited_state  = .false.
  real(kind=8)             :: threshold_excited   = 1d-2
  real(kind=8)             :: fraction_excited    = 2d-1
  ! --------------------------------------------------------------------------

  ! public functions:
  public :: gas_from_ramses_leaves_ions, gas_from_ramses_leaves,gas_from_list,gas_get_scatter_flag,gas_scatter,dump_gas,get_gas_velocity
  public :: read_gas,gas_destructor,read_gas_composition_params,print_gas_composition_params
  !--PEEL--
  public :: gas_peeloff_weight,gas_get_tau
  !--LEEP--
  !Val
  public :: gas_get_n_CD, gas_get_CD
  !laV
  public :: ion_number, mix_HI_metals

  !--CORESKIP-- push variable from module_lines_model up so that module_photon knows about it... 
  public :: HI_core_skip 
  !--PIKSEROC-- 

contains

  
  ! --------------------------------------------------------------------------
  subroutine gas_from_ramses_leaves_ions(repository,snapnum,nleaf,nvar,ramses_var,ion_density,g)

    use module_ramses

    character(2000),intent(in)                          :: repository 
    integer(kind=4),intent(in)                          :: snapnum
    integer(kind=4),intent(in)                          :: nleaf,nvar
    real(kind=8),intent(in),dimension(nvar,nleaf)       :: ramses_var
    real(kind=8),intent(in),dimension(ion_number,nleaf) :: ion_density
    type(gas),dimension(:),allocatable,intent(out)      :: g
    integer(kind=4)                                     :: ileaf
    real(kind=8),dimension(:),allocatable               :: T, nhi, metallicity, nhii
    real(kind=8),dimension(:,:),allocatable             :: v

    ! allocate gas-element array
    allocate(g(nleaf))
    do ileaf=1,nleaf
       allocate(g(ileaf)%n(ion_number))
    end do

    box_size_cm = ramses_get_box_size_cm(repository,snapnum)

    ! compute velocities in cm / s
    write(*,*) '-- module_gas_composition : extracting velocities from ramses '
    allocate(v(3,nleaf))
    call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
    do ileaf = 1,nleaf
       g(ileaf)%v = v(:,ileaf)
    end do
    deallocate(v)


    ! get nHI and temperature from ramses
    write(*,*) '-- module_gas_composition : extracting nHI and T from ramses '
    allocate(T(nleaf),nhi(nleaf))
    call ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nhi)

    ! density
    if(mix_HI_metals) then
       do ileaf=1,nleaf
          g(ileaf)%n(1) = nhi(ileaf)
       end do
       do ileaf=1,nleaf
          g(ileaf)%n(2:ion_number+1) = ion_density(:,ileaf)
       end do
    else
       do ileaf=1,nleaf
          g(ileaf)%n(:) = ion_density(:,ileaf)
       end do
    end if

 
    ! compute thermal velocity 
    g(:)%dopwidth = sqrt((2.0d0*kb/amu)*T) ! [ cm/s ]  

    ! get ndust (pseudo dust density from Laursen, Sommer-Larsen, Andersen 2009)
    write(*,*) '-- module_gas_composition : extracting ndust from ramses '
    allocate(metallicity(nleaf),nhii(nleaf))
    call ramses_get_metallicity(nleaf,nvar,ramses_var,metallicity)
    call ramses_get_nh_cgs(repository,snapnum,nleaf,nvar,ramses_var,nhii)
    nhii = nhii - nhi
    do ileaf = 1,nleaf
       g(ileaf)%ndust = metallicity(ileaf) / Zref * ( nhi(ileaf) + f_ion*nhii(ileaf) )   ! [ /cm3 ]
    end do
    deallocate(metallicity,T,nhi,nhii)

    return

  end subroutine gas_from_ramses_leaves_ions
  ! --------------------------------------------------------------------------


   !Val --
  subroutine gas_from_list(ncells, cell_pos, cell_l, gas_leaves)

    integer(kind=4), intent(out)		:: ncells
    integer(kind=4), allocatable, intent(out)	:: cell_l(:)
    real(kind=8), allocatable, intent(out)	:: cell_pos(:,:)
    type(gas), allocatable, intent(out)		:: gas_leaves(:)

    print*,'gas_from_list not available yet in module_gas_composition.f90'
    stop

  end subroutine gas_from_list
  !--laV
 
  

  subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var, g)

    ! define gas contents from ramses raw data

    character(2000),intent(in)        :: repository 
    integer(kind=4),intent(in)        :: snapnum
    integer(kind=4),intent(in)        :: nleaf,nvar
    real(kind=8),intent(in)           :: ramses_var(nvar,nleaf)
    type(gas),dimension(:),allocatable,intent(out) :: g

    print*, 'gas_from_ramses_leaves not available in module_gas_composition.f90'
    stop

    return

  end subroutine gas_from_ramses_leaves
  

  function get_gas_velocity(cell_gas)
    type(gas),intent(in)      :: cell_gas
    real(kind=8),dimension(3) :: get_gas_velocity
    get_gas_velocity(:) = cell_gas%v(:)
    return
  end function get_gas_velocity


  function gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs,iran)

    ! --------------------------------------------------------------------------
    ! Decide whether a scattering event occurs, and if so, on which element
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - cell_gas : a mix of stuff and dust
    ! - distance_to_border_cm : the maximum distance the photon may travel (before leaving the cell)
    ! - nu_cell : photon frequency in cell's frame [ Hz ]
    ! - tau_abs : optical depth at which the next scattering event will occur
    ! - iran    : random generator state of the photon 
    ! OUTPUTS:
    ! - distance_to_border_cm : comes out as the distance to scattering event (if there is an event)
    ! - tau_abs : 0 if a scatter occurs, decremented by tau_cell if photon escapes cell. 
    ! - gas_get_scatter_flag : 0 [no scatter], >0 scatter
    ! --------------------------------------------------------------------------

    ! check whether scattering occurs within cell (scatter_flag > 0) or not (scatter_flag==0)
    type(gas),intent(in)                  :: cell_gas
    real(kind=8),intent(inout)            :: distance_to_border_cm
    real(kind=8),intent(in)               :: nu_cell
    real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur.
    integer(kind=4),intent(inout)         :: iran
    integer(kind=4)                       :: gas_get_scatter_flag, i
    real(kind=8)                          :: tau, tau_ions, tau_dust, tau_cell, tirage, proba1, proba2

    ! compute optical depths for different components of the gas.
    tau_dust = get_tau_dust(cell_gas%ndust, distance_to_border_cm, nu_cell)
    tau_ions = get_tau_lines(cell_gas%n(:), cell_gas%dopwidth, distance_to_border_cm, nu_cell)
    tau_cell = tau_dust + tau_ions
    
    if (tau_abs > tau_cell) then  ! photon is due for absorption outside the cell 
       gas_get_scatter_flag = 0   ! no scatter
       tau_abs = tau_abs - tau_cell
       if (tau_abs.lt.0.0d0) then
          print*, 'tau_abs est negatif'
          stop
       endif
    else  ! the scattering happens inside the cell.
       tirage = ran3(iran)
       if(tirage < tau_dust/tau_cell) then
          gas_get_scatter_flag = 1            !Absorption by dust  (Achtung, not the same convention as normal Rascas)
       else
          tau = tau_dust
          do i=1,n_line_tot
             tau = tau + get_tau_single_line(i, cell_gas%n(ion_index(i)), cell_gas%dopwidth, distance_to_border_cm, nu_cell)
             if(tirage < tau/tau_cell) then
                gas_get_scatter_flag = i+1
                exit
             end if
          end do
       end if

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

    if(flag==1) then
       call scatter_dust(cell_gas%v, nu_cell, k, nu_ext, iran, ilost)
       if(ilost==1) flag=-1
    else
       call scatter_line(flag-1, cell_gas%v, cell_gas%dopwidth, nu_cell, k, nu_ext, iran, xcrit)
    end if

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
    real(kind=8)            :: tau_ions, tau_dust

    ! compute optical depths for different components of the gas.
    tau_ions = get_tau_lines(cell_gas%n(:), cell_gas%dopwidth, distance_cm, nu_cell)
    tau_dust = get_tau_dust(cell_gas%ndust, distance_cm, nu_cell)
    gas_get_tau = tau_ions + tau_dust

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

    if(flag==1) then
       gas_peeloff_weight = dust_peeloff_weight(cell_gas%v, nu_ext, kin, kout)
    else
       gas_peeloff_weight = line_peeloff_weight(flag-1, cell_gas%v, cell_gas%dopwidth, nu_ext, kin, kout, iran)
    end if

  end function gas_peeloff_weight
  !--LEEP--

  !Val
  function gas_get_n_CD()

    integer(kind=4)           :: gas_get_n_CD

    gas_get_n_CD = ion_number + 1

  end function gas_get_n_CD
  !laV

  !Val
  function gas_get_CD(cell_gas, distance_cm)

    real(kind=8), intent(in) :: distance_cm
    type(gas),intent(in)     :: cell_gas
    real(kind=8)             :: gas_get_CD(ion_number+1)

    gas_get_CD(1:ion_number) = distance_cm * cell_gas%n(:)
    gas_get_CD(ion_number+1) = distance_cm * cell_gas%ndust

  end function gas_get_CD
  !laV


  subroutine dump_gas(unit,g)
    type(gas),dimension(:),intent(in) :: g
    integer(kind=4),intent(in)        :: unit
    integer(kind=4)                   :: i,nleaf
    real(kind=8)                      :: mean_T, n_tot
    nleaf = size(g)
    write(unit) (g(i)%v(:), i=1,nleaf)
    write(unit) (g(i)%n(:), i=1,nleaf)
    write(unit) (g(i)%dopwidth, i=1,nleaf)
    write(unit) (g(i)%ndust, i=1,nleaf)
    write(unit) box_size_cm
    mean_T = 0d0
    n_tot = 0d0
    do i=1,nleaf
       n_tot = n_tot + g(i)%n(1)
       mean_T = mean_T + g(i)%dopwidth*g(i)%dopwidth * amu / 2 / kb * g(i)%n(1)
    end do
    mean_T = mean_T / n_tot
    print*, 'mean temperature, weighted by density of first species, = ', mean_T
  end subroutine dump_gas

  
  subroutine read_gas(unit,n,g)
    integer,intent(in)                             :: unit,n
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer                                        :: i, count_test

    allocate(g(1:n))
    do i=1,n
       allocate(g(i)%n(ion_number))
    end do

    !-----------------------Velocity----------------------------
    if(no_v) then
       do i=1,n
          g(i)%v(:) = 0d0
       end do
       read(unit)
    else
       read(unit) (g(i)%v(:),i=1,n)
    end if
    !-----------------------------------------------------------

    !-----------------------Density-----------------------------
    if(test_excited_state) then
       if(ion_number /= 2) then
          print*, 'Problem, for the moment the test_excited state is done for a case with ion_number = 2. Please use test_excited=F or ion_number=2.'
          stop
       end if
       count_test = 0
       read(unit) (g(i)%n(1),i=1,n)
       do i=1,n
          if(g(i)%n(1) > threshold_excited) then
             count_test = count_test + 1
             g(i)%n(2) = fraction_excited*g(i)%n(1)
             g(i)%n(1) = g(i)%n(1) - g(i)%n(2)
          else
             g(i)%n(2) = 1d-18
          end if
       end do
       print*, 'total number of cells, and number of cells with density above threshold : ', n, count_test
    else

       read(unit) (g(i)%n(:),i=1,n)
       
       if(file_excited) then
          
          if(ion_number /= 2) then
             print*, 'Problem, if doing the test with excited state, you must have ion_number = 2, for now. (Doesnt work with OI 1302 - SiII 1304)'
             stop
          end if

          do i=1,n
             g(i)%n(2) = g(i)%n(2) * g(i)%n(1)
             g(i)%n(1) = g(i)%n(1) - g(i)%n(2)
             if(g(i)%n(1) < 0d0) then
                print*, 'Problem, negatif n(1), i = ', i
             end if
          end do

       end if
    end if


    !-----------------------------------------------------------

    !-----------------------Doppler width-----------------------
    read(unit) (g(i)%dopwidth,i=1,n)
    
    ! Application of v_turb to the Doppler width
    if(v_turb > 0d0) then
       if((ion_number > 1) .and. (.not. test_excited_state) .and. (.not. file_excited)) then
          print*, 'Sorry, the v_turb option is not compatible with having several different ions, please complain to valentin.mauerhofer@unige.ch'
          stop
       end if
       do i=1,n
          g(i)%dopwidth = sqrt(g(i)%dopwidth**2 + (1d5*v_turb)**2)
       end do
    end if
    !-----------------------------------------------------------
    
    !-----------------------Dust--------------------------------
    if(no_dust) then
       do i=1,n
          g(i)%ndust = 0d0
       end do
       read(unit)
    else
       read(unit) (g(i)%ndust,i=1,n)
    end if
    !-----------------------------------------------------------

    !-----------------------Box Size----------------------------
    read(unit) box_size_cm
    !-----------------------------------------------------------
    
    print*, 'min/max of vgas(1) : ',minval(g(:)%v(1)),maxval(g(:)%v(1))
    print*, 'min/max of ndust : ',minval(g(:)%ndust),maxval(g(:)%ndust)

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
          case ('v_turb')
             read(value,*) v_turb
          case ('no_v')
             read(value,*) no_v
          case ('no_dust')
             read(value,*) no_dust
           case ('test_excited_state')
              read(value,*) test_excited_state
           case ('threshold_excited')
              read(value,*) threshold_excited
           case ('fraction_excited')
              read(value,*) fraction_excited
           case ('file_excited')
              read(value,*) file_excited
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
       write(unit,'(a,ES10.3)') '  v_turb          = ',v_turb
       write(unit,'(a,L1)')     '  no_v            = ',no_v
       write(unit,'(a,L1)')     '  no_dust         = ',no_dust
       write(unit,'(a,L1)')     '  test_excited_state    = ',test_excited_state
       if(test_excited_state) then
          write(unit,'(a,ES10.3)') '  threshold_excited          = ',threshold_excited
          write(unit,'(a,ES10.3)') '  fraction_excited          = ',fraction_excited
       end if
       write(unit,'(a,L1)')     '  file_excited    = ',file_excited
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
       write(*,'(a,ES10.3)') '  v_turb          = ',v_turb
       write(*,'(a,L1)')     '  no_v            = ',no_v
       write(*,'(a,L1)')     '  no_dust         = ',no_dust
       write(*,'(a,L1)')     ' test_excited_state         = ',test_excited_state
       if(test_excited_state) then
          write(*,'(a,ES10.3)') '  threshold_excited          = ',threshold_excited
          write(*,'(a,ES10.3)') '  fraction_excited          = ',fraction_excited
       end if
       write(*,'(a,L1)')     '  file_excited    = ',file_excited
       write(*,'(a)')             ' '
       call print_line_params
       write(*,'(a)')             ' '
       call print_dust_params
    end if

    return

  end subroutine print_gas_composition_params


  !FOR SEVERAL ELEMENTS
!    ! --------------------------------------------------------------------------
!   subroutine gas_from_ramses_leaves_ions(repository,snapnum,nleaf,nvar,ramses_var,ion_number,ion_density,g)

!     use module_ramses

!     character(2000),intent(in)                          :: repository 
!     integer(kind=4),intent(in)                          :: snapnum, ion_number
!     integer(kind=4),intent(in)                          :: nleaf,nvar
!     real(kind=8),intent(in),dimension(nvar,nleaf)       :: ramses_var
!     real(kind=8),intent(in),dimension(ion_number,nleaf) :: ion_density
!     type(gas),dimension(:),allocatable,intent(out)      :: g
!     integer(kind=4)                                     :: ileaf
!     real(kind=8),dimension(:),allocatable               :: T, nhi!, metallicity, nhii
!     real(kind=8),dimension(:,:),allocatable             :: v


!     if(ion_number /= n_elem) then
!        print*, 'Problem with the number of elements'
!        stop
!     end if
    
!     ! allocate gas-element array
!     allocate(g(nleaf))
!     do ileaf=1,nleaf
! !       g(ileaf)%n_elem = ion_number
!        allocate(g(ileaf)%n(ion_number))
!     end do

        
!        box_size_cm = ramses_get_box_size_cm(repository,snapnum)

!        ! compute velocities in cm / s
!        write(*,*) '-- module_gas_composition : extracting velocities from ramses '
!        allocate(v(3,nleaf))
!        call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
!        do ileaf = 1,nleaf
!           g(ileaf)%v = v(:,ileaf)
!        end do
!        deallocate(v)

!        ! get nHI and temperature from ramses
!        write(*,*) '-- module_gas_composition : extracting T from ramses '
!        allocate(T(nleaf),nhi(nleaf))
!        call ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nhi)
       
!        ! compute thermal velocity 
!        ! ++++++ TURBULENT VELOCITY >>>>> parameter to add and use here
!        g(:)%dopwidth = sqrt((2.0d0*kb)*T) ! [ cm/s ]  

!        ! get ndust (pseudo dust density from Laursen, Sommer-Larsen, Andersen 2009)
!        ! write(*,*) '-- module_gas_composition_SiII_1260_dust : extracting ndust from ramses '
!        ! allocate(metallicity(nleaf),nhii(nleaf))
!        ! call ramses_get_metallicity(nleaf,nvar,ramses_var,metallicity)
!        ! call ramses_get_nh_cgs(repository,snapnum,nleaf,nvar,ramses_var,nhii)
!        ! nhii = nhii - nhi
!        ! do ileaf = 1,nleaf
!        !    g(ileaf)%ndust = metallicity(ileaf) / Zref * ( nhi(ileaf) + f_ion*nhii(ileaf) )   ! [ /cm3 ]
!        ! end do
!        ! deallocate(metallicity,T,nhi,nhii)
!        deallocate(T,nhi)

!        do ileaf=1,nleaf
!           g(ileaf)%n(:) = ion_density(:,ileaf) 
!        end do

!     return


!   end subroutine gas_from_ramses_leaves_ions
!   ! --------------------------------------------------------------------------


end module module_gas_composition
