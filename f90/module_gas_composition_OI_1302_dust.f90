module module_gas_composition

  ! Mix of OI and dust 
  ! This modules handles one transition in absorption at 1302 A and dust

  use module_OI_1302_model
  use module_dust_model
  use module_random
  use module_ramses
  use module_constants
  use module_spectra

  implicit none

  private

  type, public :: gas
     ! fluid
     real(kind=8) :: v(3)      ! gas velocity [cm/s]
     ! OI
     ! -> density is computed with Krome
     !s
     real(kind=8) :: nOI       ! numerical density of OI  [#/cm3]
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
  character(2000)          :: input_ramses_file             !Path to the file 'input_ramses', with all the info on the simulation (ncells, box_size, positions, velocities, nHI, nHII, Z, dopwidth
  character(2000)          :: Ion_file                      !Path to the file with the OI densities
  real(kind=8)             :: f_ion           = 0.01        ! ndust = (n_HI + f_ion*n_HII) * Z/Zsun [Laursen+09]
  real(kind=8)             :: Zref            = 0.005       ! reference metallicity. Should be ~ 0.005 for SMC and ~ 0.01 for LMC. 
  ! possibility to overwrite ramses values with an ad-hoc model 
  logical                  :: gas_overwrite       = .false. ! if true, define cell values from following parameters
  real(kind=8)             :: fix_nOI             = 1.0d0   ! ad-hoc SiII density (1/cm3)
  real(kind=8)             :: fix_vth             = 1.0d5   ! ad-hoc thermal velocity (cm/s)
  real(kind=8)             :: fix_vel             = 0.0d0   ! ad-hoc cell velocity (cm/s) -> NEED BETTER PARAMETERIZATION for more than static... 
  real(kind=8)             :: fix_ndust           = 0.0d0
  real(kind=8)             :: fix_box_size_cm     = 1.0d8   ! ad-hoc box size in cm.
  ! miscelaneous
  logical                  :: verbose             = .true. ! display some run-time info on this module
  logical                  :: HI_core_skip        = .false.
  ! --------------------------------------------------------------------------

  ! public functions:
  public :: gas_from_ramses_leaves, gas_from_ramses_leaves_ions, gas_from_list, get_gas_velocity,gas_get_scatter_flag,gas_scatter,dump_gas
  public :: read_gas,gas_destructor,read_gas_composition_params,print_gas_composition_params
  !Val
  public :: gas_get_CD !, gas_get_n_CD
  !laV
  !--PEEL--
  public :: gas_peeloff_weight,gas_get_tau!,gas_get_tau_gray
  !--LEEP--

  public :: HI_core_skip


contains

    
  subroutine gas_from_list(ncells, cell_pos, cell_l, gas_leaves)

    integer(kind=4), intent(out)		:: ncells
    integer(kind=4), allocatable, intent(out)	:: cell_l(:)
    real(kind=8), allocatable, intent(out)	:: cell_pos(:,:)
    type(gas), allocatable, intent(out)		:: gas_leaves(:)

    integer(kind=4)				:: i
    real(kind=8), allocatable			:: nhi(:), nhii(:), metallicity(:)


    open(unit=20, file=input_ramses_file, status='old', form='unformatted')
    open(unit=21, file=Ion_file, status='old', form='unformatted') ! to comment for test

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
    read(20) nhi ! to comment for test
    read(20) nhii ! to comment for test
    read(20) metallicity ! to comment for test
    read(20) gas_leaves%dopwidth
    gas_leaves%dopwidth = gas_leaves%dopwidth / sqrt(15.999)  !15.999 is the atomic mass unit of Oxygen.  In the data file I printed sqrt(2*kb*T/m_u), so it's correct for Hydrogen, but other elements have to be divided by sqrt(mass of the element in atomic units). ! to comment for test
    !read(20) gas_leaves%ndust !For test
    !read(20) gas_leaves%nOI !For test
    read(21) gas_leaves%nOI ! to comment for test

    close(20) ; close(21)

    gas_leaves%ndust = metallicity / Zref * ( nhi + f_ion*nhii )   ! [ /cm3 ] ! to comment for test


    if (verbose) print*,'boxsize in cm : ', box_size_cm
    if (verbose) print*,'min/max of vth   : ',minval(gas_leaves(:)%dopwidth),maxval(gas_leaves(:)%dopwidth)
    if (verbose) print*,'min/max of nOI : ',minval(gas_leaves(:)%nOI),maxval(gas_leaves(:)%nOI)
    if (verbose) print*,'min/max of ndust : ',minval(gas_leaves(:)%ndust),maxval(gas_leaves(:)%ndust)

    deallocate(nhi, nhii, metallicity)


  end subroutine gas_from_list
  ! --------------------------------------------------------------------------
  !--laV

  !Val--
  ! --------------------------------------------------------------------------
  subroutine gas_from_ramses_leaves_ions(repository,snapnum,nleaf,nvar,ramses_var,OI_density,g)

    use module_ramses

    character(2000),intent(in)                     :: repository 
    integer(kind=4),intent(in)                     :: snapnum
    integer(kind=4),intent(in)                     :: nleaf,nvar
    real(kind=8),intent(in),dimension(nvar,nleaf)  :: ramses_var
    real(kind=8),intent(in),dimension(nleaf)       :: OI_density
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
       write(*,*) '-- module_gas_composition_OI_1302_dust : extracting velocities from ramses '
       allocate(v(3,nleaf))
       call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
       do ileaf = 1,nleaf
          g(ileaf)%v = v(:,ileaf)
       end do
       deallocate(v)

       ! get nHI and temperature from ramses
       write(*,*) '-- module_gas_composition_OI_1302_dust : extracting nHI and T from ramses '
       allocate(T(nleaf),nhi(nleaf))
       call ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nhi)
       
       ! compute thermal velocity 
       ! ++++++ TURBULENT VELOCITY >>>>> parameter to add and use here
       g(:)%dopwidth = sqrt((2.0d0*kb/mp/15.999)*T) ! [ cm/s ]   !15.999 is the atomic mass of oxygen

       ! get ndust (pseudo dust density from Laursen, Sommer-Larsen, Andersen 2009)
       write(*,*) '-- module_gas_composition_OI_1302_dust : extracting ndust from ramses '
       allocate(metallicity(nleaf),nhii(nleaf))
       call ramses_get_metallicity(nleaf,nvar,ramses_var,metallicity)
       call ramses_get_nh_cgs(repository,snapnum,nleaf,nvar,ramses_var,nhii)
       nhii = nhii - nhi
       do ileaf = 1,nleaf
          g(ileaf)%ndust = metallicity(ileaf) / Zref * ( nhi(ileaf) + f_ion*nhii(ileaf) )   ! [ /cm3 ]
       end do
       deallocate(metallicity,T,nhi,nhii)

       g(:)%nOI = OI_density
    end if

    return


  end subroutine gas_from_ramses_leaves_ions
  ! --------------------------------------------------------------------------
  
  ! --------------------------------------------------------------------------
  subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var,g,dom)

    ! define gas contents from ramses raw data, using Krome
    !NOT READY YET !!!!!!

    use module_ramses
    use module_spectra
    use module_domain
    !use krome_main
    !use krome_user
    !use krome_user_commons

    character(2000),intent(in)                     :: repository 
    integer(kind=4),intent(in)                     :: snapnum
    integer(kind=4),intent(in)                     :: nleaf,nvar
    type(domain),optional,intent(in)               :: dom
    real(kind=8),intent(in),dimension(nvar,nleaf)  :: ramses_var
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                                :: ileaf, count, nSEDgroups, i
    real(kind=8),dimension(:),allocatable          :: T, nhi, metallicity, nhii
    real(kind=8),dimension(:,:),allocatable        :: v, fractions, flux, csn
    !integer,parameter	                           :: nsp=krome_nmols, nIon=2, nPhoto=krome_nPhotoRates
   ! real(kind=8)                                   :: x(nsp),nHe,nhi_i,nhii_i,nhei,nheii,nheiii, nO, Tgas, rates(nPhoto)

    ! allocate gas-element array
    allocate(g(nleaf))

    if (gas_overwrite) then
       call overwrite_gas(g)
    else

       box_size_cm = ramses_get_box_size_cm(repository,snapnum)

       ! compute velocities in cm / s
       write(*,*) '-- module_gas_composition_OI_1302_dust : extracting velocities from ramses '
       allocate(v(3,nleaf))
       call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
       do ileaf = 1,nleaf
          g(ileaf)%v = v(:,ileaf)
       end do
       deallocate(v)

       ! get nHI and temperature from ramses
       write(*,*) '-- module_gas_composition_OI_1302_dust : extracting nHI and T from ramses '
       allocate(T(nleaf),nhi(nleaf))
       call ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nhi)

       ! compute thermal velocity 
       ! ++++++ TURBULENT VELOCITY >>>>> parameter to add and use here
       g(:)%dopwidth = sqrt((2.0d0*kb/15.999/mp)*T) ! [ cm/s ]

       ! get ndust (pseudo dust density from Laursen, Sommer-Larsen, Andersen 2009)
       write(*,*) '-- module_gas_composition_OI_1302_dust : extracting ndust from ramses '
       allocate(metallicity(nleaf),nhii(nleaf))
       call ramses_get_metallicity(nleaf,nvar,ramses_var,metallicity)
       call ramses_get_nh_cgs(repository,snapnum,nleaf,nvar,ramses_var,nhii)
       nhii = nhii - nhi
       do ileaf = 1,nleaf
          g(ileaf)%ndust = metallicity(ileaf) / Zref * ( nhi(ileaf) + f_ion*nhii(ileaf) )   ! [ /cm3 ]
       end do
   
       !get H and He ionization fractions from Ramses
       allocate(fractions(3,nleaf))
       call ramses_get_fractions(nleaf,nvar,ramses_var,fractions)
       !get flux from Ramses
       nSEDgroups = get_nSEDgroups(repository,snapnum)
       allocate(flux(nSEDgroups,nleaf))
       call ramses_get_flux(repository,snapnum,nleaf,nvar,nSEDgroups,ramses_var,flux)
       
       !call compute_csn_in_domain(repository,snapnum,dom,nIon,8,csn)

       ! !init krome (mandatory)
       ! call krome_init()

       ! !Loop on leaf cells to compute nOI with Krome
       ! count = 0 !Count the number of Krome bugs
       ! print*, 'Beginning of the loop on leaf cells !'

       ! !$OMP PARALLEL &
       ! !$OMP PRIVATE(ileaf, rates, x) &
       ! !$OMP SHARED(nleaf, nhi, nhii, fractions, metallicity, T, flux, csn, g)
       ! !$OMP DO
       
       ! do ileaf=1,nleaf

       !    rates(:) = (/ (sum(flux(:,ileaf)*csn(:,i)), i=1,nIon) /)
       !    call krome_set_photoBin_rates(rates)
          
       !    x(:) = (/ max(nhii(ileaf) + 7.895d-2*(nhi(ileaf)+nhii(ileaf))*fractions(2,ileaf) + 2*7.895d-2*(nhi(ileaf)+nhii(ileaf))*fractions(3,ileaf),1d-18), max(nhi(ileaf),1d-18), max(7.895d-2*(nhi(ileaf)+nhii(ileaf))*(1d0-fractions(2,ileaf)-fractions(3,ileaf)),1d-18), max(metallicity(ileaf)/0.0134*4.9d-4*(nhi(ileaf)+nhii(ileaf)),1d-18), max(nhii(ileaf),1d-18), max(7.895d-2*(nhi(ileaf)+nhii(ileaf))*fractions(2,ileaf),1d-18), 1d-18, max(7.895d-2*(nhi(ileaf)+nhii(ileaf))*fractions(3,ileaf),1d-18), 1d-18 /)
       !    print*,'1',x(4)
       !    call krome_equilibrium(x(:), T(ileaf))
       !    print*,'2',x(4)

       !    g(ileaf)%nOI = max(x(4),1d-18)

       ! end do
       ! !$OMP END DO
       ! !$OMP END PARALLEL


       ! ! do ileaf=1,nleaf

       ! !    if(modulo(ileaf,1000)==0 .and. verbose) print*,ileaf, 'cells computed'
          
       ! !    nhi_i = nhi(ileaf) ; nhii_i = nhii(ileaf)
       ! !    nHe = 7.895d-2*(nhi_i+nhii_i) !Assuming mass fraction of helium of 0.24
       ! !    nhei = nHe*(1d0-fractions(2,ileaf)-fractions(3,ileaf)) ; nheii = nHe*fractions(2,ileaf) ; nheiii = nHe*fractions(3,ileaf)
       ! !    nO = metallicity(ileaf)/0.0134*4.9d-4*(nhi_i+nhii_i) !Assuming solar abundances
       ! !    Tgas = T(ileaf)
          
       ! !    rates(:) = (/ (sum(flux(:,ileaf)*csn(:,i)), i=1,nIon) /)
       ! !    call krome_set_photoBin_rates(rates)
          
       ! !    !To gain time on convergence to equilibrium, distinguish cases based on temperature
       ! !    if(Tgas > 1d6) then
       ! !       g(ileaf)%nOI = 1d-18
       ! !    else
       ! !       !Here, the gas is ~dominated by OIII   
       ! !       if(Tgas > 5.5d4) then
       ! !          x(:) = (/ max(nhii_i + nheii + 2*nheiii + 2*nO,1d-18), max(nhi_i,1d-18), max(nhei,1d-18), 1d-18, max(nhii_i,1d-18), max(nheii,1d-18), 1d-18, max(nheiii,1d-18), max(nO,1d-18) /)
       ! !          !Here, the gas is ~dominated by OII
       ! !       elseif(Tgas > 1.6d4) then
       ! !          x(:) = (/ max(nhii_i + nheii + 2*nheiii + nO,1d-18), max(nhi_i,1d-18), max(nhei,1d-18), 1d-18, max(nhii_i,1d-18), max(nheii,1d-18), max(nO,1d-18), max(nheiii,1d-18), 1d-18 /)
       ! !          !Dominated by OI
       ! !       else
       ! !          x(:) = (/ max(nhii_i + nheii + 2*nheiii,1d-18), max(nhi_i,1d-18), max(nhei,1d-18), max(nO,1d-18), max(nhii_i,1d-18), max(nheii,1d-18), 1d-18, max(nheiii,1d-18), 1d-18 /)
       ! !       end if
             
       ! !       call krome_equilibrium(x(:), Tgas)
             
       ! !       !Case where there is a bug (rare !)
       ! !       if(max(x(4), x(7), x(9))/nO > 1.0001) then  !cell bugs,  < 1/100000,   no real solution
       ! !          count = count + 1
       ! !          print*, 'bug'
       ! !          if(Tgas < 1.6d4) then
       ! !             g(ileaf)%nOI =  nO
       ! !          else
       ! !             g(ileaf)%nOI =  1d-18
       ! !          end if
       ! !       else !no bug
       ! !          g(ileaf)%nOI = max(x(4),1d-18)
       ! !       end if
       ! !    end if
       ! !    print*,'test6'
          
       ! ! end do
       ! deallocate(metallicity,nhi,nhii,fractions,flux,T)
    end if

    return

  end subroutine gas_from_ramses_leaves
  ! --------------------------------------------------------------------------
  !--laV


  ! ! --------------------------------------------------------------------------
  ! function gas_get_tau_gray(cell_gas, distance_cm, nu_cell)

  !   ! --------------------------------------------------------------------------
  !   ! compute opacity due to dust accross distance_cm
  !   ! --------------------------------------------------------------------------
  !   ! INPUTS:
  !   ! - cell_gas : a mix of SiII and dust
  !   ! - distance_cm : the distance along which to compute tau [cm]
  !   ! OUTPUTS:
  !   ! - gas_get_tau_gray : the total optical depth

  !   ! --------------------------------------------------------------------------

  !   ! check whether scattering occurs within cell (scatter_flag > 0) or not (scatter_flag==0)
  !   type(gas),intent(in)    :: cell_gas
  !   real(kind=8),intent(in) :: distance_cm, nu_cell
  !   real(kind=8)            :: gas_get_tau_gray

  !   ! compute optical depths for different components of the gas.
  !   gas_get_tau_gray = get_tau_dust(cell_gas%ndust, distance_cm, nu_cell)

  !   return

  ! end function gas_get_tau_gray
  ! ! --------------------------------------------------------------------------


  ! !Val
  ! function gas_get_n_CD()

  !   integer(kind=4)           :: gas_get_n_CD

  !   gas_get_n_CD = 2

  ! end function gas_get_n_CD
  ! !laV

  !Val
  function gas_get_CD(cell_gas, distance_cm)

    real(kind=8), intent(in) :: distance_cm
    type(gas),intent(in)     :: cell_gas
    real(kind=8)             :: gas_get_CD

    !gas_get_CD = distance_cm*cell_gas%nOI
    gas_get_CD = distance_cm*cell_gas%ndust

  end function gas_get_CD
  !laV


  subroutine overwrite_gas(g)

    type(gas),dimension(:),intent(inout) :: g

    box_size_cm   = fix_box_size_cm
    g(:)%nOI      = fix_nOI
    g(:)%v(1)     = fix_vel
    g(:)%v(2)     = fix_vel
    g(:)%v(3)     = fix_vel
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
    ! - cell_gas : OI (with two absorption channels) and dust
    ! - distance_to_border_cm : the maximum distance the photon may travel (before leaving the cell)
    ! - nu_cell : photon frequency in cell's frame [ Hz ]
    ! - tau_abs : optical depth at which the next scattering event will occur
    ! - iran    : random generator state of the photon 
    ! OUTPUTS:
    ! - distance_to_border_cm : comes out as the distance to scattering event (if there is an event)
    ! - tau_abs : 0 if a scatter occurs, decremented by tau_cell if photon escapes cell. 
    ! - gas_get_scatter_flag : 0 [no scatter], 1 [OI-1302 scatter], 2 [dust] 
    ! --------------------------------------------------------------------------

    type(gas),intent(in)                  :: cell_gas
    real(kind=8),intent(inout)            :: distance_to_border_cm
    real(kind=8),intent(in)               :: nu_cell
    real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur.
    integer,intent(inout)                 :: iran 
    integer(kind=4)                       :: gas_get_scatter_flag 
    real(kind=8)                          :: tau_OI_1302, tau_dust, tau_cell, proba60, x

    ! compute optical depths for different components of the gas.
    tau_OI_1302   = get_tau_OI_1302(cell_gas%nOI, cell_gas%dopwidth, distance_to_border_cm, nu_cell)
    tau_dust      = get_tau_dust(cell_gas%ndust, distance_to_border_cm, nu_cell)
    tau_cell      = tau_OI_1302 + tau_dust

    if (tau_abs > tau_cell) then  ! photon is due for absorption outside the cell 
       gas_get_scatter_flag = 0
       tau_abs = tau_abs - tau_cell
       if (tau_abs.lt.0.d0) then
          print*, 'tau_abs est negatif'
          stop
       endif
    else  ! the scattering happens inside the cell
       ! decide if it is 1302 or dust 
       proba60 = tau_OI_1302 /  tau_cell
       x = ran3(iran)
       if (x <= proba60) then
          gas_get_scatter_flag = 1 ! absorption by OI_1302
       else
          gas_get_scatter_flag = 2 ! absorption by dust
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
    ! - cell_gas : a mix of OI and dust
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
    real(kind=8)            :: tau_OI, tau_dust

    ! compute optical depths for different components of the gas.
    tau_OI   = get_tau_OI_1302(cell_gas%nOI, cell_gas%dopwidth, distance_cm, nu_cell)
    tau_dust = get_tau_dust(cell_gas%ndust, distance_cm, nu_cell)
    gas_get_tau = tau_OI + tau_dust

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
       gas_peeloff_weight = OI_1302_peeloff_weight(cell_gas%v, cell_gas%dopwidth, nu_ext, kin, kout, iran)
    case(2)
       gas_peeloff_weight = dust_peeloff_weight(cell_gas%v, nu_ext, kin, kout)
    case default
       print*,'ERROR in module_gas_composition_OI1302_dust.f90:gas_peeloff_weight - unknown case : ',flag 
       stop
    end select

  end function gas_peeloff_weight
  !--LEEP--



  subroutine gas_scatter(flag,cell_gas,nu_cell,k,nu_ext,iran)

    integer, intent(inout)                    :: flag
    type(gas), intent(in)                     :: cell_gas
    real(kind=8), intent(inout)               :: nu_cell, nu_ext
    real(kind=8), dimension(3), intent(inout) :: k
    integer, intent(inout)                    :: iran
    integer(kind=4)                           :: ilost 

    select case(flag)
    case(1)
       call scatter_OI_1302(cell_gas%v, cell_gas%dopwidth, nu_cell, k, nu_ext, iran)
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
    write(unit) (g(i)%nOI, i=1,nleaf)
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
       read(unit) (g(i)%nOI,i=1,n)
       read(unit) (g(i)%dopwidth,i=1,n)
       read(unit) (g(i)%ndust,i=1,n)
       read(unit) box_size_cm 
    end if

    if (verbose) print*,'min/max of nOI : ',minval(g(:)%nOI),maxval(g(:)%nOI)

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
          case ('fix_vth')
             read(value,*) fix_vth
          case ('fix_vel')
             read(value,*) fix_vel
          case ('fix_ndust')
             read(value,*) fix_ndust
          case ('verbose')
             read(value,*) verbose
          case ('fix_box_size_cm')
             read(value,*) fix_box_size_cm
          end select
       end do
    end if
    close(10)

    call read_ramses_params(pfile)
    call read_dust_params(pfile)
    call read_OI_1302_params(pfile)
    call read_spectra_params(pfile)

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
       !Val---
       write(unit,'(a,a)')      'input_ramses_file      = ',trim(input_ramses_file)
       write(unit,'(a,a)')      'Ion_file               = ',trim(Ion_file)
       !--Val
       write(unit,'(a,ES10.3)') '  f_ion                = ',f_ion
       write(unit,'(a,ES10.3)') '  Zref                 = ',Zref
       write(unit,'(a)')       '# overwrite parameters'
       write(unit,'(a,L1)')    '  gas_overwrite         = ',gas_overwrite
       write(unit,'(a,ES10.3)') '  fix_vth              = ',fix_vth
       write(unit,'(a,ES10.3)') '  fix_vel              = ',fix_vel
       write(unit,'(a,ES10.3)') '  fix_box_size_cm      = ',fix_box_size_cm
       write(unit,'(a)')       '# miscelaneous parameters'
       write(unit,'(a,L1)')    '  verbose               = ',verbose
       write(unit,'(a)')             ' '
       call print_ramses_params(unit)
       write(unit,'(a)')             ' '
       call print_dust_params
       call print_OI_1302_params(unit)
       call print_spectra_params(unit)
    else
       write(*,'(a,a,a)') '[gas_composition]'
       !Val---
       write(*,'(a,a)')      'input_ramses_file      = ',trim(input_ramses_file)
       write(*,'(a,a)')      'Ion_file               = ',trim(Ion_file)
       !--Val
       write(*,'(a,ES10.3)') '  f_ion                = ',f_ion
       write(*,'(a,ES10.3)') '  Zref                 = ',Zref
       write(*,'(a)')       '# overwrite parameters'
       write(*,'(a,L1)')    '  gas_overwrite         = ',gas_overwrite
       write(*,'(a,ES10.3)') '  fix_vth              = ',fix_vth
       write(*,'(a,ES10.3)') '  fix_vel              = ',fix_vel
       write(*,'(a,ES10.3)') '  fix_box_size_cm      = ',fix_box_size_cm
       write(*,'(a)')       '# miscelaneous parameters'
       write(*,'(a,L1)')    '  verbose               = ',verbose
       write(*,'(a)')             ' '
       call print_ramses_params
       write(*,'(a)')             ' '
       call print_dust_params
       call print_OI_1302_params()
       call print_spectra_params
    end if

    return

  end subroutine print_gas_composition_params



end module module_gas_composition
