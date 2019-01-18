module module_OI_1302_model

  ! This module describes the absorption of photons by OI from level  2s^2 2p^4 3P 2 to level 2s^2 2p^3 (4S0) 3s 3S0 1.
  ! This transition is at 1302.168 A.
  ! The module also implements the three decay channels (resonant and fluorescent) at 1302.168, 1304.858 and 1306.029 A. 

  use module_constants
  use module_utils, only : isotropic_direction
  use module_uparallel
  use module_random
  use module_voigt

  implicit none

  private

  ! Atomic data, taken from Scarlata and Panagia, AjJ 801, 2015 (Table 1)
  ! In this module, we use the following convention :
  ! level 1 is 2s^2 2p^4 3P 2
  ! level 2 is 2s^2 2p^4 3P 1
  ! level 3 is 2s^2 2p^4 3P 0
  ! level 4 is 2s^2 2p^3(4S0) 3s 3S0 1

  ! transition between levels 1 and 4 
  real(kind=8),parameter :: lambda14       = 1302.168d0                   ! transition wavelength [A]
  real(kind=8),parameter :: lambda14_cm    = lambda14 / cmtoA             ! [cm]
  real(kind=8),parameter :: nu14           = clight / lambda14_cm         ! [Hz]
  real(kind=8),parameter :: f14            = 5.2d-2                       ! oscillator strength
  real(kind=8),parameter :: sigma14_factor = sqrtpi*e_ch**2*f14/me/clight ! multiply by Voigt(x,a)/nu_D to get sigma.
  real(kind=8),parameter :: A41            = 3.41d8                       ! spontaneous decay [/s]

  ! transition between levels 2 and 4
  real(kind=8),parameter :: lambda24       = 1304.858                 ! transition wavelength [A]
  real(kind=8),parameter :: lambda24_cm    = lambda24 / cmtoA         ! [cm]
  real(kind=8),parameter :: nu24           = clight / lambda24_cm     ! [Hz]
  real(kind=8),parameter :: A42            = 2.03d8                   ! spontaneous decay [/s]

  ! transition between levels 3 and 4
  real(kind=8),parameter :: lambda34       = 1306.029                 ! transition wavelength [A]
  real(kind=8),parameter :: lambda34_cm    = lambda34 / cmtoA         ! [cm]
  real(kind=8),parameter :: nu34           = clight / lambda34_cm     ! [Hz]
  real(kind=8),parameter :: A43            = 6.76d7                   ! spontaneous decay [/s]

  real(kind=8),parameter   :: A41_over_all = A41 / (A41+A42+A43)
  real(kind=8),parameter   :: A412_over_all = (A41+A42) / (A41+A42+A43)

  public :: get_tau_OI_1302, scatter_OI_1302, read_OI_1302_params, print_OI_1302_params
  !--PEEL--
  public :: OI_1302_peeloff_weight
  !--LEEP--


contains

  function get_tau_OI_1302(nOI, vth, distance_to_border_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute optical depth of OI-1302.168 over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - nOI      : number density of OI ions                              [ cm^-3 ]
    ! - vth      : thermal (+ small-scale turbulence) velocity of OI ions [ cm / s ]
    ! - distance_to_border_cm : distance over which we compute tau          [ cm ]
    ! - nu_cell  : photon's frequency in the frame of the cell              [ Hz ]
    ! OUTPUT :
    ! - get_tau_OI_1302 : optical depth of Silicon's line over distance_to_border_cm
    ! --------------------------------------------------------------------------

    real(kind=8),intent(in) :: nOI,vth,distance_to_border_cm,nu_cell
    real(kind=8)            :: delta_nu_doppler,x_cell,sigma,a,h_cell,get_tau_OI_1302, test

    ! compute Doppler width and a-parameter
    delta_nu_doppler = vth / lambda14_cm
    a = A41 / (fourpi * delta_nu_doppler)

    ! cross section of OI-1302.168
    x_cell = (nu_cell - nu14) / delta_nu_doppler
    h_cell = voigt_function(x_cell,a)
    sigma  = sigma14_factor / delta_nu_doppler * h_cell
    test = voigt_function(0.6d0, 1.6d0)
    print*, 'test ', test

    get_tau_OI_1302 = sigma * nOI * distance_to_border_cm

    return

  end function get_tau_OI_1302



  subroutine scatter_OI_1302(vcell,vth,nu_cell,k,nu_ext,iran)

    ! ---------------------------------------------------------------------------------
    ! perform scattering event on a OI ion
    ! The photon is absorbed in transition 1->4 and may decay as 4->1, 4->2 or 4->3. 
    ! ---------------------------------------------------------------------------------
    ! INPUTS :
    ! - vcell    : bulk velocity of the gas (i.e. cell velocity)       [ cm / s ] 
    ! - vth      : thermal (+turbulent) velocity dispersion of H (and Si) atoms [ cm / s ] 
    ! - nu_cell  : frequency of incoming photon in cell's rest-frame   [ Hz ] 
    ! - k        : propagaction vector (normalized) 
    ! - nu_ext   : frequency of incoming photon, in external frame     [ Hz ]
    ! - iran     : random number generator seed
    ! OUTPUTS :
    ! - nu_cell  : updated frequency in cell's frame   [ Hz ]
    ! - nu_ext   : updated frequency in external frame [ Hz ]
    ! - k        : updated propagation direction
    ! _ iran     : updated value of seed
    ! ---------------------------------------------------------------------------------

    real(kind=8),intent(inout)              :: nu_cell, nu_ext
    real(kind=8),dimension(3),intent(inout) :: k
    real(kind=8),dimension(3),intent(in)    :: vcell
    real(kind=8),intent(in)                 :: vth
    integer(kind=4),intent(inout)           :: iran
    real(kind=8)                            :: delta_nu_doppler, a, x_cell, upar, ruper
    real(kind=8)                            :: r2, uper, nu_atom, mu, bu, scalar, proba41
    real(kind=8),dimension(3)               :: knew

    ! define x_cell & a
    delta_nu_doppler = vth / lambda14_cm 
    a = A41 / fourpi / delta_nu_doppler
    x_cell = (nu_cell - nu14) / delta_nu_doppler

    ! 1/ component parallel to photon's propagation
    ! -> get velocity of interacting atom parallel to propagation
    upar = get_uparallel(x_cell,a,iran)
    upar = upar * vth    ! upar is an x -> convert to a velocity

    ! 2/ component perpendicular to photon's propagation
    ruper  = ran3(iran)
    r2     = ran3(iran)
    uper   = sqrt(-log(ruper))*cos(twopi*r2)
    uper   = uper * vth  ! from x to velocity

    ! 3/ chose de-excitation channel to determine output freq. in atom's frame
    r2 = ran3(iran)
    if (r2 <= A41_over_all) then
       ! photon goes down to level 1 -> coherent scattering
       nu_atom = nu_cell - nu_ext * upar/clight ! incoming frequency in atom's frame = outcoming freq in same frame
    elseif(r2 <= A412_over_all) then
       ! photons goes down to level two ...
       nu_atom = nu24
    else
       nu_atom = nu34
    end if

    ! 4/ determine direction of scattered photon
    call isotropic_direction(knew,iran)
    mu = k(1)*knew(1) + k(2)*knew(2) + k(3)*knew(3)
    bu = sqrt(1.0d0 - mu*mu)

    ! 5/ compute atom freq. in external frame, after scattering
    scalar = knew(1) * vcell(1) + knew(2) * vcell(2) + knew(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + bu*uper)/clight)
    nu_cell = (1.d0 - scalar/clight) * nu_ext 
    k = knew

  end subroutine scatter_OI_1302




  !--PEEL--
  function OI_1302_peeloff_weight(vcell,vth,nu_ext,kin,kout,iran)

    ! ---------------------------------------------------------------------------------
    ! Compute probability that a photon coming along kin scatters off in direction kout.
    ! Also update nu_ext to external-frame frequency along kout
    ! ---------------------------------------------------------------------------------
    ! INPUTS :
    ! - vcell    : bulk velocity of the gas (i.e. cell velocity)       [ cm / s ] 
    ! - vth      : thermal (+turbulent) velocity dispersion of H atoms [ cm / s ] 
    ! - nu_ext   : frequency of incoming photon, in external frame     [ Hz ]
    ! - kin      : propagation vector (normalized)
    ! - kout     : direction after interaction (fixed)
    ! - iran     : random number generator seed
    ! OUTPUTS :
    ! - nu_ext   : updated frequency in external frame [ Hz ]
    ! _ iran     : updated value of seed
    ! ---------------------------------------------------------------------------------

    real(kind=8),intent(inout)              :: nu_ext
    real(kind=8),dimension(3),intent(in)    :: kin, kout
    real(kind=8),dimension(3),intent(in)    :: vcell
    real(kind=8),intent(in)                 :: vth
    integer(kind=4),intent(inout)           :: iran
    real(kind=8)                            :: OI_1302_peeloff_weight
    real(kind=8)                            :: delta_nu_doppler, a, x_cell, upar, ruper
    real(kind=8)                            :: r2, uper, nu_atom, mu, bu, scalar
    real(kind=8)                            :: x_atom,nu_cell


    ! compute frequency in cell's frame 
    scalar  = kin(1) * vcell(1) + kin(2) * vcell(2) + kin(3) * vcell(3)
    nu_cell = (1.d0 - scalar/clight) * nu_ext

    ! define x_cell & a
    delta_nu_doppler = vth / lambda14_cm 
    a = A41 / fourpi / delta_nu_doppler
    x_cell = (nu_cell - nu14) / delta_nu_doppler
    
    ! 1/ component parallel to photon's propagation
    ! -> get velocity of interacting atom parallel to propagation
    upar = get_uparallel(x_cell,a,iran)
    upar = upar * vth    ! upar is an x -> convert to a velocity 

    ! 2/ component perpendicular to photon's propagation
    ruper  = ran3(iran)
    r2     = ran3(iran)
    uper   = sqrt(-log(ruper))*cos(twopi*r2)
    uper   = uper * vth  ! from x to velocity

    ! 3/ chose de-excitation channel to determine output freq. in atom's frame
    r2 = ran3(iran)
    if (r2 <= A41_over_all) then
       ! photon goes down to level 1 -> coherent scattering
       nu_atom = nu_cell - nu_ext * upar/clight ! incoming frequency in atom's frame = outcoming freq in same frame
    elseif(r2 <= A412_over_all) then
       ! photons goes down to level two ...
       nu_atom = nu24
    else
       nu_atom = nu34
    end if

    ! 4/ determine direction of scattered photon
    OI_1302_peeloff_weight = 0.5d0  ! P(mu) for isotropic phase function
    mu = kin(1)*kout(1) + kin(2)*kout(2) + kin(3)*kout(3)
    bu = sqrt(1.0d0 - mu*mu)

    ! 5/ compute atom freq. in external frame, after scattering
    scalar = kout(1) * vcell(1) + kout(2) * vcell(2) + kout(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + bu*uper)/clight)

  end function OI_1302_peeloff_weight
  !--LEEP--




  subroutine read_OI_1302_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
 
    call read_uparallel_params(pfile)
    call read_voigt_params(pfile)
    
    return

  end subroutine read_OI_1302_params


  subroutine print_OI_1302_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       call print_uparallel_params(unit)
       call print_voigt_params(unit)
    else
       call print_uparallel_params()
       call print_voigt_params()
    end if

    return

  end subroutine print_OI_1302_params


end module module_OI_1302_model
