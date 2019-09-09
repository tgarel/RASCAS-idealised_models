module module_AlII_1670_model

  ! This module describes the absorption of photons by AlII from level 2p^6 3s^2 1S 0 to level 3s 3p 1P^0 1
  ! This transition is 1670.79 A.
  ! The module also implements one decay channel (resonant) at 1670.79 A.

  use module_constants
  use module_utils, only : isotropic_direction
  use module_uparallel
  use module_random
  use module_voigt

  implicit none

  private

  ! Atomic data, from the NIST database (https://www.nist.gov)
  ! In this module, we use the following convention :
  ! level 1 is 2p^6 3s^2 1S 0
  ! level 2 is 3s 3p 1P^0 1

  ! transition between levels 2 and 1
  real(kind=8),parameter :: lambda12       = 1670.79d0                    ! transition wavelength [A]
  real(kind=8),parameter :: lambda12_cm    = lambda12 / cmtoA             ! [cm]
  real(kind=8),parameter :: nu12           = clight / lambda12_cm         ! [Hz]
  real(kind=8),parameter :: f12            = 1.77d0                       ! oscillator strength
  real(kind=8),parameter :: sigma12_factor = sqrtpi*e_ch**2*f12/me/clight ! multiply by Voigt(x,a)/delta_nu_doppler to get sigma.
  real(kind=8),parameter :: A21            = 1.41d9                       ! spontaneous decay [/s]

  public :: get_tau_AlII_1670, scatter_AlII_1670, read_AlII_1670_params, print_AlII_1670_params
  !--PEEL--
  public :: AlII_1670_peeloff_weight
  !--LEEP--

contains

  function get_tau_AlII_1670(nAlII, vth, distance_to_border_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute optical depth of AlII-1670 over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - nAlII    : number density of AlII ions                              [ cm^-3 ]
    ! - vth      : thermal (+ small-scale turbulence) velocity of AlII ions [ cm / s ]
    ! - distance_to_border_cm : distance over which we compute tau          [ cm ]
    ! - nu_cell  : photon's frequency in the frame of the cell              [ Hz ]
    ! OUTPUT :
    ! - get_tau_AlII_1670 : optical depth of Magnesium's line over distance_to_border_cm
    ! --------------------------------------------------------------------------
    
    real(kind=8),intent(in) :: nAlII,vth,distance_to_border_cm,nu_cell
    real(kind=8)            :: delta_nu_doppler,x_cell,sigma,a,h_cell,get_tau_AlII_1670

    ! compute Doppler width and a-parameter
    delta_nu_doppler = vth / lambda12_cm
    a = A21 / (fourpi * delta_nu_doppler)

    ! cross section of AlII-2796
    x_cell = (nu_cell - nu12) / delta_nu_doppler
    h_cell = voigt_function(x_cell,a)
    sigma  = sigma12_factor / delta_nu_doppler * h_cell

    get_tau_AlII_1670 = sigma * nAlII * distance_to_border_cm
   
    return

  end function get_tau_AlII_1670

  
  subroutine scatter_AlII_1670(vcell,vth,nu_cell,k,nu_ext,iran)

    ! ---------------------------------------------------------------------------------
    ! perform scattering event on a AlII ion
    ! ---------------------------------------------------------------------------------
    ! INPUTS :
    ! - vcell    : bulk velocity of the gas (i.e. cell velocity)       [ cm / s ] 
    ! - vth      : thermal (+turbulent) velocity dispersion of H atoms [ cm / s ] 
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
    real(kind=8)                            :: r2, uper, nu_atom, mu, bu, scalar
    real(kind=8),dimension(3)               :: knew

    ! define x_cell & a
    delta_nu_doppler = vth / lambda12_cm 
    a = A21 / fourpi / delta_nu_doppler
    x_cell = (nu_cell - nu12) / delta_nu_doppler

    ! 1/ component parallel to photon's propagation
    ! -> get velocity of interacting atom parallel to propagation
    upar = get_uparallel(x_cell,a,iran)
    upar = upar * vth    ! upar is an x -> convert to a velocity 

    ! 2/ component perpendicular to photon's propagation
    ruper  = ran3(iran)
    r2     = ran3(iran)
    uper   = sqrt(-log(ruper))*cos(twopi*r2)
    uper   = uper * vth  ! from x to velocity

    ! 3/ incoming frequency in atom's frame
    nu_atom = nu_cell - nu_ext * upar/clight

    ! 4/ determine direction of scattered photon
    call isotropic_direction(knew,iran)
    mu = k(1)*knew(1) + k(2)*knew(2) + k(3)*knew(3)
    bu = sqrt(1.0d0 - mu*mu)
    
    ! 5/ compute atom freq. in external frame, after scattering
    scalar = knew(1) * vcell(1) + knew(2) * vcell(2) + knew(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + bu*uper)/clight)
    nu_cell = (1.0d0 - scalar/clight) * nu_ext 
    k = knew

  end subroutine scatter_AlII_1670


    !--PEEL--
  function AlII_1670_peeloff_weight(vcell,vth,nu_ext,kin,kout,iran)

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
    real(kind=8)                            :: AlII_1670_peeloff_weight
    real(kind=8)                            :: delta_nu_doppler, a, x_cell, upar, ruper
    real(kind=8)                            :: r2,uper, nu_atom, mu, bu, scalar
    real(kind=8)                            :: x_atom,nu_cell


    ! compute frequency in cell's frame 
    scalar  = kin(1) * vcell(1) + kin(2) * vcell(2) + kin(3) * vcell(3)
    nu_cell = (1.d0 - scalar/clight) * nu_ext

    ! define x_cell & a
    delta_nu_doppler = vth / lambda12_cm 
    a = A21 / fourpi / delta_nu_doppler
    x_cell = (nu_cell - nu12) / delta_nu_doppler

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
    nu_atom = nu_cell - nu_ext * upar/clight ! incoming frequency in atom's frame = outcoming freq in same frame
    
    ! 4/ determine direction of scattered photon
    AlII_1670_peeloff_weight = 0.5d0  ! P(mu) for isotropic phase function
    mu = kin(1)*kout(1) + kin(2)*kout(2) + kin(3)*kout(3)
    bu = sqrt(1.0d0 - mu*mu)
    
    ! 5/ compute atom freq. in external frame, after scattering
    scalar = kout(1) * vcell(1) + kout(2) * vcell(2) + kout(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + bu*uper)/clight)

  end function AlII_1670_peeloff_weight
!--LEEP--



  subroutine read_AlII_1670_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------
    
    character(*),intent(in) :: pfile
    
    call read_uparallel_params(pfile)
    call read_voigt_params(pfile)
    
    return
    
  end subroutine read_AlII_1670_params



  subroutine print_AlII_1670_params(unit)
    
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
    
  end subroutine print_AlII_1670_params



end module module_AlII_1670_model
