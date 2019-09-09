module module_SiII_1526_model

  ! This module describes the absorption of photons by SiII from level  3s^2 3p 2P^0 1/2 to level 3s^2 4s 2S 1/2.
  ! This transition is at 1526.72 A.
  ! The module also implements the two decay channels (resonant and fluorescent) at 1526.72 A and 1533.43 A. 
  
  use module_constants
  use module_utils, only : isotropic_direction
  use module_uparallel
  use module_random
  use module_voigt
  
  implicit none

  private

  ! Atomic data, from the NIST database (https://www.nist.gov)
  ! In this module, we use the following convention :
  ! level 1 is 3s^2 3p 2P^0 1/2
  ! level 2 is 3s^2 3p 2P^0 3/2
  ! level 3 is 3s^2 3d 2D 3/2

  ! transition between levels 1 and 3 
  real(kind=8),parameter :: lambda13       = 1526.72d0                    ! transition wavelength [A]
  real(kind=8),parameter :: lambda13_cm    = lambda13 / cmtoA             ! [cm]
  real(kind=8),parameter :: nu13           = clight / lambda13_cm         ! [Hz]
  real(kind=8),parameter :: f13            = 1.33d-1                      ! oscillator strength
  real(kind=8),parameter :: sigma13_factor = sqrtpi*e_ch**2*f13/me/clight ! multiply by Voigt(x,a)/delta_nu_doppler to get sigma.
  real(kind=8),parameter :: A31            = 3.81d8                       ! spontaneous decay [/s]

  ! transition between levels 2 and 3
  real(kind=8),parameter :: lambda23       = 1533.45d0                    ! transition wavelength [A]
  real(kind=8),parameter :: lambda23_cm    = lambda23 / cmtoA             ! [cm]
  real(kind=8),parameter :: nu23           = clight / lambda23_cm         ! [Hz]
  real(kind=8),parameter :: A32            = 7.52d8                       ! spontaneous decay [/s]
    
  real(kind=8),parameter :: Atot = A31+A32


  public :: get_tau_SiII_1526, scatter_SiII_1526, read_SiII_1526_params, print_SiII_1526_params
  !--PEEL--
  public :: SiII_1526_peeloff_weight
  !--LEEP--

  
contains

  function get_tau_SiII_1526(nSiII, vth, distance_to_border_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute optical depth of SiII-1526.72 over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - nSiII    : number density of SiII ions                              [ cm^-3 ]
    ! - vth      : thermal (+ small-scale turbulence) velocity of SiII ions [ cm / s ]
    ! - distance_to_border_cm : distance over which we compute tau          [ cm ]
    ! - nu_cell  : photon's frequency in the frame of the cell              [ Hz ]
    ! OUTPUT :
    ! - get_tau_SiII_1526 : optical depth of Silicon's line over distance_to_border_cm
    ! --------------------------------------------------------------------------
    
    real(kind=8),intent(in) :: nSiII,vth,distance_to_border_cm,nu_cell
    real(kind=8)            :: delta_nu_doppler,x_cell,sigma,a,h_cell,get_tau_SiII_1526

    ! compute Doppler width and a-parameter
    delta_nu_doppler = vth / lambda13_cm
    a = A31 / (fourpi * delta_nu_doppler)

    ! cross section of SiII-1260.42
    x_cell = (nu_cell - nu13) / delta_nu_doppler
    h_cell = voigt_function(x_cell,a)
    sigma  = sigma13_factor / delta_nu_doppler * h_cell

    get_tau_SiII_1526 = sigma * nSiII * distance_to_border_cm
   
    return
    
  end function get_tau_SiII_1526


  subroutine scatter_SiII_1526(vcell,vth,nu_cell,k,nu_ext,iran)

    ! ---------------------------------------------------------------------------------
    ! perform scattering event on a SiII ion
    ! The photon is absorbed in transition 1->3 and may decay as 3->1 or 3->2. 
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
    real(kind=8)                            :: r2, uper, nu_atom, mu, bu, scalar, proba31
    real(kind=8),dimension(3)               :: knew

    ! define x_cell & a
    delta_nu_doppler = vth / lambda13_cm 
    a = A31 / fourpi / delta_nu_doppler
    x_cell = (nu_cell - nu13) / delta_nu_doppler

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
    proba31 = A31/Atot
    if (r2 <= proba31) then
       ! photon goes down to level 1 -> coherent scattering
       nu_atom = nu_cell - nu_ext * upar/clight ! incoming frequency in atom's frame = outcoming freq in same frame
    else
       ! photons goes down to level two ...
       nu_atom = nu23 
    end if
    
    ! 4/ determine direction of scattered photon
    call isotropic_direction(knew,iran)
    mu = k(1)*knew(1) + k(2)*knew(2) + k(3)*knew(3)
    bu = sqrt(1.0d0 - mu*mu)
    
    ! 5/ compute atom freq. in external frame, after scattering
    scalar = knew(1) * vcell(1) + knew(2) * vcell(2) + knew(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + bu*uper)/clight)
    nu_cell = (1.0d0 - scalar/clight) * nu_ext 
    k = knew

  end subroutine scatter_SiII_1526


  !--PEEL--
  function SiII_1526_peeloff_weight(vcell,vth,nu_ext,kin,kout,iran)

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
    real(kind=8)                            :: SiII_1526_peeloff_weight
    real(kind=8)                            :: delta_nu_doppler, a, x_cell, upar, ruper
    real(kind=8)                            :: r2, uper, nu_atom, mu, bu, scalar, proba31
    real(kind=8)                            :: x_atom,nu_cell


    ! compute frequency in cell's frame 
    scalar  = kin(1) * vcell(1) + kin(2) * vcell(2) + kin(3) * vcell(3)
    nu_cell = (1.d0 - scalar/clight) * nu_ext

    ! define x_cell & a
    delta_nu_doppler = vth / lambda13_cm 
    a = A31 / fourpi / delta_nu_doppler
    x_cell = (nu_cell - nu13) / delta_nu_doppler

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
    proba31 = A31/Atot
    if (r2 <= proba31) then
       ! photon goes down to level 1 -> coherent scattering
       nu_atom = nu_cell - nu_ext * upar/clight ! incoming frequency in atom's frame = outcoming freq in same frame
    else
       ! photons goes down to level two ...
       nu_atom = nu23 
    end if
    
    ! 4/ determine direction of scattered photon
    SiII_1526_peeloff_weight = 0.5d0  ! P(mu) for isotropic phase function
    mu = kin(1)*kout(1) + kin(2)*kout(2) + kin(3)*kout(3)
    bu = sqrt(1.0d0 - mu*mu)
    
    ! 5/ compute atom freq. in external frame, after scattering
    scalar = kout(1) * vcell(1) + kout(2) * vcell(2) + kout(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + bu*uper)/clight)

  end function SiII_1526_peeloff_weight
!--LEEP--

  
  subroutine read_SiII_1526_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile

    call read_uparallel_params(pfile)
    call read_voigt_params(pfile)

    return

  end subroutine read_SiII_1526_params


  subroutine print_SiII_1526_params(unit)
    
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
    
  end subroutine print_SiII_1526_params


end module module_SiII_1526_model