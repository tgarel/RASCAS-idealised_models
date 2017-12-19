module module_MgII_2796_model

  ! This module describes the absorption of photons by MgII from level 2p^6 3p 2P 3/2 to level 2p^6 3s 2S
  ! This transition is 2796.35 A.
  ! The module also implements one decay channel (resonant) at 2796.35 A.

  use module_constants
  use module_utils, only : voigt_fit, isotropic_direction
  use module_uparallel
  use module_random

  implicit none

  private

  ! Atomic data, taken from Zhu et al., ApJ 815, 2015 (table 2)
  ! In this module, we use the following convention :
  ! level 1 is 2p^6 3s 2S
  ! level 2 is 2p^6 3p 2P 3/2

  ! transition between levels 2 and 1
  real(kind=8),parameter :: lambda12       = 2796.35d0                ! transition wavelength [A]
  real(kind=8),parameter :: lambda12_cm    = lambda12 / cmtoA         ! [cm]
  real(kind=8),parameter :: nu12           = clight / lambda12_cm     ! [Hz]
  real(kind=8),parameter :: f12            = 0.608d0                  ! oscillator strength
  real(kind=8),parameter :: sigma12_factor = pi*e_ch**2*f12/me/clight ! multiply by Voigt(x,a)/nu_D to get sigma.
  real(kind=8),parameter :: A21            = 2.60d8                   ! spontaneous decay [/s]

  public :: get_tau_MgII_2796, scatter_MgII_2796, read_MgII_2796_params, print_MgII_2796_params

contains

  function get_tau_MgII_2796(nMgII, vth, distance_to_border_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute optical depth of MgII-2796 over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - nMgII    : number density of MgII ions                              [ cm^-3 ]
    ! - vth      : thermal (+ small-scale turbulence) velocity of SiII ions [ cm / s ]
    ! - distance_to_border_cm : distance over which we compute tau          [ cm ]
    ! - nu_cell  : photon's frequency in the frame of the cell              [ Hz ]
    ! OUTPUT :
    ! - get_tau_MgII_2796 : optical depth of Magnesium's line over distance_to_border_cm
    ! --------------------------------------------------------------------------
    
    real(kind=8),intent(in) :: nMgII,vth,distance_to_border_cm,nu_cell
    real(kind=8)            :: nu_D,x_cell,sigma,a,h,get_tau_MgII_2796

    ! compute Doppler width and a-parameter
    nu_D = vth / lambda12_cm
    a    = A21 / (fourpi * nu_D)

    ! cross section of MgII-2796
    x_cell = (nu_cell - nu12) / nu_D
    h      = voigt_fit(x_cell,a)
    sigma  = sigma12_factor / nu_D * h

    get_tau_MgII_2796 = sigma * nMgII * distance_to_border_cm
   
    return

  end function get_tau_MgII_2796

  
  subroutine scatter_MgII_2796(vcell,vth,nu_cell,k,nu_ext,iran)

    ! ---------------------------------------------------------------------------------
    ! perform scattering event on a MgII ion
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

  end subroutine scatter_MgII_2796



  subroutine read_MgII_2796_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------
    
    character(*),intent(in) :: pfile
    
    call read_uparallel_params(pfile)
    
    return
    
  end subroutine read_MgII_2796_params



  subroutine print_MgII_2796_params(unit)
    
    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------
    
    integer(kind=4),optional,intent(in) :: unit
    
    if (present(unit)) then 
       call print_uparallel_params(unit)
    else
       call print_uparallel_params()
    end if
    
    return
    
  end subroutine print_MgII_2796_params



end module module_MgII_2796_model
