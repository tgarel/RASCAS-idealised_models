module module_D_model

  use module_constants
  use module_utils, only : voigt_fit
  use module_uparallel
  use module_random
  use module_params, only : recoil
  
  implicit none

  private

  ! Deuterium properties
  real(kind=8),parameter   :: mdeut        = 2.d0 * mp           ! Deuterium atom's mass [ g ]
  real(kind=8),parameter   :: lambda_0_cm  = 1215.673d-8 * (1.0d0 + me/mdeut) / (1.0d0 + me/mp) ! wavelength of Lya of Deuterium [ cm ]
  real(kind=8),parameter   :: nu_0         = clight / lambda_0_cm  ! frequency of Deuterium's Lya [ Hz ]
  real(kind=8),parameter   :: gamma        = 6.265d8             ! Einstein coeff. [ s^-1 ]
  real(kind=8),parameter   :: f12          = 0.416               ! Oscillator strength for Deuterium Lya.
  real(kind=8),parameter   :: sigma_factor = pi*e_ch**2*f12/ me / clight ! cross-section factor-> multiply by Voigt(x,a)/nu_D to get sigma.
    
  public :: get_tau_D, scatter_D_isotrope 

contains
  ! PUBLIC functions : 
  ! - function get_tau(ndi, vth, distance_to_border_cm, nu_cell)
  ! - subroutine scatter_isotrope(vcell,vth, nu_cell, k, nu_ext, iran)
  ! PRIVATE function :
  ! - function voigt_fit(x,a)
  
  
  function get_tau_D(ndi, vth, distance_to_border_cm, nu_cell)
    
    ! --------------------------------------------------------------------------
    ! compute optical depth of Deuterium over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - ndi      : number density of neutral D atoms                      [ cm^-3 ]
    ! - vth      : thermal (+ small-scale turbulence) velocity of D atoms [ cm / s ]
    ! - distance_to_border_cm : distance over which we compute tau        [ cm ]
    ! - nu_cell  : photon's frequency in the frame of the cell            [ Hz ]
    ! OUTPUT :
    ! - get_tau_D : optical depth of Deuterium's Lya line over distance_to_border_cm
    ! --------------------------------------------------------------------------
    
    real(kind=8),intent(in) :: ndi,vth,distance_to_border_cm,nu_cell
    real(kind=8)            :: get_tau_D
    real(kind=8)            :: delta_nu_D, a, x, h, s
    
    ! compute Doppler width and a-parameter
    delta_nu_D = vth / lambda_0_cm ! == vth * nu_0 / clight
    a          = gamma / (4.d0 * pi * delta_nu_D)
    
    ! Cross section of Deuterium
    x = (nu_cell - nu_0)/delta_nu_D
    h = voigt_fit(x,a)
    s = sigma_factor / delta_nu_D * h  
    
    ! optical depth 
    get_tau_D = s * ndi * distance_to_border_cm
    
    return
    
  end function get_tau_D
  
  
  subroutine scatter_D_isotrope(vcell,vth, nu_cell, k, nu_ext, iran)
    
    ! --------------------------------------------------------------------------
    ! perform scattering event on a Deuterium atom
    ! --------------------------------------------------------------------------
    ! INPUTS :
    ! - vcell   : bulk velocity of the gas (i.e. cell velocity)       [ cm / s ] 
    ! - vth     : thermal (+turbulent) velocity dispersion of D atoms [ cm / s ] 
    ! - nu_cell : frequency of incoming photon in cell's rest-frame   [ Hz ] 
    ! - k       : propagaction vector (normalized) 
    ! - nu_ext  : frequency of incoming photon, in external frame     [ Hz ]
    ! - iran    : random number generator seed
    ! OUTPUTS :
    ! - nu_cell : updated frequency in cell's frame   [ Hz ]
    ! - nu_ext  : updated frequency in external frame [ Hz ]
    ! - k       : updated propagation direction
    ! _ iran    : updated value of seed
    ! --------------------------------------------------------------------------
    
    real(kind=8), intent(inout)               :: nu_cell, nu_ext
    real(kind=8), dimension(3), intent(inout) :: k
    real(kind=8), dimension(3), intent(in)    :: vcell
    real(kind=8), intent(in)                  :: vth
    integer, intent(inout)                    :: iran
    real(kind=8)               :: delta_nu_doppler, a, x_cell, blah, upar, ruper
    real(kind=8)               :: r2, uper, nu_atom, phi, theta, st, mu, scalar
    real(kind=8), dimension(3) :: knew
    
    
    ! define x_cell & a
    delta_nu_doppler = vth / lambda_0_cm  ! [ Hz ]
    a      = gamma / (4.d0 * pi * delta_nu_doppler) 
    x_cell = (nu_cell - nu_0) / delta_nu_doppler
    
    ! 1/ component parallel to photon's propagation
    ! -> get velocity of interacting atom parallel to propagation
    blah = ran3(iran)
#ifdef SWITCH_OFF_UPARALLEL
    upar = 0.5  !!!!!todo get_uparallel(a,x_cell,blah)
#else
    upar = get_uparallel(a,x_cell,blah)
#endif
    upar = upar * vth    ! upar is an x -> convert to a velocity 
    
    ! 2/ component perpendicular to photon's propagation
    ruper  = ran3(iran)
    r2     = ran3(iran)
    uper   = sqrt(-log(ruper))*cos(2.d0*pi*r2)
    uper   = uper * vth  ! from x to velocity
    
    ! 3/ determine scattering angle (in atom's frame)
    nu_atom = nu_cell - nu_ext * upar/clight 
    phi   = 2.d0*pi*ran3(iran)
    theta = acos(1d0-2d0*ran3(iran))
    st = sin(theta)
    knew(1) = st*cos(phi)   !x
    knew(2) = st*sin(phi)   !y
    knew(3) = cos(theta)    !z
    mu = k(1)*knew(1) + k(2)*knew(2) + k(3)*knew(3) 
    
    ! 4/ recoil effect
    if (recoil) then
       nu_atom = nu_atom / (1.d0 + ((planck*nu_atom)/(mdeut*clight*clight))*(1.-mu))
    endif
    
    ! 5/ compute atom freq. in external frame, after scattering
    scalar  = knew(1) * vcell(1) + knew(2) * vcell(2) + knew(3)* vcell(3)
    nu_ext  = nu_atom * (1.0d0 + scalar/clight + (upar*mu + sqrt(1-mu**2)*uper)/clight)
    nu_cell = (1.d0 - scalar/clight) * nu_ext 
    k       = knew
    
  end subroutine scatter_D_isotrope
  
end module module_D_model
