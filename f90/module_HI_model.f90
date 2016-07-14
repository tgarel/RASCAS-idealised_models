module module_HI_model

  use module_constants
  use module_uparallel
  use module_random
  use module_params, only : recoil

!#ifdef POLARIZATION
!  use polar
!#endif

  implicit none

  private

  ! definition of atomic values
  real(kind=8),parameter   :: lambda_0=1215.6701                          ![A] Lya wavelength
  real(kind=8),parameter   :: gamma=6.265d8                               ! Einstein coeff = damping constant for Voigt Function(gamma_alpha)  
  real(kind=8),parameter   :: f12=0.416                                   ! oscillator strength for Ly-alpha
  ! useful pre-computed quantities
  real(kind=8),parameter   :: lambda_0_cm = lambda_0 / cmtoA              ! cm
  real(kind=8),parameter   :: nu_0 = clight / lambda_0_cm                 ! Hz
  real(kind=8),parameter   :: sigmaH_factor = pi*e_ch**2*f12/ me / clight ! H cross-section factor-> multiply by Voigt(x,a)/nu_D to get sigma.
  real(kind=8),parameter   :: gamma_over_fourpi = gamma / fourpi
  
  private :: voigt_fit

  public :: get_tau_HI, scatter_HI_isotrope
  
contains

  function get_tau_HI(nhi, dopwidth, distance_to_border_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute optical depth of Hydrogen over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - nhi      : number density of neutral HI atoms                      [ cm^-3 ]
    ! - dopwidth : thermal (+ small-scale turbulence) velocity of HI atoms [ cm / s ]
    ! - distance_to_border_cm : distance over which we compute tau        [ cm ]
    ! - nu_cell  : photon's frequency in the frame of the cell            [ Hz ]
    ! OUTPUT :
    ! - get_tau_HI : optical depth of Hydrogen's Lya line over distance_to_border_cm
    ! --------------------------------------------------------------------------
    
    real(kind=8),intent(in) :: nhi,dopwidth,distance_to_border_cm,nu_cell
    real(kind=8)            :: nu_D,x_cell,sigmaH,a,h, get_tau_HI

    ! compute Doppler width and a-parameter, for H 
    nu_D = dopwidth / lambda_0_cm 
    a    = gamma / (fourpi * nu_D)
 
    ! Cross section of H 
    x_cell  = (nu_cell - nu_0)/nu_D
    h       = voigt_fit(x_cell,a)
    sigmaH  = sigmaH_factor / nu_D * h
 
    get_tau_HI = sigmaH * nhi * distance_to_border_cm

    return

  end function get_tau_HI


  subroutine scatter_HI_isotrope(v,nHI,dopwidth, nu_cell, k, nu_ext, iran)

    ! ---------------------------------------------------------------------------------
    ! perform scattering event on a Hydrogen atom with isotrope angular redistribution
    ! ---------------------------------------------------------------------------------
    ! INPUTS :
    ! - v        : bulk velocity of the gas (i.e. cell velocity)       [ cm / s ] 
    ! - nhi      : number density of neutral HI atoms                  [ cm^-3 ]
    ! - dopwidth : thermal (+turbulent) velocity dispersion of H atoms [ cm / s ] 
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

    real(kind=8), intent(inout)               :: nu_cell, nu_ext
    real(kind=8), dimension(3), intent(inout) :: k
    real(kind=8), dimension(3), intent(in)    :: v
    real(kind=8), intent(in)                  :: nHI,dopwidth
    integer, intent(inout)                    :: iran

    real(kind=8)               :: nu_doppler, a, x_cell, blah, upar, ruper, x
    real(kind=8)               :: r2, uper, nu_atom, phi, theta, st, mu, scalar
    real(kind=8), dimension(3) :: knew

#ifdef DEBUG
    print *,'-DEBUG- scatter routine, arguments =',v,nHI,dopwidth, nu_cell, k, nu_ext, iran
#endif

    ! define x_cell & a
    nu_doppler = dopwidth * nu_0 / clight
    a = gamma / (fourpi * nu_doppler)
    x_cell = (nu_cell - nu_0) / nu_doppler

    ! 1/ component parallel to photon's propagation
    ! -> get velocity of interacting atom parallel to propagation
    blah = ran3(iran)
#ifdef SWITCH_OFF_UPARALLEL
    upar = 0.5  !!!!!todo get_uparallel(a,x_cell,blah)
#else
    upar = get_uparallel(a,x_cell,blah)
#endif
    upar = upar * dopwidth    ! upar is an x -> convert to a velocity 

    ! 2/ component perpendicular to photon's propagation
    ruper  = ran3(iran)
    r2     = ran3(iran)
    uper   = sqrt(-log(ruper))*cos(twopi*r2)
    uper   = uper * dopwidth  ! from x to velocity

    ! 3/ determine scattering angle (in atom's frame)
    nu_atom = nu_cell - nu_ext * upar/clight

    phi   = twopi*ran3(iran)
    theta = acos(1d0-2d0*ran3(iran))
    !........director cosines
    st = sin(theta)
    knew(1) = st*cos(phi)   !x
    knew(2) = st*sin(phi)   !y
    knew(3) = cos(theta)    !z
    mu = k(1)*knew(1) + k(2)*knew(2) + k(3)*knew(3) 

    ! 4/ recoil effect
    if (recoil) then
       nu_atom = nu_atom / (1.d0 + ((planck*nu_atom)/(mp*clight*clight))*(1.-mu))
    endif

    ! 5/ compute atom freq. in external frame, after scattering
    scalar = knew(1) * v(1) + knew(2) * v(2) + knew(3)* v(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + sqrt(1-mu**2)*uper)/clight)
    nu_cell = (1.d0 - scalar/clight) * nu_ext 
    k = knew


#ifdef DEBUG
    print*,'-DEBUG- scatter_HI_isotrope'
    print*,a,x_cell
    print*,'upar=',upar
    print*,'iran=',iran
    print*,'uper=',uper
    print*,'nu_atom=',nu_atom
    print*,'nu_ext =',nu_ext
    print*,'nu_cell=',nu_cell
    print*,'dopwidth=',dopwidth
    print*,'ruper, r2 =',ruper,r2
    print*,'mu, scalar =',mu,scalar
    print*,'---------------------------'
#endif

  end subroutine scatter_HI_isotrope


  subroutine scatter_HI(vcell,vth,nu_cell,k,nu_ext,iran)

    ! ---------------------------------------------------------------------------------
    ! perform scattering event on a Hydrogen atom with an anisotropic phase function
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
    !
    ! Notes on the phase function :
    ! -----------------------------
    ! - for core photons (|x| < 0.2) we use P(mu) = 11/24 + 3/24 * mu**2
    ! - for wing photons (|x| > 0.2) we use P(mu) = 3/8 * (1 + mu**2) [this is Rayleigh]
    ! with mu = cos(theta), (and theta in [0,pi]).
    ! To draw theta values, we compute the cumulative probabilities from above :
    ! - for core photons : P(< mu) = 1/2 + 11/24 * mu + 1/24 * mu**3
    ! - for wing photons : P(< mu) = 1/2 + 3/8 * mu + 1/8 * mu**3
    ! We fit the reciprocal functions of these with polynomials :
    ! - for core photons (at better than 0.515% accuracy everywhere):
    !    mu = -0.703204 x^3 + 1.054807 x^2 + 1.643182 x^1 -0.997392  
    ! - for wing photons (at better than 0.529% accuracy everywhere):
    !    mu = -24.901267 x^7 + 87.154434 x^6 -114.220525 x^5 + 67.665227 x^4 -18.389694 x^3 + 3.496531 x^2 + 1.191722 x^1 -0.998214
    ! This allows us to get mu for x a random number between 0 and 1. 
    ! ---------------------------------------------------------------------------------

    real(kind=8), intent(inout)               :: nu_cell, nu_ext
    real(kind=8), dimension(3), intent(inout) :: k
    real(kind=8), dimension(3), intent(in)    :: vcell
    real(kind=8), intent(in)                  :: vth
    integer, intent(inout)                    :: iran

    real(kind=8)               :: delta_nu_doppler, a, x_cell, blah, upar, ruper
    real(kind=8)               :: r2, uper, nu_atom, phi, theta, st, mu, scalar
    real(kind=8)               :: cti,sti,cpi,spi,ct1,st1,cp1,sp1,x,x_atom
    real(kind=8), dimension(3) :: knew

    ! define x_cell & a
    delta_nu_doppler = vth / lambda_0_cm 
    a = gamma_over_fourpi / delta_nu_doppler
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
    uper   = sqrt(-log(ruper))*cos(twopi*r2)
    uper   = uper * vth  ! from x to velocity

    ! 3/ incoming frequency in atom's frame
    nu_atom = nu_cell - nu_ext * upar/clight
    x_atom  = (nu_atom -nu_0) / delta_nu_doppler

    ! 4/ determine direction of scattered photon
    ! 4.1/ relative to incoming photon 
    phi = twopi * ran3(iran)
    x   = ran3(iran)
    if (abs(x_atom) < 0.2) then ! core scattering
       mu = ((-0.703204*x + 1.054807)* x + 1.643182) * x - 0.997392  
    else ! wing scattering
       mu = ((((((-24.901267*x + 87.154434)*x -114.220525)*x + 67.665227)*x -18.389694)*x + 3.496531)*x + 1.191722)*x -0.998214
    end if
    ! 4.2/ in external frame (box coordinates)
    ! angular description of k
    cti = k(3)
    sti = sqrt(1.d0 - cti*cti)  ! sin(theta) is positive for theta in [0,pi]. 
    if (sti > 0) then 
       cpi = k(1)/sti
       spi = k(2)/sti
    else
       cpi = 1.
       spi = 0.
    end if
    ! angular description of knew (relative to k)
    ct1 = mu
    st1 = sqrt(1.0d0 - ct1*ct1)
    cp1 = cos(phi)
    sp1 = sin(phi)
    ! knew (such that indeed knew . k = ct1)
    knew(1) = cti*cpi*st1*cp1 + sti*cpi*ct1 - spi*st1*sp1
    knew(2) = cti*spi*st1*cp1 + sti*spi*ct1 + cpi*st1*sp1
    knew(3) = -sti*st1*cp1 + cti*ct1

    ! 5/ recoil effect 
    if (recoil) then 
       nu_atom = nu_atom / (1.d0 + ((planck*nu_atom)/(mp*clight*clight))*(1.-mu))
    end if
    
    ! 6/ compute atom freq. in external frame, after scattering
    scalar = knew(1) * vcell(1) + knew(2) * vcell(2) + knew(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + sqrt(1-mu**2)*uper)/clight)
    nu_cell = (1.d0 - scalar/clight) * nu_ext 
    k = knew

  end subroutine scatter_HI

  !==================
  ! private routines 
  !==================
  
  ! peut-etre dans un module "atomic_physics_utils.f90" ?
  function voigt_fit(x,a)
    
    !use constant, only : pi,sqrtpi
    
    implicit none
    
    real(kind=8),intent(in) :: x,a
    real(kind=8)            :: voigt_fit 
    real(kind=8)            :: q,z,x2
    
    x2 = x**2
    z  = (x2 - 0.855d0) / (x2 + 3.42d0)
    if (z > 0) then 
       q = z * (1.0d0 + 21.0d0/x2) * a / pi / (x2 + 1.0d0)
       q = q * (((5.674d0*z - 9.207d0)*z + 4.421d0)*z + 0.1117)
    else
       q = 0.0d0 
    end if
    voigt_fit = q + exp(-x2) / sqrtpi
    
    return
    
  end function voigt_fit
  
  

  end module Module_HI_model
