module module_SiII_1260_model

  ! This module describes the absorption of photons by SiII from level  3s^2 3p 2P^0 1/2 to level 3s^2 3d 2D 3/2.
  ! This transition is at 1260.42 A.
  ! The module also implements the two decay channels (resonant and fluorescent) at 1260.42 A and 1265.02 A. 
  
  use module_constants
  use module_utils, only : voigt_fit, isotropic_direction

  use module_uparallel
  use module_random

  implicit none

  private

  ! Atomic data, taken from Scarlata and Panagia, AjJ 801, 2015 (Table 1)
  ! In this module, we use the following convention :
  ! level 1 is 3s^2 3p 2P^0 1/2
  ! level 2 is 3s^2 3p 2P^0 3/2
  ! level 3 is 3s^2 3d 2D 3/2

  ! transition between levels 1 and 3 
  real(kind=8),parameter :: lambda13       = 1260.42d0                ! transition wavelength [A]
  real(kind=8),parameter :: lambda13_cm    = lambda13 / cmtoA         ! [cm]
  real(kind=8),parameter :: nu13           = clight / lambda13_cm     ! [Hz]
  real(kind=8),parameter :: f13            = 1.226d0                  ! oscillator strength
  real(kind=8),parameter :: sigma13_factor = pi*e_ch**2*f13/me/clight ! multiply by Voigt(x,a)/nu_D to get sigma.
  real(kind=8),parameter :: A31            = 2.57d9                   ! spontaneous decay [/s]

  ! transition between levels 2 and 3
  real(kind=8),parameter :: lambda23       = 1265.02d0                ! transition wavelength [A]
  real(kind=8),parameter :: lambda23_cm    = lambda23 / cmtoA         ! [cm]
  real(kind=8),parameter :: nu23           = clight / lambda23_cm     ! [Hz]
  real(kind=8),parameter :: A32            = 4.73d8                   ! spontaneous decay [/s]
  
!~   real(kind=8),parameter   :: gamma_over_fourpi = gamma / fourpi
  
  real(kind=8),parameter :: A31_over_A31_plus_A32 = A31 / (A31+A32)
  
  logical                  :: HI_core_skip    = .false.     ! if true, skip scatterings in the core of the line (as in Smith+15).

  public :: get_tau_SiII_1260, get_tau_SiI_photoionization, scatter_SiII_1260, read_SiII_1260_params, print_SiII_1260_params
  !--PEEL--
  public :: SiII_1260_peeloff_weight
  !--LEEP--
  
    public :: HI_core_skip 

contains

  function get_tau_SiII_1260(nSiII, vth, distance_to_border_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute optical depth of SiII-1260.42 over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - nSiII    : number density of SiII ions                              [ cm^-3 ]
    ! - vth      : thermal (+ small-scale turbulence) velocity of SiII ions [ cm / s ]
    ! - distance_to_border_cm : distance over which we compute tau          [ cm ]
    ! - nu_cell  : photon's frequency in the frame of the cell              [ Hz ]
    ! OUTPUT :
    ! - get_tau_SiII_1260 : optical depth of Silicon's line over distance_to_border_cm
    ! --------------------------------------------------------------------------
    
    real(kind=8),intent(in) :: nSiII,vth,distance_to_border_cm,nu_cell
    real(kind=8)            :: nu_D,x_cell,sigma,a,h,get_tau_SiII_1260

    ! compute Doppler width and a-parameter
    nu_D = vth / lambda13_cm
    a    = A31 / (fourpi * nu_D)

    ! cross section of SiII-1260.42
    x_cell = (nu_cell - nu13) / nu_D
    h      = voigt_fit(x_cell,a)
    sigma  = sigma13_factor / nu_D * h

    get_tau_SiII_1260 = sigma * nSiII * distance_to_border_cm
   
    return
    
  end function get_tau_SiII_1260



  !Not really well placed
  function get_tau_SiI_photoionization(nSiI, distance, nu_cell)

    real(kind=8)      :: nSiI, distance, nu_cell, get_tau_SiI_photoionization
    real(kind=8)      :: E, E0, cs0, P, ya, yw, y0, y1, x, y
    !------------------------------------------------------------------------
    E = planck * nu_cell / evtoerg      ! photon energy in ev
    if ( E .lt. 8.15168 ) then            ! below ionization energy
       get_tau_SiI_photoionization = 0d0
       RETURN
    endif

    E0 = 2.317d1    ; cs0 = 2.506d-17  ; P  = 3.546 
    ya = 20.57      ; yw  = 2.837d-1   ; y0 = 1.627d-5  ; y1 = 4.207d-1


    x = E/E0 - y0
    y = sqrt(x**2+y1**2)

    get_tau_SiI_photoionization = &
         cs0 * ((x-1.)**2 + yw**2) * y**(0.5*P-5.5)/(1.+sqrt(y/ya))**P

    get_tau_SiI_photoionization = get_tau_SiI_photoionization*nSiI*distance

    return

  end function get_tau_SiI_photoionization




  subroutine scatter_SiII_1260(vcell,vth,nu_cell,k,nu_ext,iran)

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
    real(kind=8)                            :: r2, uper, nu_atom, mu, bu, scalar
    real(kind=8),dimension(3)               :: knew

    ! define x_cell & a
    delta_nu_doppler = vth / lambda13_cm 
    a = A31 / fourpi / delta_nu_doppler
    x_cell = (nu_cell - nu13) / delta_nu_doppler

    ! 1/ component parallel to photon's propagation
    ! -> get velocity of interacting atom parallel to propagation
#ifdef SWITCH_OFF_UPARALLEL
    upar = 0.5
#else
    upar = get_uparallel(x_cell,a,iran)
#endif
    upar = upar * vth    ! upar is an x -> convert to a velocity 

    ! 2/ component perpendicular to photon's propagation
    ruper  = ran3(iran)
    r2     = ran3(iran)
    uper   = sqrt(-log(ruper))*cos(twopi*r2)
    uper   = uper * vth  ! from x to velocity

    ! 3/ chose de-excitation channel to determine output freq. in atom's frame
    r2 = ran3(iran)
    if (r2 <= A31_over_A31_plus_A32) then
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
    nu_cell = (1.d0 - scalar/clight) * nu_ext 
    k = knew

  end subroutine scatter_SiII_1260
  
  
  
  
  !--PEEL--
  function SiII_1260_peeloff_weight(vcell,vth,nu_ext,kin,kout,iran)

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
    real(kind=8)                            :: SiII_1260_peeloff_weight
    real(kind=8)                            :: delta_nu_doppler, a, x_cell, upar, ruper
    real(kind=8)                            :: r2, uper, nu_atom, mu, bu, scalar
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
#ifdef SWITCH_OFF_UPARALLEL
    upar = 0.5
#else
    upar = get_uparallel(x_cell,a,iran)
#endif
    upar = upar * vth    ! upar is an x -> convert to a velocity 

    ! 2/ component perpendicular to photon's propagation
    ruper  = ran3(iran)
    r2     = ran3(iran)
    uper   = sqrt(-log(ruper))*cos(twopi*r2)
    uper   = uper * vth  ! from x to velocity

    ! 3/ chose de-excitation channel to determine output freq. in atom's frame
    r2 = ran3(iran)
    if (r2 <= A31_over_A31_plus_A32) then
       ! photon goes down to level 1 -> coherent scattering
       nu_atom = nu_cell - nu_ext * upar/clight ! incoming frequency in atom's frame = outcoming freq in same frame
    else
       ! photons goes down to level two ...
       nu_atom = nu23 
    end if
    
    ! 4/ determine direction of scattered photon
    SiII_1260_peeloff_weight = 0.5d0  ! P(mu) for isotropic phase function
    mu = kin(1)*kout(1) + kin(2)*kout(2) + kin(3)*kout(3)
    bu = sqrt(1.0d0 - mu*mu)
    
    ! 5/ compute atom freq. in external frame, after scattering
    scalar = kout(1) * vcell(1) + kout(2) * vcell(2) + kout(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + bu*uper)/clight)

  end function SiII_1260_peeloff_weight
!--LEEP--
  
  

  
  subroutine read_SiII_1260_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile

    call read_uparallel_params(pfile)

    return

  end subroutine read_SiII_1260_params


    subroutine print_SiII_1260_params(unit)
    
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
    
  end subroutine print_SiII_1260_params


end module module_SiII_1260_model