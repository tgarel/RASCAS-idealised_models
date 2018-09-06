module module_SiII_1190_model

  ! This module describes the absorption of photons by SiII from level  3s^2 3p 2P^0 1/2 to level 3s 3p^2 2P 3/2.
  ! This transition is at 1190.42 A.
  ! The module also implements the two decay channels (resonant and fluorescent) at 1190.42 A and 1194.50 A. 

  use module_constants
  use module_utils, only : isotropic_direction
  use module_uparallel
  use module_random
  use module_voigt

  implicit none

  private

  ! Atomic data, taken from Scarlata and Panagia, AjJ 801, 2015 (Table 1)
  ! In this module, we use the following convention :
  ! level 1 is 3s^2 3p 2P^0 1/2
  ! level 2 is 3s^2 3p 2P^0 3/2
  ! level 4 is 3s 3p^2 2P 3/2

  ! transition between levels 1 and 4
  real(kind=8),parameter :: lambda14       = 1190.42d0                    ! transition wavelength [A]
  real(kind=8),parameter :: lambda14_cm    = lambda14 / cmtoA             ! [cm]
  real(kind=8),parameter :: nu14           = clight / lambda14_cm         ! [Hz]
  real(kind=8),parameter :: f14            = 0.277d0                      ! oscillator strength
  real(kind=8),parameter :: sigma14_factor = sqrtpi*e_ch**2*f14/me/clight ! multiply by Voigt(x,a)/delta_nu_doppler to get sigma.
  real(kind=8),parameter :: A41            = 6.53d8                       ! spontaneous decay [/s]

  ! transition between levels 2 and 4
  real(kind=8),parameter :: lambda24       = 1194.50d0                    ! transition wavelength [A]
  real(kind=8),parameter :: lambda24_cm    = lambda24 / cmtoA             ! [cm]
  real(kind=8),parameter :: nu24           = clight / lambda24_cm         ! [Hz]
  real(kind=8),parameter :: A42            = 3.45d9                       ! spontaneous decay [/s]

  real(kind=8),parameter :: A41_over_A41_plus_A42 = A41 / (A41+A42)
  
  public :: get_tau_SiII_1190, scatter_SiII_1190, read_SiII_1190_params, print_SiII_1190_params

contains

  function get_tau_SiII_1190(nSiII, vth, distance_to_border_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute optical depth of SiII-1190.42 over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - nSiII    : number density of SiII ions                              [ cm^-3 ]
    ! - vth      : thermal (+ small-scale turbulence) velocity of SiII ions [ cm / s ]
    ! - distance_to_border_cm : distance over which we compute tau          [ cm ]
    ! - nu_cell  : photon's frequency in the frame of the cell              [ Hz ]
    ! OUTPUT :
    ! - get_tau_SiII_1190 : optical depth of Silicon's line over distance_to_border_cm
    ! --------------------------------------------------------------------------
    
    real(kind=8),intent(in) :: nSiII,vth,distance_to_border_cm,nu_cell
    real(kind=8)            :: delta_nu_doppler,x_cell,sigma,a,h_cell,get_tau_SiII_1190

    ! compute Doppler width and a-parameter
    delta_nu_doppler = vth / lambda14_cm
    a = A41 / (fourpi *  delta_nu_doppler)

    ! cross section of SiII-1193.28
    x_cell = (nu_cell - nu14) / delta_nu_doppler
    h_cell = voigt_function(x_cell,a)
    sigma  = sigma14_factor / delta_nu_doppler * h_cell

    get_tau_SiII_1190 = sigma * nSiII * distance_to_border_cm
   
    return

  end function get_tau_SiII_1190

  
  subroutine scatter_SiII_1190(vcell,vth,nu_cell,k,nu_ext,iran)

    ! ---------------------------------------------------------------------------------
    ! perform scattering event on a SiII ion
    ! The photon is absorbed in transition 1->4 and may decay as 4->1 or 4->2. 
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
    if (r2 <= A41_over_A41_plus_A42) then
       ! photon goes down to level 1 -> coherent scattering
       nu_atom = nu_cell - nu_ext * upar/clight ! incoming frequency in atom's frame = outcoming freq in same frame
    else
       ! photons goes down to level two ...
       nu_atom = nu24 
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

  end subroutine scatter_SiII_1190


  subroutine read_SiII_1190_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile

    call read_uparallel_params(pfile)
    call read_voigt_params(pfile)

    return

  end subroutine read_SiII_1190_params


    subroutine print_SiII_1190_params(unit)
    
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
    
  end subroutine print_SiII_1190_params


end module module_SiII_1190_model
