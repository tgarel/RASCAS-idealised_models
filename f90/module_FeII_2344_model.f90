module module_FeII_2344_model ! Fe_II UV3

  ! This module describes the absorption of photons by FeII from level  3d^6 4s 9/2 to level 3d^6 4p 7/2.
  ! This transition is at 2344.21 A. (named 2344)
  ! The module also implements the two decay channels (resonant and fluorescent) at 2344.21 A , 2365.55 A. and 2381.49 

  use module_constants
  use module_utils, only : isotropic_direction
  use module_uparallel
  use module_random
  use module_voigt

  implicit none

  private

  ! Atomic data, taken from Zhu et al, , 2015 (Table 2)
  ! In this module, we use the following convention :
  ! level 1 is 3d^6 4s 9/2
  ! level 2 is 3d^6 4s 7/2
  ! level 3 is 3d^6 4s 5/2
  ! level 4 is 3d^6 4p 7/2

  ! transition between levels 1 and 4
  real(kind=8),parameter :: lambda14       = 2344.21d0                    ! transition wavelength [A]
  real(kind=8),parameter :: lambda14_cm    = lambda14 / cmtoA             ! [cm]
  real(kind=8),parameter :: nu14           = clight / lambda14_cm         ! [Hz]
  real(kind=8),parameter :: f14            = 1.14d-1                      ! oscillator strength
  real(kind=8),parameter :: sigma14_factor = sqrtpi*e_ch**2*f14/me/clight ! multiply by Voigt(x,a)/delta_nu_doppler to get sigma.
  real(kind=8),parameter :: A41            = 1.73d8                       ! spontaneous decay [/s]

  ! transition between levels 4 and 2
  real(kind=8),parameter :: lambda24       = 2365.55d0                    ! transition wavelength [A]
  real(kind=8),parameter :: lambda24_cm    = lambda24 / cmtoA             ! [cm]
  real(kind=8),parameter :: nu24           = clight / lambda24_cm         ! [Hz]
  real(kind=8),parameter :: A42            = 5.9d7                        ! spontaneous decay [/s]

  ! transition between levels 4 and 3
  real(kind=8),parameter :: lambda34       = 2381.49d0                    ! transition wavelength [A]
  real(kind=8),parameter :: lambda34_cm    = lambda34 / cmtoA             ! [cm]
  real(kind=8),parameter :: nu34           = clight / lambda34_cm         ! [Hz]
  real(kind=8),parameter :: A43            = 3.1d7                        ! spontaneous decay [/s]
  
  real(kind=8),parameter :: Atot = A41+A42+A43
  
  public :: get_tau_FeII_2344, scatter_FeII_2344, read_FeII_2344_params, print_FeII_2344_params

contains

  function get_tau_FeII_2344(nFeII, vth, distance_to_border_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute optical depth of FeII-2344.21 over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - nFeII    : number density of FeII ions                              [ cm^-3 ]
    ! - vth      : thermal (+ small-scale turbulence) velocity of FeII ions [ cm / s ]
    ! - distance_to_border_cm : distance over which we compute tau          [ cm ]
    ! - nu_cell  : photon's frequency in the frame of the cell              [ Hz ]
    ! OUTPUT :
    ! - get_tau_FeII_2344 : optical depth of Fer's line over distance_to_border_cm
    ! --------------------------------------------------------------------------
    
    real(kind=8),intent(in) :: nFeII,vth,distance_to_border_cm,nu_cell
    real(kind=8)            :: delta_nu_doppler,x_cell,sigma,a,h_cell,get_tau_FeII_2344

    ! compute Doppler width and a-parameter
    delta_nu_doppler = vth / lambda14_cm
    a = A41 / (fourpi * delta_nu_doppler)

    ! cross section of FeII-2344
    x_cell = (nu_cell - nu14) / delta_nu_doppler
    h_cell = voigt_function(x_cell,a)
    sigma  = sigma14_factor / delta_nu_doppler * h_cell

    get_tau_FeII_2344 = sigma * nFeII * distance_to_border_cm
   
    return

  end function get_tau_FeII_2344

  
  subroutine scatter_FeII_2344(vcell,vth,nu_cell,k,nu_ext,iran)

    ! ---------------------------------------------------------------------------------
    ! perform scattering event on a FeII ion
    ! The photon is absorbed in transition 1->4 and may decay as 4->1 or 4->2 or 4->3. 
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
    real(kind=8)                            :: r2, uper, nu_atom, mu, bu, scalar, proba41, proba42
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
    proba41 = A41/Atot
    proba42 = proba41 + A42/Atot
    if (r2 <= proba41) then
       ! photon goes down to level 1 -> resonant scattering
       nu_atom = nu_cell - nu_ext * upar/clight ! incoming frequency in atom's frame = outcoming freq in same frame
    else if (r2 <= proba42) then
       ! photons goes down to level 2 ... -> fluorecent 2
       nu_atom = nu24
    else
       ! photons goes down to level 3 ... -> fluorecent 3
       nu_atom = nu34
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

  end subroutine scatter_FeII_2344


  subroutine read_FeII_2344_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------
    
    character(*),intent(in) :: pfile
    
    call read_uparallel_params(pfile)
    call read_voigt_params(pfile)
    
    return
    
  end subroutine read_FeII_2344_params


  subroutine print_FeII_2344_params(unit)
    
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
    
  end subroutine print_FeII_2344_params
  
end module module_FeII_2344_model
