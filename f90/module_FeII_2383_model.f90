module module_FeII_2383_model ! Fe_II UV2

  ! This module describes the absorption of photons by FeII from level  3d^6 4s 9/2 to level 3d^6 4p 11/2.
  ! This transition is at 2382.76 A. (named 2383)
  ! The module implements the only decay channels (resonant) at 2382.76 A. 

  use module_constants
  use module_utils, only : isotropic_direction
  use module_uparallel
  use module_random
  use module_voigt

  implicit none

  private

  ! Atomic data, taken from Zhu et al, 2015 (Table 2)
  ! In this module, we use the following convention :
  ! level 1 is 3d^6 4s 9/2
  ! level 3 is 3d^6 4p 11/2

  ! transition between levels 1 and 3
  real(kind=8),parameter :: lambda13       = 2382.76d0                    ! transition wavelength [A]
  real(kind=8),parameter :: lambda13_cm    = lambda13 / cmtoA             ! [cm]
  real(kind=8),parameter :: nu13           = clight / lambda13_cm         ! [Hz]
  real(kind=8),parameter :: f13            = 3.2d-1                       ! oscillator strength
  real(kind=8),parameter :: sigma13_factor = sqrtpi*e_ch**2*f13/me/clight ! multiply by Voigt(x,a)/delta_nu_doppler to get sigma.
  real(kind=8),parameter :: A31            = 3.13d8                       ! spontaneous decay [/s]

  
  public :: get_tau_FeII_2383, scatter_FeII_2383, read_FeII_2383_params, print_FeII_2383_params

contains

  function get_tau_FeII_2383(nFeII, vth, distance_to_border_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute optical depth of FeII-2382.76 over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - nFeII    : number density of FeII ions                              [ cm^-3 ]
    ! - vth      : thermal (+ small-scale turbulence) velocity of FeII ions [ cm / s ]
    ! - distance_to_border_cm : distance over which we compute tau          [ cm ]
    ! - nu_cell  : photon's frequency in the frame of the cell              [ Hz ]
    ! OUTPUT :
    ! - get_tau_FeII_2383 : optical depth of Fer's line over distance_to_border_cm
    ! --------------------------------------------------------------------------
    
    real(kind=8),intent(in) :: nFeII,vth,distance_to_border_cm,nu_cell
    real(kind=8)            :: delta_nu_doppler,x_cell,sigma,a,h_cell,get_tau_FeII_2383

    ! compute Doppler width and a-parameter
    delta_nu_doppler = vth / lambda13_cm
    a = A31 / (fourpi * delta_nu_doppler)

    ! cross section of FeII-2382.76
    x_cell = (nu_cell - nu13) / delta_nu_doppler
    h_cell = voigt_function(x_cell,a)
    sigma  = sigma13_factor / delta_nu_doppler * h_cell

    get_tau_FeII_2383 = sigma * nFeII * distance_to_border_cm
   
    return

  end function get_tau_FeII_2383

  
  subroutine scatter_FeII_2383(vcell,vth,nu_cell,k,nu_ext,iran)

    ! ---------------------------------------------------------------------------------
    ! perform scattering event on a FeII ion
    ! The photon is absorbed in transition 1->3 and may decay as 3->1 
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
    ! photon goes down to level 1 -> coherent scattering
    nu_atom = nu_cell - nu_ext * upar/clight ! incoming frequency in atom's frame = outcoming freq in same frame
    
    ! 4/ determine direction of scattered photon
    call isotropic_direction(knew,iran)
    mu = k(1)*knew(1) + k(2)*knew(2) + k(3)*knew(3)
    bu = sqrt(1.0d0 - mu*mu)
    
    ! 5/ compute atom freq. in external frame, after scattering
    scalar = knew(1) * vcell(1) + knew(2) * vcell(2) + knew(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + bu*uper)/clight)
    nu_cell = (1.d0 - scalar/clight) * nu_ext 
    k = knew

  end subroutine scatter_FeII_2383


  subroutine read_FeII_2383_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------
    
    character(*),intent(in) :: pfile
    
    call read_uparallel_params(pfile)
    call read_voigt_params(pfile)
    
    return
    
  end subroutine read_FeII_2383_params


  subroutine print_FeII_2383_params(unit)
    
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
    
  end subroutine print_FeII_2383_params
  
end module module_FeII_2383_model
