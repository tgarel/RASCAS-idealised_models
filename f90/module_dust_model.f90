module module_dust_model

  use module_random
  use module_constants, only : pi, clight
  use module_utils, only : anisotropic_direction_Dust
  
  implicit none

  ! Determine dust albedo at Lya wavelength (dust_albedo) and
  ! g parameter (g_dust) of the Henyey-Greenstein phase function for dust scattering:
  ! g_dust=0.68 &  dust_albedo=0.33 from Draine 2003 for MW dust, R_V=3.1
  ! dust_albedo=0.40 from Gordon et al. 1997 (ApJ, 487, 625) for SMC dust
  real(kind=8),parameter  :: dust_albedo=0.46                         ! Empirical values from Witt & Gordon 2000 (ApJ, 528, 799) for SMC dust
  real(kind=8),parameter  :: g_dust=0.73                              ! Laursen+09
  real(kind=8),parameter  :: rd= 2.d-6                                ! radius of dust grain [cm]
  real(kind=8),parameter  :: mpmd = 5.d-8                             ! proton over dust mass
  real(kind=8),parameter  :: sigmad = (pi * rd**2)/(1.-dust_albedo)   ! total abs+scattering dust cross section for general value of albedo

  contains


    function get_tau_dust(ndust_cell, distance_to_border_cm)
      
      real(kind=8),intent(in) :: ndust_cell,distance_to_border_cm
      real(kind=8)            :: get_tau_dust
      get_tau_dust = sigmad * ndust_cell * distance_to_border_cm

      return
      
    end function get_tau_dust

    
    subroutine scatter_dust(v,nu_cell,k,nu_ext,iran)

      implicit none 
      
      ! ---------------------------------------------------------------------------------
      ! perform scattering event on a Hydrogen atom with isotrope angular redistribution
      ! ---------------------------------------------------------------------------------
      ! INPUTS :
      ! - v        : bulk velocity of the gas (i.e. cell velocity)       [ cm / s ] 
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
      
      real(kind=8), intent(inout)               :: nu_cell, nu_ext ! nu_cell in RASCAS = nu_int in MCLya
      real(kind=8), dimension(3), intent(inout) :: k
      real(kind=8), dimension(3), intent(in)    :: v
      integer, intent(inout)                    :: iran
      real(kind=8)                              :: mu, aors, scalar
      logical                                   :: ok
      real(kind=8), dimension(3)                :: knew
      integer(kind=4)                           :: iescape

#ifdef DEBUG
      print *,'-DEBUG- scatter_dust routine, arguments =',v,nu_cell,k,nu_ext,iran
#endif
    
      ! interaction with dust
      aors = ran3(iran)  ! aka "Absorption OR Scattering" ... 
      if (aors.gt.dust_albedo) then ! the photon is absorbed = lost
         iescape = 0
         return
      else
         call anisotropic_direction_Dust(k,knew,mu,iran,g_dust)                         
         ! compute atom freq. in external frame, after scattering
         scalar = knew(1) * v(1) + knew(2) * v(2) + knew(3)* v(3)
         ! nu_cell has not changed; we implicitly assume that the interacting dust grain is at rest in cell's frame
         nu_ext = nu_cell/(1. - scalar/clight)
         k = knew
      end if
      
    end subroutine scatter_dust
    
  end module module_dust_model
