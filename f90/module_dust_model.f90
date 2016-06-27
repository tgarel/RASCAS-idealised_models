module module_dust_model

  implicit none

  ! Determine dust albedo at Lya wavelength (dust_albedo) and
  ! g parameter (g_dust) of the Henyey-Greenstein phase function for dust scattering:
  !    real(kind=8),parameter :: dust_albedo=0.5 ! Adopted in Verhamme et al. (2006)
  !    real(kind=8),parameter :: g_dust=0.001    ! ~isotropic scattering (default)
  !    real(kind=8),parameter :: dust_albedo=0. ! for tests without scattering
  !    real(kind=8),parameter :: dust_albedo=0.33 ! Values from Draine 2003 for MW dust, R_V=3.1
  !    real(kind=8),parameter :: g_dust=0.68      ! Single parameter of Henyey-Greenstein function (value from Draine 2003 for MW dust, R_V=3.1)
  real(kind=8),parameter :: dust_albedo=0.46 ! Empirical values from Witt & Gordon 2000 (ApJ, 528, 799) for SMC dust
  real(kind=8),parameter :: g_dust=0.77      ! 
  !    real(kind=8),parameter :: dust_albedo=0.40 ! Values from Gordon et al. 1997 (ApJ, 487, 625) for SMC dust
  !    real(kind=8),parameter :: g_dust=0.53      !   (from model calculations)

  real(kind=8),parameter  :: rd= 2.d-6       ! radius of dust grain [cm]
  real(kind=8),parameter  :: mpmd = 5.d-8    ! proton over dust mass
  real(kind=8),parameter  :: sigmad = (pi * rd**2)/(1.-dust_albedo)  ! total abs+scattering dust cross section for general value of albedo
  real(kind=8),parameter  :: sigmaH_factor = pi*e_ch**2*f12/ me / clight ! H cross-section factor-> multiply by Voigt(x,a)/nu_D to get sigma.
  real(kind=8)            :: sigmad_gal ! "effective" dust abs cross section from Draine & Lee (1984, Fig. 7) at Ly-a wavelength
  real(kind=8)            :: albedo ! proba of interacting with H


  contains


    function get_tau_dust(ndust_cell, distance_to_border_cm)

      real(kind=8),intent(in) :: ndust_cell,distance_to_border_cm

      get_tau_dust = sigmad * ndust_cell * distance_to_border_cm

    end function get_tau_dust



    subroutine scatter_dust(p)

      type(photon),intent(inout) :: p

      ! interaction with dust

      if (i_diag.gt.0) print*, "the photon interacts with dust !!!"
      aors = ran3(iran)  ! aka "Absorption OR Scattering" ... 
      if (aors.gt.dust_albedo) then ! the photon is absorbed = lost
         iescape = 0
         return
      else                   ! the photon is scattered by dust
         !........dust angular redistribution
         xxx = 0.
         if (dipol.eq.4) then
            call emission_direction(xxx) ! indicates Henyey-Greenstein fct
         else
            dip=dipol
            dipol=2 ! indicates isotropic (for dust)
            call emission_direction(xxx)
            dipol=dip
         endif
         scalar =  a0 * vx_cell + b0 * vy_cell + c0 * vz_cell
         nu_ext = nu_int/( 1. - scalar/clight)

      end if

    end subroutine scatter_dust
    
  end module module_dust_model
