module module_dust_model

  implicit none

  ! Determine dust albedo at Lya wavelength (dust_albedo) and
  ! g parameter (g_dust) of the Henyey-Greenstein phase function for dust scattering:
  ! g_dust=0.68 &  dust_albedo=0.33 from Draine 2003 for MW dust, R_V=3.1
  ! dust_albedo=0.40 from Gordon et al. 1997 (ApJ, 487, 625) for SMC dust
  real(kind=8),parameter :: dust_albedo=0.46                         ! Empirical values from Witt & Gordon 2000 (ApJ, 528, 799) for SMC dust
  real(kind=8),parameter :: g_dust=0.73                              ! Laursen+09
  real(kind=8),parameter  :: rd= 2.d-6                               ! radius of dust grain [cm]
  real(kind=8),parameter  :: mpmd = 5.d-8                            ! proton over dust mass
  real(kind=8),parameter  :: sigmad = (pi * rd**2)/(1.-dust_albedo)  ! total abs+scattering dust cross section for general value of albedo
  real(kind=8)            :: albedo                                  ! proba of interacting with H


  contains


    function get_tau_dust(ndust_cell, distance_to_border_cm)

      real(kind=8),intent(in) :: ndust_cell,distance_to_border_cm

      get_tau_dust = sigmad * ndust_cell * distance_to_border_cm

    end function get_tau_dust


    
    subroutine scatter_dust(v,ndust,nu_cell,k,nu_ext,iran)

      !type(photon),intent(inout) :: p
      !TIBO
      ! ndust argument useless here
      real(kind=8), intent(inout)               :: nu_cell, nu_ext ! nu_cell in RASCAS = nu_int in MCLya
      real(kind=8), dimension(3), intent(inout) :: k
      integer, intent(inout)                    :: iran
      real(kind=8)                              :: phi, theta, mu, aors, f, scalar
      logical                                   :: ok
      real(kind=8), dimension(3)                :: knew
      real(kind=8)                              :: nu_doppler, a, x_cell, st
      real(kind=8), intent(in)                  :: dopwidth, ndust
      !OBIT
      
      ! interaction with dust

      if (i_diag.gt.0) print*, "the photon interacts with dust !!!"
      aors = ran3(iran)  ! aka "Absorption OR Scattering" ... 
      if (aors.gt.dust_albedo) then ! the photon is absorbed = lost
         iescape = 0
         return
         !TIBO
!!$      else                   ! the photon is scattered by dust
!!$         !........dust angular redistribution
!!$         xxx = 0.
!!$         if (dipol.eq.4) then
!!$            call emission_direction(xxx) ! indicates Henyey-Greenstein fct
!!$         else
!!$            dip=dipol
!!$            dipol=2 ! indicates isotropic (for dust)
!!$            call emission_direction(xxx)
!!$            dipol=dip
!!$         endif
!!$         scalar =  a0 * vx_cell + b0 * vy_cell + c0 * vz_cell
!!$         nu_ext = nu_int/( 1. - scalar/clight)
!!$
!!$      end if
      else
         ! 1/ determine scattering angle (in atom's frame)
         
         ra=ran3(iran)
         if (dipol.eq.4) then  ! use White 79 approximation for Henyey-Greenstein phase fct:
            mu = (1.+g_dust*g_dust-((1.-g_dust*g_dust)/(1.-g_dust+2.*g_dust*ra))**2)/(2.*g_dust)
         end if
         
!!$         ! Or rejection method...
!!$         ! Divide by (1.-g_dust**2)*(1.+g_dust**2-2.*g_dust)^(-1.5) such that max(f)=1
!!$         mu = 2.0d0*ran3(iran)-1.0d0
!!$         if (dipol.eq.4) then  ! use Henyey-Greenstein phase fct only ?        
!!$            ok = .false.
!!$            do while (.not. ok)
!!$               f = (1.-g_dust**2)*(1.+g_dust**2-2.*g_dust*mu)**(-1.5) / ((1.-g_dust**2)*(1.+g_dust**2-2.*g_dust)^(-1.5))
!!$               ra = ran3(iran) 
!!$               if (ra <= f) ok = .true.
!!$            end do
!!$         end if

         theta = acos(mu)
         phi   = 2.d0*pi*ran3(iran)

         !........director cosines
         st = sin(theta)
         knew(1) = st*cos(phi)   !x
         knew(2) = st*sin(phi)   !y
         knew(3) = cos(theta)    !z
         mu = k(1)*knew(1) + k(2)*knew(2) + k(3)*knew(3) 
                  
         ! 2/ compute atom freq. in external frame, after scattering
         scalar = knew(1) * v(1) + knew(2) * v(2) + knew(3)* v(3)
         ! nu_cell has not changed; we implicitly assume that the interacting dust grain is at rest in cell's frame
         nu_ext = nu_cell/(1. - scalar/clight)
         k = knew

      end if
      !OBIT
      
    end subroutine scatter_dust
    
  end module module_dust_model
