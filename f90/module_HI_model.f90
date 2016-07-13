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

  real(kind=8),parameter   :: lambda_0=1215.6701                          ![A] Lya wavelength
  real(kind=8),parameter   :: gamma=6.265d8                               ! Einstein coeff = damping constant for Voigt Function(gamma_alpha)  
  real(kind=8),parameter   :: f12=0.416                                   ! oscillator strength for Ly-alpha
  real(kind=8),parameter   :: lambda_0_cm = lambda_0 / cmtoA              ! cm
  real(kind=8),parameter   :: nu_0 = clight / lambda_0_cm                 ! Hz
  real(kind=8),parameter   :: sigmaH_factor = pi*e_ch**2*f12/ me / clight ! H cross-section factor-> multiply by Voigt(x,a)/nu_D to get sigma.

  private :: voigt_fit

  public :: get_tau_HI, scatter_HI_isotrope
  
contains

  function get_tau_HI(nhi, dopwidth, distance_to_border_cm, nu_cell)

    ! compute frequency in cell's moving frame
    !scalar =  a0 * vx_cell + b0 * vy_cell + c0 * vz_cell
    !nu_int = (1.d0 - scalar/clight) * nu_ext
    ! => should be done before

    real(kind=8),intent(in) :: nhi,dopwidth,distance_to_border_cm,nu_cell
    real(kind=8)            :: nu_D,x_cell,sigmaH,a,h, get_tau_HI

    ! compute Doppler width and a-parameter, for H 
    nu_D = dopwidth * nu_0 / clight
    a    = gamma / (4.d0 * pi * nu_D)
 
    ! Cross section of H 
    x_cell  = (nu_cell - nu_0)/nu_D
    h       = voigt_fit(x_cell,a)
    sigmaH  = sigmaH_factor / nu_D * h
 
    get_tau_HI = sigmaH * nhi * distance_to_border_cm

    return

  end function get_tau_HI


  subroutine scatter_HI_isotrope(v,nHI,dopwidth, nu_cell, k, nu_ext, iran)

    real(kind=8), intent(inout)               :: nu_cell, nu_ext
    real(kind=8), dimension(3), intent(inout) :: k
    real(kind=8), dimension(3), intent(in)    :: v
    real(kind=8), intent(in)                  :: nHI,dopwidth
    integer, intent(inout)                    :: iran

    real(kind=8)               :: nu_doppler, a, x_cell, blah, upar, ruper
    real(kind=8)               :: r2, uper, nu_atom, phi, theta, st, mu, scalar
    real(kind=8), dimension(3) :: knew

#ifdef DEBUG
    print *,'-DEBUG- scatter routine, arguments =',v,nHI,dopwidth, nu_cell, k, nu_ext, iran
#endif

    ! define x_cell & a
    nu_doppler = dopwidth * nu_0 / clight
    a = gamma / (4.d0 * pi * nu_doppler)
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
    uper   = sqrt(-log(ruper))*cos(2.d0*pi*r2)
    uper   = uper * dopwidth  ! from x to velocity

    ! 3/ determine scattering angle (in atom's frame)
    nu_atom = nu_cell - nu_ext * upar/clight

    phi   = 2.d0*pi*ran3(iran)
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

  


!     subroutine scatter_HI(v,nHI,dopwidth, nu_cell, k, nu_ext)
      
!       real(kind=8), intent(inout)               :: nu_cell, nu_ext
!       real(kind=8), dimension(3), intent(inout) :: k
!       real(kind=8), dimension(3), intent(in)    :: v
!       real(kind=8), intent(in)                  :: nHI,dopwidth

!       ! nu_cell = donnee entree
!       ! nu_ext aussi mais on s'en fout

!       ! define x_cell & a
!       nu_dopppler = dopwidth * nu_0 / clight
!       a = gamma / (4.d0 * pi * nu_doppler)
!       x_cell = (nu_cell - nu_0) / nu_doppler

!       ! 1/ component parallel to photon's propagation
!       call uparallel(x_cell,a,iran,upar)  ! -> get velocity of interacting atom parallel to propagation
!       upar = upar * dopwidth    ! upar is an x -> convert to a velocity 
   
!       ! 2/ component perpendicular to photon's propagation
!       ruper  = ran3(iran)
!       r2     = ran3(iran)
!       uper   = sqrt(-log(ruper))*cos(2.d0*pi*r2)
!       uper   = uper * dopwidth  ! from x to velocity

!       ! 3/ determine scattering angle (in atom's frame)
!       nu_atom = nu_cell - nu_ext * upar/clight
!       x_atom = (nu_atom -nu_0)/nu_doppler

!       !----------------------------------------------------------------------------------------------------
!       ! inline subroutine emission_direction(xint=x_atom)
!       ! dipole_model = (1 == dipole, 2 == isotropic, 3 == QM, 4 == QM for Lya + Henyey-Greenstein for dust)
!       ! if (dipol.eq.4) then do dipol=3 here since dust stuff will bw done in dust_model...
!       !----------------------------------------------------------------------------------------------------
!       select case(dipole_model)
      
!       case(1) ! dipolar angular distribution, no polarization

!          thetai=theta
!          phii=phi
!          !... step 1 : we generate the emergent angles theta1, phi1, 
!          !             but in the frame where the incident direction ki
!          !             is along the z axis.
!          !..........We generate mu=cos(theta1), from the dipolar distribution
!          !          theta1 is the angle between the incident and emergent directions 
!          ra=ran3(iran)
!          ! distribution of mu (=cos theta) is \propto (1+cos^2 theta)
!          ok = .false.
!          do while (.not. ok)
!             mu = 2.0d0*ran3(iran)-1.0d0
!             f  = 0.5d0 * (1.0d0 + mu**2)  ! 0.5 factor is a normalization so that max of f is one.
!             ra = ran3(iran)
!             if (ra <= f) ok = .true.
!          end do
!          !..........We deduce theta1 from mu and generate phi1 randomly 
!          theta1=acos(mu)
!          rphi=ran3(iran)
!          phi1=2.d0*pi*rphi
!          !.... step 2 : we transform theta1 and phi1 in theta_e and phi_e 
!          !              in the initial frame.
!          !...........We perform a rotation of phii around the z axis,
!          !           and a rotation of thetai around the y axis.
!          cti=cos(thetai)   ;   sti=sin(thetai)
!          ct1=cos(theta1)   ;   st1=sin(theta1)
!          cpi=cos(phii)     ;   spi=sin(phii)
!          cp1=cos(phi1)     ;   sp1=sin(phi1)

!          a0 = cti*cpi*st1*cp1 + sti*cpi*ct1 - spi*st1*sp1
!          b0 = cti*spi*st1*cp1 + sti*spi*ct1 + cpi*st1*sp1
!          c0 = -sti*st1*cp1 + cti*ct1

!          theta = acos(c0)
!          phi   = atan2(b0,a0)
!          if (phi .lt. 0.d0) phi = phi + 2.d0 * pi  !! change from [-pi;pi] to [0,2pi]
!          !----------------------------------------------------------------------------------------------------
      
!       case(2) ! isotropic angular distribution
!          a0i = a0
!          b0i = b0
!          c0i = c0
!          rphi   = ran3(iran)
!          rtheta = ran3(iran)     
!          phi    = 2.d0*pi*rphi
!          theta  = acos(1d0-2d0*rtheta)
!          !........director cosines
!          st = sin(theta)
!          a0 = st*cos(phi)   !x
!          b0 = st*sin(phi)   !y
!          c0 = cos(theta)    !z
!          mu = a0i*a0 + b0i*b0 + c0i*c0 
!          !----------------------------------------------------------------------------------------------------

!       case(3,4) ! QM redistribution 

!          thetai=theta
!          phii=phi
!          !...........dipolar angular distribution
!          !.... step 1 : we generate the emergent angles theta1, phi1, 
!          !              but in the frame where the incident direction ki
!          !              is along the z axis.
!          !..........We generate mu=cos(theta1), from the dipolar distribution
!          !          theta1 is the angle between the incident and emergent directions 
!          ra=ran3(iran)
!          if (abs(x_atom).lt.0.2) then !core scattering
!             ! distribution of mu is \propto (11 + 3*cos^2 theta)
!             ok = .false.
!             do while (.not. ok)
!                mu = 2.0d0*ran3(iran)-1.0d0
!                f  = (11.0d0 + 3.0d0 * mu**2)/14.0d0  ! /14 factor is a normalization so that max of f is one.
!                ra = ran3(iran) 
!                if (ra <= f) ok = .true.
!             end do
!          else
!             ! distribution of mu (=cos theta) is \propto (1+cos^2 theta)
!             ok = .false.
!             do while (.not. ok)
!                mu = 2.0d0*ran3(iran)-1.0d0
!                f  = 0.5d0 * (1.0d0 + mu**2)  ! 0.5 factor is a normalization so that max of f is one.
!                ra = ran3(iran) 
!                if (ra <= f) ok = .true.
!             end do
!          endif
!          !..........We deduce theta1 from mu and generate phi1 randomly 
!          theta1=acos(mu)
!          rphi=ran3(iran)
!          phi1=2.d0*pi*rphi
!          !.... step 2 : we transform theta1 and phi1 in theta_e and phi_e 
!          !              in the initial frame.
!          !...........We perform a rotation of phii around the z axis,
!          !           and a rotation of thetai around the y axis.
!          cti=cos(thetai)   ;   sti=sin(thetai)
!          ct1=cos(theta1)   ;   st1=sin(theta1)
!          cpi=cos(phii)     ;   spi=sin(phii)
!          cp1=cos(phi1)     ;   sp1=sin(phi1)
!          a0 = cti*cpi*st1*cp1 + sti*cpi*ct1 - spi*st1*sp1
!          b0 = cti*spi*st1*cp1 + sti*spi*ct1 + cpi*st1*sp1
!          c0 = -sti*st1*cp1 + cti*ct1
!          theta = acos(c0)
!          phi   = atan2(b0,a0)
!          if (phi .lt. 0.d0) phi = phi + 2.d0 * pi  !! change from [-pi;pi] to [0,2pi]
!          !----------------------------------------------------------------------------------------------------

!       end select


!       ! 4/ treat recoil effect (in atom's frame):
!       if (recoil) then
!          nu_atom = nu_atom / (1.d0 + ((planck*nu_atom)/(mp*clight*clight))*(1.-mu))
!       endif

!       ! 5/ compute atom freq. in external frame, after scattering
!       scalar =  a0 * vx_cell + b0 * vy_cell + c0 * vz_cell
!       nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + sqrt(1-mu**2)*uper)/clight)

!       ! return nu_cell or nu_ext???
!       ! what about a0,b0,c0,theta,phi ???
!       ! k?

!     end subroutine scatter_HI





!     subroutine scatter_HI_polarization()

!       !!!type(photon),intent(inout) :: p

!       ! 1/ component parallel to photon's propagation
!       call uparallel(x_int,a,iran,upar)  ! -> get velocity of interacting atom parallel to propagation
!       upar = upar * b_cell  ! upar is an x -> convert to a velocity 
   
!       ! 2/ component perpendicular to photon's propagation
!       ruper  = ran3(iran)
!       r2     = ran3(iran)
!       uper   = sqrt(-log(ruper))*cos(2.d0*pi*r2)
!       uper   = uper * b_cell  ! from x to velocity

!       ! 3/ determine scattering angle (in atom's frame)
!       nu_atom = nu_int - nu_ext * upar/clight
!       x_atom = (nu_atom -nu_0)/nu_D

!       !----------------------------------------------------------------------------------------------------
!       ! inline subroutine emission_direction(xint=x_atom)
!       ! dipole_model = (1 == dipole, 2 == isotropic, 3 == QM, 4 == QM for Lya + Henyey-Greenstein for dust)
!       ! if (dipol.eq.4) then do dipol=3 here since dust stuff will bw done in dust_model...
!       select case(dipole_model)
!          !----------------------------------------------------------------------------------------------------
!       case(1) ! dipolar angular distribution, no polarization
! #ifdef POLARIZATION
!          ! dipole scattering : pure Rayleigh
!          a0i = a0
!          b0i = b0
!          c0i = c0
!          ok = .false.
!          do while (.not. ok)
!             rphi   = ran3(iran)
!             rtheta = ran3(iran)     
!             phi    = 2.d0*pi*rphi
!             theta  = acos(1d0-2d0*rtheta)
!             !........director cosines
!             st = sin(theta)
!             a0 = st*cos(phi)   !x
!             b0 = st*sin(phi)   !y
!             c0 = cos(theta)    !z
!             f = 1-(a0*x_pol + b0*y_pol + c0*z_pol)**2 ! 1-(e.n')**2
!             ra = ran3(iran)
!             if (ra <=f) ok = .true.
!          end do
!          mu = a0i*a0 + b0i*b0 + c0i*c0
!          call projection_pol
! #else
!          thetai=theta
!          phii=phi
!          !... step 1 : we generate the emergent angles theta1, phi1, 
!          !             but in the frame where the incident direction ki
!          !             is along the z axis.
!          !..........We generate mu=cos(theta1), from the dipolar distribution
!          !          theta1 is the angle between the incident and emergent directions 
!          ra=ran3(iran)
!          ! distribution of mu (=cos theta) is \propto (1+cos^2 theta)
!          ok = .false.
!          do while (.not. ok)
!             mu = 2.0d0*ran3(iran)-1.0d0
!             f  = 0.5d0 * (1.0d0 + mu**2)  ! 0.5 factor is a normalization so that max of f is one.
!             ra = ran3(iran)
!             if (ra <= f) ok = .true.
!          end do
!          !..........We deduce theta1 from mu and generate phi1 randomly 
!          theta1=acos(mu)
!          rphi=ran3(iran)
!          phi1=2.d0*pi*rphi
!          !.... step 2 : we transform theta1 and phi1 in theta_e and phi_e 
!          !              in the initial frame.
!          !...........We perform a rotation of phii around the z axis,
!          !           and a rotation of thetai around the y axis.
!          cti=cos(thetai)   ;   sti=sin(thetai)
!          ct1=cos(theta1)   ;   st1=sin(theta1)
!          cpi=cos(phii)     ;   spi=sin(phii)
!          cp1=cos(phi1)     ;   sp1=sin(phi1)

!          a0 = cti*cpi*st1*cp1 + sti*cpi*ct1 - spi*st1*sp1
!          b0 = cti*spi*st1*cp1 + sti*spi*ct1 + cpi*st1*sp1
!          c0 = -sti*st1*cp1 + cti*ct1

!          theta = acos(c0)
!          phi   = atan2(b0,a0)
!          if (phi .lt. 0.d0) phi = phi + 2.d0 * pi  !! change from [-pi;pi] to [0,2pi]
! #endif
!          !----------------------------------------------------------------------------------------------------
!       case(2) ! isotropic angular distribution
!          a0i = a0
!          b0i = b0
!          c0i = c0
!          rphi   = ran3(iran)
!          rtheta = ran3(iran)     
!          phi    = 2.d0*pi*rphi
!          theta  = acos(1d0-2d0*rtheta)
!          !........director cosines
!          st = sin(theta)
!          a0 = st*cos(phi)   !x
!          b0 = st*sin(phi)   !y
!          c0 = cos(theta)    !z
!          mu = a0i*a0 + b0i*b0 + c0i*c0 
! #ifdef POLARIZATION
!          call isotropic_pol ! Random polarization
! #endif
!          !----------------------------------------------------------------------------------------------------
!       case(3,4) ! QM redistribution 
! #ifdef POLARIZATION
!          if (abs(x_atom).gt.0.2) then ! wing scattering : pure Rayleigh
!             a0i = a0
!             b0i = b0
!             c0i = c0
!             ok = .false.
!             do while (.not. ok)
!                rphi   = ran3(iran)
!                rtheta = ran3(iran)     
!                phi    = 2.d0*pi*rphi
!                theta  = acos(1d0-2d0*rtheta)
!                !........director cosines
!                st = sin(theta)
!                a0 = st*cos(phi)   !x
!                b0 = st*sin(phi)   !y
!                c0 = cos(theta)    !z
!                f = 1-(a0*x_pol + b0*y_pol + c0*z_pol)**2 ! 1-(e.n')**2
!                ra = ran3(iran)
!                if (ra <=f) ok = .true.
!             end do
!             mu = a0i*a0 + b0i*b0 + c0i*c0
!             call projection_pol
!          else ! core scattering
!             ra = ran3(iran)
!             if(ra .lt. 0.666) then ! K transition or isotropic part of H transition
!                a0i = a0
!                b0i = b0
!                c0i = c0
!                rphi   = ran3(iran)
!                rtheta = ran3(iran)     
!                phi    = 2.d0*pi*rphi
!                theta  = acos(1d0-2d0*rtheta)
!                !........director cosines
!                st = sin(theta)
!                a0 = st*cos(phi)   !x
!                b0 = st*sin(phi)   !y
!                c0 = cos(theta)    !z
!                mu = a0i*a0 + b0i*b0 + c0i*c0  
!                call isotropic_pol
!             else ! Rayleigh part of H transition
!                ok = .false.
!                a0i = a0
!                b0i = b0
!                c0i = c0
!                do while (.not. ok)
!                   rphi   = ran3(iran)
!                   rtheta = ran3(iran)     
!                   phi    = 2.d0*pi*rphi
!                   theta  = acos(1d0-2d0*rtheta)
!                   !........director cosines
!                   st = sin(theta)
!                   a0 = st*cos(phi)   !x
!                   b0 = st*sin(phi)   !y
!                   c0 = cos(theta)    !z
!                   f = 1-(a0*x_pol + b0*y_pol + c0*z_pol)**2 ! 1-(e.n')**2
!                   ra = ran3(iran)
!                   if (ra <=f) ok = .true.
!                end do
!                mu = a0i*a0 + b0i*b0 + c0i*c0
!                call projection_pol
!             endif
!          endif
! #else
!          thetai=theta
!          phii=phi
!          !...........dipolar angular distribution
!          !.... step 1 : we generate the emergent angles theta1, phi1, 
!          !              but in the frame where the incident direction ki
!          !              is along the z axis.
!          !..........We generate mu=cos(theta1), from the dipolar distribution
!          !          theta1 is the angle between the incident and emergent directions 
!          ra=ran3(iran)
!          if (abs(x_atom).lt.0.2) then !core scattering
!             ! distribution of mu is \propto (11 + 3*cos^2 theta)
!             ok = .false.
!             do while (.not. ok)
!                mu = 2.0d0*ran3(iran)-1.0d0
!                f  = (11.0d0 + 3.0d0 * mu**2)/14.0d0  ! /14 factor is a normalization so that max of f is one.
!                ra = ran3(iran) 
!                if (ra <= f) ok = .true.
!             end do
!          else
!             ! distribution of mu (=cos theta) is \propto (1+cos^2 theta)
!             ok = .false.
!             do while (.not. ok)
!                mu = 2.0d0*ran3(iran)-1.0d0
!                f  = 0.5d0 * (1.0d0 + mu**2)  ! 0.5 factor is a normalization so that max of f is one.
!                ra = ran3(iran) 
!                if (ra <= f) ok = .true.
!             end do
!          endif
!          !..........We deduce theta1 from mu and generate phi1 randomly 
!          theta1=acos(mu)
!          rphi=ran3(iran)
!          phi1=2.d0*pi*rphi
!          !.... step 2 : we transform theta1 and phi1 in theta_e and phi_e 
!          !              in the initial frame.
!          !...........We perform a rotation of phii around the z axis,
!          !           and a rotation of thetai around the y axis.
!          cti=cos(thetai)   ;   sti=sin(thetai)
!          ct1=cos(theta1)   ;   st1=sin(theta1)
!          cpi=cos(phii)     ;   spi=sin(phii)
!          cp1=cos(phi1)     ;   sp1=sin(phi1)
!          a0 = cti*cpi*st1*cp1 + sti*cpi*ct1 - spi*st1*sp1
!          b0 = cti*spi*st1*cp1 + sti*spi*ct1 + cpi*st1*sp1
!          c0 = -sti*st1*cp1 + cti*ct1
!          theta = acos(c0)
!          phi   = atan2(b0,a0)
!          if (phi .lt. 0.d0) phi = phi + 2.d0 * pi  !! change from [-pi;pi] to [0,2pi]
! #endif
!          !----------------------------------------------------------------------------------------------------
!       end select


!       ! 4/ treat recoil effect (in atom's frame):
!       if (recoil) then
!          nu_atom = nu_atom / (1.d0 + ((planck*nu_atom)/(mp*clight*clight))*(1.-mu))
!       endif

!       ! 5 compute atom freq. in external frame, after scattering
!       scalar =  a0 * vx_cell + b0 * vy_cell + c0 * vz_cell
!       nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + sqrt(1-mu**2)*uper)/clight)
!       ! return nu_cell

!     end subroutine scatter_HI_polarization




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
