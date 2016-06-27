module module_D_model

  implicit none

  use module_constants
  use module_HI_model

#ifdef POLARIZATION
  use polar
#endif

  private

  ! Deuterium properties
  real(kind=8),parameter   :: nu_0_deut       = nu_0 * (1.0d0+me/mp)/(1.0d0+me/mdeut) ! Hz
  real(kind=8),parameter   :: deut2H_nb_ratio = 3.e-5                                 ! assumed Deuterium/H abundance (in number)
  real(kind=8),parameter   :: mdeut=2.*mp                                             ![g] Deuterium mass
  real(kind=8),parameter   :: sqrt_H2Deut_mass_ratio = 0.7071067811865d0              ! == sqrt(mp/mdeut) = 1/sqrt(2) 


  public

  logical              :: deuterium
  ! => read params deuterium

  contains
    ! Routine list
    ! function get_tau
    ! subroutine scatter

    function get_tau_D(nhi, dopwidth, distance_to_border_cm, nu_cell)

      ! compute frequency in cell's moving frame
      !scalar =  a0 * vx_cell + b0 * vy_cell + c0 * vz_cell
      !nu_int = (1.d0 - scalar/clight) * nu_ext
      ! => should be done before

      real(kind=8),intent(in) :: nhi,dopwidth,distance_to_border_cm,nu_cell
      real(kind=8)            :: cross_section

      real(kind=8) :: nu_D,x_cell,sigmaH,sigma_deut,a,h
      real(kind=8) :: x_deut,nu_D_deut ,a_deut,b_cell_deut

      ! compute Doppler width and a-parameter, for Deuterium
      nu_D = dopwidth * nu_0 / clight
      a    = gamma / (4.d0 * pi * nu_D)
 
      nu_D_deut   = nu_D * sqrt_H2Deut_mass_ratio
      a_deut      = a / sqrt_H2Deut_mass_ratio
      b_cell_deut = nu_D_deut / nu_0_deut * clight
 
      ! Cross section of Deuterium
      x_cell  = (nu_cell - nu_0)/nu_D
      h       = voigt_fit(x_cell,a)
      sigmaH  = sigmaH_factor / nu_D * h
      x_deut     = (nu_cell - nu_0_deut) / nu_D_deut
      h          = voigt_fit(x_deut,a_deut)
      sigma_deut = sigmaH_factor / nu_D_deut * h

      get_tau_D = sigma_deut * deut2H_nb_ratio * nhi * distance_to_border_cm

      return

    end function get_tau_D
    


    subroutine scatter_D(p,cellProps)

      type(photon_current),intent(inout) :: p
      type(gas),intent(in) :: cellProps

      real(kind=8) :: upar
      
      input_freq = 
      therm_param = 


      ! 1.1/ component parallel to photon's propagation
      call uparallel(x_deut,a_deut,iran,upar)
      upar = upar * b_cell_deut            ! upar is an x -> convert to a velocity 
      ! 1.2/ component perpendicular to photon's propagation
      ruper  = ran3(iran)
      r2     = ran3(iran)
      uper   = sqrt(-log(ruper))*cos(2.d0*pi*r2)
      uper   = uper * b_cell_deut  ! from x to velocity
      ! 3/ determine scattering angle (in atom's frame)
      nu_atom = nu_int - nu_ext * upar/clight
      x_atom = (nu_atom - nu_0_deut)/nu_D_deut
      if (dipol.eq.4) then   ! dipol=4 signals QM scattering for Lya and Henyey-Greenstein for dust
         dip=dipol           ! --> indicate that this is Lya here (dipol=3) not dust
         dipol=3
         call emission_direction(x_atom)
         dipol=dip
      else
         call emission_direction(x_atom)
      endif
      ! 4/ treat recoil effect (in atom's frame):
      if (recoil) then
         nu_atom = nu_atom / (1.d0 + ((planck*nu_atom)/(mdeut*clight*clight))*(1.-mu))
      endif
      ! 5 compute atom freq. in external frame, after scattering
      scalar =  a0 * vx_cell + b0 * vy_cell + c0 * vz_cell
      nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + sqrt(1-mu**2)*uper) / clight)

    end subroutine scatter_D

 

  end module module_D_model
