module module_SiII_model

  use module_constants
  use module_utils, only : voigt_fit, isotropic_direction, anisotropic_direction_HIcore, anisotropic_direction_Rayleigh
  use module_uparallel
  use module_random

  implicit none

  private


  !!!redéfinir les paramètres!!!

  ! definition of atomic values
  real(kind=8),parameter   :: lambda_r1=1193.28d0                         ![A] first resonant line wavelength
  real(kind=8),parameter   :: lambda_f1=1197.39d0                         ![A] first fluorescent line wavelength
  real(kind=8),parameter   :: lambda_r2=1190.42d0                         ![A] seconde resonant line wavelength
  real(kind=8),parameter   :: lambda_f2=1194.50d0                         ![A] seconde fluorescent line wavelength
 
  real(kind=8),parameter   :: gamma_r1=2.69d9                             ! Einstein coeff = damping constant for Voigt Function(gamma_alpha)  
  real(kind=8),parameter   :: gamma_f1=1.40d9                             ! Einstein coeff = damping constant for Voigt Function(gamma_alpha)  
  real(kind=8),parameter   :: gamma_r2=0.653d9                            ! Einstein coeff = damping constant for Voigt Function(gamma_alpha)  
  real(kind=8),parameter   :: gamma_f2=3.45d9                             ! Einstein coeff = damping constant for Voigt Function(gamma_alpha)  
 
  real(kind=8),parameter   :: f13=0.575                                   ! oscillator strength for "2P0(1/2) to 2P(1/2)" line
  real(kind=8),parameter   :: f23=0.150                                   ! oscillator strength for "2P0(3/2) to 2P(1/2)" line
  real(kind=8),parameter   :: f14=0.277                                   ! oscillator strength for "2P0(1/2) to 2P(3/2)" line
  real(kind=8),parameter   :: f24=0.737                                   ! oscillator strength for "2P0(3/2) to 2P(3/2)" line

  real(kind=8),parameter   :: ratio_fluo_1 = gamma_f1/(gamma_r1+gamma_f1) !ratio of fluorescent photon over resonant photon of the first line 
  real(kind=8),parameter   :: ratio_fluo_2 = gamma_f2/(gamma_r2+gamma_f2) !ratio of fluorescent photon over resonant photon of the second line 

  real(kind=8)  :: nu_D,x_cell,sigmaH,a,h
  real(kind=8)  :: tau_r1, tau_r2, nu_0

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [SiII] in the parameter file 
  ! --------------------------------------------------------------------------
  logical                  :: recoil       = .true.     ! if set to true, recoil effect is computed [default is true]
  logical                  :: isotropic    = .true.     ! if set to true, scattering events will be isotropic [default is false]
  ! --------------------------------------------------------------------------

  public :: get_tau_SiII, scatter_SiII, read_SiII_params, print_SiII_params, resonante_or_fluo

contains
  
  subroutine usefull_quantities(lambda, gamma, f_osci, vth, nu_cell)

    !--------------------------------------!
    !Calcul parameters of the choosen line !
    !--------------------------------------!

   real(kind=8),intent(in)     :: lambda,gamma,f_osci,vth,nu_cell
   real(kind=8)                :: lambda_cm, nu_Hz, sigmaH_factor, gamma_over_fourpi
    
   lambda_cm         = lambda / cmtoA                      ! cm
   nu_Hz             = clight / lambda_cm                  ! Hz
   sigmaH_factor     = pi*e_ch**2*f_osci/ me / clight      ! SiII cross-section factor-> multiply by Voigt(x,a)/nu_D to get sigma.
   gamma_over_fourpi = gamma / fourpi


   ! compute Doppler width and a-parameter, for SiII   
   nu_D = vth / lambda_cm 
   a    = gamma_over_fourpi / nu_D


   ! Cross section of SiII
   x_cell  = (nu_cell - nu_Hz)/nu_D
   h       = voigt_fit(x_cell,a)
   sigmaH  = sigmaH_factor / nu_D * h

  end subroutine usefull_quantities


  function get_tau_SiII(nsiII, vth, distance_to_border_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute optical depth of Silicon over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - nsiII    : number density of SiII atoms                            [ cm^-3 ]
    ! - vth      : thermal (+ small-scale turbulence) velocity of HI atoms [ cm / s ]
    ! - distance_to_border_cm : distance over which we compute tau         [ cm ]
    ! - nu_cell  : photon's frequency in the frame of the cell             [ Hz ]
    ! OUTPUT :
    ! - get_tau_SiII : optical depth of Silicon's line over distance_to_border_cm
    ! --------------------------------------------------------------------------
    
    real(kind=8),intent(in) :: nsiII,vth,distance_to_border_cm,nu_cell
    real(kind=8)            :: get_tau_SiII

    call usefull_quantities(lambda_r1,gamma_r1,f13,vth,nu_cell)
    tau_r1 = sigmaH * nsiII * distance_to_border_cm
    ! print*,'Tau of the first resonante line : ',tau_r1
    

    call usefull_quantities(lambda_r2,gamma_r2,f14,vth,nu_cell)
    tau_r2 = sigmaH * nsiII * distance_to_border_cm
    ! print*,'Tau of the seconde resonante line : ',tau_r2

    get_tau_SiII = tau_r1 + tau_r2

    return

  end function get_tau_SiII



  function resonante_or_fluo(ratio_rf, iran)

    !----------------------------------------------!
    !determine if the photon is fluorescent or not !
    !----------------------------------------------!
    real(kind=8),intent(in)  :: ratio_rf 
    real(kind=8)             :: emission_fluo, prob_rf, resonante_or_fluo
    integer,intent(inout)    :: iran

    prob_rf = ran3(iran)

      if (prob_rf <= ratio_rf) then
        emission_fluo = 1
      else
        emission_fluo = 0
      end if     

    resonante_or_fluo = emission_fluo
    return

  end function resonante_or_fluo

  
  subroutine scatter_SiII(vcell,vth,nu_cell,k,nu_ext,iran)

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
    ! where mu = cos(theta), (and theta in [0,pi]).
    ! ---------------------------------------------------------------------------------

    real(kind=8), intent(inout)               :: nu_cell, nu_ext
    real(kind=8), dimension(3), intent(inout) :: k
    real(kind=8), dimension(3), intent(in)    :: vcell
    real(kind=8), intent(in)                  :: vth
    integer, intent(inout)                    :: iran
    real(kind=8)               :: blah, upar, ruper
    real(kind=8)               :: r2, uper, nu_atom, mu, bu, scalar
    real(kind=8)               :: x_atom, tirage, select_fluo
    real(kind=8), dimension(3) :: knew
    character(len=8)           :: exit_statement

    tirage = ran3(iran)
    exit_statement='no_fluo'

     if (tirage <= tau_r2/(tau_r1+tau_r2))then
      call usefull_quantities(lambda_r2,gamma_r2,f14,vth,nu_cell)
      select_fluo = resonante_or_fluo(ratio_fluo_2,iran)
      if (select_fluo == 1) then
        exit_statement = 'fluo_2'
      endif
    else
      call usefull_quantities(lambda_r1,gamma_r1,f13,vth,nu_cell)
      select_fluo = resonante_or_fluo(ratio_fluo_1,iran)
      if (select_fluo == 1) then
        exit_statement = 'fluo_1'
      endif
	  endif

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

    ! 4/ change atom frequency in case of fluorescent emission
    select case (exit_statement)
    case('fluo_1')
      nu_atom = clight*cmtoA/lambda_f1
    case('fluo_2')
      nu_atom = clight*cmtoA/lambda_f2
    end select

    ! 5/ determine direction of scattered photon
    if (isotropic) then
       call isotropic_direction(knew,iran)
       mu = k(1)*knew(1) + k(2)*knew(2) + k(3)*knew(3)
       bu = sqrt(1.0d0 - mu*mu)
    else
       x_atom  = (nu_atom - nu_0) / nu_D
       if (abs(x_atom) < 0.2) then ! core scattering 
          call anisotropic_direction_HIcore(k,knew,mu,bu,iran)
       else ! wing scattering 
          call anisotropic_direction_Rayleigh(k,knew,mu,bu,iran)
       end if
    end if

    ! 6/ recoil effect 
    if (recoil) then 
       nu_atom = nu_atom / (1.d0 + ((planck*nu_atom)/(mp*clight*clight))*(1.-mu))
    end if
    
    ! 7/ compute atom freq. in external frame, after scattering
    scalar = knew(1) * vcell(1) + knew(2) * vcell(2) + knew(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + bu*uper)/clight)
    nu_cell = (1.d0 - scalar/clight) * nu_ext 
    k = knew

  end subroutine scatter_SiII

  

  subroutine read_SiII_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present

    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:4) == '[SiII]') then
          section_present = .true.
          exit
       end if
    end do
    ! read section if present
    if (section_present) then 
       do
          read (10,'(a)',iostat=err) line
          if(err/=0) exit
          if (line(1:1) == '[') exit ! next section starting... -> leave
          i = scan(line,'=')
          if (i==0 .or. line(1:1)=='#' .or. line(1:1)=='!') cycle  ! skip blank or commented lines
          name=trim(adjustl(line(:i-1)))
          value=trim(adjustl(line(i+1:)))
          i = scan(value,'!')
          if (i /= 0) value = trim(adjustl(value(:i-1)))
          select case (trim(name))
          case ('recoil')
             read(value,*) recoil
          case ('isotropic')
             read(value,*) isotropic
          end select
       end do
    end if
    close(10)
    return
    
  end subroutine read_SiII_params


  
  subroutine print_SiII_params(unit)
    
    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)') '[HI]'
       write(unit,'(a,L1)') '  recoil    = ',recoil
       ! write(unit,'(a,L1)') '  isotropic = ',isotropic
    else
       write(*,'(a,a,a)') '[HI]'
       write(*,'(a,L1)') '  recoil    = ',recoil
       ! write(*,'(a,L1)') '  isotropic = ',isotropic
    end if

    return
    
  end subroutine print_SiII_params


end module module_SiII_model
