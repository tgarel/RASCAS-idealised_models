module module_dust_model

  ! this module implements the model for dust extinction described by Laursen, Sommer-Larsen, and Andersen (2009).
  
  use module_random
  use module_constants, only : pi, clight
  use module_utils, only : anisotropic_direction_Dust
  !--PEEL--
  use module_utils, only : anisotropic_probability_dust
  !--LEEP--

  implicit none

  private

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [dust] of the parameter file 
  ! --------------------------------------------------------------------------
  character(20) :: dust_model = 'SMC' 
  ! --------------------------------------------------------------------------
  
  public  :: get_tau_dust, scatter_dust, read_dust_params, print_dust_params, get_n_dust
  private :: sigma_d, get_albedo_g
  !--PEEL--
  public  :: dust_peeloff_weight
  !--LEEP--
  
contains

  function get_n_dust(nH,nHI,nHII,Z,Zref,f_ion,nO,T)

    real(kind=8),intent(in)     :: nH,nHI,nHII,Z,Zref,f_ion,nO,T
    real(kind=8)                :: get_n_dust
    real(kind=8)                :: x, xsun, a_rr, aH_rr, b_rr, aL_rr, xt_rr, y

    select case(trim(dust_model))
    case('SMC')
      get_n_dust = Z / Zref * ( nhI + f_ion*nHII )   ! [ /cm3 ]
    case('LMC')
      get_n_dust = Z / Zref * ( nhI + f_ion*nHII )   ! [ /cm3 ]
    case('BARE-GR-S')
      ! First check the temperature -- assume sputtering destroys dust at high T
      ! consistent with RAMSES
      if (T.gt.1.d5) then
         get_n_dust = 0.d0
      else if (nO.le.0.d0) then
         get_n_dust = 0.d0
      else
         ! XCO,Z case from remy-ruyer 2014
         x     = 12.d0 + log10(nO / nH)
         xsun  = 8.69d0
         a_rr  = 2.21d0
         aH_rr = 1.d0 
         b_rr  = 0.96d0
         aL_rr = 3.10d0
         xt_rr = 8.10d0
         if (x.gt.xt_rr) then
            y = a_rr + ( aH_rr * (xsun - x) )
         else
            y = b_rr + ( aL_rr * (xsun - x) )
         endif

         y = 10.0**y   ! This is the gas to dust mass ratio

         ! if we assume that the dust properties don't change with 
         ! metallicity, the pseudo n_dust should just be nH * (G/D)sun / (G/D)sim
         ! since the scaling between n_dust and n_h is the same as m_dust and m_H
         ! The gas to dust mass ratio in the sun for the BARE-GR-S model is 1 / 0.00619
         get_n_dust = ( (1.d0 / 0.00619d0) / y) * nH
      end if
    case default
       print*,'dust model unknown',trim(dust_model)
       stop
    end select

  end function get_n_dust

  

  function get_tau_dust(ndust, distance, nu)

    ! ---------------------------------------------------------------------------------
    ! compute wavelength-dependent optical depth of dust 
    ! ---------------------------------------------------------------------------------
    ! INPUTS :
    ! - ndust    : dust pseudo-density [/cm3]
    ! - distance : distance over which to compute tau [cm]
    ! - nu       : frequency of photon in cell's frame [Hz]
    ! OUTPUTS :
    ! - tau : optical depth of dust. 
    ! ---------------------------------------------------------------------------------
    
    real(kind=8),intent(in) :: ndust,distance,nu
    real(kind=8)            :: get_tau_dust,lbda_A
    
    lbda_A = clight/nu*1d8
    get_tau_dust = sigma_d(lbda_A,dust_model) * ndust * distance
    
    return

  end function get_tau_dust



  subroutine scatter_dust(v,nu_cell,k,nu_ext,iran,ilost)

    implicit none 
    
    ! ---------------------------------------------------------------------------------
    ! perform scattering event on a dust grain with anisotropic angular redistribution
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
    ! _ ilost    : updated status of the photon telling if photon has been absorbed
    ! ---------------------------------------------------------------------------------

    real(kind=8),intent(inout)              :: nu_cell, nu_ext ! nu_cell in RASCAS = nu_int in MCLya
    real(kind=8),dimension(3),intent(inout) :: k
    real(kind=8),dimension(3),intent(in)    :: v
    integer(kind=4),intent(inout)           :: iran
    real(kind=8)                            :: mu, aors, scalar
    real(kind=8)                            :: albedo_local, g_dust_local, lbda
    integer(kind=4)                         :: i
    real(kind=8),dimension(3)               :: knew
    integer(kind=4),intent(out)             :: ilost


    ! interaction with dust
    aors = ran3(iran)  ! aka "Absorption OR Scattering" ... 
    ilost = 0

    call get_albedo_g(nu_cell,albedo_local,g_dust_local)

    
    if (aors.gt.albedo_local) then ! the photon is absorbed = lost
       ilost = 1
       return
    else
       call anisotropic_direction_Dust(k,knew,mu,iran,g_dust_local)                         
       ! compute atom freq. in external frame, after scattering
       scalar = knew(1) * v(1) + knew(2) * v(2) + knew(3)* v(3)
       ! nu_cell has not changed; we implicitly assume that the interacting dust grain is at rest in cell's frame
       nu_ext = nu_cell/(1.0d0 - scalar/clight)
       k = knew
    end if

  end subroutine scatter_dust



  !--PEEL--
  function dust_peeloff_weight(vcell,nu_ext,kin,kout)

    ! ---------------------------------------------------------------------------------
    ! Compute probability that a photon coming along kin scatters off in direction kout.
    ! Also update nu_ext to external-frame frequency along kout
    ! ---------------------------------------------------------------------------------
    ! INPUTS :
    ! - vcell    : bulk velocity of the gas (i.e. cell velocity)       [ cm / s ] 
    ! - nu_ext   : frequency of incoming photon, in external frame     [ Hz ]
    ! - kin      : propagation vector (normalized)
    ! - kout     : direction after interaction (fixed)
    ! OUTPUTS :
    ! - nu_ext   : updated frequency in external frame [ Hz ]
    ! ---------------------------------------------------------------------------------

    real(kind=8),intent(inout)              :: nu_ext
    real(kind=8),dimension(3),intent(in)    :: kin, kout
    real(kind=8),dimension(3),intent(in)    :: vcell
    real(kind=8)                            :: dust_peeloff_weight
    real(kind=8)                            :: scalar,nu_cell,albedo_local,g_dust_local

    ! compute freq. in cell frame
    scalar  = kin(1) * vcell(1) + kin(2) * vcell(2) + kin(3) * vcell(3)
    nu_cell = (1.d0 - scalar/clight) * nu_ext


    call get_albedo_g(nu_cell,albedo_local,g_dust_local)


    dust_peeloff_weight = anisotropic_probability_Dust(kin,kout,g_dust_local) * albedo_local
    
    ! nu_cell has not changed; we implicitly assume that the interacting dust grain is at rest in cell's frame
    scalar = kout(1) * vcell(1) + kout(2) * vcell(2) + kout(3)* vcell(3)
    nu_ext = nu_cell/(1.d0 - scalar/clight)

  end function dust_peeloff_weight
  !--LEEP--

  
  function sigma_d(lbda,model)

    ! returns the effective cross section of dust (per H atom), in cgs, at wavelength lbda (in A),
    ! for a given model. NB: this is the total cross section (including abs. and scat. probs). 
    ! Models can be 'SMC' or 'LMC', which are taken from fits given by Gnedin, Kravtsov, and Chen (2008). 

    ! Harley added the BARE-GR-S model which is from Zubko 2004 and also used by Remy-Ruyer 2014
    ! Data for that model can be downloaded here: https://ipag.osug.fr/RT13/RTTRUST/opa.php 
    
    implicit none 

    real(kind=8),intent(in)  :: lbda
    character(20),intent(in) :: model
    real(kind=8)             :: sigma_d,x
    integer(kind=4)          :: i
    ! SMC parameters
    real(kind=8),parameter,dimension(7) :: smc_lbda = (/0.042d0,0.08d0,0.22d0,9.7d0,18.0d0,25.0d0,0.067d0/)*1.0d4   ! [A]  
    real(kind=8),parameter,dimension(7) :: smc_a    = (/185.0d0,27.0d0,0.005d0,0.010d0,0.012d0,0.03d0,10.0d0/)
    real(kind=8),parameter,dimension(7) :: smc_b    = (/90.0d0,15.5d0,-1.95d0,-1.95d0,-1.8d0,0.0d0,1.9d0/)
    real(kind=8),parameter,dimension(7) :: smc_p    = (/2.0d0,4.0d0,2.0d0,2.0d0,2.0d0,2.0d0,4.0d0/)
    real(kind=8),parameter,dimension(7) :: smc_q    = (/2.0d0,4.0d0,2.0d0,2.0d0,2.0d0,2.0d0,15.0d0/)
    real(kind=8),parameter              :: smc_sig0 = 1.0d-22  ! [cm^-2]
    ! LMC parameters
    real(kind=8),parameter,dimension(7) :: lmc_lbda = (/0.046d0,0.08d0,0.22d0,9.7d0,18.0d0,25.0d0,0.067d0/)*1.0d4    ! [A]
    real(kind=8),parameter,dimension(7) :: lmc_a    = (/90.0d0,19.0d0,0.023d0,0.005d0,0.006d0,0.02d0,10.0d0/)
    real(kind=8),parameter,dimension(7) :: lmc_b    = (/90.0d0,21.0d0,-1.95d0,-1.95d0,-1.8d0,0.0d0,1.9d0/)
    real(kind=8),parameter,dimension(7) :: lmc_p    = (/2.0d0,4.5d0,2.0d0,2.0d0,2.0d0,2.0d0,4.0d0/)
    real(kind=8),parameter,dimension(7) :: lmc_q    = (/2.0d0,4.5d0,2.0d0,2.0d0,2.0d0,2.0d0,15.0d0/)
    real(kind=8),parameter              :: lmc_sig0 = 3.0d-22  ! [cm^-2]

    ! Harley additions for BARE-GR-S model
    real(kind=8),parameter              :: BARE_GR_S_sig0 = 1.0d-21  ! [cm^-2]
    real(kind=8),parameter              :: BARE_GR_S_lmin = 2.6989700043360187d0  ! [log10 A]
    real(kind=8),parameter              :: BARE_GR_S_lmax = 4.0d0  ! [log10 A]
    real(kind=8),parameter,dimension(223):: BARE_GR_S = (/ 1.99464 , 2.02974 , 2.06585 , 2.10377 , 2.14249 , 2.18272 , &
                                                         & 2.22755 , 2.27368 , 2.32128 , 2.37166 , 2.42422 , 2.52427 , &
                                                         & 2.6174  , 2.71557 , 2.81412 , 2.91304 , 3.01794 , 3.12497 , &
                                                         & 3.23208 , 3.33917 , 3.45689 , 3.57147 , 3.69    , 3.80869 , &
                                                         & 3.8947  , 3.96074 , 4.00219 , 4.01394 , 3.97888 , 3.91113 , &
                                                         & 3.82775 , 3.73618 , 3.63773 , 3.51273 , 3.39593 , 3.28784 , &
                                                         & 3.18822 , 3.10239 , 3.03599 , 2.97574 , 2.92081 , 2.87048 , &
                                                         & 2.8271  , 2.79807 , 2.76805 , 2.71858 , 2.66752 , 2.61272 , &
                                                         & 2.54834 , 2.48204 , 2.4154  , 2.34225 , 2.28339 , 2.22239 , &
                                                         & 2.15469 , 2.08873 , 2.02561 , 1.96582 , 1.90941 , 1.86132 , &
                                                         & 1.82145 , 1.78    , 1.7434  , 1.70746 , 1.67348 , 1.63987 , &
                                                         & 1.60782 , 1.57673 , 1.54658 , 1.51725 , 1.48906 , 1.46846 , &
                                                         & 1.44286 , 1.41809 , 1.39371 , 1.37045 , 1.3477  , 1.32561 , &
                                                         & 1.30327 , 1.28032 , 1.25827 , 1.23799 , 1.21938 , 1.20235 , &
                                                         & 1.18758 , 1.17479 , 1.1638  , 1.15694 , 1.15109 , 1.14946 , &
                                                         & 1.15016 , 1.15549 , 1.16379 , 1.1728  , 1.18264 , 1.19672 , &
                                                         & 1.21764 , 1.24159 , 1.26811 , 1.29842 , 1.33402 , 1.3744  , &
                                                         & 1.41939 , 1.4651  , 1.50855 , 1.54798 , 1.5806  , 1.6033  , &
                                                         & 1.61161 , 1.59922 , 1.56562 , 1.51892 , 1.46754 , 1.41611 , &
                                                         & 1.3648  , 1.31241 , 1.26313 , 1.21781 , 1.17673 , 1.14017 , &
                                                         & 1.10904 , 1.08115 , 1.05559 , 1.03219 , 1.01045 , 0.989836, &
                                                         & 0.970846, 0.955278, 0.94033 , 0.92545 , 0.911688, 0.898633, &
                                                         & 0.88632 , 0.874301, 0.859604, 0.847324, 0.835483, 0.823297, &
                                                         & 0.811086, 0.799162, 0.787961, 0.777052, 0.76644 , 0.75588 , &
                                                         & 0.745454, 0.734821, 0.724276, 0.713884, 0.70391 , 0.69423 , &
                                                         & 0.684606, 0.674863, 0.664968, 0.655049, 0.645267, 0.635624, &
                                                         & 0.626174, 0.616844, 0.607489, 0.598024, 0.588505, 0.579024, &
                                                         & 0.569713, 0.560556, 0.551512, 0.542515, 0.533496, 0.524434, &
                                                         & 0.515373, 0.506454, 0.497737, 0.489132, 0.480608, 0.472122, &
                                                         & 0.463641, 0.455169, 0.446736, 0.438371, 0.430147, 0.422218, &
                                                         & 0.414361, 0.406559, 0.398804, 0.391102, 0.383464, 0.375898, &
                                                         & 0.368404, 0.360986, 0.353645, 0.346381, 0.339195, 0.332091, &
                                                         & 0.325192, 0.318368, 0.311619, 0.30495 , 0.298382, 0.291942, &
                                                         & 0.285584, 0.279303, 0.2731  , 0.266978, 0.260942, 0.254995, &
                                                         & 0.249135, 0.243357, 0.23768 , 0.232126, 0.226661, 0.221291, &
                                                         & 0.216015, 0.210831, 0.205731, 0.200711, 0.195767, 0.190924, &
                                                         & 0.18618 , 0.181521, 0.176946, 0.172455, 0.168046, 0.163719, &
                                                         & 0.159472 /)

    sigma_d = 0.0d0
    select case(trim(model))
    case('SMC')
       do i = 1,7
          x = lbda / smc_lbda(i)
          sigma_d = sigma_d + smc_a(i) / (x**smc_p(i) + x**(-smc_q(i)) + smc_b(i))
       end do
       sigma_d = sigma_d * smc_sig0
    case('LMC')
       do i = 1,7
          x = lbda / lmc_lbda(i)
          sigma_d = sigma_d + lmc_a(i) / (x**lmc_p(i) + x**(-lmc_q(i)) + lmc_b(i))
       end do
       sigma_d = sigma_d * lmc_sig0
    case('BARE-GR-S')
       ! Note that this model is only good between 500 - 10000 A
       
       ! Val, for photons below 500A, setting to the cross-section at 500A, and same above 10000A:
       i = min(max(floor(223.d0 * ((log10(lbda) - BARE_GR_S_lmin) / (BARE_GR_S_lmax-BARE_GR_S_lmin))), 0), 222)
       ! laV
       i = i + 1 ! Account for fortran 1 indexing

       sigma_d = BARE_GR_S(i) * BARE_GR_S_sig0
    case default
       print*,'dust model unknown',trim(model)
       stop
    end select

    return
    
  end function sigma_d


  subroutine get_albedo_g(nu, albedo_local, g_dust_local)

    ! Returns the albedo and g parameter from the Li & Draine 2001 table,  as a function of frequency

    implicit none

    real(kind=8),intent(in)               :: nu
    real(kind=8),intent(out)              :: albedo_local, g_dust_local
    real(kind=8)                          :: lbda
    integer(kind=4)                       :: i
    real(kind=8),parameter                :: BARE_GR_S_lmin = 2.6989700043360187d0  ! [log10 A]
    real(kind=8),parameter                :: BARE_GR_S_lmax = 4.0d0  ! [log10 A]
    real(kind=8),parameter,dimension(223) :: BARE_GR_S_a = (/ 0.2356 , 0.23265, 0.22977, 0.22688, 0.2241 , 0.22139, 0.21855, &
                                                           & 0.21582, 0.21317, 0.21069, 0.20827, 0.20216, 0.19725, 0.19263, &
                                                           & 0.18843, 0.18461, 0.18087, 0.17741, 0.17429, 0.17149, 0.16808, &
                                                           & 0.16517, 0.16241, 0.15996, 0.15747, 0.15583, 0.15532, 0.15608, &
                                                           & 0.15789, 0.16068, 0.16435, 0.16865, 0.17331, 0.17776, 0.18235, &
                                                           & 0.18701, 0.19171, 0.1958 , 0.19847, 0.20101, 0.20345, 0.20581, &
                                                           & 0.20843, 0.21246, 0.21668, 0.22151, 0.22649, 0.23167, 0.23721, &
                                                           & 0.24303, 0.24777, 0.25388, 0.2592 , 0.26508, 0.27185, 0.27819, &
                                                           & 0.28389, 0.28859, 0.29309, 0.29722, 0.29972, 0.30314, 0.30671, &
                                                           & 0.31069, 0.31473, 0.31918, 0.3236 , 0.32806, 0.33248, 0.33676, &
                                                           & 0.34049, 0.34226, 0.34474, 0.34664, 0.34815, 0.34922, 0.35012, &
                                                           & 0.35111, 0.35275, 0.35508, 0.35757, 0.36033, 0.36348, 0.36798, &
                                                           & 0.3722 , 0.37615, 0.38123, 0.38513, 0.39049, 0.39514, 0.40004, &
                                                           & 0.40352, 0.40606, 0.40916, 0.41302, 0.41541, 0.41405, 0.41221, &
                                                           & 0.41019, 0.40675, 0.40054, 0.39386, 0.38679, 0.37887, 0.37121, &
                                                           & 0.36539, 0.36157, 0.35962, 0.36099, 0.36643, 0.37542, 0.38622, &
                                                           & 0.39812, 0.41081, 0.42408, 0.43758, 0.45013, 0.46013, 0.46927, &
                                                           & 0.47752, 0.48358, 0.4891 , 0.49325, 0.49627, 0.49934, 0.50235, &
                                                           & 0.50503, 0.50566, 0.5064 , 0.50714, 0.50744, 0.50782, 0.5082 , &
                                                           & 0.50856, 0.51091, 0.5121 , 0.513  , 0.51385, 0.51479, 0.51584, &
                                                           & 0.5168 , 0.51784, 0.51889, 0.51984, 0.5205 , 0.52111, 0.52177, &
                                                           & 0.52248, 0.52336, 0.52428, 0.52514, 0.52593, 0.52669, 0.52741, &
                                                           & 0.52815, 0.52917, 0.53025, 0.53133, 0.53235, 0.53329, 0.53419, &
                                                           & 0.53509, 0.53593, 0.53681, 0.53775, 0.53869, 0.53961, 0.54051, &
                                                           & 0.54145, 0.54217, 0.54251, 0.54288, 0.54323, 0.54355, 0.54384, &
                                                           & 0.54413, 0.54441, 0.54473, 0.54493, 0.54461, 0.54426, 0.54387, &
                                                           & 0.54343, 0.54296, 0.54246, 0.54192, 0.54137, 0.54079, 0.5402 , &
                                                           & 0.53959, 0.53896, 0.53828, 0.53711, 0.5359 , 0.53463, 0.53332, &
                                                           & 0.5319 , 0.53021, 0.52847, 0.52667, 0.52482, 0.52291, 0.52094, &
                                                           & 0.51892, 0.51683, 0.51469, 0.51242, 0.50991, 0.50735, 0.50472, &
                                                           & 0.50202, 0.49925, 0.49642, 0.49352, 0.49058, 0.48746, 0.4842 , &
                                                           & 0.48088, 0.47751, 0.47408, 0.4706 , 0.46705, 0.46345 /)
    real(kind=8),parameter,dimension(223) :: BARE_GR_S_g = (/ 0.769798, 0.76558 , 0.761303, 0.757024, 0.752622, 0.748124, &
                                                           & 0.74361 , 0.738938, 0.734119, 0.729173, 0.724061, 0.718785, &
                                                           & 0.713216, 0.707259, 0.701152, 0.694905, 0.688412, 0.681739, &
                                                           & 0.674935, 0.668015, 0.66176 , 0.655385, 0.648847, 0.642169, &
                                                           & 0.638838, 0.635637, 0.632192, 0.62853 , 0.626328, 0.624827, &
                                                           & 0.623112, 0.6212  , 0.619584, 0.620912, 0.621847, 0.622399, &
                                                           & 0.622583, 0.622595, 0.622687, 0.622509, 0.622069, 0.621376, &
                                                           & 0.619794, 0.615201, 0.610477, 0.607479, 0.604671, 0.602309, &
                                                           & 0.600957, 0.599819, 0.600712, 0.601186, 0.601504, 0.602078, &
                                                           & 0.603678, 0.605876, 0.608703, 0.612383, 0.615934, 0.619207, &
                                                           & 0.623875, 0.628331, 0.631256, 0.633839, 0.636042, 0.638035, &
                                                           & 0.639849, 0.641559, 0.643134, 0.644639, 0.646041, 0.647264, &
                                                           & 0.648137, 0.648568, 0.648528, 0.647673, 0.646137, 0.643718, &
                                                           & 0.64058 , 0.637256, 0.63447 , 0.631958, 0.629929, 0.627314, &
                                                           & 0.624912, 0.622719, 0.619866, 0.617266, 0.613786, 0.610274, &
                                                           & 0.606414, 0.602866, 0.599209, 0.594682, 0.589461, 0.584777, &
                                                           & 0.581277, 0.576948, 0.571605, 0.566576, 0.562945, 0.559528, &
                                                           & 0.555225, 0.550129, 0.545419, 0.54132 , 0.537808, 0.534674, &
                                                           & 0.530069, 0.524576, 0.520134, 0.517303, 0.51513 , 0.513185, &
                                                           & 0.510474, 0.507333, 0.504559, 0.503515, 0.502768, 0.502272, &
                                                           & 0.502507, 0.501861, 0.501123, 0.500717, 0.50035 , 0.500043, &
                                                           & 0.500017, 0.501197, 0.501721, 0.501608, 0.50136 , 0.501221, &
                                                           & 0.50119 , 0.501451, 0.501864, 0.501957, 0.501557, 0.500906, &
                                                           & 0.500392, 0.500072, 0.499803, 0.499785, 0.499712, 0.499217, &
                                                           & 0.498089, 0.496909, 0.495883, 0.495056, 0.49419 , 0.493291, &
                                                           & 0.492178, 0.490675, 0.488897, 0.487125, 0.485476, 0.483933, &
                                                           & 0.482453, 0.480861, 0.47901 , 0.476895, 0.474657, 0.472515, &
                                                           & 0.470476, 0.468496, 0.466477, 0.464308, 0.461932, 0.45939 , &
                                                           & 0.456789, 0.454236, 0.451767, 0.449306, 0.44679 , 0.444166, &
                                                           & 0.441423, 0.438593, 0.435731, 0.432866, 0.429983, 0.427022, &
                                                           & 0.424013, 0.420939, 0.4178  , 0.41461 , 0.411381, 0.408119, &
                                                           & 0.404822, 0.401486, 0.39811 , 0.394694, 0.391238, 0.387734, &
                                                           & 0.384133, 0.380479, 0.376773, 0.373027, 0.369248, 0.365434, &
                                                           & 0.361582, 0.357681, 0.353729, 0.349731, 0.345694, 0.341622, &
                                                           & 0.337515, 0.333365, 0.329181, 0.324981, 0.320752, 0.3165  , &
                                                           & 0.312225, 0.307916, 0.303562, 0.299154, 0.294695, 0.29023 , &
                                                           & 0.28576 , 0.281259, 0.276724, 0.272144, 0.267514, 0.262829, &
                                                           & 0.258092 /)


    ! Note that this model is only good between 500 - 10000 A
    ! Val: add min max to make it run below 500A and above 10000A, as in sigma_d
    lbda = (clight / nu) * 1.d8 ![A]
    i = min(max(floor(223.d0 * ((log10(lbda) - BARE_GR_S_lmin) / (BARE_GR_S_lmax-BARE_GR_S_lmin))), 0), 222)
    i = i + 1 ! Account for fortran 1 indexing

    albedo_local = BARE_GR_S_a(i)
    g_dust_local = BARE_GR_S_g(i)


  end subroutine get_albedo_g



  subroutine read_dust_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000)         :: line,name,value
    integer(kind=4)         :: err,i
    logical                 :: section_present

    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:6) == '[dust]') then
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
          case ('dust_model')
              write(dust_model,'(a)') trim(value) 
          end select
       end do
    end if
    close(10)

    return

  end subroutine read_dust_params



  subroutine print_dust_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)') '[dust]'
       write(unit,'(a,a)')      '  dust_model   = ',trim(dust_model)
    else
       write(*,'(a,a,a)') '[dust]'
       write(*,'(a,a)')      '  dust_model   = ',trim(dust_model)
    end if

    return

  end subroutine print_dust_params

  
end module module_dust_model

