module module_dust_model

  ! this module implements the model for dust extinction described by Laursen, Sommer-Larsen, and Andersen (2009).
  
  use module_random
  use module_constants, only : pi, clight
  use module_utils, only : anisotropic_direction_Dust, anisotropic_probability_dust

  implicit none

  private

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [dust] of the parameter file 
  ! --------------------------------------------------------------------------
  real(kind=8)  :: albedo     = 0.32    ! dust albedo [default 0.32 for Lya, from Li & Draine 2001]. 
  real(kind=8)  :: g_dust     = 0.73    ! g parameter of the Henyey-Greenstein phase function for dust scattering [default 0.73 from Li & Draine 2001]
  character(20) :: dust_model = 'SMC' 
  ! --------------------------------------------------------------------------
  
  public get_tau_dust, scatter_dust, read_dust_params, print_dust_params, dust_peeloff_weight
  private sigma_d
  
contains


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
    real(kind=8),dimension(3)               :: knew
    integer(kind=4),intent(out)             :: ilost

    ! interaction with dust
    aors = ran3(iran)  ! aka "Absorption OR Scattering" ... 
    ilost = 0
    if (aors.gt.albedo) then ! the photon is absorbed = lost
       ilost = 1
       return
    else
       call anisotropic_direction_Dust(k,knew,mu,iran,g_dust)                         
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
    real(kind=8)                            :: scalar,nu_cell

    dust_peeloff_weight = anisotropic_probability_Dust(kin,kout,g_dust)
    ! compute freq. in cell frame
    scalar  = kin(1) * vcell(1) + kin(2) * vcell(2) + kin(3) * vcell(3)
    nu_cell = (1.d0 - scalar/clight) * nu_ext
    ! nu_cell has not changed; we implicitly assume that the interacting dust grain is at rest in cell's frame
    scalar = kout(1) * vcell(1) + kout(2) * vcell(2) + kout(3)* vcell(3)
    nu_ext = nu_cell/(1.d0 - scalar/clight)

  end function dust_peeloff_weight
  !--LEEP--

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
          case ('albedo')
             read(value,*) albedo
          case ('g_dust')
             read(value,*) g_dust
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
       write(unit,'(a,ES10.3)') '  albedo       = ',albedo
       write(unit,'(a,ES10.3)') '  g_dust       = ',g_dust
       write(unit,'(a,a)')      '  dust_model   = ',trim(dust_model)
    else
       write(*,'(a,a,a)') '[dust]'
       write(*,'(a,ES10.3)') '  albedo       = ',albedo
       write(*,'(a,ES10.3)') '  g_dust       = ',g_dust
       write(*,'(a,a)')      '  dust_model   = ',trim(dust_model)
    end if

    return

  end subroutine print_dust_params

  
  function sigma_d(lbda,model)

    ! returns the effective cross section of dust (per H atom), in cgs, at wavelength lbda (in A),
    ! for a given model. NB: this is the total cross section (including abs. and scat. probs). 
    ! Models can be 'SMC' or 'LMC', which are taken from fits given by Gnedin, Kravtsov, and Chen (2008). 
    
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
    case default
       print*,'dust model unknown',trim(model)
       stop
    end select

    return
    
  end function sigma_d

  
end module module_dust_model

