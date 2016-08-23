module module_dust_model

  use module_random
  use module_constants, only : pi, clight
  use module_utils, only : anisotropic_direction_Dust

  implicit none

  private

  real(kind=8)             :: sigmad    ! dust cross section, computed at initialisation (in read_dust_params) as a function of parameter values 

  ! user-defined parameters - read from section [dust] of the parameter file 
  real(kind=8)             :: albedo       = 0.46    ! dust albedo [default is 0.46 : empirical value from Witt & Gordon 2000 (ApJ, 528, 799) for SMC dust]
  real(kind=8)             :: g_dust       = 0.73    ! g parameter of the Henyey-Greenstein phase function for dust scattering [default 0.73 from Laursen09]
  real(kind=8)             :: grain_radius = 2.0d-6  ! radius of a dust grain [cm] used to define cross section. [default 2d-6 from Verhamme+12]

  ! ---------------------------------------------------------------------------------
  ! notes on the values of the parameters : 
  ! ---------------------------------------------------------------------------------
  ! Determine dust albedo at Lya wavelength (albedo) and
  ! g parameter (g_dust) of the Henyey-Greenstein phase function for dust scattering:
  ! g_dust=0.68 &  albedo=0.33 from Draine 2003 for MW dust, R_V=3.1
  ! albedo=0.40 from Gordon et al. 1997 (ApJ, 487, 625) for SMC dust
  ! ---------------------------------------------------------------------------------

  public get_tau_dust, scatter_dust, read_dust_params, print_dust_params

contains


  function get_tau_dust(ndust_cell, distance_to_border_cm)

    real(kind=8),intent(in) :: ndust_cell,distance_to_border_cm
    real(kind=8)            :: get_tau_dust
    get_tau_dust = sigmad * ndust_cell * distance_to_border_cm

    return

  end function get_tau_dust



  subroutine scatter_dust(v,nu_cell,k,nu_ext,iran,ilost)

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
    ! _ ilost    : updated status of the photon telling if photon has been absorbed
    ! ---------------------------------------------------------------------------------

    real(kind=8), intent(inout)               :: nu_cell, nu_ext ! nu_cell in RASCAS = nu_int in MCLya
    real(kind=8), dimension(3), intent(inout) :: k
    real(kind=8), dimension(3), intent(in)    :: v
    integer, intent(inout)                    :: iran
    real(kind=8)                              :: mu, aors, scalar
    real(kind=8), dimension(3)                :: knew
    integer(kind=4),intent(out)               :: ilost

#ifdef DEBUG
    print *,'-DEBUG- scatter_dust routine, arguments =',v,nu_cell,k,nu_ext,iran
#endif

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
       nu_ext = nu_cell/(1. - scalar/clight)
       k = knew
    end if

  end subroutine scatter_dust



  subroutine read_dust_params(pfile)

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
          case ('grain_radius')
             read(value,*) grain_radius
          end select
       end do
    end if
    close(10)

    ! compute things derived from parameters 
    sigmad = (pi * grain_radius**2)/(1.-albedo)  ! cross section [cm]

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
       write(unit,'(a,ES9.3)') '  albedo       = ',albedo
       write(unit,'(a,ES9.3)') '  g_dust       = ',g_dust
       write(unit,'(a,ES9.3)') '  grain_radius = ',grain_radius
    else
       write(*,'(a,a,a)') '[dust]'
       write(*,'(a,ES9.3)') '  albedo       = ',albedo
       write(*,'(a,ES9.3)') '  g_dust       = ',g_dust
       write(*,'(a,ES9.3)') '  grain_radius = ',grain_radius
    end if

    return

  end subroutine print_dust_params


end module module_dust_model
