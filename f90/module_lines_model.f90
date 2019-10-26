module module_lines_model

  ! This module describes the absorption of photons by a general element
  
  use module_constants
  use module_utils
  use module_uparallel
  use module_random
  use module_voigt
  
  implicit none

  private

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [CreateDomDump] of the parameter file
  ! --------------------------------------------------------------------------
  integer(kind=4)               :: ion_number = 1               ! Number of elements constituting the gas (including dust ?)
  real(kind=8),allocatable      :: m_ion(:)
  ! absorption
  integer(kind=4),allocatable   :: n_line(:)
  real(kind=8),allocatable      :: lambda(:)                    ! transition wavelength [A]
  real(kind=8),allocatable      :: f(:)                         ! oscillator strength
  real(kind=8),allocatable      :: A_line(:)                    ! spontaneous decay [/s]

  ! potential fluorescent reemissions
  integer(kind=4),allocatable   :: n_fluo(:)                    ! number of fluorescent reemissions
  real(kind=8),allocatable      :: lambda_fluo(:)               ! transition wavelength [A]
  real(kind=8),allocatable      :: A_fluo(:)                    ! spontaneous decay [/s]

  !Hydrogen option
  logical                       :: recoil       = .true.        ! if set to true, recoil effect is computed [default is true]
  logical                       :: isotropic    = .false.       ! if set to true, scattering events will be isotropic [default is false]

  !--CORESKIP--
  logical                       :: HI_core_skip    = .false.    ! if true, skip scatterings in the core of the line (as in Smith+15).
  real(kind=8)                  :: xcritmax        = 1d10       ! core-skipping will truncate at min(xcrit, xcritmax) -> set to a low value to activate. 
  !--PIKSEROC--

  !Should be in gas_composition, but impossible
  !real(kind=8)                  :: vturb           = 0d0
  ! --------------------------------------------------------------------------

  integer(kind=4)               :: n_line_tot
  integer(kind=4),allocatable   :: ion_index(:)
  real(kind=8),allocatable      :: lambda_cm(:)                 ! transition wavelength [cm]
  real(kind=8),allocatable      :: nu(:)                        ! transition frequency [Hz]
  real(kind=8),allocatable      :: sigma_factor(:)              ! multiply by Voigt(x,a)/nu_D to get sigma.
  real(kind=8),allocatable      :: Atot(:)
  real(kind=8),allocatable      :: lambda_fluo_cm(:)            ! transition wavelength [cm]
  real(kind=8),allocatable      :: nu_fluo(:)                   ! transition frequency [Hz]


  public :: get_tau_lines, get_tau_single_line, scatter_line, read_line_params, print_line_params
  !--PEEL--
  public :: line_peeloff_weight
  !--LEEP--

  public :: ion_number, n_line, n_line_tot, ion_index, HI_core_skip

  
contains

    function get_tau_single_line(line_number, n, vth, distance_to_border_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute optical depth of the lines-s over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - n        : number density of ions                                   [ cm^-3 ]
    ! - vth      : thermal (+ small-scale turbulence) velocity of ions      [ cm / s ]
    ! - distance_to_border_cm : distance over which we compute tau          [ cm ]
    ! - nu_cell  : photon's frequency in the frame of the cell              [ Hz ]
    ! OUTPUT :
    ! - get_tau_line : optical depth of Silicon's line over distance_to_border_cm
     ! --------------------------------------------------------------------------
    
     integer(kind=4),intent(in)    :: line_number 
     real(kind=8),intent(in)       :: n,vth,distance_to_border_cm,nu_cell
     real(kind=8)                  :: delta_nu_doppler,x_cell,sigma,a,h_cell,get_tau_single_line
     integer(kind=4)               :: i

     ! if(ion < 1 .or. ion > ion_number) then
     !    print*, 'Problem with the ion number in function get_tau_single_line in module_lines_model, stopping'
     !    stop
     ! end if
     
     if(line_number < 1 .or. line_number > n_line_tot) then
        print*, 'Problem with the line number in function get_tau_single_line in module_lines_model, stopping'
        stop
     end if

     ! compute Doppler width and a-parameter
     delta_nu_doppler = (vth / sqrt(m_ion(ion_index(line_number)))) / lambda_cm(line_number)
     a = A_line(line_number) / (fourpi * delta_nu_doppler)

     ! cross section of SiII-1260.42
     x_cell = (nu_cell - nu(line_number)) / delta_nu_doppler
     h_cell = voigt_function(x_cell,a)
     sigma  = sigma_factor(line_number) / delta_nu_doppler * h_cell

     get_tau_single_line = sigma * n * distance_to_border_cm

     return

   end function get_tau_single_line


  function get_tau_lines(n, vth, distance_to_border_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute optical depth of the lines-s over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - n        : number density of ions                                   [ cm^-3 ]
    ! - vth      : thermal (+ small-scale turbulence) velocity of ions      [ cm / s ]
    ! - distance_to_border_cm : distance over which we compute tau          [ cm ]
    ! - nu_cell  : photon's frequency in the frame of the cell              [ Hz ]
    ! OUTPUT :
    ! - get_tau_line : optical depth of Silicon's line over distance_to_border_cm
    ! --------------------------------------------------------------------------

    real(kind=8),intent(in),dimension(ion_number) :: n
    real(kind=8),intent(in)                       :: vth,distance_to_border_cm,nu_cell
    real(kind=8)                                  :: delta_nu_doppler,x_cell,sigma,a,h_cell,get_tau_lines
    integer(kind=4)                               :: i

    get_tau_lines = 0d0

    do i=1,n_line_tot
      ! ion = get_ion_index(i)
       ! compute Doppler width and a-parameter
       delta_nu_doppler = (vth / sqrt(m_ion(ion_index(i)))) / lambda_cm(i)
       a = A_line(i) / (fourpi * delta_nu_doppler)

       ! cross section of SiII-1260.42
       x_cell = (nu_cell - nu(i)) / delta_nu_doppler
       h_cell = voigt_function(x_cell,a)
       sigma  = sigma_factor(i) / delta_nu_doppler * h_cell

       get_tau_lines = get_tau_lines + sigma * n(ion_index(i)) * distance_to_border_cm
    end do

    return

  end function get_tau_lines


  subroutine scatter_line(line_number,vcell,vth_times_sqrtm,nu_cell,k,nu_ext,iran,xcrit)

    ! ---------------------------------------------------------------------------------
    ! perform scattering event on an ion
    ! The photon is absorbed in transition 1->3 and may decay as 3->1 or 3->2. 
    ! ---------------------------------------------------------------------------------
    ! INPUTS :
    ! - vcell    : bulk velocity of the gas (i.e. cell velocity)       [ cm / s ] 
    ! - vth      : thermal (+turbulent) velocity dispersion of H (and Si) atoms [ cm / s ] 
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

    integer(kind=4),intent(in)              :: line_number
    real(kind=8),intent(inout)              :: nu_cell, nu_ext
    real(kind=8),dimension(3),intent(inout) :: k
    real(kind=8),dimension(3),intent(in)    :: vcell
    real(kind=8),intent(in)                 :: vth_times_sqrtm
    !--CORESKIP--
    real(kind=8),intent(in)                 :: xcrit
    real(kind=8)                            :: xc
    !--PIKSEROC--
    integer(kind=4),intent(inout)           :: iran
    real(kind=8)                            :: delta_nu_doppler, a, x_cell, upar, ruper
    real(kind=8)                            :: r2, uper, nu_atom, mu, bu, scalar, proba, x_atom, vth
    real(kind=8),dimension(3)               :: knew
    integer(kind=4)                         :: i, start

    if(line_number < 1 .or. line_number > n_line_tot) then
       print*, 'Problem with the line number in subroutine scatter_line in module_line_model, stopping'
       stop
    end if

    ! !--CORESKIP--  sanity check ... 
    ! if (.not. HI_core_skip .and. xcrit .ne. 0.0d0) then
    !    print*,'ERROR: core skipping is on but xcrit is not zero ... '
    !    stop
    ! end if
    ! if (HI_core_skip)  then
    !    xc = min(xcrit,xcritmax)
    ! else
    !    xc=0
    ! endif
    ! !--PIKSEROC--


    vth = vth_times_sqrtm/sqrt(m_ion(ion_index(line_number)))
    
    ! define x_cell & a
    delta_nu_doppler = vth / lambda_cm(line_number) 
    a = A_line(line_number)  / fourpi / delta_nu_doppler
    x_cell = (nu_cell - nu(line_number) ) / delta_nu_doppler

    ! 1/ component parallel to photon's propagation
    ! -> get velocity of interacting atom parallel to propagation
    upar = get_uparallel(x_cell,a,iran)
    upar = upar * vth    ! upar is an x -> convert to a velocity 

    ! 2/ component perpendicular to photon's propagation
    ruper  = ran3(iran)
    r2     = ran3(iran)
    !--CORESKIP--
    uper   = sqrt(xc**2-log(ruper))*cos(twopi*r2)
    !--PIKSEROC--
    uper   = uper * vth  ! from x to velocity

    ! 3/ chose de-excitation channel to determine output freq. in atom's frame
    if(n_fluo(line_number) == 0) then
       nu_atom = nu_cell - nu_ext * upar/clight
    else
       r2 = ran3(iran)
       if(r2 < A_line(line_number)/Atot(line_number)) then
          nu_atom = nu_cell - nu_ext * upar/clight
       else
          proba = A_line(line_number)/Atot(line_number)
          start = 1 + sum(n_fluo(1:line_number-1))
          do i=1,n_fluo(line_number)
             if(r2 < A_fluo(start+i-1)/Atot(line_number) + proba) then
                nu_atom = nu_fluo(start+i-1)
                exit
             end if
             proba = proba +  A_fluo(start+i-1)/Atot(line_number)
          end do
       end if
    end if
    
    ! 4/ determine direction of scattered photon
    !if(isotropic) then
       call isotropic_direction(knew,iran)
       mu = k(1)*knew(1) + k(2)*knew(2) + k(3)*knew(3)
       bu = sqrt(1.0d0 - mu*mu)
    ! else
    !    !Assuming it's hydrogen...
    !    x_atom  = (nu_atom -nu(line_number)) / delta_nu_doppler
    !    if (abs(x_atom) < 0.2d0) then ! core scattering 
    !       call anisotropic_direction_HIcore(k,knew,mu,bu,iran)
    !    else ! wing scattering 
    !       call anisotropic_direction_Rayleigh(k,knew,mu,bu,iran)
    !    end if
    ! end if

    ! ! 5/ recoil effect 
    ! if (recoil) then
    !    !Assuming it's hydrogen...
    !    nu_atom = nu_atom / (1.0d0 + ((planck*nu_atom)/(mp*clight*clight))*(1.0d0-mu))
    ! end if

    ! 6/ compute atom freq. in external frame, after scattering
    scalar = knew(1) * vcell(1) + knew(2) * vcell(2) + knew(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + bu*uper)/clight)
    nu_cell = (1.0d0 - scalar/clight) * nu_ext 
    k = knew

  end subroutine scatter_line


  !--PEEL--
  function line_peeloff_weight(line_number,vcell,vth_times_sqrtm,nu_ext,kin,kout,iran)

    ! ---------------------------------------------------------------------------------
    ! Compute probability that a photon coming along kin scatters off in direction kout.
    ! Also update nu_ext to external-frame frequency along kout
    ! ---------------------------------------------------------------------------------
    ! INPUTS :
    ! - vcell    : bulk velocity of the gas (i.e. cell velocity)       [ cm / s ] 
    ! - vth      : thermal (+turbulent) velocity dispersion of H atoms [ cm / s ] 
    ! - nu_ext   : frequency of incoming photon, in external frame     [ Hz ]
    ! - kin      : propagation vector (normalized)
    ! - kout     : direction after interaction (fixed)
    ! - iran     : random number generator seed
    ! OUTPUTS :
    ! - nu_ext   : updated frequency in external frame [ Hz ]
    ! _ iran     : updated value of seed
    ! ---------------------------------------------------------------------------------

    integer(kind=4),intent(in)              :: line_number
    real(kind=8),intent(inout)              :: nu_ext
    real(kind=8),dimension(3),intent(in)    :: kin, kout
    real(kind=8),dimension(3),intent(in)    :: vcell
    real(kind=8),intent(in)                 :: vth_times_sqrtm
    integer(kind=4),intent(inout)           :: iran
    real(kind=8)                            :: line_peeloff_weight
    real(kind=8)                            :: delta_nu_doppler, a, x_cell, upar, ruper
    real(kind=8)                            :: r2, uper, nu_atom, mu, bu, scalar, proba, vth
    real(kind=8)                            :: x_atom,nu_cell
    integer(kind=4)                         :: i, start
    
    if(line_number < 1 .or. line_number > n_line_tot) then
       print*, 'Problem with the line number in subroutine line_peeloff_weight in module_line_model, stopping'
       stop
    end if

    vth = vth_times_sqrtm / sqrt(m_ion(ion_index(line_number)))

    ! compute frequency in cell's frame 
    scalar  = kin(1) * vcell(1) + kin(2) * vcell(2) + kin(3) * vcell(3)
    nu_cell = (1.d0 - scalar/clight) * nu_ext

    ! define x_cell & a
    delta_nu_doppler = vth / lambda_cm(line_number) 
    a = A_line(line_number) / fourpi / delta_nu_doppler
    x_cell = (nu_cell - nu(line_number)) / delta_nu_doppler

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
    if(n_fluo(line_number) == 0) then
       nu_atom = nu_cell - nu_ext * upar/clight
    else
       r2 = ran3(iran)
       if(r2 < A_line(line_number)/Atot(line_number)) then
          nu_atom = nu_cell - nu_ext * upar/clight
       else
          proba = A_line(line_number)/Atot(line_number)
          start = 1 + sum(n_fluo(1:line_number-1))
          do i=1,n_fluo(line_number)
             if(r2 < A_fluo(start+i-1)/Atot(line_number) + proba) then
                nu_atom = nu_fluo(start+i-1)
                exit
             end if
             proba = proba +  A_fluo(start+i-1)/Atot(line_number)
          end do
       end if
    end if
    
    ! 4/ determine direction of scattered photon
   ! if(isotropic) then
       line_peeloff_weight = 0.5d0  ! P(mu) for isotropic phase function
       mu = kin(1)*kout(1) + kin(2)*kout(2) + kin(3)*kout(3)
       bu = sqrt(1.0d0 - mu*mu)
    ! else
    !    !Assuming it's hydrogen...
    !    x_atom  = (nu_atom -nu(line_number)) / delta_nu_doppler
    !    if (abs(x_atom) < 0.2) then ! core scattering 
    !       line_peeloff_weight = anisotropic_probability_HIcore(kin,kout,mu,bu)
    !    else ! wing scattering 
    !       line_peeloff_weight = anisotropic_probability_Rayleigh(kin,kout,mu,bu)
    !    end if
    ! end if

    ! ! 5/ recoil effect 
    ! if (recoil) then
    !    !Assuming it's hydrogen...
    !    nu_atom = nu_atom / (1.d0 + ((planck*nu_atom)/(mp*clight*clight))*(1.-mu))
    ! end if
    
    ! 6/ compute atom freq. in external frame, after scattering
    scalar = kout(1) * vcell(1) + kout(2) * vcell(2) + kout(3)* vcell(3)
    nu_ext = nu_atom * (1.0d0 + scalar/clight + (upar*mu + bu*uper)/clight)

  end function line_peeloff_weight
  !--LEEP--


  subroutine set_ion_index

    integer(kind=4) :: i,k

    allocate(ion_index(n_line_tot))

    k=1
    do i=1,ion_number
       ion_index(k:k+n_line(i)-1) = i
       k = k + n_line(i)
    end do
    !print*, 'set_ion_index : ', ion_index

  end subroutine set_ion_index
 

  
  subroutine read_line_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000)         :: line,name,value
    integer(kind=4)         :: err,i,j,k
    logical                 :: section_present

    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:7) == '[lines]') then
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
          case ('ion_number')
             read(value,*) ion_number
          allocate(m_ion(ion_number), n_line(ion_number))
          case ('m_ion')
             read(value,*) m_ion(:)
          case ('n_line')
             read(value,*) n_line(:)
             n_line_tot = sum(n_line)
             call set_ion_index()
          allocate(lambda(n_line_tot),lambda_cm(n_line_tot),nu(n_line_tot),sigma_factor(n_line_tot),f(n_line_tot),A_line(n_line_tot),Atot(n_line_tot),n_fluo(n_line_tot))
          case ('lambda')
             read(value,*) lambda(:)
             lambda_cm(:) = lambda(:)/cmtoA
             nu(:) = clight/lambda_cm(:)
          case ('f')
             read(value,*) f(:)
             sigma_factor(:) = sqrtpi*e_ch**2*f(:)/me/clight
          case ('A')
             read(value,*) A_line(:)
          case ('n_fluo')
             read(value,*) n_fluo(:)
             allocate(lambda_fluo(sum(n_fluo)), lambda_fluo_cm(sum(n_fluo)), nu_fluo(sum(n_fluo)), A_fluo(sum(n_fluo)))
          case ('lambda_fluo')
             read(value,*) lambda_fluo(:)
             lambda_fluo_cm(:) = lambda_fluo(:)/cmtoA
             nu_fluo(:) = clight/lambda_fluo_cm(:)
          case ('A_fluo')
             read(value,*) A_fluo(:)
             k=1
             do j=1,n_line_tot
                Atot(j) = A_line(j) + sum(A_fluo(k:k+n_fluo(j)-1))
                k = k + n_fluo(j)
             end do
          case ('recoil')
             read(value,*) recoil
          case ('isotropic')
             read(value,*) isotropic
          !--CORESKIP--
          case ('HI_core_skip') 
             read(value,*) HI_core_skip
          case ('xcritmax')
             read(value,*) xcritmax
          !--PIKSEROC--

          end select
       end do
    end if

    close(10)

    call read_uparallel_params(pfile)
    call read_voigt_params(pfile)

    return

  end subroutine read_line_params


  subroutine print_line_params(unit)
    
    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit
    character(100) :: fmt1,fmt2,fmt3,fmt4,fmt5

    if (present(unit)) then
       write(unit,'(a,a,a)')       '[lines]'
       write(unit,'(a,i3)')        '  ion_number      = ',ion_number
       write(fmt1,'(a,i3,a)')  '(a,',ion_number,'(ES12.5,1x))'
       write(fmt2,'(a,i3,a)')   '(a,',ion_number,'(i3))'
       write(unit,fmt1)            '  m_ion           = ',m_ion(:)
       write(unit,fmt2)            '  n_line          = ',n_line(:)
       write(fmt3,'(a,i3,a)')   '(a,',n_line_tot,'(ES12.5,1x))'
       write(fmt4,'(a,i3,a)')   '(a,',n_line_tot,'(i3))'
       write(unit,fmt3)            '  lambda          = ',lambda(:)
       write(unit,fmt3)            '  f               = ',f(:)
       write(unit,fmt3)            '  A_line          = ',A_line(:)
       write(unit,fmt4)            '  n_fluo          = ',n_fluo(:)
       if(sum(n_fluo(:)) > 0) then
          write(fmt5,'(a,i3,a)')   '(a,',sum(n_fluo(:)),'(ES12.5,1x))'
          write(unit,fmt5)            '  lambda_fluo     = ',lambda_fluo(:)
          write(unit,fmt5)            '  A_fluo          = ',A_fluo(:)
       end if
       write(unit,'(a,L1)')        '  recoil          = ',recoil
       write(unit,'(a,L1)')        '  isotropic       = ',isotropic
       write(unit,'(a,L1)')        '  HI_core_skip    = ',HI_core_skip
       if(HI_core_skip) write(unit,'(a,ES12.5)')    '  xcritmax        = ',xcritmax
       write(*,'(a)')             ' '
       call print_uparallel_params(unit)
       call print_voigt_params(unit)
    else
       write(*,'(a,a,a)')       '[lines]'
       write(*,'(a,i3)')        '  ion_number      = ',ion_number
       write(fmt1,'(a,i3,a)')  '(a,',ion_number,'(ES12.5,1x))'
       write(fmt2,'(a,i3,a)')   '(a,',ion_number,'(i3))'
       write(*,fmt1)            '  m_ion           = ',m_ion(:)
       write(*,fmt2)            '  n_line          = ',n_line(:)
       write(fmt3,'(a,i3,a)')   '(a,',n_line_tot,'(ES12.5,1x))'
       write(fmt4,'(a,i3,a)')   '(a,',n_line_tot,'(i3))'
       write(*,fmt3)            '  lambda          = ',lambda(:)
       write(*,fmt3)            '  f               = ',f(:)
       write(*,fmt3)            '  A_line          = ',A_line(:)
       write(*,fmt4)            '  n_fluo          = ',n_fluo(:)
       if(sum(n_fluo(:)) > 0) then
          write(fmt5,'(a,i3,a)')   '(a,',sum(n_fluo(:)),'(ES12.5,1x))'
          write(*,fmt5)            '  lambda_fluo     = ',lambda_fluo(:)
          write(*,fmt5)            '  A_fluo          = ',A_fluo(:)
       end if
       write(*,'(a,L1)')        '  recoil          = ',recoil
       write(*,'(a,L1)')        '  isotropic       = ',isotropic
       write(*,'(a,L1)')        '  HI_core_skip    = ',HI_core_skip
       if(HI_core_skip) write(*,'(a,ES12.5)')    '  xcritmax        = ',xcritmax
       write(*,'(a)')             ' '
       call print_uparallel_params()
       call print_voigt_params()
    end if
    
    return
    
  end subroutine print_line_params


end module module_lines_model
