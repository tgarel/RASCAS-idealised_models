module module_lyc_model

  ! This module describes the absorption of photons by a general element
  
  use module_constants
  use module_utils
  use module_uparallel
  use module_random
  
  implicit none

  private

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [LyC] of the parameter file
  ! --------------------------------------------------------------------------
  logical                       :: analytic                     ! Analytical computation of HI photoionization cross-sections. Otherwise uses Verner+96

  !--CORESKIP--
  logical                       :: HI_core_skip    = .false.    ! if true, skip scatterings in the core of the line (as in Smith+15).
  real(kind=8)                  :: xcritmax        = 1d10       ! core-skipping will truncate at min(xcrit, xcritmax) -> set to a low value to activate. 
  !--PIKSEROC--
  ! --------------------------------------------------------------------------


  public :: get_tau_gas, read_line_params, print_line_params

  public :: HI_core_skip

  
contains


  function get_tau_gas(nHI, nHeI, nHeII, distance_to_border_cm, nu_cell)

    ! --------------------------------------------------------------------------
    ! compute optical depth of the ionizing photon over a given distance
    ! --------------------------------------------------------------------------
    ! INPUTS:
    ! - nHI,nHeI,nHeII        : number density of HI, HeI and HeII           [ cm^-3 ]
    ! - distance_to_border_cm : distance over which we compute tau           [ cm ]
    ! - nu_cell               : photon's frequency in the frame of the cell  [ Hz ]
    ! OUTPUT :
    ! - get_tau_gas           : optical depth of ionizing photon over distance_to_border_cm
    ! --------------------------------------------------------------------------
    
    real(kind=8),intent(in)                       :: nHI,nHeI,nHeII,distance_to_border_cm,nu_cell
    real(kind=8)                                  :: get_tau_gas, nu1, eps, E, E_erg, x, y, tau_hi, tau_hei, tau_heii
    integer(kind=4)                               :: i

    tau_hi = 0d0 ; tau_hei = 0d0 ; tau_heii = 0d0

    E = hp*nu_cell/evtoerg

    if(E>13.598) then
       if(analytic) then
          nu1 = 13.598 * evtoerg / hp
          eps = sqrt(nu_cell/nu1 - 1)
          tau_hi = 6.3e-18 * (nu1/nu_cell)**4 * exp(4 - 4*atan(eps)/eps) / (1 - exp(-2*pi/eps))
       else
          x = E/4.298e-01
          tau_hi = 5.475e-14*( (x-1)**2 ) * x**(0.5*2.963 - 5.5) / (1+sqrt(x/3.288e1))**2.963
       end if

       if(E>24.587) then
          x = E/13.61 - 4.434d-1
          y = sqrt(x**2 + 2.136**2)
          tau_hei = 9.492e-16 * ( (x-1)**2 + 2.039 ) * y**(0.5*3.188 - 5.5) / (1+sqrt(y/1.469))**3.188

          if(E>54.4178) then
             if(analytic) then
                nu1 = 54.4178 * evtoerg / hp
                eps = sqrt(nu_cell/nu1 - 1)
                tau_heii = 6.3e-18 / 4 * (nu1/nu_cell)**4 * exp(4 - 4*atan(eps)/eps) / (1 - exp(-2*pi/eps))
             else
                x = E/1.72
                tau_heii = 1.369e-14*( (x-1)**2 ) * x**(0.5*2.963 - 5.5) / (1+sqrt(x/3.288e1))**2.963
             end if
          end if
       end if
    end if

    get_tau_gas = (tau_hi*nHI + tau_hei*nHeI + tau_heii*nHeII) * distance_to_border_cm

    return

  end function get_tau_gas
 

  
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
          case('analytic')
             read(value,*) analytic
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
       write(unit,'(a,a,a)')       '[LyC]'
       write(unit,'(a,L1)')        '  analytic        = ',analytic
       write(unit,'(a,L1)')        '  HI_core_skip    = ',HI_core_skip
       if(HI_core_skip) write(unit,'(a,ES12.5)')    '  xcritmax        = ',xcritmax
       write(*,'(a)')             ' '
       call print_uparallel_params(unit)
    else
       write(*,'(a,a,a)')       '[LyC]'
       write(*,'(a,L1)')        '  analytic        = ',analytic
       write(*,'(a,L1)')        '  HI_core_skip    = ',HI_core_skip
       if(HI_core_skip) write(*,'(a,ES12.5)')    '  xcritmax        = ',xcritmax
       write(*,'(a)')             ' '
       call print_uparallel_params()
    end if
    
    return
    
  end subroutine print_line_params


end module module_lyc_model
