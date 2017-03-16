program PhotonsFromStars
  
  ! generate photons emitted by star particles within a given domain

  use module_photon
  use module_utils
  use module_domain
  use module_random
  use module_constants
  use module_ramses

  implicit none
  
  type(domain)    :: emission_domain
  integer         :: narg
  character(2000) :: parameter_file
  real(kind=8),allocatable :: star_pos(:,:),star_age(:),star_mass(:),star_vel(:,:),star_met(:)
  integer(kind=4) :: iran,i,nstars
  integer(kind=8) :: ilast,j
  real(kind=8) :: scalar,nu,r1,r2
  type(photon_init),dimension(:),allocatable :: photgrid
  ! for analysis purposes (a posteriori weighting) we want to save the emitter-frame
  ! frequency (here the freq. in the emitting stellar particle's frame)
  real(kind=8), allocatable :: nu_star(:)
  ! SED-related variables
  integer(kind=4) :: sed_nage,sed_nmet,imet,iage
  integer(kind=8) :: nflux,n
  real(kind=8),allocatable :: sed_age(:),sed_met(:),sweight(:),sed_nphot(:,:),sed_F_0(:,:),sed_beta(:,:)
  integer(kind=8), allocatable :: cum_flux_prob(:)
  real(kind=8),allocatable :: star_beta(:) 
  real(kind=8) :: total_flux,minflux,check_flux,f0,beta,betaplus2,sigma_nu,lambda
  logical:: ok
  
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [PhotonsFromStars] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: outputfile = 'PhotICs.dat' ! file to which outputs will be written
  character(2000)           :: repository = './'          ! ramses run directory (where all output_xxxxx dirs are).
  integer(kind=4)           :: snapnum = 1                ! ramses output number to use

  ! --- domain whithin which star particles will be selected (should be within computational domain used for RT). 
  character(10)             :: star_dom_type      = 'sphere'         ! shape type of domain  // default is sphere.
  real(kind=8),dimension(3) :: star_dom_pos       = (/0.5,0.5,0.5/)  ! center of domain [code units]
  real(kind=8)              :: star_dom_rsp       = 0.3              ! radius of spher [code units]
  real(kind=8)              :: star_dom_size      = 0.3              ! size of cube [code units]
  real(kind=8)              :: star_dom_rin       = 0.0              ! inner radius of shell [code units]
  real(kind=8)              :: star_dom_rout      = 0.3              ! outer radius of shell [code units]
  real(kind=8)              :: star_dom_thickness = 0.1              ! thickness of slab [code units]
  
  ! --- define star particles luminosities (in nb of photons)
  ! Two options for weight_type:
  !   - weight_type=='PowLaw' : we assume each star particle emits F_lbda = F_0 (lbda / lbda_0)**beta (in erg/s/A),
  !                             with F_0 and beta given in the weight_input_file for each age and metallicity (and
  !                             lbda_0 also read from that file). We use the wavelength range specified below (spec_powlaw_lmin_ang
  !                             and spec_powlaw_lmax_Ang) to compute the integral of F_lbda for each star particle as its weight.
  !   - weight_type=='Mono'   : Here the weight_input_file simply provides a number of photons (per sec.) emitted as a function of
  !                             age and metallicity, and we use this as a weight.
  character(30)             :: weight_type       = 'PowLaw'  ! May be 'PowLaw', 'Mono'
  character(2000)           :: weight_input_file = 'F1600.txt' ! file containing weights from SEDs
  real(kind=8)              :: weight_l0_Ang ! this is the lbda_0 above or the monochromatic wavelength, and is read from weight file 

  ! --- define how star particles emit (i.e. the star-particle-frame spectral shape)
  ! Three options here :
  ! - spec_type=='Mono'  : we emit all photons at the same wavelength (in star's frame)
  ! - spec_type=='Gauss' : we sample a Gaussian distribution ...
  ! - spec_type=='PowLaw' : we sample a power-law continuum between two wavelengths. 
  character(30)             :: spec_type = 'Gauss' ! May be 'Mono', 'Gauss', 'PowLaw' ...   
  ! parameters for spec_type == 'Mono'
  real(kind=8)              :: spec_mono_lambda0_Ang = 1215.6701  ! emission wavelength [A] -> this is the parameter read from file. 
  real(kind=8)              :: spec_mono_nu0          ! emission frequency [Hz] -> computed from spec_mono_lambda0_Ang
  ! parameters for spec_type == 'Gauss'
  real(kind=8)              :: spec_gauss_lambda0_Ang = 1215.6701 ! emission wavelength at center [A] -> read from file.
  real(kind=8)              :: spec_gauss_nu0         ! central frequency [Hz] -> computed from spec_mono_lambda0_Ang
  real(kind=8)              :: spec_gauss_sigma_kms = 10.0                  ! line width in velocity [km/s] -> read from file. 
  ! ------ spec_type == 'PowLaw' : a power-law fit to continuum of each star particle, vs. its age and met.
  real(kind=8)              :: spec_powlaw_lmin_Ang = 1120.    ! min wavelength to sample (should be in the range where fit was made ...)
  real(kind=8)              :: spec_powlaw_lmax_Ang = 1320.    ! max ...
  
  ! --- miscelaneous
  integer(kind=4)           :: nphot   = 1000000      ! number of photons to generate
  integer(kind=4)           :: ranseed = -100         ! seed for random generator
  logical                   :: verbose = .true.
  logical                   :: cosmo   = .true.         ! cosmo flag
  ! --------------------------------------------------------------------------



  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: PhotonsFromStars path/to/params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_PhotonsFromStars_params(parameter_file)
  if (verbose) call print_PhotonsFromStars_params
  ! ------------------------------------------------------------


  ! --------------------------------------------------------------------------------------
  ! define domain within which stars may shine
  ! --------------------------------------------------------------------------------------
  select case(star_dom_type)
  case('sphere')
     call domain_constructor_from_scratch(emission_domain,star_dom_type, &
          xc=star_dom_pos(1),yc=star_dom_pos(2),zc=star_dom_pos(3),r=star_dom_rsp)
  case('shell')
     call domain_constructor_from_scratch(emission_domain,star_dom_type, &
          xc=star_dom_pos(1),yc=star_dom_pos(2),zc=star_dom_pos(3),r_inbound=star_dom_rin,r_outbound=star_dom_rout)
  case('cube')
     call domain_constructor_from_scratch(emission_domain,star_dom_type, & 
          xc=star_dom_pos(1),yc=star_dom_pos(2),zc=star_dom_pos(3),size=star_dom_size)
  case('slab')
     call domain_constructor_from_scratch(emission_domain,star_dom_type, &
          xc=star_dom_pos(1),yc=star_dom_pos(2),zc=star_dom_pos(3),thickness=star_dom_thickness)
  end select
  ! --------------------------------------------------------------------------------------

  
  ! --------------------------------------------------------------------------------------
  ! read star particles within domain
  ! --------------------------------------------------------------------------------------
  if (verbose) write(*,*) '> reading star particles'
  call ramses_read_stars_in_domain(repository,snapnum,emission_domain,star_pos,star_age,star_mass,star_vel,star_met,cosmo)
  ! --------------------------------------------------------------------------------------


  ! --------------------------------------------------------------------------------------
  ! Compute star-particle weights
  ! --------------------------------------------------------------------------------------
  ! read weight files (tables generated from e.g. Bruzual & Charlot models)
  open(unit=15,file=weight_input_file,status='old',form='formatted')
  read(15,*) ! skip header
  read(15,*) ! skip header
  read(15,*) ! skip header
  read(15,*) sed_nage,sed_nmet
  allocate(sed_age(sed_nage),sed_met(sed_nmet))
  read(15,*) sed_age ! in Myr
  read(15,*) sed_met ! metallicities (absolute)
  read(15,*) weight_l0_Ang  ! ref. wavelength in A
  select case (trim(weight_type))
  case ('Mono')
     allocate(sed_nphot(sed_nage,sed_nmet))
     do imet = 1,sed_nmet
        read(15,*) sed_nphot(:,imet)  ! nb phots (at weight_l0_Ang) per sec, per solar mass. 
     end do
  case ('PowLaw')
     allocate(sed_F_0(sed_nage,sed_nmet),sed_beta(sed_nage,sed_nmet))
     do imet = 1,sed_nmet
        read(15,*) sed_F_0(:,imet)  ! normalisation of power law (in erg/s/A/Msun)
        read(15,*) sed_beta(:,imet) ! slope of power law
     end do
  case default
     write(*,'(a,a,a)') 'weight_type ',trim(weight_type),' not implemented ... '
     stop
  end select
  close(15)
  ! compute the weight of each star particle
  nstars = size(star_age)
  print*,'nstars = ',nstars
  allocate(sweight(nstars))
  if (trim(weight_type) == 'PowLaw') allocate(star_beta(nstars))
  do i = 1,nstars
     ! pick SED with closest metallicity and age (no interpolation for now)
     call locatedb(sed_met,sed_nmet,star_met(i),imet)
     if (imet < 1) imet = 1
     call locatedb(sed_age,sed_nage,star_age(i),iage)
     if (iage < 1) iage = 1
     ! correct mass for feedback (we need mass of stars formed)
     sweight(i) = star_mass(i) / 1.989d33  ! M_sun
     if (sed_age(iage) < 10.) then ! SNs go off at 10Myr ... 
        sweight(i) = sweight(i)/0.8  !! correct for recycling ... we want the mass of stars formed ...
     end if
     ! compute luminosity
     select case (trim(weight_type))
     case('Mono')
        sweight(i) = sweight(i) * sed_nphot(iage,imet)  ! nb of photons per sec. 
     case('PowLaw')
        ! integrate powerlaw from min to max wavelengths.
        ! We actually want to sample the number of photons per lambda bin.
        ! Given that F_lbda = F_0 (lbda / lbda_0)**beta (in erg/s/A),
        ! the number of photons (in /s/A) is N_lbda = F_0*lbda_0/hc * (lbda/lbda_0)**(1+beta).
        ! (NB: the first lbda_0 here has to be in cm)
        ! This integrates to (in #/s) :
        ! (F_0 lbda_0 / hc) * lbda_0/(beta+2)  * [ (lbda_max/lbda_0)**(2+beta) - (lbda_min/lbda_0)**(2+beta)]
        ! (NB: the first lbda_0 here is in cm, the second in A). 
        ! OR, if beta == -2, the integral is
        ! (F_0*lbda_0/hc) * lbda_0 * ln(lbda_max/lbda_min)     [again, first lbda_0 in cm, second in A]
        f0   = sed_F_0(iage,imet)
        beta = sed_beta(iage,imet)
        if (beta == -2.0d0) then
           sweight(i) = sweight(i) * (f0*weight_l0_Ang*1e-8/planck/clight)
           sweight(i) = sweight(i) * weight_l0_Ang * log(spec_powlaw_lmax_Ang/spec_powlaw_lmin_Ang)
        else
           sweight(i) = sweight(i) * (f0*weight_l0_Ang*1e-8*weight_l0_Ang/planck/clight/(2.+beta))
           sweight(i) = sweight(i) * ( (spec_powlaw_lmax_Ang/weight_l0_Ang)**(2.+beta) - (spec_powlaw_lmin_Ang/weight_l0_Ang)**(2.+beta) )
        end if
        ! -> sweight is the number of photons in [lbda_min;lbda_max]
        star_beta(i) = beta ! keep for later. 
     end select
  end do
  ! --------------------------------------------------------------------------------------

  
  ! --------------------------------------------------------------------------------------
  ! define linear sampling of number of photons 
  ! --------------------------------------------------------------------------------------
  ! compute the total number of photons emitted per second by the sources
  total_flux = 0.0d0
  do i=1,nstars
     total_flux = total_flux + sweight(i)
  end do
  if (verbose) write(*,*) '> Total luminosity (nb of photons per second): ',total_flux

  ! it may happen that the range of luminosities is too large for linear sampling with our method ... 
  ! In that case we need to ignore faint particles:
  minflux = minval(sweight)
  if (total_flux / minflux > 2d8) minflux = total_flux / 2d8  ! NB: dont go much higher than 1d8 (to stay below a few GB RAM). 
  ! check that we dont loose significant flux by sampling only particles with sweight > minflux
  check_flux = 0.0d0
  do i=1,nstars
     if (sweight(i)>minflux) check_flux = check_flux+sweight(i)
  end do
  if (verbose) write(*,*) '> We sample this fraction of total flux: ',check_flux / total_flux
  if ((total_flux - check_flux) / total_flux > 0.001) then
     print*,'Flux losses > 0.1 percent... change algorithm ...'
     ! debug - stop
  end if
  ! construct the cumulative flux distribution, with enough bins to have the smallest star-particle flux in a bin. 
  allocate(cum_flux_prob(int(3*total_flux / minflux,kind=8)))
  ilast = 1
  do i=1,nstars
     if (sweight(i) > minflux) then 
        n = int(3*sweight(i)/minflux,kind=8)
        cum_flux_prob(ilast:ilast+n) = i
        ilast = ilast + n
     end if
  end do
  nflux = ilast
  print*,'nflux, size(cum_fllux_prob):', nflux, size(cum_flux_prob)
  ! --------------------------------------------------------------------------------------

  
  ! --------------------------------------------------------------------------------------
  ! now we can draw integers from 1 to nflux and assign photons to stars ...
  ! --------------------------------------------------------------------------------------
  allocate(photgrid(nphot),nu_star(nphot))
  iran = ranseed
  do i = 1,nphot
     j = int(ran3(iran)*nflux,kind=8)+1
     if (j > nflux) j = nflux
     j = cum_flux_prob(j) 
     photgrid(i)%ID    = i
     photgrid(i)%x_em  = star_pos(:,j)
     photgrid(i)%iran  = -i 
     call isotropic_direction(photgrid(i)%k_em,iran)
     select case(trim(spec_type))
     case('PowLaw')
        ! sample F_lbda = F_0 (lbda / lbda_0)**beta (in erg/s/A) ...
        ! -> we actually want to sample the nb of photons : N_lbda = F_lbda * lbda / (hc) = F_0*lbda_0/(h*c) * (lbda/lbda_0)**(beta+1)
        ! FOR BETA /= 2 : 
        ! -> the probability of drawing a photon with l in [lbda_min;lbda] is:
        !      P(<lbda) = (lbda**(2+beta) - lbda_min**(2+beta))/(lbda_max**(2+beta)-lbda_min**(2+beta))
        ! -> and thus for a random number x in [0,1], we get
        !      lbda = [ lbda_min**(2+beta) + x * ( lbda_max**(2+beta) - lbda_min**(2+beta) ) ]**(1/(2+beta))
        ! FOR BETA == 2:
        ! -> the probability of drawing a photon with l in [lbda_min;lbda] is:
        !      P(<lbda) = log(lbda/lbda_min) / log(lbda_max/lbda_min)
        ! -> and thus for a random number x in [0,1], we get
        !      lbda = lbda_min * exp[ x * log(lbda_max/lbda_min)] 
        r1   = ran3(iran)
        if (star_beta(j) == -2.0d0) then
           nu = spec_powlaw_lmin_Ang * exp(r1 * log(spec_powlaw_lmax_Ang / spec_powlaw_lmin_Ang) ) ! this is lbda [A]
           nu = clight / (nu*1e-8) ! this is freq. [Hz]
        else
           betaplus2 = star_beta(j) + 2.0d0
           nu   = (spec_powlaw_lmin_Ang**betaplus2 + r1 * (spec_powlaw_lmax_Ang**betaplus2 - spec_powlaw_lmin_Ang**betaplus2))**(1./betaplus2) ! this is lbda [A]
           nu   = clight / (nu*1e-8) ! this is freq. [Hz]
        end if
     case('Gauss')
!!$        ! rejection method 
!!$        ok = .false.
!!$        sigma_nu = (spec_gauss_sigma_kms * 1d5 * spec_gauss_nu0 / clight)
!!$        do while (.not. ok)
!!$           nu = clight / (ran3(iran) * 20. - 10. + spec_gauss_lambda0_Ang) * 1e8
!!$           r1 = ran3(iran)
!!$           r2 = exp( -(nu - spec_gauss_nu0)**2 / 2./sigma_nu**2)
!!$           if (r1 <= r2) ok = .true. 
!!$        end do
        ! ohter method (not convincingly working ...) 
        r1 = ran3(iran)
        r2 = ran3(iran)
        nu = sqrt(-2.*log(r1)) * cos(2.0d0*pi*r2)
        nu = (spec_gauss_sigma_kms * 1d5 * spec_gauss_nu0 / clight) * nu + spec_gauss_nu0
     case('Mono')
        nu = spec_mono_nu0
     end select
     nu_star(i) = nu  ! star-particle-frame frequency
     ! now put in external frame using particle's velocity. 
     scalar = photgrid(i)%k_em(1)*star_vel(1,j) + photgrid(i)%k_em(2)*star_vel(2,j) + photgrid(i)%k_em(3)*star_vel(3,j)
     photgrid(i)%nu_em = nu / (1d0 - scalar/clight)
  end do
  ! --------------------------------------------------------------------------------------


  ! --------------------------------------------------------------------------------------
  ! write ICs
  ! --------------------------------------------------------------------------------------
  if (verbose) write(*,*) '> writing file'
  open(unit=14, file=trim(outputfile), status='unknown', form='unformatted', action='write')
  write(14) nphot      ! nb of MC photons 
  write(14) total_flux ! nb of real photons (per sec).
  write(14) ranseed
  write(14) (photgrid(i)%ID,i=1,nphot)
  write(14) (photgrid(i)%nu_em,i=1,nphot)
  write(14) (photgrid(i)%x_em(:),i=1,nphot)
  write(14) (photgrid(i)%k_em(:),i=1,nphot)
  write(14) (photgrid(i)%iran,i=1,nphot)
  write(14) (nu_star(i),i=1,nphot)
  close(14)
  ! --------------------------------------------------------------------------------------


  ! --------------------------------------------------------------------------------------
  ! deallocations 
  ! --------------------------------------------------------------------------------------
  deallocate(star_pos,star_vel,star_mass,star_age,star_met)
  deallocate(sed_age,sed_met)
  select case(trim(weight_type))
  case('Mono')
     deallocate(sed_nphot)
  case('PowLaw')
     deallocate(sed_F_0,sed_beta)
  end select
  deallocate(sweight)
  if (trim(weight_type) == 'PowLaw') deallocate(star_beta)
  deallocate(cum_flux_prob,photgrid,nu_star)
  ! --------------------------------------------------------------------------------------

  
contains
  
  subroutine read_PhotonsFromStars_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
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
       if (line(1:18) == '[PhotonsFromStars]') then
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
          case ('outputfile')
             write(outputfile,'(a)') trim(value)
          case ('repository')
             write(repository,'(a)') trim(value)
          case ('snapnum')
             read(value,*) snapnum
          case ('star_dom_type')
             write(star_dom_type,'(a)') trim(value)
          case ('star_dom_pos')
             read(value,*) star_dom_pos(1),star_dom_pos(2),star_dom_pos(3)
          case ('star_dom_rsp')
             read(value,*) star_dom_rsp
          case ('star_dom_size')
             read(value,*) star_dom_size
          case ('star_dom_rin')
             read(value,*) star_dom_rin
          case ('star_dom_rout')
             read(value,*) star_dom_rout
          case ('star_dom_thickness')
             read(value,*) star_dom_thickness
          case ('weight_type')
             write(weight_type,'(a)') trim(value)
          case ('weight_input_file')
             write(weight_input_file,'(a)') trim(value)
          case ('spec_type')
             write(spec_type,'(a)') trim(value)
          case ('spec_mono_lambda0_Ang')
             read(value,*) spec_mono_lambda0_Ang
             spec_mono_nu0  = clight/spec_mono_lambda0_Ang*1d8
          case ('spec_gauss_lambda0_Ang')
             read(value,*) spec_gauss_lambda0_Ang
             spec_gauss_nu0 = clight / spec_gauss_lambda0_Ang * 1d8
          case ('spec_gauss_sigma_kms')
             read(value,*) spec_gauss_sigma_kms
          case ('spec_powlaw_lmin_Ang')
             read(value,*) spec_powlaw_lmin_Ang
          case ('spec_powlaw_lmax_Ang')
             read(value,*) spec_powlaw_lmax_Ang
          case ('nphot')
             read(value,*) nphot
          case ('ranseed')
             read(value,*) ranseed
          case ('verbose')
             read(value,*) verbose
          case ('cosmo')
             read(value,*) cosmo
          case default
             write(*,'(a,a,a)') 'WARNING: parameter ',trim(name),' unknown '
          end select
       end do
    end if
    close(10)

    return

  end subroutine read_PhotonsFromStars_params

  
  subroutine print_PhotonsFromStars_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)')         '[PhotonsFromStars]'
       write(unit,'(a)')             '# input / output parameters'
       write(unit,'(a,a)')           '  outputfile      = ',trim(outputfile)
       write(unit,'(a,a)')           '  repository      = ',trim(repository)
       write(unit,'(a,i5)')          '  snapnum         = ',snapnum
       write(unit,'(a)')             '# computational domain parameters'
       write(unit,'(a,a)')           '  star_dom_type      = ',trim(star_dom_type)
       write(unit,'(a,3(ES9.3,1x))') '  star_dom_pos       = ',star_dom_pos(1),star_dom_pos(2),star_dom_pos(3)
       select case (trim(star_dom_type))
       case ('sphere')
          write(unit,'(a,ES9.3)')       '  star_dom_rsp       = ',star_dom_rsp
       case ('shell')
          write(unit,'(a,ES9.3)')       '  star_dom_rin       = ',star_dom_rin
          write(unit,'(a,ES9.3)')       '  star_dom_rout      = ',star_dom_rout
       case('cube')
          write(unit,'(a,ES9.3)')       '  star_dom_size      = ',star_dom_size
       case('slab')
          write(unit,'(a,ES9.3)')       '  star_dom_thickness = ',star_dom_thickness
       end select
       write(unit,'(a)')             '# Particle weights '
       write(unit,'(a,a)')           '  weight_type        = ',trim(weight_type)
       write(unit,'(a,a)')           '  weight_input_file  = ',trim(weight_input_file)
       write(unit,'(a)')             '# Spectral shape '
       write(unit,'(a,a)')           '  spec_type               = ',trim(spec_type)
       select case(trim(spec_type))
       case('Mono')
          write(unit,'(a,es9.3,a)')     '  spec_mono_lambda0_Ang   = ',spec_mono_lambda0_Ang, ' ! [A]' 
       case('Gauss')
          write(unit,'(a,es9.3,a)')     '  spec_gauss_lambda0_Ang  = ',spec_gauss_lambda0_Ang, ' ! [A]' 
          write(unit,'(a,es9.3,a)')     '  spec_gauss_sigma_kms = ',spec_gauss_sigma_kms, ' ! [km/s]'
       case('PowLaw')
          write(unit,'(a,es9.3,a)')     '  spec_powlaw_lmin_Ang    = ',spec_powlaw_lmin_Ang, ' ! [A]' 
          write(unit,'(a,es9.3,a)')     '  spec_powlaw_lmax_Ang    = ',spec_powlaw_lmax_Ang, ' ! [A]'
       end select
       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,i8)')          '  nphot           = ',nphot
       write(unit,'(a,i8)')          '  ranseed         = ',ranseed
       write(unit,'(a,L1)')          '  verbose         = ',verbose
       write(unit,'(a,L1)')          '  cosmo           = ',cosmo
       write(unit,'(a)')             ' '
    else
       write(*,'(a,a,a)')         '[PhotonsFromStars]'
       write(*,'(a)')             '# input / output parameters'
       write(*,'(a,a)')           '  outputfile      = ',trim(outputfile)
       write(*,'(a,a)')           '  repository      = ',trim(repository)
       write(*,'(a,i5)')          '  snapnum         = ',snapnum
       write(*,'(a)')             '# computational domain parameters'
       write(*,'(a,a)')           '  star_dom_type      = ',trim(star_dom_type)
       write(*,'(a,3(ES9.3,1x))') '  star_dom_pos       = ',star_dom_pos(1),star_dom_pos(2),star_dom_pos(3)
       select case (trim(star_dom_type))
       case ('sphere')
          write(*,'(a,ES9.3)')       '  star_dom_rsp       = ',star_dom_rsp
       case ('shell')
          write(*,'(a,ES9.3)')       '  star_dom_rin       = ',star_dom_rin
          write(*,'(a,ES9.3)')       '  star_dom_rout      = ',star_dom_rout
       case('cube')
          write(*,'(a,ES9.3)')       '  star_dom_size      = ',star_dom_size
       case('slab')
          write(*,'(a,ES9.3)')       '  star_dom_thickness = ',star_dom_thickness
       end select
       write(*,'(a)')             '# Particle weights '
       write(*,'(a,a)')           '  weight_type        = ',trim(weight_type)
       write(*,'(a,a)')           '  weight_input_file  = ',trim(weight_input_file)
       write(*,'(a)')             '# Spectral shape '
       write(*,'(a,a)')           '  spec_type               = ',trim(spec_type)
       select case(trim(spec_type))
       case('Mono')
          write(*,'(a,es9.3,a)')     '  spec_mono_lambda0_Ang   = ',spec_mono_lambda0_Ang, ' ! [A]' 
       case('Gauss')
          write(*,'(a,es9.3,a)')     '  spec_gauss_lambda0_Ang  = ',spec_gauss_lambda0_Ang, ' ! [A]' 
          write(*,'(a,es9.3,a)')     '  spec_gauss_sigma_kms = ',spec_gauss_sigma_kms, ' ! [km/s]'
       case('PowLaw')
          write(*,'(a,es9.3,a)')     '  spec_powlaw_lmin_Ang    = ',spec_powlaw_lmin_Ang, ' ! [A]' 
          write(*,'(a,es9.3,a)')     '  spec_powlaw_lmax_Ang    = ',spec_powlaw_lmax_Ang, ' ! [A]'
       end select
       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,i8)')          '  nphot           = ',nphot
       write(*,'(a,i8)')          '  ranseed         = ',ranseed
       write(*,'(a,L1)')          '  verbose         = ',verbose
       write(*,'(a,L1)')          '  cosmo           = ',cosmo
       write(*,'(a)')             ' '

    end if

    return

  end subroutine print_PhotonsFromStars_params

  
  subroutine locatedb(xx,n,x,j)

    ! subroutine which locates the position j of a value x in an array xx of n elements
    ! NB : here xx is double precision
    
    implicit none
    
    integer(kind=4) ::  n,j,jl,ju,jm
    real(kind=8)    ::  xx(n),x
    
    jl = 0
    ju = n+1
    do while (ju-jl > 1) 
       jm = (ju+jl)/2
       if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    j = jl

    return

  end subroutine locatedb


  
end program PhotonsFromStars


  
