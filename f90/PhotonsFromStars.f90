program PhotonsFromStars
  
  ! generate photons emitted by star particles within a given domain

  use module_utils
  use module_domain
  use module_random
  use module_constants
  use module_ramses
  use module_ssp_lib
  
  implicit none
  
  type(domain)             :: emission_domain
  character(2000)          :: parameter_file
  real(kind=8),allocatable :: star_pos(:,:),star_age(:),star_mass(:),star_vel(:,:),star_met(:)
  integer(kind=4)          :: iran,i,nstars,narg
  real(kind=8)             :: scalar, r2, r3
  ! for analysis purposes (a posteriori weighting) we want to save the emitter-frame
  ! frequency (here the freq. in the emitting stellar particle's frame)
  real(kind=8),allocatable :: nu_star(:)
  real(kind=8)             :: total_flux
  type(SSPgrid)            :: NdotGrid
  real(kind=8),allocatable :: low_prob(:), low_prob2(:), nu_em(:), Ndot(:), sweight(:), lbin(:)
  real(kind=8),allocatable :: x_em(:,:), k_em(:,:), NdotStar(:,:), v_em(:,:)
  integer(kind=4)          :: ilow, iphot, iseed, ilow2
  real(kind=8)             :: lambda0, k(3), lambdamin, lambdamax, nu, spec_gauss_nu0, lambda_star, weight
  
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
  
  ! --- define how star particles emit (i.e. the star-particle-frame spectral shape)
  ! Four options here :
  ! - spec_type=='Mono'   : we emit all photons at the same wavelength (in star's frame)
  ! - spec_type=='Gauss'  : we sample a Gaussian distribution ...
  ! - spec_type=='PowLaw' : we sample a power-law continuum between two wavelengths. 
  ! - spec_type=='Table'  : we sample a tabulated spectrum. 
  character(30)             :: spec_type = 'Mono'               ! May be 'Mono', 'Gauss', 'PowLaw' ...   
  character(1000)           :: spec_SSPdir = '../libs/SSPlibs/' ! the SSP lib directory
  ! parameters for spec_type == 'Mono'
  real(kind=8)              :: spec_mono_lambda0 = 1216.        ! emission wavelength [A]
  ! parameters for spec_type == 'Gauss'
  real(kind=8)              :: spec_gauss_lambda0 = 1216.       ! central wavelength [A]
  real(kind=8)              :: spec_gauss_sigma_kms = 10.0      ! line width in velocity [km/s] -> read from file. 
  ! parameters for spec_type == 'PowLaw' : a power-law fit to continuum of each star particle, vs. its age and met.
  real(kind=8)              :: spec_powlaw_lmin_Ang = 1120.     ! min wavelength to sample (should be in the range where fit was made ...)
  real(kind=8)              :: spec_powlaw_lmax_Ang = 1320.     ! max ...
  ! parameters for spec_type == 'Table'
  real(kind=8)              :: spec_table_lmin_Ang = 1120.      ! min wavelength to sample
  real(kind=8)              :: spec_table_lmax_Ang = 1320.      ! max ...
  
  ! --- miscelaneous
  integer(kind=4)           :: nphotons = 1000000      ! number of photons to generate
  integer(kind=4)           :: ranseed  = -100         ! seed for random generator
  logical                   :: verbose  = .true.
  ! --- parameters for star particles/feedback in simulation
  logical                   :: recompute_particle_initial_mass = .false.
  real(kind=8)              :: tdelay_SN = 10.        ! [Myr] SNs go off at tdelay_SN ... 
  real(kind=8)              :: recyc_frac = 0.8       ! correct for recycling ... we want the mass of stars formed ...
  
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
  if (verbose) write(*,*) 'Reading star particles'
  call ramses_read_stars_in_domain(repository,snapnum,emission_domain,star_pos,star_age,star_mass,star_vel,star_met)
  ! --------------------------------------------------------------------------------------
  
  
  print*,'Nstars read =',size(star_mass)
  print*,'minmax pos =',minval(star_pos),maxval(star_pos)
  print*,'minmax vel =',minval(star_vel),maxval(star_vel)
  print*,'minmax mass =',minval(star_mass),maxval(star_mass)
  print*,'minmax age =',minval(star_age),maxval(star_age)
  print*,'minmax met =',minval(star_met),maxval(star_met)
  
  
  ! --------------------------------------------------------------------------------------
  ! Compute luminosity/spectrum for each star-particles
  ! --------------------------------------------------------------------------------------
  
  
  call init_ssp_lib(spec_SSPdir)
  
  select case(trim(spec_type))
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
  case('PowLaw')
     print*,'Not implemented yet...'
     stop
     
     ! steps
     ! 1/ linear fit -> get grid of (Fo,Beta) = f(age,met)
     ! 2/ integrate this to get NdotGrid
     ! 3/ age-Z interpolation to get NdotStar
     ! 4/ draw emitting stars from the weights
     ! 5/ draw frequency?
     
     
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
  case('Gauss')

     print*,'Gauss spectral type'
     lambda0 = spec_gauss_lambda0
     call ssp_lib_extract_subset(lambda0, lambda0, NdotGrid)  ! charge les SEDs si pas deja fait, et extrait les flux a weight_l0_Ang
     ! -> def le nb phots (at weight_l0_Ang) per sec, per solar mass -> dans un objet class_twod_table (nAge, nZ, 1)
     print*,spec_gauss_lambda0, lambda0
     spec_gauss_nu0 = clight / (lambda0*1e-8) ! Hz
     
     ! compute the nb of photons per second emitted by each star particle
     allocate(Ndot(1))
     nstars = size(star_age)
     allocate(sweight(nstars))
     do i = 1,nstars
        !print*,i,star_age(i), star_age(i)/1.e3, log10(star_met(i))
        call ssp_lib_interpolate(NdotGrid, star_age(i)/1.e3, log10(star_met(i)), Ndot)    ! Ndot number of photons / s / A / Msun
        sweight(i) = Ndot(1) * star_mass(i) / msun  ! M_sun
        if(recompute_particle_initial_mass)then
           if (star_age(i) < tdelay_SN) then       ! SNs go off at tdelay_SN ... 
              sweight(i) = sweight(i)/recyc_frac   ! correct for recycling ... we want the mass of stars formed ...
           end if
        end if
     end do

     ! calcul pour chaque particule la luminosite inferieure de son bin dans la distribution cumulative.. 
     ! compute the total number of photons emitted per second by the sources
     call compute_cum_low_prob(nstars, sweight, low_prob, total_flux)
     if (verbose) write(*,*) 'Total luminosity (nb of photons per second): ',total_flux

     ! for each photon packet, draw the emitting star
     allocate(x_em(1:3,1:nphotons), k_em(1:3,1:nphotons), nu_em(1:nphotons), nu_star(1:nphotons))
     allocate(v_em(1:3,1:nphotons))
     
     iseed = ranseed
     do iphot = 1,nphotons

        call binary_search(iseed, nstars, low_prob, ilow)
        ! draw photon's ICs from star ilow
        ! give photon the position of the star
        x_em(:,iphot) = star_pos(:,ilow) 
        ! draw propagation direction
        call isotropic_direction(k,iseed)
        k_em(:,iphot) = k
        ! compute frequency in star frame using the Box-Muller method
        r2 = ran3(iran)
        r3 = ran3(iran)
        nu = sqrt(-2.*log(r2)) * cos(2.0d0*pi*r3)
        nu = (spec_gauss_sigma_kms * 1d5 * spec_gauss_nu0 / clight) * nu + spec_gauss_nu0
        nu_star(iphot) =  nu    ! clight / (lambda0*1e-8)  ! Hz
        ! compute frequency in external frame 
        scalar = k(1)*star_vel(1,ilow) + k(2)*star_vel(2,ilow) + k(3)*star_vel(3,ilow)
        nu_em(iphot)  = nu_star(iphot) / (1d0 - scalar/clight)
        ! store velocity of the star, for peeling-off only
        v_em(:,iphot) = star_vel(:,ilow)
     end do

     
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
  case('Table')
     !print*,'Not implemented yet...'
     !stop

     print*,'Table spectral type'
     lambdamin = spec_table_lmin_Ang
     lambdamax = spec_table_lmax_Ang
     call ssp_lib_extract_subset(lambdamin,lambdamax,NdotGrid) ! extract the SSP age-met grid of Nphotons in lambda range [lambdamin;lambdamax]
     allocate(Ndot(NdotGrid%nlambda))
     
     ! compute the nb of photons per second emitted by each star particle
     nstars = size(star_age)
     allocate(sweight(nstars))
     ! need to store Ndot for each star into NdotStar
     allocate(NdotStar(nstars,NdotGrid%nlambda))
     do i = 1,nstars
        ! interpolate NdotGrid(lambda)
        call ssp_lib_interpolate(NdotGrid, star_age(i)/1.e3, log10(star_met(i)), Ndot)
        Ndot = Ndot * star_mass(i) / msun  ! nb of photons / s / A
        if(recompute_particle_initial_mass)then
           if (star_age(i) < tdelay_SN) then       ! SNs go off at tdelay_SN ... 
              Ndot = Ndot/recyc_frac               ! correct for recycling ... we want the mass of stars formed ...
           end if
        end if
        NdotStar(i,:) = Ndot(:)
        ! integrate nphotPerSecPerMsun(nlambda)
        call ssp_lib_integrate(NdotGrid%lambda, Ndot, NdotGrid%nlambda, weight)
        sweight(i) = weight
     enddo
     
     ! calcul pour chaque particule la luminosite inferieure de son bin dans la distribution cumulative.. 
     ! compute the total number of photons emitted per second by the sources
     call compute_cum_low_prob(nstars, sweight, low_prob, total_flux)
     if (verbose) write(*,*) 'Total luminosity (nb of photons per second): ',total_flux
     
     ! for each photon packet, draw the emitting star
     allocate(x_em(1:3,1:nphotons), k_em(1:3,1:nphotons), nu_em(1:nphotons), nu_star(1:nphotons))
     allocate(v_em(1:3,1:nphotons))
     allocate(low_prob2(NdotGrid%nlambda+1))
     iseed = ranseed

     ! compute lbin
     allocate(lbin(NdotGrid%nlambda))
     lbin(1) = NdotGrid%lambda(1)
     do i = 2,NdotGrid%nlambda-1
        lbin(i) = (NdotGrid%lambda(i+1) + NdotGrid%lambda(i))/2.0d0
     enddo
     lbin(NdotGrid%nlambda) = NdotGrid%lambda(NdotGrid%nlambda)

     do iphot = 1,nphotons

        call binary_search(iseed, nstars, low_prob, ilow)

        ! photon emitted from star ilow
        ! give photon the position of the star
        !print*,iphot,ilow
        x_em(:,iphot) = star_pos(:,ilow) 
        ! draw propagation direction
        call isotropic_direction(k,iseed)
        k_em(:,iphot) = k

        ! compute frequency in star frame... here is the difficulty....

        ! 1/ find the frequency/lambda bin
        low_prob2(1) = 0.0d0
        low_prob2(2) = NdotStar(ilow,1) * (NdotGrid%lambda(2) - NdotGrid%lambda(1))/2.  ! => P(lambda(1)
        do i = 3,NdotGrid%nlambda
           low_prob2(i) = low_prob2(i-1) + NdotStar(ilow,i-1) * (NdotGrid%lambda(i) - NdotGrid%lambda(i-2))/2.
        end do
        low_prob2 = low_prob2 / low_prob2(NdotGrid%nlambda)
        low_prob2(NdotGrid%nlambda+1) = 1.1  ! higher than upper limit 
        call binary_search(iseed, NdotGrid%nlambda, low_prob2, ilow2)
        ! photon emitted in the bin ilow2;ilow2+1
        ! => photon emitted at Ndot(ilow2), which means lambda between lambda(ilow2)+lambda(ilow2-1)/2. and  lambda(ilow2+1)+lambda(ilow2)/2. 
        ! if ilow2 = 1 => lambda between lambda(1) and  lambda(2)+lambda(1)/2.
        ! if ilow2 = nlambda => lambda between lambda(ilow2)+lambda(ilow2-1)/2. and lambda(ilow2)
        
        ! 2/ get lamba_em
        ! 0th order solution is a flat distribution in this bin
        r2 = ran3(iseed)
        lambda_star = lbin(ilow2) + r2 * (lbin(ilow2+1)-lbin(ilow2))
        
        nu_star(iphot) = clight / (lambda_star*1e-8) ! Hz
        ! compute frequency in external frame 
        scalar = k(1)*star_vel(1,ilow) + k(2)*star_vel(2,ilow) + k(3)*star_vel(3,ilow)
        nu_em(iphot)  = nu_star(iphot) / (1d0 - scalar/clight)

        ! store velocity of the star, for peeling-off only
        v_em(:,iphot) = star_vel(:,ilow)
     enddo
     

! --------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------
  case('Mono')
     
     print*,'Monochromatic spectral type'
     lambda0 = spec_mono_lambda0
     call ssp_lib_extract_subset(lambda0, lambda0, NdotGrid)  ! charge les SEDs si pas deja fait, et extrait les flux a weight_l0_Ang
     ! -> def le nb phots (at weight_l0_Ang) per sec, per solar mass -> dans un objet class_twod_table (nAge, nZ, 1)
     print*,spec_mono_lambda0, lambda0
     allocate(Ndot(1))
     
     ! compute the nb of photons per second emitted by each star particle
     nstars = size(star_age)
     allocate(sweight(nstars))
     do i = 1,nstars
        !print*,i,star_age(i), star_age(i)/1.e3, log10(star_met(i))
        call ssp_lib_interpolate(NdotGrid, star_age(i)/1.e3, log10(star_met(i)), Ndot) ! Ndot number of photons / s / A / Msun
        sweight(i) = Ndot(1) * star_mass(i) / msun  ! number of photons / s / A
        if(recompute_particle_initial_mass)then
           if (star_age(i) < tdelay_SN) then       ! SNs go off at tdelay_SN ... 
              sweight(i) = sweight(i)/recyc_frac   ! correct for recycling ... we want the mass of stars formed ...
           end if
        end if
      end do
     
     ! calcul pour chaque particule la luminosite inferieure de son bin dans la distribution cumulative.. 
     ! compute the total number of photons emitted per second by the sources
     call compute_cum_low_prob(nstars, sweight, low_prob, total_flux)
     if (verbose) write(*,*) 'Total luminosity (nb of photons per second): ',total_flux
     
     ! for each photon packet, draw the emitting star
     allocate(x_em(1:3,1:nphotons), k_em(1:3,1:nphotons), nu_em(1:nphotons), nu_star(1:nphotons))
     allocate(v_em(1:3,1:nphotons))
     
     iseed = ranseed
     do iphot = 1,nphotons

        call binary_search(iseed, nstars, low_prob, ilow)
        ! draw photon's ICs from star ilow

        ! give photon the position of the star
        !print*,iphot,ilow
        x_em(:,iphot) = star_pos(:,ilow) 
        ! draw propagation direction
        call isotropic_direction(k,iseed)
        k_em(:,iphot) = k
        ! compute frequency in star frame
        nu_star(iphot) =  clight / (lambda0*1e-8)  ! Hz
        ! compute frequency in external frame 
        scalar = k(1)*star_vel(1,ilow) + k(2)*star_vel(2,ilow) + k(3)*star_vel(3,ilow)
        nu_em(iphot)  = nu_star(iphot) / (1d0 - scalar/clight)
        ! store velocity of the star, for peeling-off only
        v_em(:,iphot) = star_vel(:,ilow)
     end do

  end select
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------


  ! --------------------------------------------------------------------------------------
  ! write ICs
  ! --------------------------------------------------------------------------------------
  if (verbose) then
     write(*,*) 'Writing file: ',trim(outputfile)
     write(*,*) ' '
  endif
  open(unit=14, file=trim(outputfile), status='unknown', form='unformatted', action='write')
  write(14) nphotons    ! nb of MC photons = nb of photon packets
  write(14) total_flux  ! nb of real photons (per sec).
  write(14) ranseed
  write(14) (i,i=1,nphotons) ! ID
  write(14) (nu_em(i),i=1,nphotons)
  write(14) (x_em(:,i),i=1,nphotons)
  write(14) (k_em(:,i),i=1,nphotons)
  write(14) (-i,i=1,nphotons) ! seeds
  write(14) (nu_star(i),i=1,nphotons)
  write(14) (v_em(:,i),i=1,nphotons)
  close(14)
  ! --------------------------------------------------------------------------------------


  ! --------------------------------------------------------------------------------------
  ! deallocations 
  ! --------------------------------------------------------------------------------------
  deallocate(star_pos,star_vel,star_mass,star_age,star_met)
  deallocate(sweight, Ndot)
  deallocate(nu_star, nu_em, x_em, k_em, v_em)
  ! --------------------------------------------------------------------------------------

  
contains

  subroutine compute_cum_low_prob(n, weight, low_prob, cumtot)
    
    implicit none
    integer(kind=4),intent(in)                        :: n
    real(kind=8),dimension(n),intent(in)              :: weight
    real(kind=8),dimension(:),allocatable,intent(out) :: low_prob
    real(kind=8),intent(out)                          :: cumtot
    integer(kind=4)                                   :: i
    
    ! compute the total number of photons emitted per second by the sources
    cumtot = 0.0d0
    do i=1,n
       cumtot = cumtot + weight(i)
    end do
    
    ! calcul pour chaque particule la luminosite inferieure de son bin dans la distribution cumulative.. 
    allocate(low_prob(n+1))
    low_prob(1) = 0.0d0
    do i = 2,n
       low_prob(i) = low_prob(i-1) + weight(i-1)
    end do
    low_prob = low_prob / (low_prob(n)+weight(n))
    low_prob(n+1) = 1.1d0  ! higher than upper limit 
    
    return
  end subroutine compute_cum_low_prob
  
  
  subroutine binary_search(iseed, nbin, low_prob, ilow)
    
    implicit none
    integer(kind=4),intent(inout)             :: iseed
    integer(kind=4),intent(in)                :: nbin
    real(kind=8),intent(in),dimension(nbin+1) :: low_prob
    integer(kind=4),intent(out)               :: ilow
    integer(kind=4)                           :: iup, imid
    real(kind=8)                              :: mid, r1
    
    r1 = ran3(iseed)
    ! binary search
    iup = nbin+1
    ilow = 1
    do while (iup - ilow > 1)
       imid = (iup+ilow)/2
       mid  = low_prob(imid)
       if (r1 >= mid) then 
          ilow = imid
       else
          iup = imid
       end if
    end do
    ! check
    if (.not. (r1 >= low_prob(ilow) .and. r1 < low_prob(iup) )) then
       print*,'ERROR'
       stop
    end if
    return
  end subroutine binary_search
  
  
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
          case ('spec_type')
             write(spec_type,'(a)') trim(value)
          case ('spec_SSPdir')
             write(spec_SSPdir,'(a)') trim(value)
          case ('spec_mono_lambda0')
             read(value,*) spec_mono_lambda0
          case ('spec_gauss_lambda0')
             read(value,*) spec_gauss_lambda0
          case ('spec_gauss_sigma_kms')
             read(value,*) spec_gauss_sigma_kms
          case ('spec_powlaw_lmin_Ang')
             read(value,*) spec_powlaw_lmin_Ang
          case ('spec_powlaw_lmax_Ang')
             read(value,*) spec_powlaw_lmax_Ang
          case ('spec_table_lmin_Ang')
             read(value,*) spec_table_lmin_Ang
          case ('spec_table_lmax_Ang')
             read(value,*) spec_table_lmax_Ang
          case ('nPhotonPackets')
             read(value,*) nphotons
          case ('ranseed')
             read(value,*) ranseed
          case ('verbose')
             read(value,*) verbose
          case ('recompute_particle_initial_mass')
             read(value,*) recompute_particle_initial_mass
          case ('tdelay_SN')
             read(value,*) tdelay_SN
          case ('recyc_frac')
             read(value,*) recyc_frac
          case default
             write(*,'(a,a,a)') '> WARNING: parameter ',trim(name),' unknown '
          end select
       end do
    end if
    close(10)
    
    call read_ramses_params(pfile)
    
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
       write(unit,'(a,3(ES10.3,1x))') '  star_dom_pos       = ',star_dom_pos(1),star_dom_pos(2),star_dom_pos(3)
       select case (trim(star_dom_type))
       case ('sphere')
          write(unit,'(a,ES10.3)')       '  star_dom_rsp       = ',star_dom_rsp
       case ('shell')
          write(unit,'(a,ES10.3)')       '  star_dom_rin       = ',star_dom_rin
          write(unit,'(a,ES10.3)')       '  star_dom_rout      = ',star_dom_rout
       case('cube')
          write(unit,'(a,ES10.3)')       '  star_dom_size      = ',star_dom_size
       case('slab')
          write(unit,'(a,ES10.3)')       '  star_dom_thickness = ',star_dom_thickness
       end select
       write(unit,'(a)')             '# Spectral shape '
       write(unit,'(a,a)')           '  spec_type               = ',trim(spec_type)
       write(unit,'(a,a)')           '  spec_SSPdir             = ',trim(spec_SSPdir)
       select case(trim(spec_type))
       case('Mono')
          write(unit,'(a,ES10.3,a)')     '  spec_mono_lambda0   = ',spec_mono_lambda0, ' ! [A]'
       case('Gauss')
          write(unit,'(a,ES10.3,a)')     '  spec_gauss_lambda0    = ',spec_gauss_lambda0,   ' ! [A]'
          write(unit,'(a,ES10.3,a)')     '  spec_gauss_sigma_kms  = ',spec_gauss_sigma_kms, ' ! [km/s]'
       case('PowLaw')
          write(unit,'(a,es10.3,a)')     '  spec_powlaw_lmin_Ang    = ',spec_powlaw_lmin_Ang, ' ! [A]' 
          write(unit,'(a,es10.3,a)')     '  spec_powlaw_lmax_Ang    = ',spec_powlaw_lmax_Ang, ' ! [A]'
       case('Table')
          write(unit,'(a,es10.3,a)')     '  spec_table_lmin_Ang    = ',spec_table_lmin_Ang, ' ! [A]' 
          write(unit,'(a,es10.3,a)')     '  spec_table_lmax_Ang    = ',spec_table_lmax_Ang, ' ! [A]'
       end select
       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,i8)')          '  nPhotonPackets  = ',nphotons
       write(unit,'(a,i8)')          '  ranseed         = ',ranseed
       write(unit,'(a,L1)')          '  verbose         = ',verbose
       write(unit,'(a,L1)')          '  recompute_particle_initial_mass = ',recompute_particle_initial_mass
       if(recompute_particle_initial_mass)then
          write(unit,'(a,L1)')          '  tdelay_SN       = ',tdelay_SN
          write(unit,'(a,L1)')          '  recyc_frac      = ',recyc_frac
       endif
       write(unit,'(a)')             ' '
       call print_ramses_params(unit)
    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ' '
       write(*,'(a,a,a)')         '[PhotonsFromStars]'
       write(*,'(a)')             '# input / output parameters'
       write(*,'(a,a)')           '  outputfile      = ',trim(outputfile)
       write(*,'(a,a)')           '  repository      = ',trim(repository)
       write(*,'(a,i5)')          '  snapnum         = ',snapnum
       write(*,'(a)')             '# computational domain parameters'
       write(*,'(a,a)')           '  star_dom_type      = ',trim(star_dom_type)
       write(*,'(a,3(ES10.3,1x))') '  star_dom_pos       = ',star_dom_pos(1),star_dom_pos(2),star_dom_pos(3)
       select case (trim(star_dom_type))
       case ('sphere')
          write(*,'(a,ES10.3)')       '  star_dom_rsp       = ',star_dom_rsp
       case ('shell')
          write(*,'(a,ES10.3)')       '  star_dom_rin       = ',star_dom_rin
          write(*,'(a,ES10.3)')       '  star_dom_rout      = ',star_dom_rout
       case('cube')
          write(*,'(a,ES10.3)')       '  star_dom_size      = ',star_dom_size
       case('slab')
          write(*,'(a,ES10.3)')       '  star_dom_thickness = ',star_dom_thickness
       end select
       write(*,'(a)')             '# Spectral shape '
       write(*,'(a,a)')           '  spec_type               = ',trim(spec_type)
       write(*,'(a,a)')           '  spec_SSPdir             = ',trim(spec_SSPdir)
       select case(trim(spec_type))
       case('Mono')
          write(*,'(a,ES10.3,a)')     '  spec_mono_lambda0       = ',spec_mono_lambda0, ' ! [A]'
       case('Gauss')
          write(*,'(a,ES10.3,a)')     '  spec_gauss_lambda0      = ',spec_gauss_lambda0, ' ! [A]'
          write(*,'(a,es10.3,a)')     '  spec_gauss_sigma_kms    = ',spec_gauss_sigma_kms, ' ! [km/s]'
       case('PowLaw')
          write(*,'(a,es10.3,a)')     '  spec_powlaw_lmin_Ang    = ',spec_powlaw_lmin_Ang, ' ! [A]' 
          write(*,'(a,es10.3,a)')     '  spec_powlaw_lmax_Ang    = ',spec_powlaw_lmax_Ang, ' ! [A]'
       case('Table')
          write(*,'(a,es10.3,a)')     '  spec_table_lmin_Ang     = ',spec_table_lmin_Ang, ' ! [A]' 
          write(*,'(a,es10.3,a)')     '  spec_table_lmax_Ang     = ',spec_table_lmax_Ang, ' ! [A]'
       end select
       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,i8)')          '  nPhotonPackets  = ',nphotons
       write(*,'(a,i8)')          '  ranseed         = ',ranseed
       write(*,'(a,L1)')          '  verbose         = ',verbose
       write(*,'(a,L1)')          '  recompute_particle_initial_mass = ',recompute_particle_initial_mass
       if(recompute_particle_initial_mass)then
          write(*,'(a,L1)')          '  tdelay_SN       = ',tdelay_SN
          write(*,'(a,L1)')          '  recyc_frac      = ',recyc_frac
       endif
       write(*,'(a)')             ' '
       call print_ramses_params
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ' '
    end if

    return

  end subroutine print_PhotonsFromStars_params

  
end program PhotonsFromStars


  
