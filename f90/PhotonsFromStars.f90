program PhotonsFromStars
  
  ! generate photons emitted by star particles within a given domain

  use module_utils
  use module_domain
  use module_random
  use module_constants
  use module_ramses

  use module_ssp_lib
  
  implicit none
  
  type(domain)    :: emission_domain
  character(2000) :: parameter_file
  real(kind=8),allocatable :: star_pos(:,:),star_age(:),star_mass(:),star_vel(:,:),star_met(:)
  integer(kind=4) :: iran,i,nstars,narg
  !!integer(kind=8) :: ilast,j
  real(kind=8)    :: scalar,r1,r2, r3
  ! for analysis purposes (a posteriori weighting) we want to save the emitter-frame
  ! frequency (here the freq. in the emitting stellar particle's frame)
  real(kind=8),allocatable    :: nu_star(:)
  ! SED-related variables
  !!!integer(kind=4)             :: sed_nage,sed_nmet,imet,iage
  !!integer(kind=8)             :: nflux,n,ix
  !!!real(kind=8),allocatable    :: sed_age(:),sed_met(:),sweight(:),sed_nphot(:,:),sed_F_0(:,:),sed_beta(:,:)
  !integer(kind=4),allocatable :: star_iage(:),star_imet(:)
  !real(kind=8),allocatable    :: star_beta(:) 
  real(kind=8)                :: total_flux   !!!,minflux,check_flux,f0 !!,beta,betaplus2,lambda,x,dx1,dx2,dx
  !real(kind=8) :: dxage1,dxmet1,w,dxage2,dxmet2

  type(SSPgrid) :: NdotGrid
  real(kind=8),allocatable :: low_prob(:), low_prob2(:), x_em(:,:), k_em(:,:), nu_em(:), Ndot(:), NdotStar(:,:), sweight(:)
  integer(kind=4) :: ilow, iup, imid, iphot, iseed, ilow2, iup2, imid2
  real(kind=8) :: lambda0, mid, k(3), lambdamin, lambdamax, mid2, nu, spec_gauss_nu0

  
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
  character(30)             :: spec_type = 'Mono'           ! May be 'Mono', 'Gauss', 'PowLaw' ...   
  character(1000)           :: spec_SSPdir = '/Users/leo/ASTROPHYSICS/ASTRO/SSPlibs/'
  ! parameters for spec_type == 'Mono'
  real(kind=8)              :: spec_mono_lambda0 = 1216.                ! emission frequency [Hz] -> computed from weight_l0_Ang
  ! parameters for spec_type == 'Gauss'
  real(kind=8)              :: spec_gauss_lambda0 = 1216.                ! central frequency [Hz] -> computed from weight_l0_Ang
  real(kind=8)              :: spec_gauss_sigma_kms = 10.0   ! line width in velocity [km/s] -> read from file. 
  ! parameters for spec_type == 'PowLaw' : a power-law fit to continuum of each star particle, vs. its age and met.
  real(kind=8)              :: spec_powlaw_lmin_Ang = 1120.  ! min wavelength to sample (should be in the range where fit was made ...)
  real(kind=8)              :: spec_powlaw_lmax_Ang = 1320.  ! max ...
  ! parameters for spec_type == 'Table'
  real(kind=8)              :: spec_table_lmin_Ang = 1120.  ! min wavelength to sample
  real(kind=8)              :: spec_table_lmax_Ang = 1320.  ! max ...

  ! --- parameters for star particles/feedback in simulation
  real(kind=8)              :: tdelay_SN = 10.        ![Myr] SNs go off at tdelay_SN ... 
  real(kind=8)              :: recyc_frac = 0.8       ! correct for recycling ... we want the mass of stars formed ...
  
  ! --- miscelaneous
  integer(kind=4)           :: nphot   = 1000000      ! number of photons to generate
  integer(kind=4)           :: ranseed = -100         ! seed for random generator
  logical                   :: verbose = .true.
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
  call ramses_read_stars_in_domain(repository,snapnum,emission_domain,star_pos,star_age,star_mass,star_vel,star_met)
  ! --------------------------------------------------------------------------------------
  

  print*,size(star_mass)
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
        !if (sed_age(iage) < tdelay_SN) then ! SNs go off at 10Myr ... 
        !   sweight(i) = sweight(i)/recyc_frac  !! correct for recycling ... we want the mass of stars formed ...
        !end if
     end do
     ! compute the total number of photons emitted per second by the sources
     total_flux = 0.0d0
     do i=1,nstars
        total_flux = total_flux + sweight(i)
     end do
     if (verbose) write(*,*) '> Total luminosity (nb of photons per second): ',total_flux
  
     ! calcul pour chaque particule la luminosite inferieure de son bin dans la distribution cumulative.. 
     allocate(low_prob(nstars+1))
     low_prob(1) = 0.0d0
     do i = 2,nstars
        low_prob(i) = low_prob(i-1) + sweight(i-1)
     end do
     low_prob = low_prob / low_prob(nstars)
     low_prob(nstars+1) = 1.1  ! higher than upper limit 

     ! for each photon packet, draw the emitting star
     allocate(x_em(1:3,1:nphot), k_em(1:3,1:nphot), nu_em(1:nphot), nu_star(1:nphot))

     iseed = ranseed
     do iphot = 1,nphot
        r1 = ran3(iseed)
        ! binary search
        iup = nstars
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
           print*,'hi Harley ;) '
        end if
        ! draw photon's ICs from star ilow
        ! give photon the position of the star
        !print*,iphot,ilow
        x_em(:,iphot) = star_pos(:,ilow) 
        ! draw propagation direction
        call isotropic_direction(k,iseed)
        k_em(:,iphot) = k
        ! compute frequency in star frame
        r2 = ran3(iran)
        r3 = ran3(iran)
        nu = sqrt(-2.*log(r3)) * cos(2.0d0*pi*r3)
        nu = (spec_gauss_sigma_kms * 1d5 * spec_gauss_nu0 / clight) * nu + spec_gauss_nu0
        nu_star(iphot) =  nu    ! clight / (lambda0*1e-8)  ! Hz
        ! compute frequency in external frame 
        scalar = k(1)*star_vel(1,ilow) + k(2)*star_vel(2,ilow) + k(3)*star_vel(3,ilow)
        nu_em(iphot)  = nu_star(iphot) / (1d0 - scalar/clight)
     end do

     
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
  case('Table')
     !print*,'Not implemented yet...'
     !stop

     lambdamin = spec_table_lmin_Ang
     lambdamax = spec_table_lmax_Ang
     call ssp_lib_extract_subset(lambdamin,lambdamax,NdotGrid) ! extract the SSP age-met grid of Nphot in lambda range [lambdamin;lambdamax]
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
        NdotStar(i,:) = Ndot(:)
        ! integrate nphotPerSecPerMsun(nlambda)
        sweight(i) = trapz1(NdotGrid%lambda, Ndot, NdotGrid%nlambda)
     enddo
     ! compute the total number of photons emitted per second by the sources
     total_flux = 0.0d0
     do i=1,nstars
        total_flux = total_flux + sweight(i)
     end do
     if (verbose) write(*,*) '> Total luminosity (nb of photons per second): ',total_flux
     
     ! calcul pour chaque particule la luminosite inferieure de son bin dans la distribution cumulative.. 
     allocate(low_prob(nstars+1))
     low_prob(1) = 0.0d0
     do i = 2,nstars
        low_prob(i) = low_prob(i-1) + sweight(i-1)
     end do
     low_prob = low_prob / low_prob(nstars)
     low_prob(nstars+1) = 1.1  ! higher than upper limit 

     ! for each photon packet, draw the emitting star
     allocate(x_em(1:3,1:nphot), k_em(1:3,1:nphot), nu_em(1:nphot), nu_star(1:nphot))
     allocate(low_prob2(NdotGrid%nlambda+1))
     iseed = ranseed
     do iphot = 1,nphot
        r1 = ran3(iseed)
        ! binary search
        iup = nstars
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
           print*,'hi Harley ;) '
        end if
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
        do i = 2,NdotGrid%nlambda
           low_prob2(i) = low_prob2(i-1) + NdotStar(ilow,i-1)
        end do
        low_prob2 = low_prob2 / low_prob2(NdotGrid%nlambda)
        low_prob2(NdotGrid%nlambda+1) = 1.1  ! higher than upper limit 
        r2 = ran3(iseed)
        ! binary search
        iup2 = NdotGrid%nlambda
        ilow2 = 1
        do while (iup2 - ilow2 > 1)
           imid2 = (iup2+ilow2)/2
           mid2  = low_prob2(imid2)
           if (r2 >= mid2) then 
              ilow2 = imid2
           else
              iup2 = imid2
           end if
        end do
        ! check
        if (.not. (r2 >= low_prob2(ilow2) .and. r2 < low_prob2(iup2) )) then
           print*,'hi Harley ;) '
        end if
        ! photon emitted in the bin ilow2;ilow2+1
        ! get lamba_em
        ! TO DO
        nu_star(iphot) = clight / (NdotGrid%lambda(ilow2)*1e-8) ! Hz
        ! compute frequency in external frame 
        scalar = k(1)*star_vel(1,ilow) + k(2)*star_vel(2,ilow) + k(3)*star_vel(3,ilow)
        nu_em(iphot)  = nu_star(iphot) / (1d0 - scalar/clight)
     enddo
     

! --------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------
  case('Mono')
     
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
        call ssp_lib_interpolate(NdotGrid, star_age(i)/1.e3, log10(star_met(i)), Ndot)    ! Ndot number of photons / s / A / Msun
        sweight(i) = Ndot(1) * star_mass(i) / msun  ! M_sun
        !if (sed_age(iage) < tdelay_SN) then ! SNs go off at 10Myr ... 
        !   sweight(i) = sweight(i)/recyc_frac  !! correct for recycling ... we want the mass of stars formed ...
        !end if
     end do
     ! compute the total number of photons emitted per second by the sources
     total_flux = 0.0d0
     do i=1,nstars
        total_flux = total_flux + sweight(i)
     end do
     if (verbose) write(*,*) '> Total luminosity (nb of photons per second): ',total_flux
     
     ! calcul pour chaque particule la luminosite inferieure de son bin dans la distribution cumulative.. 
     allocate(low_prob(nstars+1))
     low_prob(1) = 0.0d0
     do i = 2,nstars
        low_prob(i) = low_prob(i-1) + sweight(i-1)
     end do
     low_prob = low_prob / low_prob(nstars)
     low_prob(nstars+1) = 1.1  ! higher than upper limit 

     ! for each photon packet, draw the emitting star
     allocate(x_em(1:3,1:nphot), k_em(1:3,1:nphot), nu_em(1:nphot), nu_star(1:nphot))

     iseed = ranseed
     do iphot = 1,nphot
        r1 = ran3(iseed)
        ! binary search
        iup = nstars
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
           print*,'hi Harley ;) '
        end if
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
     end do

  end select
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------


  ! --------------------------------------------------------------------------------------
  ! write ICs
  ! --------------------------------------------------------------------------------------
  if (verbose) write(*,*) '> writing file'
  open(unit=14, file=trim(outputfile), status='unknown', form='unformatted', action='write')
  write(14) nphot      ! nb of MC photons 
  write(14) total_flux ! nb of real photons (per sec).
  write(14) ranseed
  write(14) (i,i=1,nphot) ! ID
  write(14) (nu_em(i),i=1,nphot)
  write(14) (x_em(:,i),i=1,nphot)
  write(14) (k_em(:,i),i=1,nphot)
  write(14) (-i,i=1,nphot) ! seeds
  write(14) (nu_star(i),i=1,nphot)
  close(14)
  ! --------------------------------------------------------------------------------------


  ! --------------------------------------------------------------------------------------
  ! deallocations 
  ! --------------------------------------------------------------------------------------
  deallocate(star_pos,star_vel,star_mass,star_age,star_met)
  !select case(trim(weight_type))
  !case('Mono')
  !   deallocate(sed_nphot)
  !case('PowLaw')
  !   deallocate(sed_F_0,sed_beta,star_beta)
  !case ('Table')
  !   deallocate(sed_nphot,spec_table_lofx,star_iage,star_imet)
  !end select
  deallocate(sweight)
  deallocate(nu_star, nu_em, x_em, k_em)
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
    logical         :: section_present,ok
    
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
          case default
             write(*,'(a,a,a)') '> WARNING: parameter ',trim(name),' unknown '
          end select
       end do
    end if
    close(10)

    ! test for compatibility of parameters
    !ok = (weight_type=='Mono' .and. spec_type=='Mono') .or. (weight_type=='Mono' .and. spec_type=='Gauss')
    !ok = ok .or. (weight_type=='PowLaw' .and. spec_type=='PowLaw') .or. (weight_type=='Table' .and. spec_type=='Table')
    !if (.not. ok) then
    !   write(*,'(a,a,a,a)') '> ERROR: incompatible options : weight_type==',trim(weight_type),' and spec_type==',trim(spec_type)
    !   stop
    !end if

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
       select case(trim(spec_type))
       case('Gauss')
          write(unit,'(a,ES10.3,a)')     '  spec_gauss_sigma_kms = ',spec_gauss_sigma_kms, ' ! [km/s]'
       case('PowLaw')
          write(unit,'(a,es10.3,a)')     '  spec_powlaw_lmin_Ang    = ',spec_powlaw_lmin_Ang, ' ! [A]' 
          write(unit,'(a,es10.3,a)')     '  spec_powlaw_lmax_Ang    = ',spec_powlaw_lmax_Ang, ' ! [A]'
       end select
       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,i8)')          '  nphot           = ',nphot
       write(unit,'(a,i8)')          '  ranseed         = ',ranseed
       write(unit,'(a,L1)')          '  verbose         = ',verbose
       write(unit,'(a)')             ' '
       call print_ramses_params(unit)
    else
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
       select case(trim(spec_type))
       case('Gauss')
          write(*,'(a,es10.3,a)')     '  spec_gauss_sigma_kms = ',spec_gauss_sigma_kms, ' ! [km/s]'
       case('PowLaw')
          write(*,'(a,es10.3,a)')     '  spec_powlaw_lmin_Ang    = ',spec_powlaw_lmin_Ang, ' ! [A]' 
          write(*,'(a,es10.3,a)')     '  spec_powlaw_lmax_Ang    = ',spec_powlaw_lmax_Ang, ' ! [A]'
       end select
       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,i8)')          '  nphot           = ',nphot
       write(*,'(a,i8)')          '  ranseed         = ',ranseed
       write(*,'(a,L1)')          '  verbose         = ',verbose
       write(*,'(a)')             ' '
       call print_ramses_params
    end if

    return

  end subroutine print_PhotonsFromStars_params

  
end program PhotonsFromStars


  
