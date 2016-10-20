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
  real(kind=8),allocatable :: star_pos(:,:),star_age(:),star_mass(:),star_vel(:,:)
  real(kind=8),allocatable :: star_pos2(:,:),star_vel2(:,:)
  integer(kind=4) :: i,nstars,nyoung,ilast,j,iran
  real(kind=8) :: minmass,scalar,nu
  type(photon_init),dimension(:),allocatable :: photgrid

  
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [PhotonsFromStars] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: outputfile = 'PhotICs.dat' ! file to which outputs will be written
  character(2000)           :: repository = './'      ! ramses run directory (where all output_xxxxx dirs are).
  integer(kind=4)           :: snapnum = 1            ! ramses output number to use
  ! --- domain whithin which star particles will be selected (should be within computational domain used for RT). 
  character(10)             :: star_dom_type      = 'sphere'         ! shape type of domain  // default is sphere.
  real(kind=8),dimension(3) :: star_dom_pos       = (/0.5,0.5,0.5/)  ! center of domain [code units]
  real(kind=8)              :: star_dom_rsp       = 0.3              ! radius of spher [code units]
  real(kind=8)              :: star_dom_size      = 0.3              ! size of cube [code units]
  real(kind=8)              :: star_dom_rin       = 0.0              ! inner radius of shell [code units]
  real(kind=8)              :: star_dom_rout      = 0.3              ! outer radius of shell [code units]
  real(kind=8)              :: star_dom_thickness = 0.1              ! thickness of slab [code units]
  ! --- which stars shine
  integer(kind=4)           :: nphot   = 1000000      ! number of photons to generate
  real(kind=8)              :: max_age = 10.0d0       ! stars older than this don't shine [Myr]
  ! --- how stars shine
  character(30)             :: spec_type = 'monochromatic' ! how to draw frequencies
  ! ------ spec_type == 'monochromatic' : all photons at the same frequency nu_0
  real(kind=8)              :: nu_0    = clight/1215.6701d-8 ! emission frequency [Hz] 
  ! ------ spec_type == 'flat_fnu' : photons have a flat distribution in nu, between nu_min and nu_max
  real(kind=8)              :: nu_min  = clight/1221d-8    ! min frequency [Hz]
  real(kind=8)              :: nu_max  = clight/1210d-8    ! max frequency [Hz]
  ! --- miscelaneous
  integer(kind=4)           :: ranseed = -100         ! seed for random generator
  logical                   :: verbose = .true.
  logical                   :: cosmo = .true.         ! cosmo flag
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
  call ramses_read_stars_in_domain(repository,snapnum,emission_domain,star_pos,star_age,star_mass,star_vel,cosmo)
  ! --------------------------------------------------------------------------------------

  
  ! --------------------------------------------------------------------------------------
  ! keep only stars younger than max_age and oversample according to mass of particles
  ! --------------------------------------------------------------------------------------
  if (verbose) write(*,*) '> selecting young star particles '
  nstars = size(star_age)
  nyoung = 0
  minmass = minval(star_mass)
  ! count stars in unit of min-mass-star-particles so that stars will emit proportionally to their mass ... 
  do i = 1,nstars
     if (star_age(i) <= max_age) nyoung = nyoung + nint(star_mass(i)/minmass)
  end do
  if (nyoung == 0) then
     print*,'No young stars ... '
     stop
  end if
  ! copy arrays
  allocate(star_pos2(3,nyoung),star_vel2(3,nyoung))
  ilast = 1
  do i = 1,nstars
     if (star_age(i) <= max_age) then
        do j = ilast, ilast + nint(star_mass(i)/minmass)-1
           star_pos2(:,j) = star_pos(:,i)
           star_vel2(:,j) = star_pos(:,i)
        end do
        ilast = ilast + nint(star_mass(i)/minmass)
     end if
  end do
  ! --------------------------------------------------------------------------------------

  
  ! --------------------------------------------------------------------------------------
  ! make particles shine
  ! --------------------------------------------------------------------------------------
  if (verbose) write(*,*) '> generating photons'
  allocate(photgrid(nphot))
  iran = ranseed
  do i = 1,nphot
     ! pick a star particle
     j = int(ran3(iran)*nyoung) + 1
     if (j > nyoung) j = nyoung
     ! define photon accordingly
     photgrid(i)%ID    = i
     photgrid(i)%x_em  = star_pos2(:,j)
     photgrid(i)%iran  = iran
     call isotropic_direction(photgrid(i)%k_em,iran)
     ! define star-particle-frame emission frequency
     select case(trim(spec_type))
     case('monochromatic')
        nu = nu_0
     case('flat_fnu')
        nu = ran3(iran) * (nu_max-nu_min)
     case default
        print*,'ERROR: unknown spec_type :',trim(spec_type)
     end select
     ! knowing the direction of emission and the velocity of the source (star particle), we
     ! compute the external-frame frequency :
     scalar = photgrid(i)%k_em(1)*star_vel2(1,j) + photgrid(i)%k_em(2)*star_vel2(2,j) + photgrid(i)%k_em(3)*star_vel2(3,j)
     photgrid(i)%nu_em = nu / (1d0 - scalar/clight)
  end do
  ! --------------------------------------------------------------------------------------

  ! --------------------------------------------------------------------------------------
  ! write ICs
  ! --------------------------------------------------------------------------------------
  if (verbose) write(*,*) '> writing file'
  open(unit=14, file=trim(outputfile), status='unknown', form='unformatted', action='write')
  write(14) nphot
  write(14) ranseed
  write(14) (photgrid(i)%ID,i=1,nphot)
  write(14) (photgrid(i)%nu_em,i=1,nphot)
  write(14) (photgrid(i)%x_em(:),i=1,nphot)
  write(14) (photgrid(i)%k_em(:),i=1,nphot)
  write(14) (photgrid(i)%iran,i=1,nphot)
  close(14)
  ! --------------------------------------------------------------------------------------

  deallocate(star_pos2,star_vel2,star_pos,star_vel,star_mass,star_age,photgrid)
  
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
          case ('star_dom_pos')
             read(value,*) star_dom_pos(1),star_dom_pos(2),star_dom_pos(3)
          case ('star_dom_rsp')
             read(value,*) star_dom_rsp
          case ('star_dom_type')
             write(star_dom_type,'(a)') trim(value)
          case ('verbose')
             read(value,*) verbose
          case ('ranseed')
             read(value,*) ranseed
          case ('cosmo')
             read(value,*) cosmo
          case ('outputfile')
             write(outputfile,'(a)') trim(value)
          case ('repository')
             write(repository,'(a)') trim(value)
          case ('snapnum')
             read(value,*) snapnum
          case ('max_age')
             read(value,*) max_age
          case ('spec_type')
             write(spec_type,'(a)') trim(value)
          case ('nu_0')
             read(value,*) nu_0
          case ('nu_min')
             read(value,*) nu_min
          case ('nu_max')
             read(value,*) nu_max
          case ('nphot')
             read(value,*) nphot
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
       write(unit,'(a,a)')           '  star_dom_type   = ',trim(star_dom_type)
       write(unit,'(a,3(ES9.3,1x))') '  star_dom_pos    = ',star_dom_pos(1),star_dom_pos(2),star_dom_pos(3)
       write(unit,'(a,ES9.3)')       '  star_dom_rsp    = ',star_dom_rsp
       write(unit,'(a)')             '# how stars shine'
       write(unit,'(a,i8)')          '  nphot           = ',nphot
       write(unit,'(a,es9.3,a)')     '  max_age         = ',max_age, ' ! [Myr]' 
       write(unit,'(a,a)')           '  spec_type       = ',trim(spec_type)
       select case(trim(spec_type))
       case('monochromatic')
          write(unit,'(a,es9.3,a)')     '  nu_0            = ',nu_0, ' ! [Hz]'
       case('flat_fnu')
          write(unit,'(a,es9.3,a)')     '  nu_min          = ',nu_min, ' ! [Hz]'
          write(unit,'(a,es9.3,a)')     '  nu_max          = ',nu_max, ' ! [Hz]'
       case default
          print*,'ERROR: unknown spec_type :',trim(spec_type)
       end select
       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,i8)')          '  ranseed         = ',ranseed
       write(unit,'(a,L1)')          '  verbose         = ',verbose
       write(unit,'(a,L1)')          '  cosmo           = ',cosmo
       write(unit,'(a)')             ' '
    else
       write(*,'(a,a,a)')         '[PhotonsFromStars]'
       write(*,'(a)')             '# input / output parameters'
       write(*,'(a,a)')           '  outputfile    = ',trim(outputfile)
       write(*,'(a,a)')           '  repository    = ',trim(repository)
       write(*,'(a,i5)')          '  snapnum       = ',snapnum
       write(*,'(a)')             '# computational domain parameters'
       write(*,'(a,a)')           '  star_dom_type = ',trim(star_dom_type)
       write(*,'(a,3(ES9.3,1x))') '  star_dom_pos  = ',star_dom_pos(1),star_dom_pos(2),star_dom_pos(3)
       write(*,'(a,ES9.3)')       '  star_dom_rsp  = ',star_dom_rsp
       write(*,'(a)')             '# how stars shine'
       write(*,'(a,i8)')          '  nphot         = ',nphot
       write(*,'(a,es9.3,a)')     '  max_age       = ',max_age, ' ! [Myr]'
       write(*,'(a,a)')           '  spec_type     = ',trim(spec_type)
       select case(trim(spec_type))
       case('monochromatic')
          write(*,'(a,es9.3,a)')     '  nu_0            = ',nu_0, ' ! [Hz]'
       case('flat_fnu')
          write(*,'(a,es9.3,a)')     '  nu_min          = ',nu_min, ' ! [Hz]'
          write(*,'(a,es9.3,a)')     '  nu_max          = ',nu_max, ' ! [Hz]'
       case default
          print*,'ERROR: unknown spec_type :',trim(spec_type)
       end select
       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,i8)')          '  ranseed     = ',ranseed
       write(*,'(a,L1)')          '  verbose    = ',verbose
       write(*,'(a,L1)')          '  cosmo      = ',cosmo
       write(*,'(a)')             ' '       
    end if

    return

  end subroutine print_PhotonsFromStars_params


  
end program PhotonsFromStars


  
