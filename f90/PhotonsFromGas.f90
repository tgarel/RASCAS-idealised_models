program PhotonsFromGas
  
  ! generate photons emitted by star particles within a given domain

  use module_photon
  use module_utils
  use module_domain
  use module_random
  use module_constants
  use module_ramses

  implicit none
  
  type(domain)    :: emission_domain
  character(2000) :: parameter_file
  integer(kind=4) :: iran,i,narg
  integer(kind=8) :: ilast,j
  real(kind=8)    :: scalar,nu,r1,r2
  type(photon_init),dimension(:),allocatable :: photgrid
  ! for analysis purposes (a posteriori weighting) we want to save the emitter-frame
  ! frequency (here the freq. in the emitting stellar particle's frame)
  real(kind=8),allocatable    :: nu_gas(:)
  ! SED-related variables
  integer(kind=8)             :: nflux,n
  integer(kind=8),allocatable :: cum_flux_prob(:)
  real(kind=8)                :: total_flux,minflux,check_flux
  
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [PhotonsFromStars] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: outputfile = 'PhotICs.dat' ! file to which outputs will be written
  character(2000)           :: repository = './'          ! ramses run directory (where all output_xxxxx dirs are).
  integer(kind=4)           :: snapnum = 1                ! ramses output number to use

  ! --- domain whithin which gas cells will be selected (should be within computational domain used for RT). 
  character(10)             :: gas_dom_type      = 'sphere'         ! shape type of domain  // default is sphere.
  real(kind=8),dimension(3) :: gas_dom_pos       = (/0.5,0.5,0.5/)  ! center of domain [code units]
  real(kind=8)              :: gas_dom_rsp       = 0.3              ! radius of spher [code units]
  real(kind=8)              :: gas_dom_size      = 0.3              ! size of cube [code units]
  real(kind=8)              :: gas_dom_rin       = 0.0              ! inner radius of shell [code units]
  real(kind=8)              :: gas_dom_rout      = 0.3              ! outer radius of shell [code units]
  real(kind=8)              :: gas_dom_thickness = 0.1              ! thickness of slab [code units]
  
  ! --- define how star particles emit (i.e. the star-particle-frame spectral shape)
  ! Four options here :
  ! - spec_type=='Mono'   : we emit all photons at the same wavelength (in star's frame)
  ! - spec_type=='Gauss'  : we sample a Gaussian distribution ...
  ! - spec_type=='PowLaw' : we sample a power-law continuum between two wavelengths. 
  ! - spec_type=='Table'  : we sample a tabulated spectrum. 
  character(30)             :: spec_type = 'Gauss'           ! May be 'Mono', 'Gauss', 'PowLaw' ...   
  ! parameters for spec_type == 'Mono'
  real(kind=8)              :: spec_mono_nu0                 ! emission frequency [Hz] -> computed from weight_l0_Ang
  ! parameters for spec_type == 'Gauss'
  real(kind=8)              :: spec_gauss_nu0                ! central frequency [Hz] -> computed from weight_l0_Ang
  real(kind=8)              :: spec_gauss_sigma_kms = -1.0   ! line width in velocity [km/s] -> read from file. 
  
  real(kind=8)              :: weight_l0_Ang=1215.67 ! this is the lbda_0 above or the monochromatic

  ! --- miscelaneous
  integer(kind=4)           :: nphot   = 1000000      ! number of photons to generate
  integer(kind=4)           :: ranseed = -100         ! seed for random generator
  logical                   :: verbose = .true.
  logical                   :: cosmo   = .true.       ! cosmo flag
  ! --------------------------------------------------------------------------
  character(2000)           :: DataDir = './'                                   ! where input files below are
  character(2000)           :: DomDumpFile  = 'domain_decomposition_params.dat' ! the file describing the outputs
  character(2000)           :: line, file_compute_dom
  character(2000),allocatable,dimension(:) :: mesh_file_list, domain_file_list
  integer(kind=4)           :: ndomain,cell_level,ind
  integer(kind=8)           :: icell,ileaf,ioct
  real(kind=8),dimension(3) :: posoct,cell_corner,cell_pos,cell_vel
  real(kind=8)              :: cell_size,cell_vol,emiss,eLya_erg
  ! --------------------------------------------------------------------------
  type(mesh)                :: meshdom
  real(kind=8)              :: sigma_kms = -1.0   ! line width in velocity [km/s]  
  logical                   :: is_inside=.false.  !
  ! --------------------------------------------------------------------------

  eLya_erg=10.16 * 1.60218d-12

  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: PhotonsFromGas path/to/params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_PhotonsFromGas_params(parameter_file)
  if (verbose) call print_PhotonsFromGas_params
  ! ------------------------------------------------------------



  ! --------------------------------------------------------------------------------------
  ! define domain within which stars may shine
  ! --------------------------------------------------------------------------------------
  select case(gas_dom_type)
  case('sphere')
     call domain_constructor_from_scratch(emission_domain,gas_dom_type, &
          xc=gas_dom_pos(1),yc=gas_dom_pos(2),zc=gas_dom_pos(3),r=gas_dom_rsp)
  case('shell')
     call domain_constructor_from_scratch(emission_domain,gas_dom_type, &
          xc=gas_dom_pos(1),yc=gas_dom_pos(2),zc=gas_dom_pos(3),r_inbound=gas_dom_rin,r_outbound=gas_dom_rout)
  case('cube')
     call domain_constructor_from_scratch(emission_domain,gas_dom_type, & 
          xc=gas_dom_pos(1),yc=gas_dom_pos(2),zc=gas_dom_pos(3),size=gas_dom_size)
  case('slab')
     call domain_constructor_from_scratch(emission_domain,gas_dom_type, &
          xc=gas_dom_pos(1),yc=gas_dom_pos(2),zc=gas_dom_pos(3),thickness=gas_dom_thickness)
  end select
  ! --------------------------------------------------------------------------------------



  ! -------------------- Read domain list from CreateDomDump param file --------------------
  if (verbose) print *,'--> reading domain and mesh...'
  open(unit=18,file=DomDumpFile,status='old',form='formatted')
  read(18,'(a)') line ; i = scan(line,'=') ; file_compute_dom = trim(DataDir)//trim(adjustl(line(i+1:)))
  read(18,'(a)') line ; i = scan(line,'=') ; read(line(i+1:),*) ndomain
  allocate(mesh_file_list(ndomain),domain_file_list(ndomain))
  do j = 1, ndomain
     read(18,'(a)') line ; i = scan(line,'=') ; domain_file_list(j) = trim(DataDir)//trim(adjustl(line(i+1:)))
     read(18,'(a)') line ; i = scan(line,'=') ; mesh_file_list(j) = trim(DataDir)//trim(adjustl(line(i+1:)))
  end do
  close(18)
  ! ------------------------------------------------------------------------------------------
  if(ndomain>1)then
    write(*,*)'Domain decomposition is not ready yet for PhotonsFromGas'
    stop
  endif
  call mesh_from_file(mesh_file_list(1),meshdom)


  box_size_cm = ramses_get_box_size_cm(repository,snapnum)








  ! --------------------------------------------------------------------------------------
  ! check gas info within ALL domain; TODO domain selection
  ! --------------------------------------------------------------------------------------
  ! sum up the emissivity
  total_flux = 0.0d0
  minflux = 1d99
  do icell=1,meshdom%nCell
     ileaf = meshdom%son(icell)
     if(ileaf<0)then ! if leaf
        ind   = int((icell - meshdom%nCoarse - 1) / meshdom%nOct + 1,kind=4)
        ioct  = icell - meshdom%nCoarse - (ind - 1) * meshdom%nOct
        posoct(:) = meshdom%xoct(ioct,:)
        cell_level   = meshdom%octlevel(ioct)
        cell_corner  = get_cell_corner(posoct,ind,cell_level)
        cell_size    = 0.5d0**cell_level 

        cell_pos     = cell_corner + 0.5*cell_size
        cell_vol     = (cell_size*box_size_cm)**3.0

        ! check if the point is inside the domain we want
        is_inside = domain_contains_point(cell_pos,emission_domain)

        if(is_inside)then
           emiss = meshdom%gas(-ileaf)%emiss ! [erg/cm^3/s]
           emiss = emiss * cell_vol / eLya_erg ! [erg/s] -> [#/s]
           total_flux   = total_flux + emiss
        else
           emiss = 0.0
        endif

        if(emiss.lt.minflux.and.emiss.gt.0)minflux=emiss

!        if(ind>5)then
!        print *, 'cell_level =', cell_level
!        print *, 'cell_corner =', cell_corner 
!        print *, 'nHI =', meshdom%gas(-ileaf)%nHI
!        stop
!        endif
     endif
  end do 
 
  ! --------------------------------------------------------------------------------------
  ! define linear sampling of number of photons 
  ! --------------------------------------------------------------------------------------
  if (verbose) write(*,*) '> Total emissivity (nb of photons per second): ',total_flux

  ! it may happen that the range of luminosities is too large for linear sampling with our method ... 
  ! In that case we need to ignore faint particles:
  if (total_flux / minflux > 2d8) minflux = total_flux / 2d8  ! NB: dont go much higher than 1d8 (to stay below a few GB RAM). 

  ! check that we dont loose significant flux by sampling only particles with sweight > minflux
  check_flux = 0.0d0
  do icell=1,meshdom%nCell
     ileaf = meshdom%son(icell)
     if(ileaf<0)then ! if leaf

        ind         = int((icell - meshdom%nCoarse - 1) / meshdom%nOct + 1,kind=4)
        ioct        = icell - meshdom%nCoarse - (ind - 1) * meshdom%nOct
        posoct(:)   = meshdom%xoct(ioct,:)
        cell_level  = meshdom%octlevel(ioct)
        cell_corner = get_cell_corner(posoct,ind,cell_level)
        cell_size   = 0.5d0**cell_level 

        cell_pos    = cell_corner + 0.5*cell_size
        cell_vol    = (cell_size*box_size_cm)**3.0
        ! check if the point is inside the domain we want
        is_inside = domain_contains_point(cell_pos,emission_domain)

        if(is_inside)then
           emiss = meshdom%gas(-ileaf)%emiss ! [erg/cm^3/s]
           emiss = emiss * cell_vol / eLya_erg ! [erg/s] -> [#/s]
        else 
           emiss = 0.0
        endif

        if(emiss.ge.minflux)  check_flux=check_flux+emiss
     endif
  end do 
 
  if (verbose) write(*,*) '> We sample this fraction of total flux: ',check_flux / total_flux
  if ((total_flux - check_flux) / total_flux > 0.001) then
     print*,'> Flux losses > 0.1 percent... change algorithm ...',(total_flux - check_flux) / total_flux
     ! debug - stop
  end if
  ! construct the cumulative flux distribution, with enough bins to have the smallest star-particle flux in a bin. 
  allocate(cum_flux_prob(int(3*total_flux / minflux,kind=8)))
  ilast = 1
  cum_flux_prob = -99

  do icell=1,meshdom%nCell
     ileaf = meshdom%son(icell)
     if(ileaf<0)then ! if leaf

        ind         = int((icell - meshdom%nCoarse - 1) / meshdom%nOct + 1,kind=4)
        ioct        = icell - meshdom%nCoarse - (ind - 1) * meshdom%nOct
        posoct(:)   = meshdom%xoct(ioct,:)
        cell_level  = meshdom%octlevel(ioct)
        cell_corner = get_cell_corner(posoct,ind,cell_level)
        cell_size   = 0.5d0**cell_level 

        cell_pos    = cell_corner + 0.5*cell_size
        cell_vol    = (cell_size*box_size_cm)**3.0

        ! check if the point is inside the domain we want
        is_inside = domain_contains_point(cell_pos,emission_domain)

        if(is_inside)then
           emiss = meshdom%gas(-ileaf)%emiss ! [erg/cm^3/s]
           emiss = emiss * cell_vol / eLya_erg ! [erg/s] -> [#/s]
        else 
           emiss = 0.0
        endif

        if(emiss.ge.minflux)then
           n = int(3*emiss/minflux,kind=8)
           cum_flux_prob(ilast:ilast+n-1) = icell ! be careful
           ilast = ilast + n
        endif
     endif
  end do 

  nflux = ilast-1
  print*,'> nflux, size(cum_flux_prob):', nflux, size(cum_flux_prob)
  ! --------------------------------------------------------------------------------------

  
  ! --------------------------------------------------------------------------------------
  ! now we can draw integers from 1 to nflux and assign photons to stars according to spectral shape.
  ! --------------------------------------------------------------------------------------
  if (trim(spec_type) == 'Mono')  spec_mono_nu0  = clight/weight_l0_Ang*1d8
  if (trim(spec_type) == 'Gauss') spec_gauss_nu0 = clight/weight_l0_Ang*1d8
  allocate(photgrid(nphot),nu_gas(nphot))
  iran = -abs(ranseed)
  do i = 1,nphot
     j = int(ran3(iran)*nflux,kind=8)+1
     if (j > nflux) j = nflux
     icell = cum_flux_prob(j)
     photgrid(i)%ID    = i

     ileaf = meshdom%son(icell) ! has to be a leaf cell, so son(icell) < 0
     ind   = int((icell - meshdom%nCoarse - 1) / meshdom%nOct + 1,kind=4)
     ioct  = icell - meshdom%nCoarse - (ind - 1) * meshdom%nOct
     posoct(:) = meshdom%xoct(ioct,:)
     cell_level   = meshdom%octlevel(ioct)
     cell_corner  = get_cell_corner(posoct,ind,cell_level)
     cell_size    = 0.5d0**cell_level 
     cell_pos     = cell_corner + 0.5*cell_size
     photgrid(i)%x_em  = cell_pos
     photgrid(i)%iran  = -i 
     call isotropic_direction(photgrid(i)%k_em,iran)
     if(spec_gauss_sigma_kms>0) then
        sigma_kms = spec_gauss_sigma_kms
     else
        sigma_kms = meshdom%gas(-ileaf)%dopwidth / 1d5 
     endif
     select case(trim(spec_type))
     case('Gauss')
        r1 = ran3(iran)
        r2 = ran3(iran)
        nu = sqrt(-2.*log(r1)) * cos(2.0d0*pi*r2)
        nu = (sigma_kms * 1d5 * spec_gauss_nu0 / clight) * nu + spec_gauss_nu0
     case('Mono')
        nu = spec_mono_nu0
     end select

     nu_gas(i) = nu  ! gas-particle-frame frequency
     ! now put in external frame using particle's velocity.

     cell_vel(:) = meshdom%gas(-ileaf)%v(:) ! [cm/s] 
     scalar = photgrid(i)%k_em(1)*cell_vel(1) + photgrid(i)%k_em(2)*cell_vel(2) + photgrid(i)%k_em(3)*cell_vel(3)
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
  write(14) (nu_gas(i),i=1,nphot)
  close(14)
  ! --------------------------------------------------------------------------------------


  ! --------------------------------------------------------------------------------------
  ! deallocations 
  ! --------------------------------------------------------------------------------------
  deallocate(cum_flux_prob,photgrid,nu_gas)
  ! --------------------------------------------------------------------------------------

  
contains
  
  subroutine read_PhotonsFromGas_params(pfile)

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
       if (line(1:18) == '[PhotonsFromGas]') then
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
          case ('DataDir')
             write(DataDir,'(a)') trim(value)
          case ('snapnum')
             read(value,*) snapnum
          case ('gas_dom_type')
             write(gas_dom_type,'(a)') trim(value)
          case ('gas_dom_pos')
             read(value,*) gas_dom_pos(1),gas_dom_pos(2),gas_dom_pos(3)
          case ('gas_dom_rsp')
             read(value,*) gas_dom_rsp
          case ('gas_dom_size')
             read(value,*) gas_dom_size
          case ('gas_dom_rin')
             read(value,*) gas_dom_rin
          case ('gas_dom_rout')
             read(value,*) gas_dom_rout
          case ('gas_dom_thickness')
             read(value,*) gas_dom_thickness
          case ('spec_type')
             write(spec_type,'(a)') trim(value)
          case ('spec_gauss_sigma_kms')
             read(value,*) spec_gauss_sigma_kms
          case ('nphot')
             read(value,*) nphot
          case ('ranseed')
             read(value,*) ranseed
          case ('verbose')
             read(value,*) verbose
          case ('cosmo')
             read(value,*) cosmo
          case ('ramses_simple_binary') ! TK
             read(value,*) ramses_simple_binary
          case ('DomDumpFile') ! TK
             write(DomDumpFile,'(a)') trim(value)
          case default
             write(*,'(a,a,a)') '> WARNING: parameter ',trim(name),' unknown '
          end select
       end do
    end if
    close(10)

    return

  end subroutine read_PhotonsFromGas_params

  
  subroutine print_PhotonsFromGas_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)')         '[PhotonsFromGas]'
       write(unit,'(a)')             '# input / output parameters'
       write(unit,'(a,a)')           '  outputfile      = ',trim(outputfile)
       write(unit,'(a,a)')           '  repository      = ',trim(repository)
       write(unit,'(a,i5)')          '  snapnum         = ',snapnum
       write(unit,'(a,a)')           '  DomDumpFile     = ',trim(DomDumpFile)
       write(unit,'(a)')             '# computational domain parameters'
       write(unit,'(a,a)')           '  gas_dom_type      = ',trim(gas_dom_type)
       write(unit,'(a,3(ES10.3,1x))') '  gas_dom_pos       = ',gas_dom_pos(1),gas_dom_pos(2),gas_dom_pos(3)
       select case (trim(gas_dom_type))
       case ('sphere')
          write(unit,'(a,ES10.3)')       '  gas_dom_rsp       = ',gas_dom_rsp
       case ('shell')
          write(unit,'(a,ES10.3)')       '  gas_dom_rin       = ',gas_dom_rin
          write(unit,'(a,ES10.3)')       '  gas_dom_rout      = ',gas_dom_rout
       case('cube')
          write(unit,'(a,ES10.3)')       '  gas_dom_size      = ',gas_dom_size
       case('slab')
          write(unit,'(a,ES10.3)')       '  gas_dom_thickness = ',gas_dom_thickness
       end select
       write(unit,'(a)')             '# Spectral shape '
       write(unit,'(a,a)')           '  spec_type               = ',trim(spec_type)
       select case(trim(spec_type))
       case('Gauss')
          write(unit,'(a,ES10.3,a)')     '  spec_gauss_sigma_kms = ',spec_gauss_sigma_kms, ' ! [km/s]'
       end select
       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,i8)')          '  nphot           = ',nphot
       write(unit,'(a,i8)')          '  ranseed         = ',ranseed
       write(unit,'(a,L1)')          '  verbose         = ',verbose
       write(unit,'(a,L1)')          '  cosmo           = ',cosmo
       write(unit,'(a)')             ' '
    else
       write(*,'(a,a,a)')         '[PhotonsFromGas]'
       write(*,'(a)')             '# input / output parameters'
       write(*,'(a,a)')           '  outputfile      = ',trim(outputfile)
       write(*,'(a,a)')           '  repository      = ',trim(repository)
       write(*,'(a,i5)')          '  snapnum         = ',snapnum
       write(*,'(a,a)')           '  DomDumpFile     = ',trim(DomDumpFile)
       write(*,'(a)')             '# computational domain parameters'
       write(*,'(a,a)')           '  gas_dom_type      = ',trim(gas_dom_type)
       write(*,'(a,3(ES10.3,1x))') '  gas_dom_pos       = ',gas_dom_pos(1),gas_dom_pos(2),gas_dom_pos(3)
       select case (trim(gas_dom_type))
       case ('sphere')
          write(*,'(a,ES10.3)')       '  gas_dom_rsp       = ',gas_dom_rsp
       case ('shell')
          write(*,'(a,ES10.3)')       '  gas_dom_rin       = ',gas_dom_rin
          write(*,'(a,ES10.3)')       '  gas_dom_rout      = ',gas_dom_rout
       case('cube')
          write(*,'(a,ES10.3)')       '  gas_dom_size      = ',gas_dom_size
       case('slab')
          write(*,'(a,ES10.3)')       '  gas_dom_thickness = ',gas_dom_thickness
       end select
       write(*,'(a)')             '# Spectral shape '
       write(*,'(a,a)')           '  spec_type               = ',trim(spec_type)
       select case(trim(spec_type))
       case('Gauss')
          write(*,'(a,es10.3,a)')     '  spec_gauss_sigma_kms = ',spec_gauss_sigma_kms, ' ! [km/s]'
       end select
       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,i8)')          '  nphot           = ',nphot
       write(*,'(a,i8)')          '  ranseed         = ',ranseed
       write(*,'(a,L1)')          '  verbose         = ',verbose
       write(*,'(a,L1)')          '  cosmo           = ',cosmo
       write(*,'(a,L1)')          '  ramses_simple_binary = ',ramses_simple_binary
       write(*,'(a)')             ' '

    end if

    return

  end subroutine print_PhotonsFromGas_params

  
end program PhotonsFromGas


  
