program main

  ! generate photons emitted by star particles within a given domain,  simplified case to compute column densities. No frequencies for the photons

  use module_photon
  use module_utils
  use module_domain
  use module_random
  use module_constants
  use module_ramses
  use module_spectra

  implicit none

  type(domain)    :: emission_domain
  character(2000) :: parameter_file
  real(kind=8),allocatable :: star_pos(:,:),star_age(:),star_mass(:),star_vel(:,:),star_met(:), star_L(:,:), star_L_tot(:)
  integer(kind=4) :: iran,i,nstars,narg, nSEDgroups
  integer(kind=8) :: ilast,j

  real(kind=8), allocatable :: x_em(:,:), k_em(:,:), weight(:)

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

  ! --- miscelaneous
  integer(kind=4)           :: nphot   = 1000000      ! number of photons to generate
  real(kind=8),dimension(3) :: direction = (/1d0,0d0,0d0/)
  character(2000)           :: method  = 'fromstars' 
  logical                   :: verbose = .true.
  ! --------------------------------------------------------------------------



  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: IC_Columndensity path/to/params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_IC_ColumnDensity_params(parameter_file)
  if (verbose) call print_IC_ColumnDensity_params
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



  if(method == 'fromstars') then

     ! --------------------------------------------------------------------------------------
     ! read star particles within domain
     ! --------------------------------------------------------------------------------------
     if (verbose) write(*,*) '> reading star particles'
     call ramses_read_stars_in_domain(repository,snapnum,emission_domain,star_pos,star_age,star_mass,star_vel,star_met)
     deallocate(star_vel)
     ! --------------------------------------------------------------------------------------
     nphot = size(star_age)
     nSEDgroups = get_nOptBins()   !From module_spectra

     !initialize SED properties
     call init_SED_table()

     !compute SED luminosities,  in #photons/s
     allocate(star_L(nSEDgroups,nphot), star_L_tot(nSEDgroups)) ; star_L_tot = 0d0
     do i=1,nphot
        call inp_sed_table(star_age(i)/1d3, star_met(i), 1, .false., star_L(:,i))    !Third variable :  1 for L[#photons/s],  3 for mean energy in bin,  2+2*Iion for mean cross-section,  Iion: 1:HI,2:HeI,3:HeII,4:SiI, etc
        star_L(:,i) = star_L(:,i)*star_mass(i)/msun
        star_L_tot(:) = star_L_tot(:) + star_L(:,i)
     end do

     print*, 'number of stars in domain = ', nphot
     allocate(x_em(3,nphot), k_em(3,nphot), weight(nphot))

     do i=1,nphot
        x_em(:,i) = star_pos(:,i)
        k_em(:,i) = direction(:)
        weight(i) = star_L(1,i)/star_L_tot(1)    !Have to use nSEDgroups = 1 !
     end do
     deallocate(star_pos, star_mass, star_met, star_age)
  else
     print*, 'method ', method, ' not known for the moment'
  end if


  ! --------------------------------------------------------------------------------------
  ! write ICs
  ! --------------------------------------------------------------------------------------
  if (verbose) write(*,*) '> writing file'
  open(unit=14, file=trim(outputfile), status='unknown', form='unformatted', action='write')
  write(14) nphot      ! nb of MC photons
  do i=1,nphot
     write(14) x_em(1,i), x_em(2,i), x_em(3,i), k_em(1,i), k_em(2,i), k_em(3,i), weight(i)
  end do
  close(14)
  ! --------------------------------------------------------------------------------------


contains

  subroutine read_IC_ColumnDensity_params(pfile)

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
       if (line(1:18) == '[IC_ColumnDensity]') then
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
          case ('nphot')
             read(value,*) nphot
          case ('direction')
             read(value,*) direction(1), direction(2), direction(3)
          case('method')
             write(method,'(a)') trim(value)
          case ('verbose')
             read(value,*) verbose
          case default
             write(*,'(a,a,a)') '> WARNING: parameter ',trim(name),' unknown '
          end select
       end do
    end if
    close(10)

    call read_ramses_params(pfile)
    call read_spectra_params(pfile)


    return

  end subroutine read_IC_ColumnDensity_params


  subroutine print_IC_ColumnDensity_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)')         '[IC_ColumnDensity]'
       write(unit,'(a)')             '# input / output parameters'
       write(unit,'(a,a)')           '  outputfile      = ',trim(outputfile)
       write(unit,'(a,a)')           '  repository      = ',trim(repository)
       write(unit,'(a,i5)')          '  snapnum         = ',snapnum
       write(unit,'(a)')             '# computational domain parameters'
       write(unit,'(a,a)')           '  star_dom_type      = ',trim(star_dom_type)
       write(unit,'(a,3(ES10.3,1x))') '  star_dom_pos       = ',star_dom_pos(1),star_dom_pos(2),star_dom_pos(3)

       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,i8)')          '  nphot           = ',nphot
       write(unit,'(a,a)')           '  method          = ',trim(method)
       write(unit,'(a,3(ES10.3,1x))')'  direction       = ',direction(1), direction(2), direction(3)
       write(unit,'(a,L1)')          '  verbose         = ',verbose
       write(unit,'(a)')             ' '
       call print_ramses_params(unit)
       call print_spectra_params(unit)
    else
       write(*,'(a,a,a)')         '[IC_ColumnDensity]'
       write(*,'(a)')             '# input / output parameters'
       write(*,'(a,a)')           '  outputfile      = ',trim(outputfile)
       write(*,'(a,a)')           '  repository      = ',trim(repository)
       write(*,'(a,i5)')          '  snapnum         = ',snapnum
       write(*,'(a)')             '# computational domain parameters'
       write(*,'(a,a)')           '  star_dom_type      = ',trim(star_dom_type)
       write(*,'(a,3(ES10.3,1x))') '  star_dom_pos       = ',star_dom_pos(1),star_dom_pos(2),star_dom_pos(3)

       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,i8)')          '  nphot           = ',nphot
       write(*,'(a,3(ES10.3,1x))')'  direction       = ',direction(1), direction(2), direction(3)
       write(*,'(a,a)')           '  method          = ',trim(method)
       write(*,'(a,L1)')          '  verbose         = ',verbose
       write(*,'(a)')             ' '
       call print_ramses_params
       call print_spectra_params
    end if

    return

  end subroutine print_IC_ColumnDensity_params

end program main

