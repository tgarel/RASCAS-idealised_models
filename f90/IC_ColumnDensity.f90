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
  real(kind=8),allocatable :: star_pos(:,:),star_age(:),star_mass(:),star_vel(:,:),star_met(:), star_L(:,:)
  real(kind=8)    :: star_L_tot, L_faint
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
  character(2000)           :: method  = 'fromstars'  !Either fromstar, or onepoint
  !real(kind=8)              :: frac_lum = 1d0
  logical                   :: verbose = .true.

  real(kind=8),dimension(3) :: x_phot = (/ 0.5, 0.5, 0.5 /)
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
     nSEDgroups = get_nOptBins()   !From module_spectra, Have to use nSEDgroups = 1 !

     !initialize SED properties
     call init_SED_table()

     !compute SED luminosities,  in #photons/s
     allocate(star_L(nSEDgroups,nphot)) ; star_L_tot = 0d0
     do i=1,nphot
        call inp_sed_table(star_age(i)/1d3, star_met(i), 1, .false., star_L(:,i))    !Third variable :  1 for L[#photons/s],  3 for mean energy in bin,  2+2*Iion for mean cross-section,  Iion: 1:HI,2:HeI,3:HeII,4:SiI, etc
        star_L(:,i) = star_L(:,i)*star_mass(i)/msun
        star_L_tot = star_L_tot + star_L(1,i)
     end do

     ! call hpsort_real(nphot,star_L(1,:))
     ! print*,star_L(1,nphot)/star_L_tot, star_L(1,nphot-1)/star_L_tot
     ! i=0 ; L_faint = 0d0
     ! do while(L_faint/star_L_tot < 1d0-frac_lum)
     !    i=i+1
     !    L_faint = L_faint + star_L(1,i)
     ! end do
     ! print*,i,L_faint,star_L_tot
     

     print*, 'number of stars in domain = ', nphot
     allocate(x_em(3,nphot), k_em(3,nphot), weight(nphot))


     do i=1,nphot
        x_em(:,i) = star_pos(:,i)
        weight(i) = star_L(1,i)
     end do
     
     deallocate(star_pos, star_mass, star_met, star_age)
  elseif(method == 'onepoint') then
     nphot = 1
     allocate(x_em(3,nphot), k_em(3,nphot), weight(nphot))
     x_em(:,1) = x_phot(:)
     weight(1) = 1
  else
     
     print*, 'method ', method, ' not known for the moment'
  end if


  ! --------------------------------------------------------------------------------------
  ! write ICs
  ! --------------------------------------------------------------------------------------
  if (verbose) write(*,*) '> writing file'
  open(unit=14, file=trim(outputfile), status='unknown', form='unformatted', action='write')
  write(14) nphot      ! nb of MC photons
  write(14) (x_em(:,i), i=1,nphot)
  write(14) (weight(i), i=1,nphot)
  close(14)
  ! --------------------------------------------------------------------------------------


contains

  ! subroutine hpsort_real(n,ra) 

  !   implicit none

  !   integer(kind=4),intent(in)    :: n
  !   real(kind=8),intent(inout) :: ra(n)
  !   real(kind=8) :: rra
  !   integer(kind=4)               :: l, ir, i, j 

  !   if (n < 2) return
  !   l=n/2+1
  !   ir=n
  !   do
  !      if(l.gt.1) then 
  !         l=l-1
  !         rra=ra(l) 
  !      else
  !         rra=ra(ir) 
  !         ra(ir)=ra(1) 
  !         ir=ir-1 
  !         if(ir.eq.1)then
  !            ra(1)=rra
  !            return 
  !         endif
  !      endif
  !      i=l
  !      j=l+l
  !      do while (j.le.ir)
  !         if(j.lt.ir)then
  !            if(ra(j).lt.ra(j+1)) j=j+1
  !         end if
  !         if(rra.lt.ra(j))then 
  !            ra(i)=ra(j)
  !            i=j
  !            j=j+j 
  !         else
  !            j=ir+1 
  !         endif
  !      end do
  !      ra(i)=rra
  !   end do

  !   return

  ! end subroutine hpsort_real

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
          case ('x_phot')
             read(value,*) x_phot(1),x_phot(2),x_phot(3)
          case ('nphot')
             read(value,*) nphot
          case('method')
             write(method,'(a)') trim(value)
          case ('verbose')
             read(value,*) verbose
         ! case ('frac_lum')
           !  read(value,*) frac_lum
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
       write(unit,'(a,3(ES10.3,1x))')'  star_dom_pos       = ',star_dom_pos(1),star_dom_pos(2),star_dom_pos(3)

       write(unit,'(a,3(ES10.3,1x))')'  x_phot             = ',x_phot(1),x_phot(2),x_phot(3)

       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,i8)')          '  nphot           = ',nphot
       write(unit,'(a,a)')           '  method          = ',trim(method)
       !write(unit,'(a,ES10.3,1x)')   '  frac_lum        = ',frac_lum
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

       write(*,'(a,3(ES10.3,1x))') '  x_phot             = ',x_phot(1),x_phot(2),x_phot(3)

       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,i8)')          '  nphot           = ',nphot
       write(*,'(a,a)')           '  method          = ',trim(method)
       !write(*,'(a,ES10.3,1x)')   '  frac_lum        = ',frac_lum
       write(*,'(a,L1)')          '  verbose         = ',verbose
       write(*,'(a)')             ' '
       call print_ramses_params
       call print_spectra_params
    end if

    return

  end subroutine print_IC_ColumnDensity_params

end program main

