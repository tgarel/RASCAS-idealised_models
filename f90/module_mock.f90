module module_mock

  use module_constants, only:clight
  
  private

  type mockObs 
     ! direction of observation
     real(kind=8) :: kobs(3),kobs_perp_1(3),kobs_perp_2(3)
     ! center of observation [code units]
     real(kind=8) :: center(3)
     ! --- FLUX --- 
     real(kind=8)             :: flux_aperture         ! collect only photons within a circle of radius flux_aperture
     real(kind=8)             :: flux                  ! flux 
     ! --- SPECTRUM --- 
     integer(kind=4)          :: spec_npix = 0         ! nb of pixels
     real(kind=8)             :: spec_aperture         ! collect only photons within a circle of radius spec_aperture
     real(kind=8)             :: spec_lmin, spec_lmax  ! min/max of spectrum [Angstrom]
     real(kind=8),allocatable :: spectrum(:)           ! actual spectrum. 
     ! --- IMAGE --- 
     integer(kind=4)          :: image_npix = 0        ! nb of pixels across image
     real(kind=8)             :: image_side            ! extent of observation [box units]
     real(kind=8),allocatable :: image(:,:)            ! actual image. 
     ! --- CUBE ---
     ! ... 
     ! useful variables
     real(kind=8)    :: flux_aperture2, spec_aperture2 
     logical         :: compute_image
     logical         :: compute_spectrum
  end type mockObs
  type(mockObs),allocatable :: mock(:)

  ! parameters in the [mock] section :
  integer(kind=4) :: nDirections = 0
  character(2000) :: mock_parameter_file
  character(2000) :: mock_outputfilename  ! Prefix for output files (including absolute path) -> will be followed by "_image.xxxxx" or "_spectrum.xxxxx", with xxxxx the cpu number.


  ! global parameter setting peeling-off on or off.
  logical         :: peeling_off

  
  ! public variables 
  public :: peeling_off, nDirections, mock
  ! public functions 
  public :: read_mock_params, mock_line_of_sight, mock_point_in_spectral_aperture, mock_point_in_flux_aperture, mock_point_in_image
  public :: mock_projected_pos, peel_to_flux, peel_to_map, peel_to_spec, dump_mocks
contains

  
  subroutine mock_init()
    
    implicit none
    integer(kind=4) :: idir,unit=33
    
    ! initialise mocks
    if (nDirections > 0) then
       peeling_off = .true.
       
       allocate(mock(nDirections))
       open(unit=unit,file=mock_parameter_file,status='old',action='read',form='formatted')
       do idir = 1,nDirections
          call read_a_mock_param_set(unit,idir)
          
          ! initialise flux
          mock(idir)%flux = 0.0d0
          mock(idir)%flux_aperture2 = mock(idir)%flux_aperture*mock(idir)%flux_aperture

          ! initialise spectrum 
          mock(idir)%compute_spectrum = .false.
          if (mock(idir)%spec_npix > 0) then 
             allocate(mock(idir)%spectrum(mock(idir)%spec_npix))
             mock(idir)%spectrum = 0.0d0
             mock(idir)%compute_spectrum = .true.
             mock(idir)%spec_aperture2 = mock(idir)%spec_aperture*mock(idir)%spec_aperture
          end if

          ! initialise image 
          mock(idir)%compute_image = .false.
          if (mock(idir)%image_npix > 0) then 
             allocate(mock(idir)%image(mock(idir)%image_npix,mock(idir)%image_npix))
             mock(idir)%image = 0.0d0
             mock(idir)%compute_image = .true.
          end if
          
          ! define direction of observation (normalise vector)
          mock(idir)%kobs = mock(idir)%kobs / sqrt(mock(idir)%kobs(1)*mock(idir)%kobs(1)+&
               & mock(idir)%kobs(2)*mock(idir)%kobs(2)+mock(idir)%kobs(3)*mock(idir)%kobs(3))
          ! define basis for sky plane
          if (mock(idir)%kobs(1) < 1.d0) then
             ! kobs_perp_1 is cross prod. of kobs and x 
             mock(idir)%kobs_perp_1(1) = 0.
             mock(idir)%kobs_perp_1(2) = mock(idir)%kobs(3)
             mock(idir)%kobs_perp_1(3) = -mock(idir)%kobs(2)
             mock(idir)%kobs_perp_1 = mock(idir)%kobs_perp_1 / sqrt(mock(idir)%kobs_perp_1(1)*mock(idir)%kobs_perp_1(1)+&
                  & mock(idir)%kobs_perp_1(2)*mock(idir)%kobs_perp_1(2)+mock(idir)%kobs_perp_1(3)*mock(idir)%kobs_perp_1(3))
             ! kobs_perp_2 = kobs x kobs_perp_1
             mock(idir)%kobs_perp_2(1) = mock(idir)%kobs(2)*mock(idir)%kobs_perp_1(3) - mock(idir)%kobs(3)*mock(idir)%kobs_perp_1(2)
             mock(idir)%kobs_perp_2(2) = -mock(idir)%kobs(1)*mock(idir)%kobs_perp_1(3)
             mock(idir)%kobs_perp_2(3) = mock(idir)%kobs(1)*mock(idir)%kobs_perp_1(2)
             mock(idir)%kobs_perp_2 = mock(idir)%kobs_perp_2 / sqrt(mock(idir)%kobs_perp_2(1)*mock(idir)%kobs_perp_2(1)+ &
                  & mock(idir)%kobs_perp_2(2)*mock(idir)%kobs_perp_2(2)+mock(idir)%kobs_perp_2(3)*mock(idir)%kobs_perp_2(3))       
          else
             ! kobs is aligned with x -> use kobs_perp_1 = y, kobs_perp_2 = z
             mock(idir)%kobs_perp_1 = (/0.0d0,1.0d0,0.0d0/)
             mock(idir)%kobs_perp_2 = (/0.0d0,0.0d0,1.0d0/)
          end if
       end do
       close(unit)
    else
       peeling_off = .false. 
    end if
    
    return

  contains

    subroutine read_a_mock_param_set(unit,idir)
      implicit none
      integer(kind=4),intent(in) :: unit,idir
      read(unit,*) mock(idir)%kobs
      read(unit,*) mock(idir)%center
      read(unit,*) mock(idir)%flux_aperture
      read(unit,*) mock(idir)%spec_npix, mock(idir)%spec_aperture, mock(idir)%spec_lmin, mock(idir)%spec_lmax
      read(unit,*) mock(idir)%image_npix, mock(idir)%image_side
      return 
    end subroutine read_a_mock_param_set
    
  end subroutine mock_init

  
  function mock_line_of_sight(idir)
    implicit none
    integer(kind=4),intent(in)  :: idir
    real(kind=8) :: mock_line_of_sight(3)
    mock_line_of_sight = mock(idir)%kobs
  end function mock_line_of_sight

  
  function mock_point_in_spectral_aperture(pp,idir)
    implicit none
    logical  :: mock_point_in_spectral_aperture
    real(kind=8),intent(in) :: pp(2)
    integer(kind=4),intent(in) :: idir
    mock_point_in_spectral_aperture = (pp(1)*pp(1) + pp(2)*pp(2) < mock(idir)%spec_aperture2)
    return
  end function mock_point_in_spectral_aperture

  
  function mock_point_in_flux_aperture(pp,idir)
    implicit none
    logical  :: mock_point_in_flux_aperture
    real(kind=8),intent(in) :: pp(2)
    integer(kind=4),intent(in) :: idir 
    mock_point_in_flux_aperture = (pp(1)*pp(1) + pp(2)*pp(2) < mock(idir)%flux_aperture2)
    return
  end function mock_point_in_flux_aperture

  
  function mock_point_in_image(pp,idir)
    implicit none
    logical  :: mock_point_in_image
    real(kind=8),intent(in) :: pp(2)
    integer(kind=4),intent(in) :: idir
    real(kind=8) :: dx
    dx = 0.5d0*mock(idir)%image_side
    mock_point_in_image = ((pp(1) < dx) .and. (pp(1) > -dx) .and. (pp(2) < dx) .and. (pp(2) > -dx))
    return
  end function mock_point_in_image

  
  subroutine mock_projected_pos(pos,pp,idir)
    implicit none
    real(kind=8),intent(in)    :: pos(3)
    real(kind=8),intent(inout) :: pp(2)
    integer(kind=4),intent(in) :: idir
    pp(1) = (pos(1)-mock(idir)%center(1))*mock(idir)%kobs_perp_1(1) + (pos(2)-mock(idir)%center(2))*mock(idir)%kobs_perp_1(2) + (pos(3)-mock(idir)%center(3))*mock(idir)%kobs_perp_1(3)
    pp(2) = (pos(1)-mock(idir)%center(1))*mock(idir)%kobs_perp_2(1) + (pos(2)-mock(idir)%center(2))*mock(idir)%kobs_perp_2(2) + (pos(3)-mock(idir)%center(3))*mock(idir)%kobs_perp_2(3)
    return 
  end subroutine mock_projected_pos
    

  subroutine peel_to_flux(peel_contribution,idir)
    implicit none
    real(kind=8),intent(in) :: peel_contribution
    integer(kind=4),intent(in) :: idir
    mock(idir)%flux = mock(idir)%flux + peel_contribution
  end subroutine peel_to_flux
    
  
  subroutine peel_to_map(pp,peel_contribution,idir)

    implicit none

    real(kind=8),intent(in) :: pp(2),peel_contribution
    integer(kind=4),intent(in) :: idir
    integer(kind=4) :: ix,iy,n
    real(kind=8) :: dx
    dx = mock(idir)%image_side
    n  = mock(idir)%image_npix
    ix = int((pp(1) + 0.5d0 * dx) /dx * n)
    iy = int((pp(2) + 0.5d0 * dx) /dx * n)
    if (ix>0 .and. ix<=n .and. iy>0 .and. iy<=n) mock(idir)%image(ix,iy) = mock(idir)%image(ix,iy) + peel_contribution
    
    return

  end subroutine peel_to_map

  
  subroutine peel_to_spec(peel_nu,peel_contribution,idir)
    implicit none

    real(kind=8),intent(in) :: peel_nu,peel_contribution
    integer(kind=4),intent(in) :: idir
    real(kind=8)            :: lambda 
    integer(kind=4) :: i

    lambda = clight / peel_nu * 1d8 ! [Angstrom]
    i = int( (lambda - mock(idir)%spec_lmin) / (mock(idir)%spec_lmax - mock(idir)%spec_lmin) * mock(idir)%spec_npix)
    if ((i > 0) .and. (i<=mock(idir)%spec_npix)) then
       mock(idir)%spectrum(i) = mock(idir)%spectrum(i) + peel_contribution
    end if
    
    return
  end subroutine peel_to_spec

  subroutine dump_mocks(rank)
    implicit none
    integer(kind=4),intent(in) :: rank
    character(2000)            :: filename
    integer(kind=4)            :: i,j,iunit=133,sunit=134,funit=135,idir
    logical :: iopen=.false.,sopen=.false.,fopen=.false.

    
    do idir = 1,nDirections
       ! save flux
       if (.not. fopen) then 
          write(filename,'(a,a,i5.5)') trim(mock_outputfilename),'_flux.',rank
          open(unit=funit,file=filename,form='unformatted',status='unknown')
          fopen = .true.
       end if
       write(funit) mock(idir)%flux_aperture, mock(idir)%flux
       ! save spectrum
       if (mock(idir)%compute_spectrum) then 
          if (.not. sopen) then 
             write(filename,'(a,a,i5.5)') trim(mock_outputfilename),'_spectrum.',rank
             open(unit=sunit,file=filename,form='unformatted',status='unknown')
             sopen = .true.
          end if
          write(sunit) mock(idir)%spec_npix
          write(sunit) mock(idir)%spec_aperture,mock(idir)%spec_lmin,mock(idir)%spec_lmax
          write(sunit) (mock(idir)%spectrum(i),i=1,mock(idir)%spec_npix)
       end if
       ! save image
       if (mock(idir)%compute_image) then 
          if (.not. iopen) then 
             write(filename,'(a,a,i5.5)') trim(mock_outputfilename),'_image.',rank
             open(unit=iunit,file=filename,form='unformatted',status='unknown')
             iopen = .true.
          end if
          write(iunit) mock(idir)%image_npix
          write(iunit) mock(idir)%image_side
          write(iunit) (mock(idir)%center(i),i=1,3)
          write(iunit) ((mock(idir)%image(i,j),i=1,mock(idir)%image_npix),j=1,mock(idir)%image_npix)
       end if
    end do
    
    if (iopen) close(iunit)
    if (sopen) close(sunit)
    if (fopen) close(funit)

    return
  end subroutine dump_mocks
    

  subroutine read_mock_params(pfile)
    
    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    !
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000)         :: line,name,value
    integer(kind=4)         :: err,i
    logical                 :: section_present

    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:6) == '[mock]') then
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
          case ('nDirections')
             read(value,*) nDirections
          case ('mock_parameter_file')
             write(mock_parameter_file,'(a)') trim(value)
          case ('mock_outputfilename')
             write(mock_outputfilename,'(a)') trim(value)
          end select
       end do
    end if
    close(10)

    call mock_init
       
    return

  end subroutine read_mock_params


end module module_mock
