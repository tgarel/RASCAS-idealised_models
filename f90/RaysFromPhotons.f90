program RaysFromPhotons

  use module_ray
  use module_photon

  implicit none

  type(photon_current),allocatable :: pgrid(:)
  type(ray_type),allocatable       :: rays(:)
  integer(kind=4)                  :: nrays,i,narg,n
  character(2000)                  :: parameter_file
  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [RaysFromSourceModel] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: PhotonFile = 'RascasResult.dat' ! file containing results of a rascas run ... 
  character(2000)           :: outputfile = 'RaysIC.dat'       ! file to which outputs will be written
  real(kind=8)              :: kobs(3)    = (/0.,0.,1./)             ! direction of observation . 
  logical                   :: verbose = .true.
  ! --------------------------------------------------------------------------


  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: RaysFromPhotons path/to/params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_RaysFromPhotons_params(parameter_file)
  if (verbose) call print_RaysFromPhotons_params
  ! ------------------------------------------------------------

  ! Read photons (result of a rascas run)
  if (verbose) write(*,*) 'reading photons ...'
  call read_photon_dump(PhotonFile,pgrid)

  ! spawn rays from these photons in a chose direction
  if (verbose) write(*,*) 'defining rays ... '
  nrays = size(pgrid)
  allocate(rays(nrays))
  n = 0
  do i=1,nrays
     if (pgrid(i)%status >0) then !== 1) then 
        n = n + 1
        rays(n)%ID     = pgrid(i)%ID
        rays(n)%x_em   = pgrid(i)%xlast
        rays(n)%k_em   = kobs
        rays(n)%nu_ext = pgrid(i)%nu_ext
     end if
  end do
  
  ! --------------------------------------------------------------------------------------
  ! write ICs
  ! --------------------------------------------------------------------------------------
  if (verbose) write(*,*) '--> writing file'
  open(unit=14, file=trim(outputfile), status='unknown', form='unformatted', action='write')
  write(14) n
  write(14) (rays(i)%ID,i=1,n)
  write(14) (rays(i)%nu_ext,i=1,n)
  write(14) (rays(i)%x_em(:),i=1,n)
  write(14) (rays(i)%k_em(:),i=1,n)
  close(14)
  ! --------------------------------------------------------------------------------------

contains

  subroutine read_RaysFromPhotons_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present
    real(kind=8)    :: norm 
    
    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:17) == '[RaysFromPhotons]') then
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
          case ('PhotonFile')
             write(PhotonFile,'(a)') trim(value)
          case ('kobs')
             read(value,*) kobs(1:3)
             ! force normalisation of kobs
             norm = kobs(1)*kobs(1)+kobs(2)*kobs(2)+kobs(3)*kobs(3)
             if (norm > 0) then
                kobs = kobs / sqrt(norm)
             else
                print*,'kobs has to be non zero... '
                stop
             end if
          case ('verbose')
             read(value,*) verbose
          end select
       end do
    end if
    close(10)
    return

  end subroutine read_RaysFromPhotons_params


  subroutine print_RaysFromPhotons_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)')         '[RaysFromPhotons]'
       write(unit,'(a)')             '# input / output parameters'
       write(unit,'(a,a)')           '  outputfile      = ',trim(outputfile)
       write(unit,'(a,a)')           '  PhotonFile      = ',trim(PhotonFile)
       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,L1)')          '  verbose         = ',verbose
       write(unit,'(a)')             ' '
    else
       write(*,'(a,a,a)')         '[RaysFromPhotons]'
       write(*,'(a)')             '# input / output parameters'
       write(*,'(a,a)')           '  outputfile      = ',trim(outputfile)
       write(*,'(a,a)')           '  PhotonFile      = ',trim(PhotonFile)
       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,L1)')          '  verbose         = ',verbose
       write(*,'(a)')             ' '
    end if

    return

  end subroutine print_RaysFromPhotons_params

end program RaysFromPhotons
  
