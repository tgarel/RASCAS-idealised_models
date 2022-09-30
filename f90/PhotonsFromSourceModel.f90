program PhotonsFromSourceModel

  use module_photon
  use module_utils, only: isotropic_direction, binary_search
  use module_random
  use module_constants
  use module_ssp_lib

  implicit none

  type(photon_init),dimension(:),allocatable :: photgrid
  integer(kind=4)                            :: iran, i, narg, iseed, ilow2, j
  real(kind=8)                               :: nu, r1, r2, one, betaplus2
  character(2000)                            :: parameter_file
  type(SSPgrid)                              :: NdotGrid
  real(kind=8),allocatable                   :: low_prob2(:), Ndot(:), lbin(:)
  real(kind=8)                               :: lambdamin, lambdamax, lambda_source, scalar, total_flux

  ! --------------------------------------------------------------------------
  ! user-defined parameters - read from section [PhotonsFromSourceModel] of the parameter file
  ! --------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: outputfile = 'ppic.dat'               ! file to which outputs will be written
  ! --- source type 
  character(10)             :: source_type = 'pointlike'             ! type of source model
  real(kind=8),dimension(3) :: source_pos  = (/0.5d0,0.5d0,0.5d0/)   ! position of the source [code units]
  real(kind=8),dimension(3) :: source_vel  = (/0.d0,0.d0,0.d0/)      ! velocity of the source [cm/s]
  integer(kind=4)           :: nphotons    = 1000                    ! number of photons to generate
  ! --- how source shines
  ! Four options here :
  ! - spec_type=='Mono'   : we emit all photons at the same wavelength (in source's frame)
  ! - spec_type=='Gauss'  : we sample a Gaussian distribution ...
  ! - spec_type=='PowLaw' : we sample a power-law continuum between two wavelengths.
  ! - spec_type=='Table'  : we sample a tabulated spectrum, e.g. from Bruzual and Charlot. 
  character(30)             :: spec_type = 'Gauss'              ! May be 'Mono', 'Gauss', 'PowLaw', 'Table' 
  ! parameters for spec_type == 'Mono'
  real(kind=8)              :: spec_mono_l0_Ang = 1215.67       ! emission wavelength [A] -> read from parameter file                    
  ! parameters for spec_type == 'Gauss'
  real(kind=8)              :: spec_gauss_l0_Ang = 1215.67      ! emission wavelength [A] -> read from parameter file
  real(kind=8)              :: spec_gauss_sigma_kms = 10.0      ! line width in velocity [km/s] -> read from parameter file. 
  ! parameters for spec_type == 'PowLaw' : a power-law F_lambda = F_0 * (lambda/lambda_0)**beta  (with F_0 == 1)
  real(kind=8)              :: spec_powlaw_lmin_Ang = 1120.     ! min wavelength to sample (should be in the range where fit was made ...)
  real(kind=8)              :: spec_powlaw_lmax_Ang = 1320.     ! max ...
  real(kind=8)              :: spec_powlaw_l0_Ang   = 1200      ! lambda_0 in the expression above [A]
  real(kind=8)              :: spec_powlaw_beta     = -2.3      ! beta in the expression above. 
  ! parameters for spec_type == 'Table'
  character(1000)           :: spec_SSPdir = '../libs/SSPlibs/' ! the SSP lib directory
  real(kind=8)              :: spec_table_lmin_Ang = 1120.      ! min wavelength to sample
  real(kind=8)              :: spec_table_lmax_Ang = 1320.      ! max ...
  real(kind=8)              :: spec_table_age   = 10.0          ! age of the stellar population to use [Myr]
  real(kind=8)              :: spec_table_met   = 0.02          ! metallicity of the stellar population to use
  real(kind=8)              :: spec_table_mass  = 1.e6          ! mass of the source [Msun]
  ! --- miscelaneous
  integer(kind=4)           :: ranseed = 1234                   ! seed for random generator
  logical                   :: verbose = .true.
  ! --------------------------------------------------------------------------
  ! computed from user-defined parameters
  real(kind=8)              :: spec_mono_nu0                    ! emission frequency [Hz] -> computed from spec_mono_l0_Ang
  real(kind=8)              :: spec_gauss_nu0                   ! central frequency [Hz] -> computed from spec_gauss_l0_Ang

  ! -------------------- read parameters --------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: PhotonsFromSourceModel path/to/params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_PhotonsFromSourceModel_params(parameter_file)
  if (verbose) call print_PhotonsFromSourceModel_params
  ! ------------------------------------------------------------


  if (trim(spec_type) == 'Table') then ! read SSP lib

     call init_ssp_lib(spec_SSPdir)
     lambdamin = spec_table_lmin_Ang
     lambdamax = spec_table_lmax_Ang
     call ssp_lib_extract_subset(lambdamin,lambdamax,NdotGrid) ! extract the SSP age-met grid of Nphotons in lambda range [lambdamin;lambdamax]
     allocate(Ndot(NdotGrid%nlambda))

     if (verbose) then
        write(*,*) 'Age and metallicity of the source: ',spec_table_age, spec_table_met
        write(*,*) 'Mass of the source: ',spec_table_mass
     endif
     
     ! interpolate NdotGrid(lambda)
     call ssp_lib_interpolate(NdotGrid, spec_table_age/1.e3, log10(spec_table_met), Ndot)
     Ndot = Ndot * spec_table_mass   ! nb of photons / s / A

     ! integrate nphotPerSecPerMsun(nlambda)
     ! to compute the total number of photons emitted per second by the source
     call ssp_lib_integrate(NdotGrid%lambda, Ndot, NdotGrid%nlambda, total_flux)
     if (verbose) write(*,*) 'Total luminosity (nb of photons per second): ',total_flux
     
     allocate(low_prob2(NdotGrid%nlambda+1))
     ! compute lbin
     allocate(lbin(NdotGrid%nlambda))
     lbin(1) = NdotGrid%lambda(1)
     do i = 2,NdotGrid%nlambda-1
        lbin(i) = (NdotGrid%lambda(i+1) + NdotGrid%lambda(i))/2.0d0
     enddo
     lbin(NdotGrid%nlambda) = NdotGrid%lambda(NdotGrid%nlambda)

  end if

  
  ! --------------------------------------------------------------------------------------
  ! make source shines
  ! --------------------------------------------------------------------------------------
  one = 1.0d0
  allocate(photgrid(nphotons))
  iran = -abs(ranseed)
  do i=1,nphotons
     photgrid(i)%ID    = i
     ! compute nu
     select case(trim(spec_type))
     case('Mono')
        nu = spec_mono_nu0
        total_flux = one    ! nb of real photons emitted per sec ... This is irrelevant for spec Mono but needed for spec Table
     case('Gauss')
        r1 = ran3(iran)
        r2 = ran3(iran)
        nu = sqrt(-2.*log(r1)) * cos(2.0d0*pi*r2)
        nu = (spec_gauss_sigma_kms * 1d5 * spec_gauss_nu0 / clight) * nu + spec_gauss_nu0
        total_flux = one    ! nb of real photons emitted per sec ... This is irrelevant for spec Gauss but needed for spec Table
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
        if (spec_powlaw_beta == -2.0d0) then
           nu = spec_powlaw_lmin_Ang * exp(r1 * log(spec_powlaw_lmax_Ang / spec_powlaw_lmin_Ang) ) ! this is lbda [A]
           nu = clight / (nu*1e-8) ! this is freq. [Hz]
        else
           betaplus2 = spec_powlaw_beta + 2.0d0
           nu   = (spec_powlaw_lmin_Ang**betaplus2 + r1 * (spec_powlaw_lmax_Ang**betaplus2 - spec_powlaw_lmin_Ang**betaplus2))**(1./betaplus2) ! this is lbda [A]
           nu   = clight / (nu*1e-8) ! this is freq. [Hz]
        end if
        total_flux = one    ! nb of real photons emitted per sec ... This is irrelevant for spec PowLaw but needed for spec Table
        
     case('Table')

        ! 1/ find the frequency/lambda bin
        low_prob2(1) = 0.0d0
        low_prob2(2) = Ndot(1) * (NdotGrid%lambda(2) - NdotGrid%lambda(1))/2.  ! => P(lambda(1)
        do j = 3,NdotGrid%nlambda
           low_prob2(j) = low_prob2(j-1) + Ndot(j-1) * (NdotGrid%lambda(j) - NdotGrid%lambda(j-2))/2.
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
        lambda_source = lbin(ilow2) + r2 * (lbin(ilow2+1)-lbin(ilow2))
        
        nu = clight / (lambda_source*1e-8) ! Hz

     case default
        print*,'ERROR: unknown spec_type :',trim(spec_type)
     end select

     select case(trim(source_type))
     case('pointlike')
        photgrid(i)%x_em = source_pos
        photgrid(i)%v_em = source_vel
     case default
        print*,'ERROR: unknown source_type :',trim(source_type)
     end select

     photgrid(i)%iran  = -i 
     call isotropic_direction(photgrid(i)%k_em,iran)

     ! compute frequency in external frame 
     !photgrid(i)%nu_em = nu
     scalar = photgrid(i)%k_em(1)*photgrid(i)%v_em(1) + photgrid(i)%k_em(2)*photgrid(i)%v_em(2) + photgrid(i)%k_em(3)*photgrid(i)%v_em(3)
     photgrid(i)%nu_em = nu / (1d0 - scalar/clight)

  enddo

  if (trim(spec_type) == 'Table') then ! deallocate stuff
     deallocate(Ndot,lbin,low_prob2)
  end if
  ! --------------------------------------------------------------------------------------
  ! write ICs
  ! --------------------------------------------------------------------------------------
  if (verbose) write(*,*) '--> writing file: ',trim(outputfile)
  open(unit=14, file=trim(outputfile), status='unknown', form='unformatted', action='write')
  write(14) nphotons
  write(14) total_flux ! nb of real photons emitted per sec ... This is relevant only for spec Table
  write(14) ranseed
  write(14) (photgrid(i)%ID,i=1,nphotons)
  write(14) (photgrid(i)%nu_em,i=1,nphotons)
  write(14) (photgrid(i)%x_em(:),i=1,nphotons)
  write(14) (photgrid(i)%k_em(:),i=1,nphotons)
  write(14) (photgrid(i)%iran,i=1,nphotons)
  write(14) (photgrid(i)%v_em(:),i=1,nphotons)
  close(14)
  ! --------------------------------------------------------------------------------------

  if (verbose) then
     write(*,*) '--> work done'
     write(*,*) ' '
  endif
  
  deallocate(photgrid)

contains


  subroutine read_PhotonsFromSourceModel_params(pfile)

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
       if (line(1:24) == '[PhotonsFromSourceModel]') then
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
          case ('source_pos')
             read(value,*) source_pos(1),source_pos(2),source_pos(3)
          case ('source_vel')
             read(value,*) source_vel(1),source_vel(2),source_vel(3)
          case ('source_type')
             write(source_type,'(a)') trim(value)
          case ('verbose')
             read(value,*) verbose
          case ('ranseed')
             read(value,*) ranseed
          case ('spec_type')
             write(spec_type,'(a)') trim(value)
          case ('spec_mono_l0_Ang')
             read(value,*) spec_mono_l0_Ang  ! [A]
             spec_mono_nu0 = clight / spec_mono_l0_Ang * 1d8  ! [Hz]
          case ('spec_gauss_l0_Ang')
             read(value,*) spec_gauss_l0_Ang ! [A]
             spec_gauss_nu0 = clight / spec_gauss_l0_Ang * 1d8  ! [Hz]
         case ('spec_gauss_sigma_kms')
             read(value,*) spec_gauss_sigma_kms
          case ('spec_powlaw_lmin_Ang')
             read(value,*) spec_powlaw_lmin_Ang
          case ('spec_powlaw_lmax_Ang')
             read(value,*) spec_powlaw_lmax_Ang
          case ('spec_powlaw_l0_Ang')
             read(value,*) spec_powlaw_l0_Ang
          case ('spec_powlaw_beta')
             read(value,*) spec_powlaw_beta
          case ('spec_SSPdir')
             write(spec_SSPdir,'(a)') trim(value)
          case ('spec_table_lmin_Ang')
             read(value,*) spec_table_lmin_Ang
          case ('spec_table_lmax_Ang')
             read(value,*) spec_table_lmax_Ang
          case ('spec_table_age')
             read(value,*) spec_table_age
          case ('spec_table_met')
             read(value,*) spec_table_met
          case ('spec_table_mass')
             read(value,*) spec_table_mass
          case ('nPhotonPackets')
             read(value,*) nphotons
          end select
       end do
    end if
    close(10)

    return

  end subroutine read_PhotonsFromSourceModel_params


  subroutine print_PhotonsFromSourceModel_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)')          '[PhotonsFromSourceModel]'
       write(unit,'(a)')              '# input / output parameters'
       write(unit,'(a,a)')            '  outputfile      = ',trim(outputfile)
       write(unit,'(a)')              '# source type parameters'
       write(unit,'(a,a)')            '  source_type     = ',trim(source_type)
       write(unit,'(a,3(ES10.3,1x))') '  source_pos      = ',source_pos(1),source_pos(2),source_pos(3)
       write(unit,'(a,3(ES10.3,1x))') '  source_vel      = ',source_vel(1),source_vel(2),source_vel(3)
       write(unit,'(a)')              '# how source shines'
       write(unit,'(a,i8)')           '  nPhotonPackets  = ',nphotons
       write(unit,'(a,a)')            '  spec_type       = ',trim(spec_type)
       select case(trim(spec_type))
       case('Mono')
          write(unit,'(a,es10.3,a)')     '  spec_mono_l0_Ang     = ',spec_mono_l0_Ang, ' ! [A]'
       case('Gauss')
          write(unit,'(a,es10.3,a)')     '  spec_gauss_l0_Ang    = ',spec_gauss_l0_Ang, ' ! [A]'
          write(unit,'(a,es10.3,a)')     '  spec_gauss_sigma_kms = ',spec_gauss_sigma_kms, ' ! [km/s]'
       case('PowLaw')
          write(unit,'(a,es10.3,a)')     '  spec_powlaw_lmin_Ang = ',spec_powlaw_lmin_Ang, ' ! [A]'
          write(unit,'(a,es10.3,a)')     '  spec_powlaw_lmax_Ang = ',spec_powlaw_lmax_Ang, ' ! [A]'
          write(unit,'(a,es10.3,a)')     '  spec_powlaw_l0_Ang   = ',spec_powlaw_l0_Ang, ' ! [A]'
          write(unit,'(a,es10.3)')       '  spec_powlaw_beta     = ',spec_powlaw_beta
       case('Table')
          write(unit,'(a,a)')            '  spec_SSPdir          = ',trim(spec_SSPdir)
          write(unit,'(a,es10.3,a)')     '  spec_table_lmin_Ang  = ',spec_table_lmin_Ang, ' ! [A]'
          write(unit,'(a,es10.3,a)')     '  spec_table_lmax_Ang  = ',spec_table_lmax_Ang, ' ! [A]'
          write(unit,'(a,es10.3,a)')     '  spec_table_age       = ',spec_table_age, ' ! [Myr]'
          write(unit,'(a,es10.3)')       '  spec_table_met       = ',spec_table_met
          write(unit,'(a,es10.3,a)')     '  spec_table_mass      = ',spec_table_mass, ' ! [Msun]'
       case default
          print*,'ERROR: unknown spec_type :',trim(spec_type)
       end select
       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,i8)')          '  ranseed         = ',ranseed
       write(unit,'(a,L1)')          '  verbose         = ',verbose
       write(unit,'(a)')             ' '
    else
       write(*,'(a)')              '--------------------------------------------------------------------------------'
       write(*,'(a)')              ' '
       write(*,'(a,a,a)')          '[PhotonsFromSourceModel]'
       write(*,'(a)')              '# input / output parameters'
       write(*,'(a,a)')            '  outputfile    = ',trim(outputfile)
       write(*,'(a)')              '# source type parameters'
       write(*,'(a,a)')            '  source_type   = ',trim(source_type)
       write(*,'(a,3(ES10.3,1x))') '  source_pos    = ',source_pos(1),source_pos(2),source_pos(3)
       write(*,'(a,3(ES10.3,1x))') '  source_vel    = ',source_vel(1),source_vel(2),source_vel(3)
       write(*,'(a)')              '# how source shines'
       write(*,'(a,i8)')           '  nPhotonPackets = ',nphotons
       write(*,'(a,a)')            '  spec_type      = ',trim(spec_type)
       select case(trim(spec_type))
       case('Mono')
          write(*,'(a,es10.3,a)')     '  spec_mono_l0_Ang     = ',spec_mono_l0_Ang, ' ! [A]'
       case('Gauss')
          write(*,'(a,es10.3,a)')     '  spec_gauss_l0_Ang    = ',spec_gauss_l0_Ang, ' ! [A]'
          write(*,'(a,es10.3,a)')     '  spec_gauss_sigma_kms = ',spec_gauss_sigma_kms, ' ! [km/s]'
       case('PowLaw')
          write(*,'(a,es10.3,a)')     '  spec_powlaw_lmin_Ang = ',spec_powlaw_lmin_Ang, ' ! [A]'
          write(*,'(a,es10.3,a)')     '  spec_powlaw_lmax_Ang = ',spec_powlaw_lmax_Ang, ' ! [A]'
          write(*,'(a,es10.3,a)')     '  spec_powlaw_l0_Ang   = ',spec_powlaw_l0_Ang, ' ! [A]'
          write(*,'(a,es10.3)')       '  spec_powlaw_beta     = ',spec_powlaw_beta
       case('Table')
          write(*,'(a,a)')            '  spec_SSPdir          = ',trim(spec_SSPdir)
          write(*,'(a,es10.3,a)')     '  spec_table_lmin_Ang  = ',spec_table_lmin_Ang, ' ! [A]'
          write(*,'(a,es10.3,a)')     '  spec_table_lmax_Ang  = ',spec_table_lmax_Ang, ' ! [A]'
          write(*,'(a,es10.3,a)')     '  spec_table_age       = ',spec_table_age, ' ! [Myr]'
          write(*,'(a,es10.3)')       '  spec_table_met       = ',spec_table_met
          write(*,'(a,es10.3,a)')     '  spec_table_mass      = ',spec_table_mass, ' ! [Msun]'
       case default
          print*,'ERROR: unknown spec_type :',trim(spec_type)
       end select
       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,i8)')          '  ranseed       = ',ranseed
       write(*,'(a,L1)')          '  verbose       = ',verbose
       write(*,'(a)')             ' '       
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ' '
    end if

    return

  end subroutine print_PhotonsFromSourceModel_params


end program PhotonsFromSourceModel
