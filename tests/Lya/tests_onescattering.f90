program tests_onescattering

  ! -----------------------------------------------------------------------------------------------------------------
  ! Test Lya radiative transfer
  ! -----------------------------------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------------------------------
  ! - tests.f90 runs one-scattering interactions for either HI, dust or D using ../../f90/ modules
  ! - 'T', 'x_in', 'nphotons', 'param_file', 'outfile', 'outputdir' values must be specified here 
  ! - Other parameters to edit in 'params_test.dat' (recoil, angular redist. etc)
  ! - x_out and k_out saved to file
  ! -----------------------------------------------------------------------------------------------------------------
  
  use module_HI_model
  use module_dust_model
  use module_utils
  use module_constants
  use module_D_model

  implicit none

  integer                    :: iran,i,j
  real(kind=8)               :: nu_ext,x_out
  real(kind=8)               :: nu_cell
  real(kind=8), dimension(3) :: k
  real(kind=8), dimension(3) :: vcell
  character(1000)            :: outfile,param_file
  character(100)             :: temp_st,x_input_st
  real(kind=8)               :: delta_nu_doppler, vth
  integer(kind=4)            :: ilost
  character(10)              :: aaa
  real(kind=8)               :: x_input
  integer(kind=4)            :: narg
  
  !real(kind=8)               :: x_atom
  !! HI  
  real(kind=8),parameter     :: lambda_0      = 1215.67d0        ! [A] Lya wavelength
  !! D
  real(kind=8),parameter     :: lambda_0_Deut = 1215.34d0   ! wavelength of Lya of Deuterium [A]
  
  real(kind=8),parameter     :: lambda_0_cm      = lambda_0 / cmtoA                   ! cm
  real(kind=8),parameter     :: lambda_0_Deut_cm = lambda_0_Deut / cmtoA              ! cm
  real(kind=8),parameter     :: nu_0             = clight / lambda_0_cm               ! Hz
  real(kind=8),parameter     :: nu_0_Deut        = clight / lambda_0_Deut_cm          ! Hz
    
  !! Test parameters
  character(10)                         :: species   = 'H'
  real(kind=8)                          :: T         = 1.e4
 ! real(kind=8)               :: x_in     = 0.0
  real(kind=8),dimension(:),allocatable :: x_in
  integer(kind=4)                       :: n_xin     = 1
  integer(kind=4)                       :: nphotons  = 1000000     ! number of photons to cast
  character(200)                        :: outputdir = 'LyaTests_output_data/'
  logical                               :: recoil    = .false.
  logical                               :: isotropic = .false.     ! if set to true, scattering events will be isotropic [default is false]
  real(kind=8)                          :: albedo    = 0.32
  real(kind=8)                          :: g_dust    = 0.73

  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: PhotonsFromSourceModel path/to/params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, param_file)
  call read_test_parameters(param_file)


  
  !----------------------------------
  !   test HI 
  !----------------------------------

  if (species .eq. 'H') then
     call read_HI_params(param_file)
     call print_HI_params
     
     !! Convert T into well-formatted string
     write(temp_st,'(es14.0)') T
!!$     print*,T
!!$     print*,trim(temp_st)
     write(temp_st,'(a6)') trim(adjustl(temp_st))
     print*,temp_st
     
     vcell = 0.                     ! cm/s
     vth   = sqrt((2.0d0*kb/mp)*T)  ! cm/s
     delta_nu_doppler  = vth / lambda_0_cm 

     do j=1,n_xin
        x_input = x_in(j)
        
        if (x_input < 10.0) then
           !! Convert x_in into well-formatted string
           write(x_input_st,'(f6.2)') x_input
           write(x_input_st,'(a1)') trim(adjustl(x_input_st))
        else
           !! Convert x_in into well-formatted string
           write(x_input_st,'(f6.2)') x_input
           write(x_input_st,'(a2)') trim(adjustl(x_input_st))
        end if
        if (isotropic) then
           if (recoil) then
              write(outfile,'(a,a,a,a,a,a)'),trim(outputdir),'tests_hi_xin',trim(x_input_st),'_T',trim(temp_st),'_Isotropic_Recoil.dat'
           else
              write(outfile,'(a,a,a,a,a,a)'),trim(outputdir),'tests_hi_xin',trim(x_input_st),'_T',trim(temp_st),'_Isotropic_NoRecoil.dat'
           end if
        else
           if (recoil) then
              !write(outfile,'(a,a,a,a,a,a)'),trim(outputdir),'tests_hi_xin',trim(x_input_st),'_T',trim(temp_st),'_Dipolar_Recoil.dat'
              write(outfile,'(a,a,a,a,a,a)'),trim(outputdir),'tests_hi_xin',trim(x_input_st),'_T',trim(temp_st),'_Rayleigh_Recoil.dat'
           else
             ! write(outfile,'(a,a,a,a,a,a)'),trim(outputdir),'tests_hi_xin',trim(x_input_st),'_T',trim(temp_st),'_Dipolar_NoRecoil.dat'
              write(outfile,'(a,a,a,a,a,a)'),trim(outputdir),'tests_hi_xin',trim(x_input_st),'_T',trim(temp_st),'_Rayleigh_NoRecoil.dat'
           end if
        end if
        outfile = trim(outfile)//char(0)
        
        !! Perform nphotons realizations of one HI scattering and dump nu_ext and kout to file
        !! Then one can check for angular and freq. redistribution etc
        open(unit=11, file=trim(outfile), status='replace', form='formatted', action='write')
        write(11,'(i8,1x,es14.7,1x,es14.6,1x,es14.6)') nphotons, x_input, T, delta_nu_doppler
        do i=1,nphotons
           k(1) = 0.
           k(2) = 0.
           k(3) = 1.
           nu_cell = x_input * delta_nu_doppler + nu_0
           call scatter_HI(vcell, vth, nu_cell, k, nu_ext, iran)
           x_out  = (nu_cell - nu_0)/delta_nu_doppler
           write(11,'(es14.7,1x,es14.7,1x,es14.7,1x,es14.7)') x_out,  k(1),  k(2),  k(3)
        end do
        close(11)
        
     end do
     
  end if
     
     
     
  !----------------------------------
  !   test dust
  !----------------------------------

  if (species .eq. 'dust') then
     
     call read_dust_params(param_file)
     call print_dust_params
     
     vcell = 0.                     ! cm/s
     vth   = sqrt((2.0d0*kb/mp)*T)  ! cm/s
     delta_nu_doppler  = vth / lambda_0_cm 
        
     do j=1,n_xin
        x_input = x_in(j)
        
        if (x_input < 10.0) then
           !! Convert x_in into well-formatted string
           write(x_input_st,'(f6.2)') x_input
           write(x_input_st,'(a1)') trim(adjustl(x_input_st))
        else
           !! Convert x_in into well-formatted string
           write(x_input_st,'(f6.2)') x_input
           write(x_input_st,'(a2)') trim(adjustl(x_input_st))
         
        end if

          
        write(outfile,'(a,a,a,a,a)'),trim(outputdir),'tests_dust_xin',trim(x_input_st),'_HG41_alb032.dat'
        outfile = trim(outfile)//char(0)
        
        !! Perform nphotons realizations of one Dust scattering and dump nu_ext and kout to file
        !! Then one can check for angular and freq. redistribution etc
        open(unit=11, file=trim(outfile), status='replace', form='formatted', action='write')
        write(11,'(i8,1x,es14.7,1x,es14.6,1x,es14.6)') nphotons, x_input, T, delta_nu_doppler
        do i=1,nphotons
           k(1) = 0.
           k(2) = 0.
           k(3) = 1.
           ilost = -99
           !ilost = 1 => photon absorbed
           
           !! NB: first, we need to pass extra argument 'iescape' (or 'status') to scatter_dust function in module_dust_model.f90 to know if photon abs. or scattered
           nu_cell = x_input * delta_nu_doppler + nu_0
           call scatter_dust(vcell,nu_cell,k,nu_ext,iran,ilost)
           x_out  = (nu_cell - nu_0)/delta_nu_doppler
           write(11,'(es14.7,1x,es14.7,1x,es14.7,1x,es14.7,1x,i1)') x_out,  k(1),  k(2),  k(3), ilost
        end do
        close(11)
     end do
  end if

  !----------------------------------
  !   test Deuterium
  !----------------------------------

  if (species .eq. 'D') then
     
     call read_D_params(param_file)
     call print_D_params

     !! Convert T into well-formatted string
     write(temp_st,'(es14.0)') T
!!$     print*,T
!!$     print*,trim(temp_st)
     write(temp_st,'(a6)') trim(adjustl(temp_st))
     print*,temp_st
     
     vcell = 0.                     ! cm/s
     vth   = sqrt((2.0d0*kb/mp)*T) * sqrt_H2Deut_mass_ratio ! cm/s
     delta_nu_doppler  = vth / lambda_0_Deut_cm 
     
     do j=1,n_xin
        x_input = x_in(j)

        if (x_input < 10.0) then
           !! Convert x_in into well-formatted string
           write(x_input_st,'(f6.2)') x_input
           write(x_input_st,'(a1)') trim(adjustl(x_input_st))
        else
           !! Convert x_in into well-formatted string
           write(x_input_st,'(f6.2)') x_input
           write(x_input_st,'(a2)') trim(adjustl(x_input_st))
        end if
        if (isotropic) then
           if (recoil) then
              write(outfile,'(a,a,a,a,a,a)'),trim(outputdir),'tests_D_xin',trim(x_input_st),'_T',trim(temp_st),'_Isotropic_Recoil.dat'
           else
              write(outfile,'(a,a,a,a,a,a)'),trim(outputdir),'tests_D_xin',trim(x_input_st),'_T',trim(temp_st),'_Isotropic_NoRecoil.dat'
           end if
        else
           if (recoil) then
              write(outfile,'(a,a,a,a,a,a)'),trim(outputdir),'tests_D_xin',trim(x_input_st),'_T',trim(temp_st),'_Dipolar_Recoil.dat'
           else
              write(outfile,'(a,a,a,a,a,a)'),trim(outputdir),'tests_D_xin',trim(x_input_st),'_T',trim(temp_st),'_Dipolar_NoRecoil.dat'
           end if
        end if
        outfile = trim(outfile)//char(0)
        
        
        !! Perform nphotons realizations of one D scattering and dump nu_ext and kout to file
        !! Then one can check for angular and freq. redistribution etc
        open(unit=11, file=trim(outfile), status='replace', form='formatted', action='write')
        write(11,'(i8,1x,es14.7,1x,es14.6,1x,es14.6)') nphotons, x_input, T, delta_nu_doppler
        do i=1,nphotons
           k(1) = 0.
           k(2) = 0.
           k(3) = 1.
           nu_cell = x_input * delta_nu_doppler + nu_0_Deut
           call scatter_D(vcell,vth,nu_cell,k,nu_ext,iran)
           x_out  = (nu_cell - nu_0_Deut)/delta_nu_doppler
           write(11,'(es14.7,1x,es14.7,1x,es14.7,1x,es14.7)') x_out,  k(1),  k(2),  k(3)
        end do
        close(11)
        
     end do
  end if
  
  deallocate(x_in)
  
contains
  
  subroutine read_test_parameters(pfile)


    ! -----------------------------------------------------------------------------------------------------------------
    !! This subroutine reads the Lya_tests parameters first (species, T, x_in etc).
    !! The isotropic, recoil, dust etc parameters are read using the built-in RASCAS subroutines (read_HI_params etc)
    !! but they are also read at the end of this suroutine because we need to know their value here
    ! -----------------------------------------------------------------------------------------------------------------

    
    character(*),intent(in)   :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present
    
    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:24) == '[Lya_tests]') then
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
          case ('species')
             write(species,'(a)') trim(value)
          case ('T')
             read(value,*) T
          case ('n_xin')
             read(value,*) n_xin
             print*,n_xin
             allocate(x_in(n_xin))
          case ('x_in')
             read(value,*) x_in(:)
          case ('nphotons')
             read(value,*) nphotons
          case ('outputdir')
             write(outputdir,'(a)') trim(value)
!!$          case ('recoil')
!!$             read(value,*) recoil
!!$          case ('isotropic')
!!$             read(value,*) isotropic
!!$          case ('albedo')
!!$             read(value,*) albedo
!!$          case ('g_dust')
!!$             read(value,*) g_dust
          end select
       end do
    end if
    close(10)


    !! I now read again isotropic etc to know them in the test.f90 routine
    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:24) == '[HI]') then
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
          case ('recoil')
             read(value,*) recoil
          case ('isotropic')
             read(value,*) isotropic
          end select
       end do
    end if
    close(10)


    !! I now read again isotropic etc to know them in the test.f90 routine
    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:24) == '[D]') then
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
          case ('recoil')
             read(value,*) recoil
          case ('isotropic')
             read(value,*) isotropic
          end select
       end do
    end if
    close(10)

    !! I now read again isotropic etc to know them in the test.f90 routine
    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:24) == '[dust]') then
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
          case ('albedo')
             read(value,*) albedo
          case ('g_dust')
             read(value,*) g_dust
          end select
       end do
    end if
    close(10)
    
    
    return

  end subroutine read_test_parameters

  
end program tests_onescattering
  



