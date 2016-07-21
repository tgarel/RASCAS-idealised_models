module module_params

  ! parameters for MCLya (will be at some stage read in a namelist/parameterfile)

  public

  ! MPI - load-balancing parameters
  integer                                  :: nbuffer
  logical                                  :: verbose

  ! Repository/Data
  character(2000)                          :: datadir,file_ICs,fileout

  ! Domain decomposition
  integer                                  :: ndomain
  character(2000)                          :: file_compute_dom,root_domain_file
  character(2000),dimension(:),allocatable :: domain_file_list, mesh_file_list 

  ! Special purpose parameters - Idealized cases
  logical                                  :: overwritegas
  real(kind=8)                             :: temp_fix, tau0_fix
  real(kind=8)                             :: nhi_new, vth_new

  ! Atomic physics
  logical                                  :: recoil ! if true, compute recoil effect for scatterings on H and D... 

  ! Deuterium parameters 
  real(kind=8)                             :: deut2H_nb_ratio ! nb ratio of Deuterium to Hydrogen atoms 

  ! Dust parameters for the model of Verhamme+2012 
  real(kind=8)                             :: dust_to_metal_ratio ! mass ratio of dust to metals (= 0.3, Inoue 2003)
  real(kind=8)                             :: mH_over_mdust       ! mH / mdust (= 5e-8, Draine & Lee 1984)

  
contains

  subroutine read_params
    
    implicit none 

    character(2000) :: line,name,value,infile
    integer         :: i,j,narg,err
    character(2)    :: nchar

    ! define default values for all parameters 
    nbuffer             = 1000
    verbose             = .true.
    ndomain             = 1
    datadir             = '../data/'
    file_compute_dom    = 'compute_domain.dom'
    root_domain_file    = 'domain_'
    file_ICs            = 'ICs_photons_n1e6.dat'
    fileout             = 'photons_done.dat'
    overwritegas        = .true.
    temp_fix            = 1.d4
    tau0_fix            = 1.d5
    recoil              = .true.
    deut2H_nb_ratio     = 3.d-5
    dust_to_metal_ratio = 0.3d0
    mH_over_mdust       = 5e-8
    
    ! Read namelist filename from command line argument
    narg = command_argument_count()
    IF(narg .LT. 1)THEN
       write(*,*)'You should type: MCLya params.dat [restart]'
       write(*,*)'File params.dat should contain a parameter namelist'
       write(*,*)'restart is optional'
       stop
       !!call clean_stop
    END IF
    call get_command_argument(1, infile)

!#ifndef WITHOUTMPI
!  call MPI_BCAST(infile,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
!#endif

    ! read parameters from params.dat file
    open(unit=10,file=trim(infile),status='old',form='formatted')
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       i = scan(line,'=')
       if (i==0 .or. line(1:1)=='#' .or. line(1:1)=='!') cycle  ! skip blank or commented lines
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))
       i = scan(value,'!')
       if (i /= 0) value = trim(adjustl(value(:i-1)))
       !write(*,'(a,a15,a3,a)') '> ',trim(name),' : ',trim(value)
       select case (trim(name))
       case ('verbose')
          read(value,*) verbose
       case ('datadir')
          datadir = trim(value)
       case ('nbuffer')
          read(value,*) nbuffer
       case ('ndomain')
          read(value,*) ndomain
       case ('file_compute_dom')
          file_compute_dom=trim(value)
       case ('root_domain_file')
          root_domain_file=trim(value)
       case ('file_ICs')
          file_ICs=trim(value)   
       case ('fileout')
          fileout=trim(value)
       case ('overwritegas')
          read(value,*) overwritegas
       case ('temp_fix')
          read(value,*) temp_fix
       case ('tau0_fix')
          read(value,*) tau0_fix
       case ('recoil')
          read(value,*) recoil
       case ('deut2H_nb_ratio')
          read(value,*) deut2H_nb_ratio
       case ('dust_to_metal_ratio')
          read(value,*) dust_to_metal_ratio
       case ('mH_over_mdust')
          read(value,*) mH_over_mdust
       end select
    end do
    close(10)

    file_compute_dom    = trim(datadir)//trim(file_compute_dom)
    root_domain_file    = trim(datadir)//trim(root_domain_file)
    file_ICs            = trim(datadir)//trim(file_ICs)

    if(ndomain>=1)then
       allocate(domain_file_list(ndomain),mesh_file_list(ndomain))
       do j=1,ndomain
          write(nchar,'(i2.2)') j
          domain_file_list(j) = trim(root_domain_file)//nchar//'.dom'
          mesh_file_list(j)   = trim(root_domain_file)//nchar//'.mesh'
       enddo
    else
       print*,"Cannot run without mesh domain..."
       stop
    endif
    
    return

  end subroutine read_params



  subroutine print_params
    ! print params for log

    integer :: i

    write(*,*) '--> Parameters:'
    write(*,'(a,L1)')    '                verbose : ',verbose
    write(*,'(a,i2)')    '                ndomain : ',ndomain
    write(*,'(a,i7)')    '                nbuffer : ',nbuffer
    write(*,'(a,a)')     '         data directory : ',trim(datadir)
    write(*,'(a,a)')     '   computational domain : ',trim(file_compute_dom)
    do i=1,ndomain
       write(*,'(a,a)')  '            mesh domain : ',trim(domain_file_list(i))
       write(*,'(a,a)')  '                        : ',trim(mesh_file_list(i))
    enddo
    write(*,'(a,a)')     '         IC source file : ',trim(file_ICs)
    write(*,'(a,a)')     '           Results file : ',trim(fileout)

    write(*,'(a,L1)')    '    overwrite gas props : ',overwritegas
    write(*,'(a,ES9.3)') '        temperature fix : ',temp_fix
    write(*,'(a,ES9.3)') '               tau0 fix : ',tau0_fix

    write(*,'(a,L1)')    '                 recoil : ',recoil
    
    write(*,'(a,ES9.3)') '        deut2H_nb_ratio : ',deut2H_nb_ratio
    write(*,'(a,ES9.3)') '    dust_to_metal_ratio : ',dust_to_metal_ratio
    write(*,'(a,ES9.3)') '          mH_over_mdust : ',mH_over_mdust

  end subroutine print_params

end module module_params

