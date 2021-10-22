module module_ssp_lib
  
  ! Module which contains declarations and routines/functions relative to the SSPs library

  use module_utils
  use module_constants
  
  private
  ! SSP/SED parameters (read from files)
  integer(kind=4)           :: nagebins
  real(kind=4),allocatable  :: agebins(:)        ! values of these ages (in Gyr in the code, yr in the file ...
  real(kind=4),allocatable  :: mid_agebins(:)    ! age boundaries for SED use
  integer(kind=4)           :: num_mets          ! number of metallicities produced by the SSP model
  real(kind=4),allocatable  :: metallicities(:)  ! array of model metallicities (in absolute unit), in log
  integer(kind=4)           :: nlambda           ! number of wavelengths
  real(kind=8),allocatable  :: lambda(:)         ! wavelength
  real(kind=8),allocatable  :: sed_list(:,:,:)   ! seds of SSPs (lambda,age,met)


  type, public :: SSPgrid
     real(kind=8),allocatable :: grid(:,:,:)
     real(kind=8),allocatable :: age(:),met(:),lambda(:)
     integer(kind=4)          :: nage,nmet,nlambda
  end type SSPgrid
  
  type(SSPgrid) :: fullLib
  logical       :: isLoaded = .false.
  real(kind=8)  :: Lsun_in_cgs = 3.826d33          ! erg/s using the same convention as in BC03
  
  public :: init_ssp_lib, ssp_lib_extract_subset, ssp_lib_interpolate, ssp_lib_integrate


contains

!*******************************************************************************************************

  subroutine init_ssp_lib(SSPdir,lambdamin,lambdamax)

    implicit none

    character(1000),intent(in)       :: SSPdir
    real(kind=8),intent(in),optional :: lambdamin, lambdamax
    
    ! read data
    call read_met_bins(SSPdir)
    call read_age_bins(SSPdir)
    call read_all_seds(SSPdir,lambdamin,lambdamax)

    print*,'reading SSPs...'
    print*,minval(agebins),maxval(agebins)
    print*,minval(metallicities),maxval(metallicities)
    print*,minval(lambda), maxval(lambda)
    print*, nagebins,num_mets,nlambda

    fullLib%nage = nagebins
    fullLib%nmet = num_mets
    fullLib%nlambda = nlambda    
    allocate(fullLib%age(nagebins)) ; fullLib%age(:) = agebins        ! convention x = ages
    allocate(fullLib%met(num_mets)) ; fullLib%met(:) = metallicities  ! convention y = metallicities
    allocate(fullLib%lambda(nlambda)) ; fullLib%lambda(:) = lambda    ! convention z = SEDs
    allocate(fullLib%grid(nagebins,num_mets,nlambda))
    fullLib%grid(:,:,:) =  reshape(source=sed_list,shape=(/nagebins,num_mets,nlambda/),order=(/3,1,2/)) ! convention z = SEDs

    isLoaded=.true.
    
    
    ! deallocate
    deallocate(metallicities, agebins, lambda, sed_list)
    
    
    ! spectra are in [Lsun / A / Msun]
    ! with Lsun = 3.826e33 erg/s and Msun = 2e33g
    ! lambda are in Angstrom
    ! age are in yr
  end subroutine init_ssp_lib
    
!*******************************************************************************************************
  subroutine ssp_lib_extract_subset(lmin,lmax,subLib)

    implicit none

    real(kind=8),intent(inout) :: lmin,lmax
    type(SSPgrid), intent(out) :: subLib
    integer(kind=4)            :: i0,i1,nl,i,j
    real(kind=8)               :: lambda0
    real(kind=8),allocatable   :: nu(:)

    if (lmin==lmax) then
       call locatedb(fullLib%lambda, fullLib%nlambda, lmin, i0)
       if( (lmin-fullLib%lambda(i0)) > (fullLib%lambda(i0+1)-lmin))then
          i0 = i0+1
       endif
       lambda0 = fullLib%lambda(i0)
       nl = 1
       i1 = i0
    else
       call locatedb(fullLib%lambda, fullLib%nlambda, lmin, i0)
       if( (lmin-fullLib%lambda(i0)) > (fullLib%lambda(i0+1)-lmin))then
          i0 = i0+1
       endif
       lmin = fullLib%lambda(i0)
       call locatedb(fullLib%lambda, fullLib%nlambda, lmax, i1)
       if( (lmax-fullLib%lambda(i1)) > (fullLib%lambda(i1+1)-lmax))then
          i1 = i1+1
       endif
       lmax = fullLib%lambda(i1)
       nl = i1-i0+1
    endif
    
    subLib%nage = fullLib%nage
    subLib%nmet = fullLib%nmet
    subLib%nlambda = nl
    allocate(subLib%age(subLib%nage))
    allocate(subLib%met(subLib%nmet))
    allocate(subLib%lambda(subLib%nlambda))
    allocate(subLib%grid(subLib%nage,subLib%nmet,subLib%nlambda))
    
    subLib%age = fullLib%age
    subLib%met = fullLib%met
    subLib%lambda = fullLib%lambda(i0:i1)
    subLib%grid(:,:,1:nl) = fullLib%grid(:,:,i0:i1) * Lsun_in_cgs   ! in erg/s/A/Msun
    allocate(nu(nl))
    nu(:) = clight / (subLib%lambda(:)*1e-8) ! [Hz]
    do i=1,subLib%nage
       do j=1,subLib%nmet
          subLib%grid(i,j,:) = subLib%grid(i,j,:) / planck / nu(:)  ! nb of photons / s / A / Msun 
       enddo
    enddo
    deallocate(nu)
    
    return
  end subroutine ssp_lib_extract_subset

!*******************************************************************************************************
  subroutine ssp_lib_interpolate(lib, x, y, res)

    implicit none
    type(SSPgrid),intent(in)              :: lib
    real(kind=8), intent(in)              :: x,y
    real(kind=8),dimension(:),intent(out) :: res
    integer(kind=4)                       :: j1,j2,ny,i1,i2,nx
    real(kind=8)                          :: wy1,wy2,norm,wx1,wx2
    
    ny = size(lib%met)
    ! locate bins and compute weights
    if (y <= lib%met(1)) then 
       j1 = 1
       wy1  = 1.0d0
       j2 = 2
       wy2  = 0.0d0
    else if (y >= lib%met(ny)) then 
       j1  = ny-1
       wy1 = 0.0d0
       j2  = ny
       wy2 = 1.0d0
    else
       do j2 = 1,ny
          if (y <= lib%met(j2)) exit
       end do
       j1 = j2 - 1
       wy1  = lib%met(j2) - y 
       wy2  = y - lib%met(j1)
       norm = wy1 + wy2
       wy1  = wy1 / norm
       wy2  = wy2 / norm
    end if
    

    nx = size(lib%age)
    ! locate bins and compute weights
    if (x <= lib%age(1)) then 
       i1 = 1
       wx1 = 1.0d0
       i2 = 2
       wx2 = 0.0d0
    else if (x >= lib%age(nx)) then 
       i1  = nx-1
       wx1 = 0.0d0
       i2  = nx
       wx2 = 1.0d0
    else
       do i2 = 1,nx
          if (x <= lib%age(i2)) exit
       end do
       i1 = i2 - 1
       wx1  = lib%age(i2) - x
       wx2  = x - lib%age(i1)
       norm = wx1 + wx2
       wx1  = wx1 / norm
       wx2  = wx2 / norm
    end if
    
    res(:) =  wx2 * wy2 * lib%grid(i2, j2, :)  &
            + wx1 * wy2 * lib%grid(i1, j2, :)  &
            + wx2 * wy1 * lib%grid(i2, j1, :)  &
            + wx1 * wy1 * lib%grid(i1, j1, :)
    return
  end subroutine ssp_lib_interpolate
    
!*******************************************************************************************************
  subroutine ssp_lib_integrate(x,y,n,integral)
    real(kind=8),intent(in)  :: x(n), y(n)
    integer,intent(in)       :: n
    real(kind=8),intent(out) :: integral
    integral = trapz1(x,y,n)
    return
  end subroutine ssp_lib_integrate
  
!*******************************************************************************************************
  !*************************************************************************
  FUNCTION trapz1(X,Y,N,cum)
    
    ! Integrates function Y(X) along the whole interval 1..N, using a very 
    ! simple staircase method and returns the result.
    ! Optionally, the culumative integral is returned in the cum argument.
    !-------------------------------------------------------------------------
    integer :: N,i
    real(kind=8):: trapz1
    real(kind=8):: X(N),Y(N)
    real(kind=8),optional::cum(N)
    real(kind=8),allocatable::cumInt(:)
    !-------------------------------------------------------------------------
    allocate(cumInt(N))
    cumInt(:)=0.d0
    trapz1=0.d0
    if (N.le.1) RETURN
    do i=2,N
       cumInt(i)= cumInt(i-1) + abs(X(i)-X(i-1)) * (Y(i)+Y(i-1)) / 2.d0
    end do
    trapz1 = cumInt(N)
    if(present(cum)) cum=cumInt
    deallocate(cumInt)
    return    
  END FUNCTION trapz1
!*************************************************************************

  subroutine read_met_bins(SSPdir)
    
    implicit none
    
    character(1000),intent(in) :: SSPdir
    character(1000) :: metsbin_file
    integer(kind=4) :: imet

    write(metsbin_file,'(a,a)') trim(SSPdir),"/metallicity_bins.dat"
    call test_stop(metsbin_file)
    open(unit=12,file=metsbin_file,status='old',form='formatted')
    read(12,'(i8)') num_mets
    allocate(metallicities(num_mets))
    do imet = 1,num_mets
       read(12,'(e14.6)') metallicities(imet)
    end do
    close(12)

    ! take log to make interpolations in log
    metallicities  = log10(metallicities)

    return

  end subroutine read_met_bins

!*******************************************************************************************************

  subroutine read_all_seds(SSPdir,lambdamin,lambdamax)
    
    implicit none
    
    character(1000),intent(in) :: SSPdir
    real(kind=8),intent(in),optional :: lambdamin,lambdamax
    character(1000)           :: sed_file
    integer(kind=4)           :: dum,iage,imet
    integer(kind=4)           :: ilambdamin,ilambdamax
    real(kind=8),allocatable  :: sed(:)

    write(sed_file,'(a,a)') trim(SSPdir),"/all_seds.dat"
    call test_stop(sed_file)
    open(unit=12,file=sed_file,status='old',form='unformatted')
    read(12) nlambda,dum
    allocate(sed(nlambda))
    read(12) sed(:)

    if(present(lambdamin).and.present(lambdamax))then
       ! extract only desired wavelength range 
       call locatedb(sed,nlambda,lambdamin-10.,ilambdamin)  
       call locatedb(sed,nlambda,lambdamax+10.,ilambdamax)  
       nlambda = ilambdamax - ilambdamin + 1
       allocate(lambda(nlambda))
       lambda = sed(ilambdamin:ilambdamax)
    else
       allocate(lambda(nlambda))
       lambda(:) = sed(:)
       ilambdamin = 1
       ilambdamax = nlambda
    endif
    allocate(sed_list(nlambda,nagebins,num_mets))
    
    do imet = 1,num_mets
       do iage = 1,nagebins
          read(12) sed
          sed_list(:,iage,imet) = sed(ilambdamin:ilambdamax)
       end do
    end do
    close(12)

    deallocate(sed)
    
    return

  end subroutine read_all_seds

!*******************************************************************************************************

  subroutine read_age_bins(SSPdir)
    
    implicit none
    
    character(1000),intent(in) :: SSPdir
    integer(kind=4) :: iage
    character(1000) :: filename
    
    write(filename,'(a,a)') trim(SSPdir),"/age_bins.dat"
    call test_stop(filename)
    open(unit=12,file=filename,status='old',form='formatted')
    read(12,'(i8)') nagebins
    allocate(agebins(nagebins))
    do iage = 1,nagebins
       read(12,'(e14.6)') agebins(iage)
       agebins(iage) = agebins(iage) * 1.e-9  ! convert from yr to Gyr
    end do
    close(12)
    
    ! define age boundaries where to use each SSP-SED
    allocate(mid_agebins(nagebins+1)) 
    mid_agebins(1)    = agebins(1) ! == 0 yr
    do iage = 2,nagebins
       mid_agebins(iage) = 0.5*(agebins(iage)+agebins(iage-1))
    end do
    mid_agebins(nagebins+1) = agebins(nagebins) ! == 20 Gyr

    return

  end subroutine read_age_bins

!*******************************************************************************************************

  subroutine test_stop(file)

    implicit none

    character(*) :: file
    logical      :: exist

    inquire(file=file,exist=exist)
    if (.not.exist) then
       write(*,*) '> Error: ',trim(file),' does not exist.'  
       stop
    endif

    return

  end subroutine test_stop

!*******************************************************************************************************

end module module_ssp_lib
