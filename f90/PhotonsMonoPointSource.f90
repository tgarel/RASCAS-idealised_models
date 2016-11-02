program PhotonsMonoPointSource

  use module_photon
  !!use module_HI_model
  use module_random
  use module_constants

  implicit none

  type(photon_init),dimension(:),allocatable :: photgrid
  integer                                    :: n_photon, iseed, iran, i
  real(kind=8),dimension(3)                  :: xsource
  real(kind=8)                               :: nu_source, k1, k2, k3
  character(2000)                            :: file
  real(kind=8)                               :: lambda_0, lambda_0_cm, nu_0

  ! should improve since this is already defined in module_HI_model...
  lambda_0    = 1215.67d0
  lambda_0_cm = lambda_0 / cmtoA
  nu_0        = clight/lambda_0_cm

  !----parameters-----------------
  xsource   = (/0.5d0,0.5d0,0.5d0/)
  n_photon  = 5000000
  nu_source = nu_0
  iseed     = 1234
  file      = 'ICs_photons_n5e6.dat'
  !-------------------------------

  allocate(photgrid(n_photon))
  iran = -iseed

  do i=1,n_photon
     photgrid(i)%ID    = i
     photgrid(i)%nu_em = nu_0
     photgrid(i)%x_em  = xsource
     photgrid(i)%iran  = -1*i
     !print*,iran
     call emission_source(iran,k1,k2,k3)
     photgrid(i)%k_em  = (/k1,k2,k3/)
     !print*,i,iran,k1,k2,k3
  enddo

  ! write ICs
  open(unit=14, file=trim(file), status='unknown', form='unformatted', action='write')
  write(14) n_photon
  write(14) iseed
  write(14) (photgrid(i)%ID,i=1,n_photon)
  write(14) (photgrid(i)%nu_em,i=1,n_photon)
  write(14) (photgrid(i)%x_em(:),i=1,n_photon)
  write(14) (photgrid(i)%k_em(:),i=1,n_photon)
  write(14) (photgrid(i)%iran,i=1,n_photon)
  close(14)

  deallocate(photgrid)

contains


  !================================================================
  SUBROUTINE emission_source(iran,a0,b0,c0)
    !--------------------------------------------------------------------
    !
    !  DETERMINES THE *INITIAL* EMISSION DIRECTION OF PHOTONS (FROM THE SOURCE) isotropically
    !
    implicit none

    integer,intent(inout)    :: iran
    real(kind=8),intent(out) :: a0,b0,c0
    real(kind=8)             ::rtheta,rphi,st,phi,theta     ! 0<rtheta,rphi<1, for ran3 

    ! ........generating the emission direction
    rphi   = ran3(iran)
    rtheta = ran3(iran)     

    phi   = 2.d0*pi*rphi
    theta = acos(1d0-2d0*rtheta)  ! isotropic emission

    !........director cosines
    st = sin(theta)
    a0 = st*cos(phi)   !x
    b0 = st*sin(phi)   !y
    c0 = cos(theta)    !z

  END SUBROUTINE emission_source
  !================================================================

end program PhotonsMonoPointSource
