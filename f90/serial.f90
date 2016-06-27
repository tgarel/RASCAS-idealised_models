program main

  ! serial version of MCLya 

  use module_params
  use module_photon
  use module_mesh
  use module_domain
  use module_uparallel
  use module_constants

  implicit none

  type(photon_current),dimension(:),allocatable :: photgrid
  type(mesh)                                    :: meshdom
  type(domain)                                  :: compute_dom
  integer                                       :: nphot
  real(kind=8)                                  :: tau_sphere, R_cm, density,sigma_0, temp, start, tmptime, finish


  call cpu_time(start)

  ! read params file
  call read_params
  if(verbose) call print_params

  ! read ICs photons
  print *,'--> reading ICs photons in file: ',trim(file_ICs)
  call init_photons_from_file(file_ICs,photgrid)
  nphot = size(photgrid)
  print *,'--> Nphoton =',nphot


  print *,'--> reading domain and mesh...'
  ! Get domain properties
  call domain_constructor_from_file(file_compute_dom,compute_dom)

  if(ndomain/=1)then
     print *,'ndomain /= 1 use the MPI version'
     stop
  endif
  call mesh_from_file(mesh_file_list(1),meshdom)

  print *,'--> Ndomain =',ndomain
  if(verbose)then
     print *,'    |_ ',trim(file_compute_dom)
     print *,'    |_ ',trim(mesh_file_list(1))
  endif

#ifdef DEBUG
  print *,'--> check mesh dom'
  print *,meshdom%domain
  print *,meshdom%nCoarse,meshdom%nOct,meshdom%nLeaf,meshdom%nCell
  print *,minval(meshdom%xoct(:,:)),maxval(meshdom%xoct(:,:))
  print *,minval(meshdom%nbor(:,:)),maxval(meshdom%nbor(:,:))
  print *,minval(meshdom%octlevel(:)),maxval(meshdom%octlevel(:))
  print *,minval(meshdom%son(:)),maxval(meshdom%son(:))
  print *,minval(meshdom%father(:)),maxval(meshdom%father(:))
  !!!print *,'level of leaves =',minval(meshdom%octlevel(:), mask=(meshdom%son(:)<0)),maxval(meshdom%octlevel(:), mask=(meshdom%son(:)<0))
#endif

  ! Pour les tests de la sphere uniquement 
  box_size_cm = 1.d18
  R_cm        = 0.4d0 * box_size_cm
  sigma_0     = 5.88d-14 * (temp_fix/1.d4)**(-0.5) !cm^2 from Dijkstra14
  nhi_new     = tau0_fix/sigma_0/R_cm
  vth_new     = 12.9d0*sqrt(temp_fix/1.d4)*1.d5

  if(verbose)then
     print *,'--> box size in cm ',R_cm,box_size_cm,box_size_cm/mpc*1.e6
     print *,'--> tau =',tau0_fix
     print *,'--> density =',nhi_new
     print*,'--> precision stuff',precision(density),digits(density),epsilon(density)
  endif

  ! for idealized experiments
  if(overwritegas)then
     print*,'--> overwrite gas with nHI =',nhi_new
     print*,'--> overwrite gas with vth =',vth_new
     call overwrite_mesh(meshdom,nhi_new,vth_new)
  endif


  ! load uparallel table
#ifndef SWITCH_OFF_UPARALLEL
  print*,'--> loading uparallel tables...'
  call init_uparallel_tables
#endif

  call cpu_time(tmptime)
  print '(" --> Time = ",f12.3," seconds.")',tmptime-start



  print *,'--> starting RT...'
  ! do the RT stuff
  call MCRT(nphot,photgrid,meshdom,compute_dom)

  print *,'--> RT done'

  ! write results
  ! some check
  if(verbose)then
     print *,' '
     print *,'--> Some diagnostics...'
     print *,'min max status      =',minval(photgrid%status),maxval(photgrid%status)
     print *,'min max pos x       =',minval(photgrid%xcurr(1)),maxval(photgrid%xcurr(1))
     print *,'min max pos y       =',minval(photgrid%xcurr(2)),maxval(photgrid%xcurr(2))
     print *,'min max pos z       =',minval(photgrid%xcurr(3)),maxval(photgrid%xcurr(3))
     print *,'min max nb scatt    =',minval(photgrid%nb_abs),maxval(photgrid%nb_abs)
     print *,'min max nu          =',minval(photgrid%nu_ext),maxval(photgrid%nu_ext)
     print *,'min max lambda      =',clight/minval(photgrid%nu_ext)*cmtoA,clight/maxval(photgrid%nu_ext)*cmtoA

     print *,'min max travel time =',minval(photgrid%time),maxval(photgrid%time)
     print*, 'Sphere crossing time =',R_cm / clight

     print *,'Last scattering'
     print *,'min max pos x       =',minval(photgrid%xlast(1)),maxval(photgrid%xlast(1))
     print *,'min max pos y       =',minval(photgrid%xlast(2)),maxval(photgrid%xlast(2))
     print *,'min max pos z       =',minval(photgrid%xlast(3)),maxval(photgrid%xlast(3))
  endif

  ! check order for python reading
  !print *,'=============================='
  !i = 1
  !print *,i,photgrid(i)%xlast(1),photgrid(i)%xlast(2),photgrid(i)%xlast(3)
  !i = 2
  !print *,i,photgrid(i)%xlast(1),photgrid(i)%xlast(2),photgrid(i)%xlast(3)
  !i = 999
  !print *,i,photgrid(i)%xlast(1),photgrid(i)%xlast(2),photgrid(i)%xlast(3)

  print *,' '
  print*,'--> writing results in file: ',trim(fileout)
  call dump_photons(fileout,photgrid)

  call cpu_time(finish)
  print '(" --> Time = ",f12.3," seconds.")',finish-start

end program main
