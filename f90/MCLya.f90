program main

  use module_parallel_mpi
  use module_master
  use module_worker
  use module_params
  use module_uparallel 

  implicit none

  real(kind=8) :: start,finish,intermed
  real(kind=8) :: tau_sphere, R_cm, density, sigma_0, temperature

  call cpu_time(start)

  call initialisation_mpi
  call type_derive

  nslave=nb_cpus-1

  if (rank == 0) then
     if(verbose) print*,'--> Nworker =',nslave
  end if

  !!call init_params
  ! or get run params
  call read_params
  if(rank==1) call print_params

  call MPI_BARRIER(MPI_COMM_WORLD,code)

  !! some settings
  box_size_cm = 1.d18
  R_cm        = 0.4d0 * box_size_cm
  sigma_0     = 5.88d-14 * (temp_fix/1.d4)**(-0.5) !cm^2 from Dijkstra14
  nhi_new     = tau0_fix/sigma_0/R_cm
  vth_new     = 12.9d0*sqrt(temp_fix/1.d4)*1.d5
  
  ! load uparallel table
#ifndef SWITCH_OFF_UPARALLEL
  if(rank==0) print*,'--> loading uparallel tables...'
  call init_uparallel_tables
  call cpu_time(intermed)
  if(rank==0) print '(" --> time to compute uparallel tables = ",f12.3," seconds.")',intermed-start
#endif

  ! Master - Worker separation
  if (rank == 0) then
     ! Master section, will dispatch the jobs.
     call master
  else
     ! Worker section, will mostly do radiative transfer (MCLya)
     call worker
  end if

  ! write results


  call finalisation_mpi
  call cpu_time(finish)
  if(rank==0)then
     print*,'--> work done, MPI finalized'
     print '(" --> Time = ",f12.3," seconds.")',finish-start
  endif

end program main
