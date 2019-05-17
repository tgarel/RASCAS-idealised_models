program LyaPhotonsFromGas

  use module_ramses
  use module_domain
  use module_random
  use module_constants
  use module_utils
  use module_HI_1216_model, only: lambda_0, lambda_0_cm, nu_0
  
  implicit none

  character(2000)              :: parameter_file
  type(domain)                 :: emission_domain
  integer(kind=4)              :: nleaftot,nvar,narg,nsel,i,iphot,j,iseed,lmax
  integer(kind=4), allocatable :: emitting_cells(:),leaf_level(:)
  real(kind=8),allocatable     :: x_leaf(:,:),ramses_var(:,:),recomb_em(:),coll_em(:),v_leaf(:,:),HIDopWidth(:),cell_volume_vs_level(:)
  real(kind=8),allocatable     :: coolingTime(:)
  real(kind=8)                 :: r1, r2, dx, dv, nu, scalar, recomb_total,coll_total,k(3), boxsize,maxrec,maxcol
  real(kind=8)                 :: start_photpacket,end_photpacket,x(3),dt,xmin,xmax,ymin,ymax,zmin,zmax
  logical                      :: ok
  real(kind=8),allocatable     :: nu_em(:),x_em(:,:),k_em(:,:),nu_cell(:)
  integer(kind=4),allocatable  :: cpu_list(:)
  integer(kind=4)              :: ncpu_read
  real(kind=8),allocatable     :: low_prob_rec(:),low_prob_col(:)
  integer(kind=4)              :: ilow,iup,imid
  real(kind=8)                 :: mid

  ! ---------------------------------------------------------------------------
  ! user-defined parameters - read from section [LyaPhotonsFromGas] of the parameter file
  ! ---------------------------------------------------------------------------
  ! --- input / outputs
  character(2000)           :: outputfileRec = 'LyaPhotIC.recLya' ! file to which recombination photons will be written
  character(2000)           :: outputfileCol = 'LyaPhotIC.colLya' ! file to which collisional photons will be written
  character(2000)           :: repository = './'      ! ramses run directory (where all output_xxxxx dirs are).
  integer(kind=4)           :: snapnum = 1            ! ramses output number to use
  ! --- emission domain  
  character(10)             :: emission_dom_type      = 'sphere'         ! shape type of domain  // default is a shpere.
  real(kind=8),dimension(3) :: emission_dom_pos       = (/0.5,0.5,0.5/)  ! center of domain [code units]
  real(kind=8)              :: emission_dom_rsp       = 0.3              ! radius of sphere [code units]
  real(kind=8)              :: emission_dom_size      = 0.3              ! size of cube [code units]
  real(kind=8)              :: emission_dom_rin       = 0.0              ! inner radius of shell [code units]
  real(kind=8)              :: emission_dom_rout      = 0.3              ! outer radius of shell [code units]
  real(kind=8)              :: emission_dom_thickness = 0.1              ! thickness of slab [code units]
  ! --- parameters
  integer(kind=4)           :: nphotons = 10000
  integer(kind=4)           :: ranseed  = -100
  logical                   :: doRecombs = .false.  ! set to true to sample emission from recombinations
  logical                   :: doColls   = .false.  ! set to true to sample emission from coll. excitations.
  real(kind=8)              :: tcool_resolution = 3.0d0 ! ignore collisional emission from cells with Dt > tcool / tcool_resolution 
  ! --- miscelaneous
  logical                   :: verbose = .true.
  ! ---------------------------------------------------------------------------


  ! -------------------- read parameters --------------------------------------
  narg = command_argument_count()
  if(narg .lt. 1)then
     write(*,*)'You should type: LyaPhotonsFromGas params.dat'
     write(*,*)'File params.dat should contain a parameter namelist'
     stop
  end if
  call get_command_argument(1, parameter_file)
  call read_LyaPhotonsFromGas_params(parameter_file)
  if (verbose) call print_LyaPhotonsFromGas_params
  ! ----------------------------------------------------------------------------

  if (.not. doRecombs .and. .not. doColls) then
     print*,'Nothing to do ... '
     stop
  end if

  ! ------- Define the emission domain (region from which cells emit). ---------
  select case(emission_dom_type)
  case('sphere')
     call domain_constructor_from_scratch(emission_domain,emission_dom_type, &
          xc=emission_dom_pos(1),yc=emission_dom_pos(2),zc=emission_dom_pos(3),r=emission_dom_rsp)
  case('shell')
     call domain_constructor_from_scratch(emission_domain,emission_dom_type, &
          xc=emission_dom_pos(1),yc=emission_dom_pos(2),zc=emission_dom_pos(3),r_inbound=emission_dom_rin,r_outbound=emission_dom_rout)
  case('cube')
     call domain_constructor_from_scratch(emission_domain,emission_dom_type, & 
          xc=emission_dom_pos(1),yc=emission_dom_pos(2),zc=emission_dom_pos(3),size=emission_dom_size)
  case('slab')
     call domain_constructor_from_scratch(emission_domain,emission_dom_type, &
          xc=emission_dom_pos(1),yc=emission_dom_pos(2),zc=emission_dom_pos(3),thickness=emission_dom_thickness)
  end select
  ! ----------------------------------------------------------------------------
  
  
  ! ---- Read/select leaf cells and compute their luminositites ----------------
  if (verbose) print*,'Reading cells ... '
  ! define max extent of emission domain 
  call domain_get_bounding_box(emission_domain,xmin,xmax,ymin,ymax,zmin,zmax)
  call get_cpu_list_periodic(repository, snapnum, xmin,xmax,ymin,ymax,zmin,zmax, ncpu_read, cpu_list)
  call read_leaf_cells_omp(repository, snapnum, ncpu_read, cpu_list, nleaftot, nvar, x_leaf, ramses_var, leaf_level)
  !call read_leaf_cells(repository, snapnum, nleaftot, nvar, x_leaf, ramses_var, leaf_level)
  if (verbose) print*,'done reading'
  call select_cells_in_domain(emission_domain,nleaftot,x_leaf,leaf_level,emitting_cells)
  nsel = size(emitting_cells)
  
  if (verbose) print*,'done selecting cells in domain'
  allocate(recomb_em(nsel),coll_em(nsel),HIDopWidth(nsel))
  call ramses_get_LyaEmiss_HIDopwidth(repository,snapnum,nleaftot,nvar,ramses_var,recomb_em,coll_em,HIDopWidth,sample=emitting_cells)
  if (verbose) print*,'done computing emissivities'
  allocate(coolingTime(nsel))
  call ramses_get_cooling_time(repository,snapnum,nleaftot,nvar,ramses_var,coolingTime,sample=emitting_cells)
  if (verbose) print*,'done computing cooling times (min/max=',minval(coolingTime)/(365.25d0*24.*3600.),maxval(coolingTime)/(365.25d0*24.*3600.),' yr)'
  lmax = 30
  allocate(cell_volume_vs_level(lmax))
  boxsize = ramses_get_box_size_cm(repository,snapnum)
  do i=1,lmax
     cell_volume_vs_level(i) = (boxsize * 0.5**i)**3 ! cm3
  end do
  maxrec = -1.
  maxcol = -1. 
  do i=1,nsel
     j = emitting_cells(i)
     dv = cell_volume_vs_level(leaf_level(j))
     recomb_em(i) = recomb_em(i) * dv   ! erg/s 
     if (recomb_em(i) > maxrec) maxrec = recomb_em(i)
     ! if cooling time is not resolved, set collisional emission to zero
     dt = (boxsize*0.5**leaf_level(j))/0.0125d0/clight/3.d0  ! courant condition [s]
     if (coolingTime(i) > tcool_resolution * dt) then 
        coll_em(i)   = coll_em(i) * dv
     else
        coll_em(i) = 0.0d0
     end if
     if (coll_em(i) > maxcol) maxcol = coll_em(i) 
  end do
  
  recomb_total = sum(recomb_em) / (planck*nu_0)  ! nb of photons per second
  coll_total   = sum(coll_em) / (planck*nu_0)  ! nb of photons per second
  
  print*,'coll_total,recomb_total [erg/s] = ',coll_total*(planck*nu_0),recomb_total*(planck*nu_0)
  
  recomb_em = recomb_em / maxrec
  coll_em   = coll_em / maxcol
  
  ! NewScheme -
  if (doRecombs) then 
     allocate(low_prob_rec(nsel+1))
     low_prob_rec(1) = 0.0d0
     do i = 2,nsel
        low_prob_rec(i) = low_prob_rec(i-1) + recomb_em(i-1) 
     end do
     low_prob_rec = low_prob_rec / (low_prob_rec(nsel)+recomb_em(nsel))
     low_prob_rec(nsel+1) = 1.1  ! higher than upper limit 
  end if
  if (doColls) then 
     allocate(low_prob_col(nsel+1))
     low_prob_col(1) = 0.0d0
     do i = 2,nsel
        low_prob_col(i) = low_prob_col(i-1) + coll_em(i-1) 
     end do
     low_prob_col = low_prob_col / (low_prob_col(nsel)+coll_em(nsel)) 
     low_prob_col(nsel+1) = 1.1  ! higher than upper limit 
  end if
  ! - NewScheme
  
  allocate(v_leaf(3,nleaftot))
  call ramses_get_velocity_cgs(repository,snapnum,nleaftot,nvar,ramses_var,v_leaf)
  deallocate(ramses_var)
  ! ----------------------------------------------------------------------------
  
  allocate(nu_em(nphotons),x_em(3,nphotons),k_em(3,nphotons),nu_cell(nphotons))
  
  ! --------------------------------------------------------------------------------------
  if (doRecombs) then 
     if (verbose) then
        write(*,*) "Starting to sample recombination emissivity" 
        call cpu_time(start_photpacket)
     end if
     ! --------------------------------------------------------------------------------------
     
     ! NewScheme -
     iseed = ranseed
     do iphot = 1,nphotons
        ok = .false.
        do while (.not. ok) 
           r1 = ran3(iseed)
           ! binary search
           iup = nsel
           ilow = 1
           do while (iup - ilow > 1)
              imid = (iup+ilow)/2
              mid  = low_prob_rec(imid)
              if (r1 >= mid) then 
                 ilow = imid
              else
                 iup = imid
              end if
           end do
           ! check
           if (.not. (r1 >= low_prob_rec(ilow) .and. r1 < low_prob_rec(iup) )) then
              print*,'hi Harley ;) '
           end if
           ! draw photon's ICs from cell ilow
           j  = emitting_cells(ilow)
           dx = 1.d0 / 2**leaf_level(j)
           ! draw photon position in cell
           x(1)  = (ran3(iseed)-0.5d0) * dx + x_leaf(j,1)
           x(2)  = (ran3(iseed)-0.5d0) * dx + x_leaf(j,2)
           x(3)  = (ran3(iseed)-0.5d0) * dx + x_leaf(j,3)
           ! check that photon is inside emission domain and if so keep it. 
           ! (it may happen that a cell is only partly included ...)
           if (domain_contains_point(x,emission_domain)) then 
              x_em(:,iphot)  = x(:)
              ! draw propagation direction
              call isotropic_direction(k,iseed)
              k_em(:,iphot) = k
              ! compute frequency in cell frame
              r1 = ran3(iseed)
              r2 = ran3(iseed)
              nu = sqrt(-2.*log(r1)) * cos(2.0d0*pi*r2)
              nu_cell(iphot) = (HIDopWidth(ilow) * nu_0 / clight) * nu + nu_0
              ! compute frequency in exteral frame 
              scalar = k(1)*v_leaf(1,j) + k(2)*v_leaf(2,j) + k(3)*v_leaf(3,j)
              nu_em(iphot)  = nu_cell(iphot) / (1d0 - scalar/clight)
              ok = .true.
           end if
        end do
     end do
     ! - NewScheme
     ! --------------------------------------------------------------------------------------
     ! write ICs
     ! --------------------------------------------------------------------------------------
     open(unit=14, file=outputfileRec, status='unknown', form='unformatted', action='write')
     write(14) nphotons      ! nb of MC photons 
     write(14) recomb_total  ! nb of real photons (per sec).
     write(14) ranseed
     write(14) (i,i=1,nphotons) ! ID
     write(14) (nu_em(i),i=1,nphotons)
     write(14) (x_em(:,i),i=1,nphotons)
     write(14) (k_em(:,i),i=1,nphotons)
     write(14) (-i,i=1,nphotons) ! seeds
     write(14) (nu_cell(i),i=1,nphotons)
     close(14)
     ! --------------------------------------------------------------------------------------
     if (verbose) then 
        call cpu_time(end_photpacket)
        print*, 'time to draw recombinations = ',end_photpacket-start_photpacket,' seconds.'
     end if
  end if
  ! --------------------------------------------------------------------------------------
  
  
  
  ! --------------------------------------------------------------------------------------
  if (doColls) then 
     if (verbose) then
        write(*,*) "Starting to sample collisional emissivity" 
        call cpu_time(start_photpacket)
     end if
     ! NewScheme -
     iseed = ranseed
     do iphot = 1,nphotons
        ok = .false.
        do while (.not. ok) 
           r1 = ran3(iseed)
           ! binary search
           iup = nsel
           ilow = 1
           do while (iup - ilow > 1)
              imid = (iup+ilow)/2
              mid  = low_prob_col(imid)
              if (r1 >= mid) then 
                 ilow = imid
              else
                 iup = imid
              end if
           end do
           ! check
           if (.not. (r1 >= low_prob_col(ilow) .and. r1 < low_prob_col(iup) )) then
              print*,'hi Harley ;) '
           end if
           ! draw photon's ICs from cell ilow
           j  = emitting_cells(ilow)
           dx = 1.d0 / 2**leaf_level(j)
           ! draw photon position in cell
           x(1)  = (ran3(iseed)-0.5d0) * dx + x_leaf(j,1)
           x(2)  = (ran3(iseed)-0.5d0) * dx + x_leaf(j,2)
           x(3)  = (ran3(iseed)-0.5d0) * dx + x_leaf(j,3)
           ! check that photon is inside emission domain and if so keep it. 
           ! (it may happen that a cell is only partly included ...)
           if (domain_contains_point(x,emission_domain)) then 
              x_em(:,iphot)  = x(:)
              ! draw propagation direction
              call isotropic_direction(k,iseed)
              k_em(:,iphot) = k
              ! compute frequency in cell frame
              r1 = ran3(iseed)
              r2 = ran3(iseed)
              nu = sqrt(-2.*log(r1)) * cos(2.0d0*pi*r2)
              nu_cell(iphot) = (HIDopWidth(ilow) * nu_0 / clight) * nu + nu_0
              ! compute frequency in exteral frame 
              scalar = k(1)*v_leaf(1,j) + k(2)*v_leaf(2,j) + k(3)*v_leaf(3,j)
              nu_em(iphot)  = nu_cell(iphot) / (1d0 - scalar/clight)
              ok = .true.
           end if
        end do
     end do
     ! - NewScheme
     ! --------------------------------------------------------------------------------------
     ! write ICs
     ! --------------------------------------------------------------------------------------
     open(unit=14, file=outputfileCol, status='unknown', form='unformatted', action='write')
     write(14) nphotons      ! nb of MC photons 
     write(14) coll_total    ! nb of real photons (per sec).
     write(14) ranseed
     write(14) (i,i=1,nphotons) ! ID
     write(14) (nu_em(i),i=1,nphotons)
     write(14) (x_em(:,i),i=1,nphotons)
     write(14) (k_em(:,i),i=1,nphotons)
     write(14) (-i,i=1,nphotons) ! seeds
     write(14) (nu_cell(i),i=1,nphotons)
     close(14)
     ! --------------------------------------------------------------------------------------
     if (verbose) then 
        call cpu_time(end_photpacket)
        print*, 'time to draw collisions = ',end_photpacket-start_photpacket,' seconds.'
        print*,' '
     end if
  end if
  ! --------------------------------------------------------------------------------------
  
  deallocate(nu_em,x_em,k_em,nu_cell)
  
contains


    subroutine read_LyaPhotonsFromGas_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! subroutine which reads parameters of current module in the parameter file pfile
    ! default parameter values are set at declaration (head of module)
    !
    ! ALSO read parameter form used modules (mesh)
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present
    logical         :: ndomain_present 
    
    section_present = .false.
    ndomain_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:19) == '[LyaPhotonsFromGas]') then
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
          case ('outputfileRec')
             write(outputfileRec,'(a)') trim(value)
          case ('outputfileCol')
             write(outputfileCol,'(a)') trim(value)
          case ('repository')
             write(repository,'(a)') trim(value)
          case ('snapnum')
             read(value,*) snapnum
          case('emission_dom_type')
             write(emission_dom_type,'(a)') trim(value)
          case ('emission_dom_pos')
             read(value,*) emission_dom_pos(1),emission_dom_pos(2),emission_dom_pos(3)
          case ('emission_dom_rsp')
             read(value,*) emission_dom_rsp
          case ('emission_dom_rin')
             read(value,*) emission_dom_rin
          case ('emission_dom_rout')
             read(value,*) emission_dom_rout
          case ('emission_dom_size')
             read(value,*) emission_dom_size
          case ('emission_dom_thickness')
             read(value,*) emission_dom_thickness
          case('nphotons')
             read(value,*) nphotons
          case('ranseed')
             read(value,*) ranseed 
          case('doRecombs')
             read(value,*) doRecombs
          case('doColls')
             read(value,*) doColls
          case('tcool_resolution')
             read(value,*) tcool_resolution
          case ('verbose')
             read(value,*) verbose

          end select
       end do
    end if
    close(10)

    call read_ramses_params(pfile)
    
    return

  end subroutine read_LyaPhotonsFromGas_params

  
  subroutine print_LyaPhotonsFromGas_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)')         '[LyaPhotonsFromGas]'
       write(unit,'(a)')             '# input / output parameters'
       write(unit,'(a,a)')           '  outputfileRec   = ',trim(outputfileRec)
       write(unit,'(a,a)')           '  outputfileCol   = ',trim(outputfileCol)
       write(unit,'(a,a)')           '  repository      = ',trim(repository)
       write(unit,'(a,i5)')          '  snapnum         = ',snapnum
       write(unit,'(a)')             '# emissionational domain parameters'
       write(unit,'(a,a)')           '  emission_dom_type      = ',trim(emission_dom_type)
       write(unit,'(a,3(ES10.3,1x))')'  emission_dom_pos       = ',emission_dom_pos(1),emission_dom_pos(2),emission_dom_pos(3)
       select case(emission_dom_type)
       case ('sphere')
          write(unit,'(a,ES10.3)')      '  emission_dom_rsp       = ',emission_dom_rsp
       case ('shell')
          write(unit,'(a,ES10.3)')      '  emission_dom_rin       = ',emission_dom_rin
          write(unit,'(a,ES10.3)')      '  emission_dom_rout      = ',emission_dom_rout
       case ('cube')
          write(unit,'(a,ES10.3)')      '  emission_dom_size      = ',emission_dom_size
       case ('slab')
          write(unit,'(a,ES10.3)')      '  emission_dom_thickness = ',emission_dom_thickness
       end select
       write(unit,'(a)')             '# parameters'
       write(unit,'(a,i10)')           '  nphotons         = ',nphotons
       write(unit,'(a,i10)')           '  ranseed          = ',ranseed
       write(unit,'(a,L1)')            '  doRecombs        = ',doRecombs
       write(unit,'(a,L1)')            '  doColls          = ',doColls  
       write(unit,'(a,ES10.3)')        '  tcool_resolution = ',tcool_resolution
       write(unit,'(a)')             '# miscelaneous parameters'
       write(unit,'(a,L1)')          '  verbose            = ',verbose
       write(unit,'(a)')             ' '
       call print_ramses_params(unit)

    else
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ' '
       write(*,'(a,a,a)')         '[LyaPhotonsFromGas]'
       write(*,'(a)')             '# input / output parameters'
       write(*,'(a,a)')           '  outputfileRec = ',trim(outputfileRec)
       write(*,'(a,a)')           '  outputfileCol = ',trim(outputfileCol)
       write(*,'(a,a)')           '  repository = ',trim(repository)
       write(*,'(a,i5)')          '  snapnum    = ',snapnum
       write(*,'(a)')             '# emissionational domain parameters'
       write(*,'(a,a)')           '  emission_dom_type      = ',trim(emission_dom_type)
       write(*,'(a,3(ES10.3,1x))')'  emission_dom_pos       = ',emission_dom_pos(1),emission_dom_pos(2),emission_dom_pos(3)
       select case(emission_dom_type)
       case ('sphere')
          write(*,'(a,ES10.3)')      '  emission_dom_rsp       = ',emission_dom_rsp
       case ('shell')
          write(*,'(a,ES10.3)')      '  emission_dom_rin       = ',emission_dom_rin
          write(*,'(a,ES10.3)')      '  emission_dom_rout      = ',emission_dom_rout
       case ('cube')
          write(*,'(a,ES10.3)')      '  emission_dom_size      = ',emission_dom_size
       case ('slab')
          write(*,'(a,ES10.3)')      '  emission_dom_thickness = ',emission_dom_thickness
       end select
       write(*,'(a)')             '# parameters'
       write(*,'(a,i10)')           '  nphotons      = ',nphotons
       write(*,'(a,i10)')           '  ranseed       = ',ranseed
       write(*,'(a,L1)')            '  doRecombs     = ',doRecombs
       write(*,'(a,L1)')            '  doColls       = ',doColls
       write(*,'(a,ES10.3)')        '  tcool_resolution = ',tcool_resolution
       write(*,'(a)')             '# miscelaneous parameters'
       write(*,'(a,L1)')          '  verbose         = ',verbose
       write(*,'(a)')             ' '
       call print_ramses_params()
       write(*,'(a)')             ' '
       write(*,'(a)')             '--------------------------------------------------------------------------------'
       write(*,'(a)')             ' '
    end if

    return

  end subroutine print_LyaPhotonsFromGas_params


  
end program LyaPhotonsFromGas

