module module_uparallel

  private

  !
  ! contains tables (and functions to fill them up) which are used to 
  ! generate atom parallel velocoties by interpolation ... 
  ! 
  ! JB - 02/2011
  ! 
  ! 
  ! HOW TO 
  ! 
  ! 1/ call init_uparallel_tables once for a whole run.
  ! 2/ upar = get_uparallel(a,x,c), where a is a, x is x_int, and c is a random number between 0 and 1. 
  !

  ! two entries to the tables : x_int and a
  ! -> binning of a in log
  ! -> binning of x_int linear
  ! NB: THE VALUES OF THE FOLLOWING PARAMETERS HAVE BEEN TESTED AND OPTIMIZED ... 
  !     -> DO NOT CHANGE THEM WITHOUT A GOOD REASON AND LOTS OF TESTING. 

  ! sampling of a (in a reasonable range -> which depends on vturb ... )
  real(kind=8),parameter    :: amin     = -3.8d0  ! in log 
  real(kind=8),parameter    :: amax     = -1.7d0
  integer(kind=4),parameter :: nbins_a  = 50
  real(kind=8),parameter    :: a_step   = (amax-amin)/real(nbins_a-1,8)
  real(kind=8)              :: a_bins(nbins_a)

  ! sampling of x_int
  real(kind=8),parameter    :: xintmin    = -50.0d0
  real(kind=8),parameter    :: xintmax    = 0.0d0
  integer(kind=4),parameter :: nbins_xint = 500
  real(kind=8)              :: xint_bins(nbins_xint)
  real(kind=8),parameter    :: xint_step  = (xintmax-xintmin)/real(nbins_xint-1,8)

  ! tabulation of cumulative probability function 
  !       C(u) = a/(pi H) \int_-infinity^u e^-u2 / ((xint-u)2 + a2)
  ! --> in fact, we tabulate the inverse function u(c), on regular bins of c, which is what we need. 
  ! -> use three binnings for better accuracy ... 
  ! + cminimin < c < cmin    -> regular bins of log(c)
  ! + cmin < c < cmax -> regular bins of c
  ! + cmax < c < 1    -> regular bins of log(1-c)
  real(kind=8),parameter    :: cminimin = 1.d-10
  real(kind=8),parameter    :: cmin = 0.1d0
  real(kind=8),parameter    :: cmax = 0.9d0
  real(kind=8),parameter    :: cmaximax = 1.0d0 - 1.d-10
  integer(kind=4),parameter :: nbins_c1 = 5000
  integer(kind=4),parameter :: nbins_c2 = 2000
  integer(kind=4),parameter :: nbins_c3 = 2000
  real(kind=8)              :: c_bins1(nbins_c1)
  real(kind=8)              :: c_step1
  real(kind=8)              :: c_bins2(nbins_c2)
  real(kind=8)              :: c_step2
  real(kind=8)              :: c_bins3(nbins_c3)
  real(kind=8)              :: c_step3
  ! tables u(c) for each range of c's. First dim is c, then x_int, then a.
  real(kind=8)              :: table1(nbins_c1,nbins_xint,nbins_a)  ! low values of c
  real(kind=8)              :: table2(nbins_c2,nbins_xint,nbins_a)  ! mid values of c
  real(kind=8)              :: table3(nbins_c3,nbins_xint,nbins_a)  ! high values of c

  logical :: tables_initialized = .false.
  
  public :: get_uparallel
  
contains

  subroutine init_uparallel_tables
    
    implicit none 

    integer(kind=4) :: i,j
    !real(kind=8)    :: a,xin
    !integer(kind=4) :: xbin

    
    ! 1/ DEFINE BINNING OF TABLES 
    ! ---------------------------
    ! fill up a-bin values
    do i = 1,nbins_a
       a_bins(i) = amin + (i-1)*a_step
    end do
    a_bins = 10.0d0**a_bins

    ! fill up xint bins 
    do i = 1,nbins_xint
       xint_bins(i) = xintmin + (i-1)*xint_step
    end do
    
    ! fill up c bins
    c_step1 = (log10(cmin) - log10(cminimin))/real(nbins_c1-1,8)
    c_step2 = (cmax - cmin)/real(nbins_c2-1,8)
    c_step3 = (log10(1.0d0 - cmax) - log10(1.0d0-cmaximax))/real(nbins_c3-1,8)
    do i = 1,nbins_c1 
       c_bins1(i) = log10(cminimin) + (i-1) * c_step1
    end do
    c_bins1 = 10.0d0**c_bins1
    do i = 1,nbins_c2
       c_bins2(i) = cmin + (i-1) * c_step2
    end do
    do i = 1,nbins_c3
       c_bins3(nbins_c3-i+1) = log10(1.0d0-cmaximax) + (i-1)*c_step3
    end do
    c_bins3 = 1.0d0 - 10.0d0**(c_bins3)


    ! 2/ DEFINE TABLES (3 in total)
    ! -----------------------------
    do i = 1,nbins_a 
       do j = 1,nbins_xint
          call get_c_of_u(i,j)
       end do
    end do

    return
    
  end subroutine init_uparallel_tables


  subroutine get_c_of_u(ia,ix)
    
    ! 1/ compute C(u) = a / (pi H(a,x)) * \int_{umin}^u exp(-u^2) / [(x_int - u)^2 + a^2] du 
    ! with u ranging from umin to umax (linearly sampled)
    ! 
    ! 2/ invert into u(c), with special binning for c... 
    ! 
    ! 3/ save into tables (table1, table2, table2) appropriately

    implicit none 
    
    integer(kind=4),intent(in) :: ia,ix
    integer(kind=4),parameter  :: n = 200000  ! nb of steps for the integration. This should NOT be changed ! high res _is_ needed
    real(kind=8)               :: f(0:n),c_of_u(n),utable(0:n)
    real(kind=8),parameter     :: umin = -10.0d0
    real(kind=8),parameter     :: umax = 5.0d0
    real(kind=8),parameter     :: du = (umax-umin)/real(n-1,8)
    integer(kind=4)            :: i,j
    real(kind=8)               :: u,g,a2,c,cb,xx,yy,ub,x_int,a
    real(kind=8)               :: t1(nbins_c1),t2(nbins_c2),t3(nbins_c3)

    a     = a_bins(ia)
    x_int = xint_bins(ix)

    ! 1/ Integral of cumulative probability function 
    ! ----------------------------------------------
    ! compute f
    a2 = a*a
    do i = 0,n
       u = umin + (i-1)*du
       utable(i) = u
       g = (x_int - u)**2 + a2
       f(i) = exp(-u*u) / g
    end do
    ! integrate f into c_of_u (with trapezium rule)
    c_of_u(1) = 0.5d0*(f(1)+f(0))*du
    do i = 2,n
       c_of_u(i) = 0.5d0*(f(i)+f(i-1))*du + c_of_u(i-1)
    end do
    ! normalize (without call to Voigt function -> just force numerical norm)
    c_of_u = c_of_u / c_of_u(n)

    
    ! 2/ inversion of cumulative probability function
    ! -----------------------------------------------
    ! -> build u(c) for 3 different binnings of c (log(c), c, log(1-c)) depending 
    ! on c values. Save u(c) directly into global variables tableX. 

    ! 2.1 / first table for low values of c (sampled in log)
    t1 = 0.0d0
    j  = 1
    cb = log10(c_bins1(j))
    i  = 1
    do while (i <= n)
       c = log10(c_of_u(i))
       if (c > cb) then 
          ! interpolate value of u at cb
          xx = c - cb
          yy = cb - log10(c_of_u(i-1))
          ub = utable(i-1) * xx + utable(i) * yy 
          t1(j) = ub / (xx + yy)
          ! increment j 
          j  = j + 1
          if (j > nbins_c1) exit
          cb = log10(c_bins1(j))
       end if
       if (c <= cb) i = i + 1
    end do
    ! 2.2/ second table for intermediate values of c
    t2 = 0.0d0
    j = 1
    cb = c_bins2(j)
    i  = 1
    do while (i <= n)
       c = c_of_u(i)
       if (c > cb) then 
          ! interpolate value of u at cb
          xx = c - cb
          yy = cb - c_of_u(i-1)
          ub = utable(i-1) * xx + utable(i) * yy 
          t2(j) = ub / (xx + yy)
          ! increment j 
          j  = j + 1
          if (j > nbins_c2) exit
          cb = c_bins2(j)
       end if
       if (c <= cb) i = i + 1
    end do
    ! 2.3/ third table for values of c close to one
    t3 = 0.0d0
    j = 1
    cb = c_bins3(j)
    i  = 1
    do while (i <= n)
       c = c_of_u(i)
       if (c > cb) then 
          ! interpolate value of u at cb
          xx = c - cb
          yy = cb - c_of_u(i-1)
          ub = utable(i-1) * xx + utable(i) * yy 
          t3(j) = ub / (xx + yy)
          ! increment j 
          j  = j + 1
          if (j > nbins_c3) exit
          cb = c_bins3(j)
       end if
       if (c <= cb) i = i + 1
    end do

    ! 3/ save tables into global variables
    table1(:,ix,ia) = t1
    table2(:,ix,ia) = t2
    table3(:,ix,ia) = t3

    return
    
  end subroutine get_c_of_u

  function get_uparallel(a,xxin,c)

    ! this function returns an atom's parallel velocity drawn from probability function 
    ! P(u) =  a / (pi H(a,x)) * exp(-u^2) / [(x_int - u)^2 + a^2]
    ! The function uses tables precomputed for a fine (enough) grid of x_int and a values

    implicit none
    
    real(kind=8),intent(in) :: a,xxin,c
    real(kind=8)            :: xin
    real(kind=8)            :: get_uparallel
    integer(kind=4)         :: ia,ix,cbin,ic
    real(kind=8)            :: la,ccc,lc
    real(kind=8)            :: dc1,dc2

    ! initialise on first call
    if (.not. tables_initialized) then
       print*,'initialising u_parallel tables'
       call init_uparallel_tables
       print*,'--done initialising u_parallel tables'
    end if
    
    xin   = sign(xxin,-1.0d0) ! -> force x_int to be negative
    get_uparallel = 0.0d0
    
    ! get nearest tabulated a 
    la = log10(a)
    if (la <= amin) then 
       ia = 1
    else if (la >= amax) then 
       ia = nbins_a
    else
       ia = nint( (la-amin)/a_step ) + 1
    end if

    ! Get nearest tabulated x_int
    if (xin <= xintmin) then 
       ix = 1
    else if (xin > xintmax) then 
       print*,' this is actually a problem ... '
       stop
    else
       ix = nint( (xin-xintmin)/xint_step ) + 1
    end if

    ! get bin of c (and table index)
    if (c < cmin) then 
       ! bins in log(c)
       cbin = 1
       lc   = log10(c)
       ic   = int((lc-log10(cminimin))/c_step1) + 1
       if (ic < 1) then
          ic  = 1
          dc1 = 0.0d0
          dc2 = 1.0d0
       else if (ic >= nbins_c1) then 
          ic  = nbins_c1 -1
          dc1 = 1.0d0
          dc2 = 0.0d0
       else
          dc1 = lc - log10(c_bins1(ic))
          dc2 = log10(c_bins1(ic+1)) - lc
       end if
    else if (c <= cmax) then
       ! linear c bins
       cbin = 2
       ic   = int((c - cmin)/c_step2)+1
       if (ic < 1) then 
          ic  = 1
          dc1 = 0.0d0
          dc2 = 1.0d0
       else if (ic >= nbins_c2) then 
          ic = nbins_c2 - 1
          dc1 = 1.0d0
          dc2 = 0.0d0
       else 
          dc1 = c - c_bins2(ic)
          dc2 = c_bins2(ic+1) - c
       end if
    else ! c > cmax
       ! bins in log(1-c) ... and clearly not a log of a headache ... 
       ! -> linear interpolation, even though sampling is log(1-c)
       cbin = 3
       ccc = log10(1.0d0-c)
       ic   = nbins_c3 - int((ccc - log10(1.0d0-cmaximax))/c_step3) - 1
       if (ic < 1) then 
          ic = 1
          dc1 = 0.0d0
          dc2 = 1.0d0
       else if (ic >= nbins_c3) then 
          ic = nbins_c3-1
          dc1 = 1.0d0
          dc2 = 0.0d0
       else
          dc1 = ccc - log10(1.0d0 - c_bins3(ic))
          dc2 = log10(1.0d0 - c_bins3(ic+1)) - ccc
       end if
    end if

    select case(cbin) 
    case(1)
       get_uparallel = table1(ic,ix,ia)*dc2 + table1(ic+1,ix,ia)*dc1
    case(2)
       get_uparallel = table2(ic,ix,ia)*dc2 + table2(ic+1,ix,ia)*dc1
    case(3)
       get_uparallel = table3(ic,ix,ia)*dc2 + table3(ic+1,ix,ia)*dc1
    end select
    get_uparallel = get_uparallel / (dc1+dc2)

    ! if we changed sign of xin (i.e. xxin > 0), then apply symetry to uparallel
    if (xxin > 0.0) get_uparallel = sign(get_uparallel,-get_uparallel) ! -> assign sign of -uparallel to uparallel...

    return

  end function get_uparallel

end module module_uparallel



