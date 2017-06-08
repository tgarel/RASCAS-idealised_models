module module_idealised_models

  use module_random
  use module_constants, only : kb, mp, pi
  
  implicit none

  private

  !! common params to all models
  real(kind=8)             :: fix_vel             = -999.0d0         ! ad-hoc cell velocity (cm/s) -> NEED BETTER PARAMETERIZATION for more than static... 
  real(kind=8)             :: fix_box_size_cm     = -999.0d0         ! ad-hoc box size in cm.

  !! uniform model
  real(kind=8)             :: fix_nhi             = -999.0d0         ! ad-hoc HI density (H/cm3)
  real(kind=8)             :: fix_vth             = -999.0d0         ! ad-hoc thermal velocity (cm/s)
  real(kind=8)             :: fix_ndust           = -999.0d0         ! ad-hoc dust number density (/cm3)

  !! sphere models
  real(kind=8)             :: fix_tauH            = -999.0d0         ! total HI opacity
  real(kind=8)             :: fix_temp            = -999.0d0         ! temperature (K)
  real(kind=8)             :: fix_taudust         = -999.0d0         ! total dust opacity

  !! sphere_gradient model only
  character(100)           :: fix_dens_grad       = 'isothermal'  ! form of the density profile: isothermal, powerlaw
  real(kind=8)             :: rho_0               = -999.0d0
  real(kind=8)             :: vel_alpha           = -999.0d0
  real(kind=8)             :: dens_alpha          = -999.0d0

  !! disc models
  real(kind=8)             :: fix_mtot_hi         = -999.0d0         ! total HI mass [Msun]
  real(kind=8)             :: fix_rd_disc         = -999.0d0         ! Disc scalelength [box units]
  real(kind=8)             :: fix_zd_disc         = -999.0d0         ! Disc scale height [box units]
  
  !! Scaralata models
  real(kind=8)             :: rsf       = -0.1
  real(kind=8)             :: rw        = -0.4
  real(kind=8)             :: vel_gamma = -999.0d0
  real(kind=8)             :: rho_rsf   = -999.0d0
  
  ! miscelaneous
  logical                  :: verbose             = .false.       ! display some run-time info on this module
  
  ! public functions:
  public :: read_overwrite_params, print_overwrite_params, disc_thin, disc_thick, sphere_homogen_velfix, sphere_homogen_velgrad, sphere_densgrad_velfix, sphere_densgrad_velgrad, shell_homogen_velfix, sphere_homogen_steidel, sphere_homogen_velgrad_ct_outflow_rate, sphere_homogen_velgrad_rad_pressure, sphere_scarlata15
  
  
contains

    !!+++++++++++++++++++++++++++++ SPHERE a la Scarlata+15 : Velocity gradient and ct outflow rate (=>rho(r)) ++++++++++++++++++++++++++++++++++++++

  subroutine sphere_scarlata15(vx_ideal,vy_ideal,vz_ideal,nh_ideal,dopwidth_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)

    implicit none

    ! point source at center and medium transparent at R < rsf and at R > rw
    
    integer(kind=4)                :: nMC
    real(kind=8)                   :: dx_cell     ! Cell size
    real(kind=8)                   :: volfrac,volfrac2
    real(kind=8)                   :: dist_cell,dist_cell_min,dist_cell_max,dist2,dist_mc2,dist_mc
    real(kind=8)                   :: vx_mc,vy_mc,vz_mc
    real(kind=8)                   :: xcell_ideal,ycell_ideal,zcell_ideal
    real(kind=8)                   :: xmin,xmax,xmc,ymin,ymax,ymc,zmin,zmax,zmc
    real(kind=8),intent(inout)     :: nh_ideal,vx_ideal,vy_ideal,vz_ideal,dopwidth_ideal
    real(kind=8)                   :: r
    integer(kind=4)                :: localseed
    real(kind=8)                   :: rho
    real(kind=8)                   :: max_dist, dist_cell_mc
    
    ! define sphere radius : from param file now
!!$    sphere_radius = 0.45d0          ! Maximum extent of the gas  
!!$    max_dist      = sphere_radius
!!$    rw            = ! as in S+15: Radius at which gas speed reaches V_infinity 
!!$    min_dist      = rsf  !  in S+15
    
    vx_mc = 0.0d0
    vy_mc = 0.0d0
    vz_mc = 0.0d0
    dist2 = 0.0d0
    dist_cell_mc = 0.0d0

    ! xcell, ycell and zcell are in frame with origin at bottom-left corner of box
    dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 + (zcell_ideal-0.5d0)**2  ! in frame with origin at center of box

    dist_cell = sqrt(dist2)
    dist_cell_max = dist_cell + dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    dist_cell_min = dist_cell - dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
!    if (dist_cell_min > rw .or. dist_cell_max < rsf) then  ! cell completely out of sphere or completely within Rsf
!       vx_ideal  = 0.0d0
!       vy_ideal  = 0.0d0
!       vz_ideal  = 0.0d0
!       nh_ideal  = 0.0d0
!    end if
    if ((dist_cell_max > rw .and. dist_cell_min < rw) .or. (dist_cell_max > rsf .and. dist_cell_min < rsf)) then ! cell partially within sphere
       ! bounds of cell
       xmin = xcell_ideal - 0.5d0*dx_cell 
       xmax = xcell_ideal + 0.5d0*dx_cell
       ymin = ycell_ideal - 0.5d0*dx_cell 
       ymax = ycell_ideal + 0.5d0*dx_cell
       zmin = zcell_ideal - 0.5d0*dx_cell 
       zmax = zcell_ideal + 0.5d0*dx_cell
       
       ! use Monte Carlo to compute density of cell (i.e. fraction of its volume within the sphere)
       volfrac = 0.0d0
       nMC = 0
       do while(volfrac .lt. 1000.d0)
          if (volfrac .lt. 1.d0 .and. nMC .eq. 100) then
             exit
          end if
          r = ran3(localseed)
          xmc = r * dx_cell + xmin
          r = ran3(localseed)
          ymc = r * dx_cell + ymin
          r = ran3(localseed)
          zmc = r * dx_cell + zmin
          dist_mc2 = 0.0d0
          dist_mc2 = (xmc-0.5d0)**2 + (ymc-0.5d0)**2 + (zmc-0.5d0)**2
          dist_mc = sqrt(dist_mc2)
          if (dist_mc < rw .and. dist_mc > rsf) then ! point within sphere
             volfrac = volfrac + 1.0d0
             ! fix_vel = V0 from S+15
             vx_mc = vx_mc + fix_vel * (xmc - 0.5d0) / rw 
             vy_mc = vy_mc + fix_vel * (ymc - 0.5d0) / rw
             vz_mc = vz_mc + fix_vel * (zmc - 0.5d0) / rw
             dist_cell_mc = dist_cell_mc + dist_mc
          end if
          nMC = nMC + 1
       end do
       
       if(volfrac .ne. 0) then
          ! Velocity
          vx_ideal = vx_mc / volfrac
          vy_ideal = vy_mc / volfrac
          vz_ideal = vz_mc / volfrac
          ! redefine dist_cell as mean of dist_mc
          dist_cell = dist_cell_mc / volfrac
       end if
       volfrac2 = volfrac / dble(nMC)
       
    end if
    if (dist_cell_max < rw .and. dist_cell_min > rsf) then ! cell completely within sphere
       vx_ideal = fix_vel * (xcell_ideal-0.5d0) / rw
       vy_ideal = fix_vel * (ycell_ideal-0.5d0) / rw
       vz_ideal = fix_vel * (zcell_ideal-0.5d0) / rw
       volfrac2 = 1.0
    end if
    
    ! assign density, ndust
    rho      = rho_rsf *  (dist_cell / rsf)**(-vel_gamma-2.)
    nh_ideal = rho * volfrac2
 
    if (dist_cell_min > rw .or. dist_cell_max < rsf) then  ! cell completely out of sphere or completely within Rsf                                                       
       vx_ideal  = 0.0d0
       vy_ideal  = 0.0d0
       vz_ideal  = 0.0d0
       nh_ideal  = 0.0d0
    end if
 
 ! Give all cells shell_temp (cells not in sphere won't matter for RT as dnh = 0...)
 dopwidth_ideal = sqrt((2.0d0*kb/mp)*fix_temp) 
    
 return
 
  end subroutine sphere_scarlata15


  
  !!+++++++++++++++++++++++++++++ THIN DISC ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
   subroutine disc_thin(vx_ideal,vy_ideal,vz_ideal,nh_ideal,dopwidth_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)

    implicit none
    
    integer(kind=4)                :: nMC
    real(kind=8)                   :: dx_cell     ! Cell size
    real(kind=8)                   :: volfrac,volfrac2
    real(kind=8)                   :: dist_cell,dist_cell_min,dist_cell_max,dist2,dist_mc2,dist_mc
    real(kind=8)                   :: vx_mc,vy_mc,vz_mc
    real(kind=8)                   :: xcell_ideal,ycell_ideal,zcell_ideal
    real(kind=8)                   :: xmin,xmax,xmc,ymin,ymax,ymc,zmin,zmax,zmc
    real(kind=8),intent(inout)     :: nh_ideal,vx_ideal,vy_ideal,vz_ideal,dopwidth_ideal
    real(kind=8)                   :: r
    integer(kind=4)                :: localseed
    real(kind=8)                   :: max_dist, dist_cell_mc
    real(kind=8),parameter         :: msun_to_gram  = 1.9891d33 ! g
    real(kind=8),parameter         :: gram_to_atoms = 1.673534d-24 ! Hydrogen mass in g
    real(kind=8)                   :: Sigma0, dz, rdisc_max , logSigma0, log_mtot


    ! fix_mtot_hi = fix_mtot_hi * msun_to_gram
    ! fix_mtot_hi = fix_mtot_hi / gram_to_atoms    ! Number of HI atoms
    ! in log10
    log_mtot = log10(fix_mtot_hi) + log10(msun_to_gram) - log10(gram_to_atoms)    ! Msun -> g -> atoms
   
    dz          = dx_cell                        ! thickness of disc = size of cell
    rdisc_max   = 0.45d0                         ! rdisc_max*fix_box_size_cm should be equal to the HI radius measured for LARS object
    max_dist    = rdisc_max
    vx_mc       = 0.0d0
    vy_mc       = 0.0d0
    vz_mc       = 0.0d0
    nh_ideal    = 0.d0
    dist_cell_mc = 0.0d0

    !Sigma0    = fix_mtot_hi / 2.d0 / pi / (fix_rd_disc * fix_box_size_cm)**2 ! atoms/cm^2 - Central HI surface density
    !! in log10
    logSigma0  = log_mtot - log10(2.d0) - log10(pi) - 2.d0 * log10(fix_rd_disc * fix_box_size_cm)
    Sigma0    = 10.d0**logSigma0
        
    ! Select only cells which are exactly in the z=0 plane. Others cells have nh_ideal=0.d0 (from initialization)
    if (abs(zcell_ideal-0.5d0) .lt. 0.5*dz) then !  if (abs(zcell_ideal-0.5d0) .eq. 0.d0) should also work
       
       ! Then, I should now have only cells with (zcell_ideal-0.5d0)=0
       dist2 = 0.0d0
       dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 ! squared distance of cell from center in xy plane
       
       dist_cell = sqrt(dist2)
       dist_cell_max = dist_cell + dx_cell / sqrt(2.0d0) ! dx_cell / sqrt(2.0) is half the longest length in a square
       dist_cell_min = dist_cell - dx_cell / sqrt(2.0d0) 
       if (dist_cell_min > max_dist) then  ! cell completely out of disc 
          vx_ideal  = 0.0d0
          vy_ideal  = 0.0d0
          vz_ideal  = 0.0d0
          nh_ideal  = 0.0d0
       else
          ! bounds of cell
          xmin = xcell_ideal - 0.5d0*dx_cell 
          xmax = xcell_ideal + 0.5d0*dx_cell
          ymin = ycell_ideal - 0.5d0*dx_cell 
          ymax = ycell_ideal + 0.5d0*dx_cell
          zmin = zcell_ideal - 0.5d0*dx_cell 
          zmax = zcell_ideal + 0.5d0*dx_cell
          
          ! use Monte Carlo to compute density of cell
          volfrac = 0.0d0
          nMC = 0
          do while(volfrac .lt. 1000.d0)
             if (volfrac .lt. 1.d0 .and. nMC .eq. 100) then
                exit
             end if
             r   = ran3(localseed)
             xmc = r * dx_cell + xmin
             r   = ran3(localseed)
             ymc = r * dx_cell + ymin
             r   = ran3(localseed)
             zmc = r * dx_cell + zmin
             ! test if point is in the disc
             dist_mc2 = 0.0d0
             dist_mc2 = (xmc-0.5d0)**2 + (ymc-0.5d0)**2 + (zmc-0.5d0)**2
             dist_mc  = sqrt(dist_mc2)
             if (dist_mc < max_dist) then ! point within sphere
                volfrac = volfrac + 1.0d0   
!!$             vx_mc = vx_mc + fix_vel * (xmc - 0.5d0) / dist_mc
!!$             vy_mc = vy_mc + fix_vel * (ymc - 0.5d0) / dist_mc
!!$             vz_mc = vz_mc + fix_vel * (zmc - 0.5d0) / dist_mc
                dist_cell_mc = dist_cell_mc + dist_mc
             end if
             nMC = nMC + 1
          end do
          
          volfrac2 = volfrac / dble(nMC)
          
       if(volfrac .ne. 0) then
          ! Velocity
          vx_ideal = vx_mc / volfrac
          vy_ideal = vy_mc / volfrac
          vz_ideal = vz_mc / volfrac
          ! redefine dist_cell as mean of dist_mc
          dist_cell = dist_cell_mc / volfrac
       end if
          
          ! Static disc for now....
          vx_ideal  = 0.0d0
          vy_ideal  = 0.0d0
          vz_ideal  = 0.0d0
          
       ! assign density, ndust
          nh_ideal = Sigma0 / (dz * fix_box_size_cm) * exp(-dist_cell/fix_rd_disc) * volfrac2 ! atoms/cm^3 -  here i take center of cell... could use MC points (weighting?)
          
          if(volfrac2 .gt. 1.0d0) then
             print*,volfrac2,volfrac,nMC
             stop
          endif
       end if
    else
       vx_ideal  = 0.0d0
       vy_ideal  = 0.0d0
       vz_ideal  = 0.0d0
       nh_ideal  = 0.0d0
    end if
    
    ! Give all cells shell_temp (cells not in sphere won't matter for RT as dnh = 0...)
    dopwidth_ideal = sqrt((2.0d0*kb/mp)*fix_temp) 

!!$    if (xcell_ideal < 0.5 .and. ycell_ideal < 0.5 .and. nh_ideal > 0.d0) then
!!$       print*,'XY = ',xcell_ideal,ycell_ideal,sqrt((xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2)-dx_cell / sqrt(2.0d0),nh_ideal
!!$    end if
       
    return
     
  end subroutine disc_thin

  !!+++++++++++++++++++++++++++++ THICK DISC ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  subroutine disc_thick(vx_ideal,vy_ideal,vz_ideal,nh_ideal,dopwidth_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)

    implicit none
    
    integer(kind=4)                :: nMC
    real(kind=8)                   :: dx_cell     ! Cell size
    real(kind=8)                   :: volfrac,volfrac2
    real(kind=8)                   :: dist_cell,dist_cell_min,dist_cell_max,dist2,dist_mc2,dist_mc
    real(kind=8)                   :: vx_mc,vy_mc,vz_mc
    real(kind=8)                   :: xcell_ideal,ycell_ideal,zcell_ideal
    real(kind=8)                   :: xmin,xmax,xmc,ymin,ymax,ymc,zmin,zmax,zmc
    real(kind=8),intent(inout)     :: nh_ideal,vx_ideal,vy_ideal,vz_ideal,dopwidth_ideal
    real(kind=8)                   :: r
    integer(kind=4)                :: localseed
    real(kind=8)                   :: max_dist, dist_cell_mc
    real(kind=8),parameter         :: msun_to_gram  = 1.9891d33 ! g
    real(kind=8),parameter         :: gram_to_atoms = 1.673534d-24 ! Hydrogen mass in g
    real(kind=8)                   :: Sigma0, dz, r_max , logSigma0, log_mtot, mtot_hi_tilde
    !real(kind=8)                   :: fix_zd_disc

    !! Disc in XY plane, with rho(z) profile too
    
    ! fix_mtot_hi = fix_mtot_hi * msun_to_gram
    ! fix_mtot_hi = fix_mtot_hi / gram_to_atoms    ! Number of HI atoms

    dz          = dx_cell                        ! thickness of disc = size of cell
    r_max       = 0.45d0                         ! rdisc_max*fix_box_size_cm should be equal to the HI radius measured for LARS object
    max_dist    = r_max
    vx_mc       = 0.0d0
    vy_mc       = 0.0d0
    vz_mc       = 0.0d0
    nh_ideal    = 0.0d0
    dist_cell_mc = 0.0d0

    ! I use mtot_hi_tilde instead of fix_mtot_hi such that fix_mtot_hi is encompssed into r_max
    mtot_hi_tilde = fix_mtot_hi / (1.0-(r_max/fix_rd_disc+1.0)*exp(-r_max/fix_rd_disc))
    ! in log10
    ! log_mtot = log10(fix_mtot_hi) + log10(msun_to_gram) - log10(gram_to_atoms)    ! Msun -> g -> atoms
    log_mtot = log10(mtot_hi_tilde) + log10(msun_to_gram) - log10(gram_to_atoms)    ! Msun -> g -> atoms
    
    ! Cut off cells at R(x,y,z) > r_max, i.e. out of sphere domain of R = r_max
    dist2 = 0.0d0
    dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 + (zcell_ideal-0.5d0)**2 ! squared distance of cell from center in xy plane
    
    dist_cell = sqrt(dist2)
    dist_cell_max = dist_cell + dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    dist_cell_min = dist_cell - dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    if (dist_cell_min > max_dist) then  ! cell completely out of disc
       vx_ideal  = 0.0d0
       vy_ideal  = 0.0d0
       vz_ideal  = 0.0d0
       nh_ideal  = 0.0d0
    else ! cell partially or completely within disc
       if (dist_cell_max > max_dist) then ! cell partially within disc
          
          ! bounds of cell
          xmin = xcell_ideal - 0.5d0*dx_cell 
          xmax = xcell_ideal + 0.5d0*dx_cell
          ymin = ycell_ideal - 0.5d0*dx_cell 
          ymax = ycell_ideal + 0.5d0*dx_cell
          zmin = zcell_ideal - 0.5d0*dx_cell 
          zmax = zcell_ideal + 0.5d0*dx_cell
          
          ! use Monte Carlo to compute density of cell
          volfrac = 0.0d0
          nMC = 0
          do while(volfrac .lt. 1000.d0)
             if (volfrac .lt. 1.d0 .and. nMC .eq. 100) then
                exit
             end if
             r   = ran3(localseed)
             xmc = r * dx_cell + xmin
             r   = ran3(localseed)
             ymc = r * dx_cell + ymin
             r   = ran3(localseed)
             zmc = r * dx_cell + zmin
             dist_mc2 = 0.0d0
             dist_mc2 = (xmc-0.5d0)**2 + (ymc-0.5d0)**2 + (zmc-0.5d0)**2
             dist_mc  = sqrt(dist_mc2)
             if (dist_mc < max_dist) then ! point within sphere domain
                volfrac = volfrac + 1.0d0   
!!$             vx_mc = vx_mc + fix_vel * (xmc - 0.5d0) / dist_mc
!!$             vy_mc = vy_mc + fix_vel * (ymc - 0.5d0) / dist_mc
!!$             vz_mc = vz_mc + fix_vel * (zmc - 0.5d0) / dist_mc
                dist_cell_mc = dist_cell_mc + dist_mc
             end if
             nMC = nMC + 1
          end do

          if(volfrac .ne. 0) then
             ! redefine dist_cell as mean of dist_mc
             dist_cell = dist_cell_mc / volfrac
          end if
          
          volfrac2 = volfrac / dble(nMC)
          
       else ! cell completely within disc
          volfrac2 = 1.0
       end if
       
       ! Static disc for now....
       vx_ideal  = 0.0d0
       vy_ideal  = 0.0d0
       vz_ideal  = 0.0d0
       
       ! assign density, ndust
       ! log_mtot in atoms here
       nh_ideal = 10.d0**log_mtot / 4.d0 / pi / (fix_zd_disc*fix_box_size_cm) / (fix_rd_disc*fix_box_size_cm)**2. * exp(-dist_cell/fix_rd_disc) * 1.d0 / dcosh(-abs(zcell_ideal-0.5d0)/fix_zd_disc)**2.
       nh_ideal = nh_ideal * volfrac2
       
       !! nh_ideal = Sigma0 / (dz * fix_box_size_cm) * exp(-dist_cell/fix_rd_disc) * exp(-abs(zcell_ideal-0.5d0)/fix_zd_disc) * volfrac2 ! atoms/cm^3 -  here i take center of cell... could use MC points (weighting?)
       
       if(volfrac2 .gt. 1.0d0) then
          print*,volfrac2,volfrac,nMC
          stop
       endif
    end if
    
    ! Give all cells shell_temp (cells not in sphere won't matter for RT as dnh = 0...)
    dopwidth_ideal = sqrt((2.0d0*kb/mp)*fix_temp) 
    
    return
    
  end subroutine disc_thick
  
 !!+++++++++++++++++++++++++++++ HOMOGEN SHELL with Fixed Velocity ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  subroutine shell_homogen_velfix(vx_ideal,vy_ideal,vz_ideal,nh_ideal,dopwidth_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)
    
    implicit none
    
    integer(kind=4)                :: nMC
    real(kind=8)                   :: dx_cell     ! Cell size
    real(kind=8)                   :: volfrac,volfrac2
    real(kind=8)                   :: dist_cell,dist_cell_min,dist_cell_max,dist2,dist_mc2,dist_mc
    real(kind=8)                   :: vx_mc,vy_mc,vz_mc
    real(kind=8)                   :: xcell_ideal,ycell_ideal,zcell_ideal
    real(kind=8)                   :: xmin,xmax,xmc,ymin,ymax,ymc,zmin,zmax,zmc
    real(kind=8),intent(inout)     :: nh_ideal,vx_ideal,vy_ideal,vz_ideal,dopwidth_ideal
    real(kind=8)                   :: r
    integer(kind=4)                :: localseed
    real(kind=8)                   :: sphere_radius,rho
    real(kind=8)                   :: max_dist,min_dist, dist_cell_mc
    real(kind=8)                   :: shell_dr
   ! real(kind=8)                   :: max_dist2,min_dist2,dist,scale_rho

  
    ! define sphere radius
    sphere_radius = 0.45d0          ! code units
    max_dist      = sphere_radius

    shell_dr      = max_dist / 5.0d0
    min_dist      = max_dist - shell_dr
  
    rho  = rho_shell_from_tau0_T_rs_Lbox(fix_tauH,fix_temp,shell_dr)
    
    vx_mc = 0.0d0
    vy_mc = 0.0d0
    vz_mc = 0.0d0
    dist2 = 0.0d0
    dist_cell_mc = 0.0d0

    ! xcell, ycell and zcell are in frame with origin at bottom-left corner of box
    dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 + (zcell_ideal-0.5d0)**2  ! in frame with origin at center of box

    dist_cell = sqrt(dist2)
    dist_cell_max = dist_cell + dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    dist_cell_min = dist_cell - dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    if (dist_cell_min > max_dist .or. dist_cell_max < min_dist) then  ! cell not in sphere region
       vx_ideal  = 0.0d0
       vy_ideal  = 0.0d0
       vz_ideal  = 0.0d0
       nh_ideal  = 0.0d0
    else
       ! bounds of cell
       xmin = xcell_ideal - 0.5d0*dx_cell 
       xmax = xcell_ideal + 0.5d0*dx_cell
       ymin = ycell_ideal - 0.5d0*dx_cell 
       ymax = ycell_ideal + 0.5d0*dx_cell
       zmin = zcell_ideal - 0.5d0*dx_cell 
       zmax = zcell_ideal + 0.5d0*dx_cell
       
       ! use Monte Carlo to compute density of cell (i.e. fraction of its volume within the sphere)
       volfrac = 0.0d0
       nMC = 0
       do while(volfrac .lt. 1000.d0)
          if (volfrac .lt. 1.d0 .and. nMC .eq. 100) then
             exit
          end if
          r = ran3(localseed)
          xmc = r * dx_cell + xmin
          r = ran3(localseed)
          ymc = r * dx_cell + ymin
          r = ran3(localseed)
          zmc = r * dx_cell + zmin
          ! test if point is in the sphere 
          dist_mc2 = 0.0d0
          dist_mc2 = (xmc-0.5d0)**2 + (ymc-0.5d0)**2 + (zmc-0.5d0)**2
          dist_mc = sqrt(dist_mc2)
          if (dist_mc < max_dist .and. dist_mc > min_dist) then ! point within sphere
             volfrac = volfrac + 1.0d0
             vx_mc = vx_mc + fix_vel * (xmc - 0.5d0) / dist_mc
             vy_mc = vy_mc + fix_vel * (ymc - 0.5d0) / dist_mc
             vz_mc = vz_mc + fix_vel * (zmc - 0.5d0) / dist_mc
             dist_cell_mc = dist_cell_mc + dist_mc
          end if
          nMC = nMC + 1
       end do
    
       volfrac2 = volfrac / dble(nMC)
       
       if(volfrac .ne. 0) then
          vx_ideal = vx_mc / volfrac
          vy_ideal = vy_mc / volfrac
          vz_ideal = vz_mc / volfrac
          ! redefine dist_cell as mean of dist_mc
          dist_cell = dist_cell_mc / volfrac
       end if

       ! assign density, ndust
       nh_ideal = rho * volfrac2
       
    end if

    dopwidth_ideal = sqrt((2.0d0*kb/mp)*fix_temp) 
    
    return
     
  end subroutine shell_homogen_velfix

  
  !!+++++++++++++++++++++++++++++ HOMOGEN SPHERE with Fixed Velocity ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  subroutine sphere_homogen_velfix(vx_ideal,vy_ideal,vz_ideal,nh_ideal,dopwidth_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)
    
    implicit none
    
    integer(kind=4)                :: nMC
    real(kind=8)                   :: dx_cell     ! Cell size
    real(kind=8)                   :: volfrac,volfrac2
    real(kind=8)                   :: dist_cell,dist_cell_min,dist_cell_max,dist2,dist_mc2,dist_mc
    real(kind=8)                   :: vx_mc,vy_mc,vz_mc
    real(kind=8)                   :: xcell_ideal,ycell_ideal,zcell_ideal
    real(kind=8)                   :: xmin,xmax,xmc,ymin,ymax,ymc,zmin,zmax,zmc
    real(kind=8),intent(inout)     :: nh_ideal,vx_ideal,vy_ideal,vz_ideal,dopwidth_ideal
    real(kind=8)                   :: r
    integer(kind=4)                :: localseed
    real(kind=8)                   :: sphere_radius,rho
    real(kind=8)                   :: max_dist,min_dist, dist_cell_mc

    ! define sphere radius
    sphere_radius = 0.45d0          ! code units
    max_dist      = sphere_radius
    min_dist      = 0.d0

    rho  = rho_sphere_from_tau0_T_rs_Lbox(fix_tauH,fix_temp,sphere_radius)
    
    vx_mc = 0.0d0
    vy_mc = 0.0d0
    vz_mc = 0.0d0
    dist2 = 0.0d0
    dist_cell_mc = 0.0d0

    ! xcell, ycell and zcell are in frame with origin at bottom-left corner of box
    dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 + (zcell_ideal-0.5d0)**2  ! in frame with origin at center of box

    dist_cell = sqrt(dist2)
    dist_cell_max = dist_cell + dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    dist_cell_min = dist_cell - dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    if (dist_cell_min > max_dist) then  ! cell completely out of sphere
       vx_ideal  = 0.0d0
       vy_ideal  = 0.0d0
       vz_ideal  = 0.0d0
    else  ! cell partially or completely within sphere
       if (dist_cell_max > max_dist) then ! cell partially within sphere
          ! bounds of cell
          xmin = xcell_ideal - 0.5d0*dx_cell 
          xmax = xcell_ideal + 0.5d0*dx_cell
          ymin = ycell_ideal - 0.5d0*dx_cell 
          ymax = ycell_ideal + 0.5d0*dx_cell
          zmin = zcell_ideal - 0.5d0*dx_cell 
          zmax = zcell_ideal + 0.5d0*dx_cell
          
          ! use Monte Carlo to compute density of cell (i.e. fraction of its volume within the sphere)
          volfrac = 0.0d0
          nMC = 0
          do while(volfrac .lt. 1000.d0)
             if (volfrac .lt. 1.d0 .and. nMC .eq. 100) then
                exit
             end if
             r = ran3(localseed)
             xmc = r * dx_cell + xmin
             r = ran3(localseed)
             ymc = r * dx_cell + ymin
             r = ran3(localseed)
             zmc = r * dx_cell + zmin
             dist_mc2 = 0.0d0
             dist_mc2 = (xmc-0.5d0)**2 + (ymc-0.5d0)**2 + (zmc-0.5d0)**2
             dist_mc = sqrt(dist_mc2)
             if (dist_mc < max_dist .and. dist_mc > min_dist) then ! point within sphere
                volfrac = volfrac + 1.0d0
                vx_mc = vx_mc + fix_vel * (xmc - 0.5d0) / dist_mc
                vy_mc = vy_mc + fix_vel * (ymc - 0.5d0) / dist_mc
                vz_mc = vz_mc + fix_vel * (zmc - 0.5d0) / dist_mc
                dist_cell_mc = dist_cell_mc + dist_mc
             end if
             nMC = nMC + 1
          end do

          if(volfrac .ne. 0) then
             ! Velocity
             vx_ideal = vx_mc / volfrac
             vy_ideal = vy_mc / volfrac
             vz_ideal = vz_mc / volfrac
             ! redefine dist_cell as mean of dist_mc
             dist_cell = dist_cell_mc / volfrac
          end if
          
          volfrac2 = volfrac / dble(nMC)
          
       else ! cell completely within sphere
          vx_ideal = fix_vel * (xcell_ideal-0.5d0) / dist_cell
          vy_ideal = fix_vel * (ycell_ideal-0.5d0) / dist_cell
          vz_ideal = fix_vel * (zcell_ideal-0.5d0) / dist_cell
          volfrac2 = 1.0
       end if
       
       ! assign density, ndust
       nh_ideal = rho * volfrac2
       
    end if
    
    ! Give all cells shell_temp (cells not in sphere won't matter for RT as dnh = 0...)
    dopwidth_ideal = sqrt((2.0d0*kb/mp)*fix_temp) 
    
    return
    
  end subroutine sphere_homogen_velfix


  !!+++++++++++++++++++++++++++++ HOMOGEN SPHERE with Velocity gradient ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sphere_homogen_velgrad(vx_ideal,vy_ideal,vz_ideal,nh_ideal,dopwidth_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)

    implicit none
    
    integer(kind=4)                :: nMC
    real(kind=8)                   :: dx_cell     ! Cell size
    real(kind=8)                   :: volfrac,volfrac2
    real(kind=8)                   :: dist_cell,dist_cell_min,dist_cell_max,dist2,dist_mc2,dist_mc
    real(kind=8)                   :: vx_mc,vy_mc,vz_mc
    real(kind=8)                   :: xcell_ideal,ycell_ideal,zcell_ideal
    real(kind=8)                   :: xmin,xmax,xmc,ymin,ymax,ymc,zmin,zmax,zmc
    real(kind=8),intent(inout)     :: nh_ideal,vx_ideal,vy_ideal,vz_ideal,dopwidth_ideal
    real(kind=8)                   :: r
    integer(kind=4)                :: localseed
    real(kind=8)                   :: sphere_radius,rho
    real(kind=8)                   :: max_dist,min_dist, dist_cell_mc

    ! define sphere radius
    sphere_radius = 0.45d0          ! code units
    max_dist      = sphere_radius
    min_dist      = 0.d0

    rho  = rho_sphere_from_tau0_T_rs_Lbox(fix_tauH,fix_temp,sphere_radius)
    
    vx_mc        = 0.0d0
    vy_mc        = 0.0d0
    vz_mc        = 0.0d0
    dist2        = 0.0d0
    dist_cell_mc = 0.0d0
    
    ! xcell, ycell and zcell are in frame with origin at bottom-left corner of box
    dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 + (zcell_ideal-0.5d0)**2  ! in frame with origin at center of box

    dist_cell = sqrt(dist2)
    dist_cell_max = dist_cell + dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    dist_cell_min = dist_cell - dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    if (dist_cell_min > max_dist .or. dist_cell_max < min_dist) then   ! cell completely out of sphere
       vx_ideal  = 0.0d0
       vy_ideal  = 0.0d0
       vz_ideal  = 0.0d0
    else ! cell partially or completely within sphere
       if (dist_cell_max > max_dist) then ! cell partially within sphere
          ! bounds of cell
          xmin = xcell_ideal - 0.5d0*dx_cell 
          xmax = xcell_ideal + 0.5d0*dx_cell
          ymin = ycell_ideal - 0.5d0*dx_cell 
          ymax = ycell_ideal + 0.5d0*dx_cell
          zmin = zcell_ideal - 0.5d0*dx_cell 
          zmax = zcell_ideal + 0.5d0*dx_cell
          
          ! use Monte Carlo to compute density of cell (i.e. fraction of its volume within the sphere)
          volfrac = 0.0d0
          nMC = 0
          do while(volfrac .lt. 1000.d0)
             if (volfrac .lt. 1.d0 .and. nMC .eq. 100) then
                exit
             end if
             r = ran3(localseed)
             xmc = r * dx_cell + xmin
             r = ran3(localseed)
             ymc = r * dx_cell + ymin
             r = ran3(localseed)
             zmc = r * dx_cell + zmin
             dist_mc2 = 0.0d0
             dist_mc2 = (xmc-0.5d0)**2 + (ymc-0.5d0)**2 + (zmc-0.5d0)**2
             dist_mc = sqrt(dist_mc2)
             if (dist_mc < max_dist .and. dist_mc > min_dist) then ! point within sphere
                volfrac = volfrac + 1.0d0
                vx_mc = vx_mc + fix_vel * (xmc - 0.5d0) / max_dist
                vy_mc = vy_mc + fix_vel * (ymc - 0.5d0) / max_dist
                vz_mc = vz_mc + fix_vel * (zmc - 0.5d0) / max_dist
                dist_cell_mc = dist_cell_mc + dist_mc
             end if
             nMC = nMC + 1
          end do
          
          if(volfrac .ne. 0) then
             ! Velocity
             vx_ideal = vx_mc / volfrac
             vy_ideal = vy_mc / volfrac
             vz_ideal = vz_mc / volfrac
             ! redefine dist_cell as mean of dist_mc
             dist_cell = dist_cell_mc / volfrac
          end if
          
          volfrac2 = volfrac / dble(nMC)
          
       else ! cell completely within sphere
          vx_ideal = fix_vel * (xcell_ideal-0.5d0) / max_dist
          vy_ideal = fix_vel * (ycell_ideal-0.5d0) / max_dist
          vz_ideal = fix_vel * (zcell_ideal-0.5d0) / max_dist
          volfrac2 = 1.0
       end if
       
       ! assign density, ndust
       nh_ideal = rho * volfrac2
    end if
       
    dopwidth_ideal = sqrt((2.0d0*kb/mp)*fix_temp) 
    
    return
     
  end subroutine sphere_homogen_velgrad

 !!+++++++++++++++++++++++++++++ SPHERE with Velocity gradient and ct outflow rate (=>rho(r)) ++++++++++++++++++++++++++++++++++++++

  subroutine sphere_homogen_velgrad_ct_outflow_rate(vx_ideal,vy_ideal,vz_ideal,nh_ideal,dopwidth_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)

    implicit none
    
    integer(kind=4)                :: nMC
    real(kind=8)                   :: dx_cell     ! Cell size
    real(kind=8)                   :: volfrac,volfrac2
    real(kind=8)                   :: dist_cell,dist_cell_min,dist_cell_max,dist2,dist_mc2,dist_mc
    real(kind=8)                   :: vx_mc,vy_mc,vz_mc
    real(kind=8)                   :: xcell_ideal,ycell_ideal,zcell_ideal
    real(kind=8)                   :: xmin,xmax,xmc,ymin,ymax,ymc,zmin,zmax,zmc
    real(kind=8),intent(inout)     :: nh_ideal,vx_ideal,vy_ideal,vz_ideal,dopwidth_ideal
    real(kind=8)                   :: r
    integer(kind=4)                :: localseed
    real(kind=8)                   :: sphere_radius,rho
    real(kind=8)                   :: max_dist, min_dist, rho_0, rmin, rmax, dist_cell_mc

    !! Here V=0 and rho=0 at r < min_dist
    
    ! define sphere radius
    sphere_radius = 0.45d0          ! code units
    max_dist      = sphere_radius
    min_dist      = 0.02d0

    rmin          = min_dist
    rmax          = max_dist
        
    vx_mc        = 0.0d0
    vy_mc        = 0.0d0
    vz_mc        = 0.0d0
    dist2        = 0.0d0
    dist_cell_mc = 0.0d0
    
    ! xcell, ycell and zcell are in frame with origin at bottom-left corner of box
    dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 + (zcell_ideal-0.5d0)**2  ! in frame with origin at center of box

    dist_cell = sqrt(dist2)
    dist_cell_max = dist_cell + dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    dist_cell_min = dist_cell - dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    if (dist_cell_min > max_dist .or. dist_cell_max < min_dist) then  ! cell not in sphere region
       vx_ideal  = 0.0d0
       vy_ideal  = 0.0d0
       vz_ideal  = 0.0d0
       nh_ideal  = 0.0d0
    else
       ! bounds of cell
       xmin = xcell_ideal - 0.5d0*dx_cell 
       xmax = xcell_ideal + 0.5d0*dx_cell
       ymin = ycell_ideal - 0.5d0*dx_cell 
       ymax = ycell_ideal + 0.5d0*dx_cell
       zmin = zcell_ideal - 0.5d0*dx_cell 
       zmax = zcell_ideal + 0.5d0*dx_cell
       
       ! use Monte Carlo to compute density of cell (i.e. fraction of its volume within the sphere)
       volfrac = 0.0d0
       nMC = 0
       do while(volfrac .lt. 1000.d0)
          if (volfrac .lt. 1.d0 .and. nMC .eq. 100) then
             exit
          end if
          r = ran3(localseed)
          xmc = r * dx_cell + xmin
          r = ran3(localseed)
          ymc = r * dx_cell + ymin
          r = ran3(localseed)
          zmc = r * dx_cell + zmin
          ! test if point is in the sphere 
          dist_mc2 = 0.0d0
          dist_mc2 = (xmc-0.5d0)**2 + (ymc-0.5d0)**2 + (zmc-0.5d0)**2
          dist_mc = sqrt(dist_mc2)
          if (dist_mc < max_dist .and. dist_mc > min_dist) then ! point within sphere
             volfrac = volfrac + 1.0d0
             vx_mc = vx_mc + fix_vel * (xmc - 0.5d0) / max_dist
             vy_mc = vy_mc + fix_vel * (ymc - 0.5d0) / max_dist
             vz_mc = vz_mc + fix_vel * (zmc - 0.5d0) / max_dist
             dist_cell_mc = dist_cell_mc + dist_mc
          end if
          nMC = nMC + 1
       end do
    
       volfrac2 = volfrac / dble(nMC)
       
       if(volfrac .ne. 0) then
       ! Velocity
          vx_ideal = vx_mc / volfrac
          vy_ideal = vy_mc / volfrac
          vz_ideal = vz_mc / volfrac
          ! redefine dist_cell as mean of dist_mc
          dist_cell = dist_cell_mc / volfrac
       end if
       ! assign density, ndust
       !! 5.898d-14 = sigma_lya_0 (cm2)
       !! rho(r) for ct outflow rate assuming V(r) = V0 * (r/rmax)**vel_alpha
       rho_0    = (vel_alpha+1.d0) * fix_tauH / (5.898d-14 * (rmax*fix_box_size_cm)**(vel_alpha+2.) * ((rmin*fix_box_size_cm)**(-1.-vel_alpha)-(rmax*fix_box_size_cm)**(-1.-vel_alpha))) ! at / cm3
       rho      = rho_0 *  (rmax / dist_cell)**(vel_alpha+2.)
       nh_ideal = rho * volfrac2
       
       if(volfrac2 .gt. 1.0) then
          print*,volfrac2,volfrac,nMC
       endif
    end if
       
    ! Give all cells shell_temp (cells not in sphere won't matter for RT as dnh = 0...)
    dopwidth_ideal = sqrt((2.0d0*kb/mp)*fix_temp) 
    
    return
     
  end subroutine sphere_homogen_velgrad_ct_outflow_rate

!!+++++++++++++++++++++++++++++ SPHERE with Velocity gradient from Steidel10 ++++++++++++++++++++++++++++++++++++++

  subroutine sphere_homogen_steidel(vx_ideal,vy_ideal,vz_ideal,nh_ideal,dopwidth_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)

    implicit none
    
    integer(kind=4)                :: nMC
    real(kind=8)                   :: dx_cell     ! Cell size
    real(kind=8)                   :: volfrac,volfrac2
    real(kind=8)                   :: dist_cell,dist_cell_min,dist_cell_max,dist2,dist_mc2,dist_mc
    real(kind=8)                   :: vx_mc,vy_mc,vz_mc
    real(kind=8)                   :: xcell_ideal,ycell_ideal,zcell_ideal
    real(kind=8)                   :: xmin,xmax,xmc,ymin,ymax,ymc,zmin,zmax,zmc
    real(kind=8),intent(inout)     :: nh_ideal,vx_ideal,vy_ideal,vz_ideal,dopwidth_ideal
    real(kind=8)                   :: r
    integer(kind=4)                :: localseed
    real(kind=8)                   :: sphere_radius,rho
    real(kind=8)                   :: max_dist, min_dist, rmin, rmax, dist_cell_mc

    !! Here V=0 and rho=0 at r < min_dist
    
    ! define sphere radius
    sphere_radius = 0.45d0          ! code units
    max_dist      = sphere_radius
    min_dist      = 0.02d0

    rmin          = min_dist
    rmax          = max_dist
    
    rho  = rho_sphere_from_tau0_T_rs_Lbox(fix_tauH,fix_temp,sphere_radius)
    
    vx_mc = 0.0d0
    vy_mc = 0.0d0
    vz_mc = 0.0d0
    dist2 = 0.0d0
    dist_cell_mc = 0.0d0
    
    ! xcell, ycell and zcell are in frame with origin at bottom-left corner of box
    dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 + (zcell_ideal-0.5d0)**2  ! in frame with origin at center of box

    dist_cell = sqrt(dist2)
    dist_cell_max = dist_cell + dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    dist_cell_min = dist_cell - dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    if (dist_cell_min > max_dist .or. dist_cell_max < min_dist) then   ! cell completely out of sphere
       vx_ideal  = 0.0d0
       vy_ideal  = 0.0d0
       vz_ideal  = 0.0d0
       nh_ideal  = 0.0d0
    else ! cell partially or completely within sphere
       if (dist_cell_max > max_dist) then ! cell partially within sphere
          ! bounds of cell
          xmin = xcell_ideal - 0.5d0*dx_cell 
          xmax = xcell_ideal + 0.5d0*dx_cell
          ymin = ycell_ideal - 0.5d0*dx_cell 
          ymax = ycell_ideal + 0.5d0*dx_cell
          zmin = zcell_ideal - 0.5d0*dx_cell 
          zmax = zcell_ideal + 0.5d0*dx_cell
          
          ! use Monte Carlo to compute density of cell (i.e. fraction of its volume within the sphere)
          volfrac = 0.0d0
          nMC = 0
          do while(volfrac .lt. 1000.d0)
             if (volfrac .lt. 1.d0 .and. nMC .eq. 100) then
                exit
             end if
             r = ran3(localseed)
             xmc = r * dx_cell + xmin
             r = ran3(localseed)
             ymc = r * dx_cell + ymin
             r = ran3(localseed)
             zmc = r * dx_cell + zmin
             dist_mc2 = 0.0d0
             dist_mc2 = (xmc-0.5d0)**2 + (ymc-0.5d0)**2 + (zmc-0.5d0)**2
             dist_mc = sqrt(dist_mc2)
             if (dist_mc < max_dist .and. dist_mc > min_dist) then ! point within sphere
                volfrac = volfrac + 1.0d0
                vx_mc = vx_mc + fix_vel * sqrt((rmin**(1.-vel_alpha)-dist_mc**(1.-vel_alpha)) / (rmin**(1.-vel_alpha)-rmax**(1.-vel_alpha))) * (xmc - 0.5d0) / dist_mc
                vy_mc = vy_mc + fix_vel * sqrt((rmin**(1.-vel_alpha)-dist_mc**(1.-vel_alpha)) / (rmin**(1.-vel_alpha)-rmax**(1.-vel_alpha))) * (ymc - 0.5d0) / dist_mc
                vz_mc = vz_mc + fix_vel * sqrt((rmin**(1.-vel_alpha)-dist_mc**(1.-vel_alpha)) / (rmin**(1.-vel_alpha)-rmax**(1.-vel_alpha))) * (zmc - 0.5d0) / dist_mc
                dist_cell_mc = dist_cell_mc + dist_mc
             end if
             nMC = nMC + 1
          end do

          if(volfrac .ne. 0) then
             ! Velocity
             vx_ideal = vx_mc / volfrac
             vy_ideal = vy_mc / volfrac
             vz_ideal = vz_mc / volfrac
             ! redefine dist_cell as mean of dist_mc
             dist_cell = dist_cell_mc / volfrac
          end if
          
          volfrac2 = volfrac / dble(nMC)
          
       else ! cell completely within sphere
          vx_ideal = fix_vel * sqrt((rmin**(1.-vel_alpha)-dist_cell**(1.-vel_alpha)) / (rmin**(1.-vel_alpha)-rmax**(1.-vel_alpha))) * (xcell_ideal - 0.5d0) / dist_cell
          vy_ideal = fix_vel * sqrt((rmin**(1.-vel_alpha)-dist_cell**(1.-vel_alpha)) / (rmin**(1.-vel_alpha)-rmax**(1.-vel_alpha))) * (ycell_ideal - 0.5d0) / dist_cell
          vz_ideal = fix_vel * sqrt((rmin**(1.-vel_alpha)-dist_cell**(1.-vel_alpha)) / (rmin**(1.-vel_alpha)-rmax**(1.-vel_alpha))) * (zcell_ideal - 0.5d0) / dist_cell
          volfrac2 = 1.0
       end if
       
       ! assign density, ndust
       nh_ideal = rho * volfrac2
    end if
    
    dopwidth_ideal = sqrt((2.0d0*kb/mp)*fix_temp) 
    
    return
    
  end subroutine sphere_homogen_steidel

  !!+++++++++++++++++++++++++++++ HOMOGEN SPHERE + radiation pressure feedback ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sphere_homogen_velgrad_rad_pressure(vx_ideal,vy_ideal,vz_ideal,nh_ideal,dopwidth_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)

    implicit none
    
    integer(kind=4)                :: nMC
    real(kind=8)                   :: dx_cell     ! Cell size
    real(kind=8)                   :: volfrac,volfrac2
    real(kind=8)                   :: dist_cell,dist_cell_min,dist_cell_max,dist2,dist_mc2,dist_mc
    real(kind=8)                   :: vx_mc,vy_mc,vz_mc
    real(kind=8)                   :: xcell_ideal,ycell_ideal,zcell_ideal
    real(kind=8)                   :: xmin,xmax,xmc,ymin,ymax,ymc,zmin,zmax,zmc
    real(kind=8),intent(inout)     :: nh_ideal,vx_ideal,vy_ideal,vz_ideal,dopwidth_ideal
    real(kind=8)                   :: r
    integer(kind=4)                :: localseed
    real(kind=8)                   :: sphere_radius,rho
    real(kind=8)                   :: max_dist, min_dist, rmin, rmax, rg, test_vneg1, test_vneg2, dist_cell_mc

    ! define sphere radius
    sphere_radius = 0.45d0          ! code units
    max_dist      = sphere_radius
    min_dist      = 0.02d0

    rmin          = min_dist
    rmax          = max_dist
    rg            = 0.0651d0   ! = where V peaks
    rho  = rho_sphere_from_tau0_T_rs_Lbox(fix_tauH,fix_temp,sphere_radius)
    
    vx_mc = 0.0d0
    vy_mc = 0.0d0
    vz_mc = 0.0d0
    dist2 = 0.0d0
    dist_cell_mc = 0.0d0
    
    ! xcell, ycell and zcell are in frame with origin at bottom-left corner of box
    dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 + (zcell_ideal-0.5d0)**2  ! in frame with origin at center of box

    dist_cell = sqrt(dist2)
    dist_cell_max = dist_cell + dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    dist_cell_min = dist_cell - dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube

    test_vneg1 = rg * (1.d0 / rmin - 1.d0 / dist_cell_min) + log(rmin / dist_cell_min)
    test_vneg2 = rg * (1.d0 / rmin - 1.d0 / dist_cell_max) + log(rmin / dist_cell_max)

    !!if (test_vneg .le. 0.d0) then
    if (dist_cell_min > max_dist .or. dist_cell_max < min_dist .or. test_vneg1 .le. 0.d0 .or. test_vneg2 .le. 0.d0) then  ! cell not in sphere region
       vx_ideal  = 0.0d0
       vy_ideal  = 0.0d0
       vz_ideal  = 0.0d0
       nh_ideal  = 0.0d0
    else
       ! bounds of cell
       xmin = xcell_ideal - 0.5d0*dx_cell 
       xmax = xcell_ideal + 0.5d0*dx_cell
       ymin = ycell_ideal - 0.5d0*dx_cell 
       ymax = ycell_ideal + 0.5d0*dx_cell
       zmin = zcell_ideal - 0.5d0*dx_cell 
       zmax = zcell_ideal + 0.5d0*dx_cell
       
       ! use Monte Carlo to compute density of cell (i.e. fraction of its volume within the sphere)
       volfrac = 0.0d0
       nMC = 0
       do while(volfrac .lt. 1000.d0)
          if (volfrac .lt. 1.d0 .and. nMC .eq. 100) then
             exit
          end if
          r = ran3(localseed)
          xmc = r * dx_cell + xmin
          r = ran3(localseed)
          ymc = r * dx_cell + ymin
          r = ran3(localseed)
          zmc = r * dx_cell + zmin
          ! test if point is in the sphere 
          dist_mc2 = 0.0d0
          dist_mc2 = (xmc-0.5d0)**2 + (ymc-0.5d0)**2 + (zmc-0.5d0)**2
          dist_mc = sqrt(dist_mc2)
          if (dist_mc < max_dist .and. dist_mc > min_dist) then ! point within sphere
             volfrac = volfrac + 1.0d0
             !! V(r) from Prochaska+11 - a la Murray
             !print*,dist_mc,(1.d0 / rmin - 1.d0 / dist_mc),log(rmin / dist_mc)
             !test_vneg = rg * (1.d0 / rmin - 1.d0 / dist_mc) + log(rmin / dist_mc)
             !if (test_vneg .le. 0.d0) then
             !   vx_mc = vx_mc + 0.d0
             !   vy_mc = vy_mc + 0.d0
             !   vz_mc = vz_mc + 0.d0
             !else
                !print*,fix_vel,rg,rmin
                !print*,dist_mc,(1.d0 / rmin - 1.d0 / dist_mc),log(rmin / dist_mc),test_vneg,(xmc - 0.5d0)
             vx_mc = vx_mc + fix_vel * sqrt(rg * (1.d0 / rmin - 1.d0 / dist_mc) + log(rmin / dist_mc)) * (xmc - 0.5d0) / dist_mc
             vy_mc = vy_mc + fix_vel * sqrt(rg * (1.d0 / rmin - 1.d0 / dist_mc) + log(rmin / dist_mc)) * (ymc - 0.5d0) / dist_mc
             vz_mc = vz_mc + fix_vel * sqrt(rg * (1.d0 / rmin - 1.d0 / dist_mc) + log(rmin / dist_mc)) * (zmc - 0.5d0) / dist_mc
             dist_cell_mc = dist_cell_mc + dist_mc
             !end if
          end if
          nMC = nMC + 1
       end do
    
       volfrac2 = volfrac / dble(nMC)
       
       if(volfrac .ne. 0) then
       ! Velocity
          vx_ideal = vx_mc / volfrac
          vy_ideal = vy_mc / volfrac
          vz_ideal = vz_mc / volfrac
          ! redefine dist_cell as mean of dist_mc
          dist_cell = dist_cell_mc / volfrac
       end if
       ! assign density, ndust
       nh_ideal = rho * volfrac2
       if(volfrac2 .gt. 1.0) then
          print*,volfrac2,volfrac,nMC
       endif
    end if
       
    ! Give all cells shell_temp (cells not in sphere won't matter for RT as dnh = 0...)
    dopwidth_ideal = sqrt((2.0d0*kb/mp)*fix_temp) 
    
    return
     
  end subroutine sphere_homogen_velgrad_rad_pressure

  
  !!+++++++++++++++++++++++++++++ SPHERE with density gradient & fixed velocity ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sphere_densgrad_velfix(vx_ideal,vy_ideal,vz_ideal,nh_ideal,dopwidth_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)

    implicit none

    !! Use only for fix_vel = 0 for now 
    
    integer(kind=4)                :: nMC
    real(kind=8)                   :: dx_cell     ! Cell size
    real(kind=8)                   :: volfrac,volfrac2
    real(kind=8)                   :: dist_cell,dist_cell_min,dist_cell_max,dist2,dist_mc2,dist_mc
    real(kind=8)                   :: vx_mc,vy_mc,vz_mc
    real(kind=8)                   :: xcell_ideal,ycell_ideal,zcell_ideal
    real(kind=8)                   :: xmin,xmax,xmc,ymin,ymax,ymc,zmin,zmax,zmc
    real(kind=8),intent(inout)     :: nh_ideal,vx_ideal,vy_ideal,vz_ideal,dopwidth_ideal
    real(kind=8)                   :: r
    integer(kind=4)                :: localseed
    real(kind=8)                   :: sphere_radius,rho
    real(kind=8)                   :: max_dist, min_dist, rmin, rmax, dist_cell_mc
   ! real(kind=8)                   :: dens_alpha

    ! define sphere radius
    sphere_radius = 0.45d0          ! code units
    max_dist      = sphere_radius
    min_dist      = 0.02d0

    vx_mc = 0.0d0
    vy_mc = 0.0d0
    vz_mc = 0.0d0
    dist2 = 0.0d0
    dist_cell_mc = 0.0d0
    
    ! xcell, ycell and zcell are in frame with origin at bottom-left corner of box
    dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 + (zcell_ideal-0.5d0)**2  ! in frame with origin at center of box

    dist_cell = sqrt(dist2)

    dist_cell_max = dist_cell + dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    dist_cell_min = dist_cell - dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    if (dist_cell_min > max_dist .or. dist_cell_max < min_dist) then  ! cell completely out of sphere
       vx_ideal  = 0.0d0
       vy_ideal  = 0.0d0
       vz_ideal  = 0.0d0
       nh_ideal  = 0.0d0
    else  ! cell partially or completely within sphere
       if (dist_cell_max > max_dist .and. dist_cell_min < max_dist) then ! cell partially within sphere
          ! bounds of cell
          xmin = xcell_ideal - 0.5d0*dx_cell 
          xmax = xcell_ideal + 0.5d0*dx_cell
          ymin = ycell_ideal - 0.5d0*dx_cell 
          ymax = ycell_ideal + 0.5d0*dx_cell
          zmin = zcell_ideal - 0.5d0*dx_cell 
          zmax = zcell_ideal + 0.5d0*dx_cell
          
          ! use Monte Carlo to compute density of cell (i.e. fraction of its volume within the sphere)
          volfrac = 0.0d0
          nMC = 0
          do while(volfrac .lt. 1000.d0)
             if (volfrac .lt. 1.d0 .and. nMC .eq. 100) then
                exit
             end if
             r = ran3(localseed)
             xmc = r * dx_cell + xmin
             r = ran3(localseed)
             ymc = r * dx_cell + ymin
             r = ran3(localseed)
             zmc = r * dx_cell + zmin
             dist_mc2 = 0.0d0
             dist_mc2 = (xmc-0.5d0)**2 + (ymc-0.5d0)**2 + (zmc-0.5d0)**2
             dist_mc = sqrt(dist_mc2)
             if (dist_mc < max_dist .and. dist_mc > min_dist) then ! point within sphere
                volfrac = volfrac + 1.0d0
                vx_mc = vx_mc + fix_vel * (xmc - 0.5d0) / dist_mc
                vy_mc = vy_mc + fix_vel * (ymc - 0.5d0) / dist_mc
                vz_mc = vz_mc + fix_vel * (zmc - 0.5d0) / dist_mc
                dist_cell_mc = dist_cell_mc + dist_mc
             end if
             nMC = nMC + 1
          end do

          if(volfrac .ne. 0) then
             vx_ideal = vx_mc / volfrac
             vy_ideal = vy_mc / volfrac
             vz_ideal = vz_mc / volfrac
             ! redefine dist_cell as mean of dist_mc
             dist_cell = dist_cell_mc / volfrac
          end if
          
          volfrac2 = volfrac / dble(nMC)
          
       else if (dist_cell_max < max_dist .and. dist_cell_min > min_dist) then ! cell completely within sphere
          vx_ideal = fix_vel * (xcell_ideal-0.5d0) / dist_cell
          vy_ideal = fix_vel * (ycell_ideal-0.5d0) / dist_cell
          vz_ideal = fix_vel * (zcell_ideal-0.5d0) / dist_cell
          volfrac2 = 1.0
       end if
       
       !! assign density, ndust
       !! rho(r) = rho0 * (rmax/r)**2: comes from ct_outflow_rate for vel_alpha = 0 and fix_vel = 0
       rho_0    = fix_tauH / (5.898d-14 * (rmax*fix_box_size_cm)**2. * ((rmin*fix_box_size_cm)**(-1.)-(rmax*fix_box_size_cm)**(-1.))) ! at / cm3
       rho      = rho_0 *  (rmax / dist_cell)**(2.)
       nh_ideal = rho * volfrac2
       
    end if
    
    ! Give all cells shell_temp (cells not in sphere won't matter for RT as dnh = 0...)
    dopwidth_ideal = sqrt((2.0d0*kb/mp)*fix_temp) 
    
    return
    
  end subroutine sphere_densgrad_velfix
  
   !!+++++++++++++++++++++++++++++ SPHERE with Velocity & density gradient ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sphere_densgrad_velgrad(vx_ideal,vy_ideal,vz_ideal,nh_ideal,dopwidth_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)

    implicit none
    
    integer(kind=4)                :: nMC
    real(kind=8)                   :: dx_cell     ! Cell size
    real(kind=8)                   :: volfrac,volfrac2
    real(kind=8)                   :: dist_cell,dist_cell_min,dist_cell_max,dist2,dist_mc2,dist_mc
    real(kind=8)                   :: vx_mc,vy_mc,vz_mc
    real(kind=8)                   :: xcell_ideal,ycell_ideal,zcell_ideal
    real(kind=8)                   :: xmin,xmax,xmc,ymin,ymax,ymc,zmin,zmax,zmc
    real(kind=8),intent(inout)     :: nh_ideal,vx_ideal,vy_ideal,vz_ideal,dopwidth_ideal
    real(kind=8)                   :: r
    integer(kind=4)                :: localseed
    real(kind=8)                   :: sphere_radius
    real(kind=8)                   :: max_dist,min_dist, dist_cell_mc
    real(kind=8)                   :: vel_prof!, vel_alpha, dens_alpha

    ! define sphere radius
    sphere_radius = 0.45d0          ! code units
    max_dist      = sphere_radius
    min_dist      = 0.d0
    
    vx_mc = 0.0d0
    vy_mc = 0.0d0
    vz_mc = 0.0d0
    dist2 = 0.0d0
    dist_cell_mc = 0.0d0
    ! xcell, ycell and zcell are in frame with origin at bottom-left corner of box
    dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 + (zcell_ideal-0.5d0)**2  ! in frame with origin at center of box

    dist_cell = sqrt(dist2)
    if (dist_cell .le. 1.01d0*max_dist/10.d0) then
       vel_prof = 0.d0
    else
       vel_prof = (1.d0 - max_dist / (10.d0 * dist_cell))**vel_alpha
    end if

    dist_cell_max = dist_cell + dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    dist_cell_min = dist_cell - dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
    if (dist_cell_min > max_dist .or. dist_cell_max < min_dist) then  ! cell not in sphere region
       vx_ideal  = 0.0d0
       vy_ideal  = 0.0d0
       vz_ideal  = 0.0d0
    else
       ! bounds of cell
       xmin = xcell_ideal - 0.5d0*dx_cell 
       xmax = xcell_ideal + 0.5d0*dx_cell
       ymin = ycell_ideal - 0.5d0*dx_cell 
       ymax = ycell_ideal + 0.5d0*dx_cell
       zmin = zcell_ideal - 0.5d0*dx_cell 
       zmax = zcell_ideal + 0.5d0*dx_cell
       
       ! use Monte Carlo to compute density of cell (i.e. fraction of its volume within the sphere)
       volfrac = 0.0d0
       nMC = 0
       do while(volfrac .lt. 1000.d0)
          if (volfrac .lt. 1.d0 .and. nMC .eq. 100) then
             exit
          end if
          r = ran3(localseed)
          xmc = r * dx_cell + xmin
          r = ran3(localseed)
          ymc = r * dx_cell + ymin
          r = ran3(localseed)
          zmc = r * dx_cell + zmin
          ! test if point is in the sphere 
          dist_mc2 = 0.0d0
          dist_mc2 = (xmc-0.5d0)**2 + (ymc-0.5d0)**2 + (zmc-0.5d0)**2
          dist_mc = sqrt(dist_mc2)
          if (dist_mc < max_dist .and. dist_mc > min_dist) then ! point within sphere
             volfrac = volfrac + 1.0d0
             vx_mc = vx_mc + fix_vel * vel_prof * (xmc - 0.5d0) / dist_mc
             vy_mc = vy_mc + fix_vel * vel_prof * (ymc - 0.5d0) / dist_mc
             vz_mc = vz_mc + fix_vel * vel_prof * (zmc - 0.5d0) / dist_mc
             dist_cell_mc = dist_cell_mc + dist_mc
             ! Actually, I could do the same for nh instead of computing nh at mid-cell... nh = nh + nh(dist_mc) then divide by volfrac and multiply by volfrac2 (?)
          end if
          nMC = nMC + 1
       end do
    
       volfrac2 = volfrac / dble(nMC)
       
       if(volfrac .ne. 0) then
       ! Velocity
          vx_ideal = vx_mc / volfrac
          vy_ideal = vy_mc / volfrac
          vz_ideal = vz_mc / volfrac
          ! redefine dist_cell as mean of dist_mc
          dist_cell = dist_cell_mc / volfrac
       end if

       ! assign density, ndust

       if (fix_dens_grad .eq. 'powerlaw') then
          ! Avoids nh to diverge at r=0 => if R < 5 pc => nh(r) = nh(5pc)
          ! nh(5pc)=354317.0d0 cm^-3
          ! Rmax=0.45 <=> 630 pc   => 5pc <=> 0.00357143
          ! Use 20pc actually:  nh(20pc)=127.529 cm^-3 - 20pc <=> 0.0142857
          ! Use 63pc actually:  nh(63pc)=0.18 cm^-3 - 63pc <=> 0.045                                                                          
          if (dist_cell .le. 0.045) then
             nh_ideal = 0.18  ! cm^-3
          else
             nh_ideal = rho_0 * (10.d0 / max_dist * dist_cell)**dens_alpha ! cm^-3
          end if
       end if
       !if (nh_ideal .gt. 1.) then
       !   print*,'nh = ',nh_ideal
       !end if

       nh_ideal = nh_ideal * volfrac2
       
       if(volfrac2 .gt. 1.0) then
          print*,volfrac2,volfrac,nMC
       endif
    end if
       
    ! Give all cells shell_temp (cells not in sphere won't matter for RT as dnh = 0...)
    dopwidth_ideal = sqrt((2.0d0*kb/mp)*fix_temp) 
    
    return
     
  end subroutine sphere_densgrad_velgrad


  
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function rho_sphere_from_tau0_T_rs_Lbox(tauH,temp,rs)
    
    implicit none 
    
    real(kind=8),intent(in)     :: tauH,temp,rs
    real(kind=8)                :: rho_sphere_from_tau0_T_rs_Lbox!,fix_box_size_cmBIS
    
   ! fix_box_size_cmBIS = 1.d22 ! cm -> later: read this from pram file
    
    ! fix_box_size_cm is fixed by user in param file
    rho_sphere_from_tau0_T_rs_Lbox = 1.d20 * tauH * (temp/1.d4)**0.5 / (5.898d6 * rs * fix_box_size_cm) ! cm^-3
    
    return
    
  end function rho_sphere_from_tau0_T_rs_Lbox
  

  function rho_shell_from_tau0_T_rs_Lbox(tauH,temp,dr)
    
    implicit none 
    
    real(kind=8),intent(in)     :: tauH,temp,dr
    real(kind=8)                :: rho_shell_from_tau0_T_rs_Lbox!,fix_box_size_cmBIS

    rho_shell_from_tau0_T_rs_Lbox = 1.d20 * tauH * (temp/1.d4)**0.5 / (5.898d6 * dr * fix_box_size_cm) ! cm^-3

    return
    
  end function rho_shell_from_tau0_T_rs_Lbox

  
  subroutine read_overwrite_params(overwrite_model)

    implicit none
    
    character(2000) :: param_file
    character(1000) :: line,name,value
    integer(kind=4) :: err,i
    logical         :: section_present
    character(50)   :: overwrite_model != 'uniform'

    call get_command_argument(1, param_file)

    section_present = .false.
    open(unit=10,file=trim(param_file),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:50) == '[gas_composition]') then
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
          case ('overwrite_model')
             read(value,*) overwrite_model
          case ('verbose')
             read(value,*) verbose
          end select

          !! uniform model
          if (overwrite_model .eq. 'uniform') then
             select case (trim(name))
             case ('fix_nhi')
                read(value,*) fix_nhi
             case ('fix_vth')
                read(value,*) fix_vth
             case ('fix_vel')
                read(value,*) fix_vel
             case ('fix_box_size_cm')
                read(value,*) fix_box_size_cm
             end select
          end if
          !! sphere_homogen_velfix model
          if (overwrite_model .eq. 'sphere_homogen_velfix') then
             select case (trim(name))
             case ('fix_tauH')
                read(value,*) fix_tauH
             case ('fix_temp')
                read(value,*) fix_temp
             case ('fix_taudust')
                read(value,*) fix_taudust
             case ('fix_vel')
                read(value,*) fix_vel
             case ('fix_box_size_cm')
                read(value,*) fix_box_size_cm
             end select
          end if
           !! shell_homogen_velfix model
          if (overwrite_model .eq. 'shell_homogen_velfix') then
             select case (trim(name))
             case ('fix_tauH')
                read(value,*) fix_tauH
             case ('fix_temp')
                read(value,*) fix_temp
             case ('fix_taudust')
                read(value,*) fix_taudust
             case ('fix_vel')
                read(value,*) fix_vel
             case ('fix_box_size_cm')
                read(value,*) fix_box_size_cm
             end select
          end if
          !! sphere_gradient model
          if (overwrite_model .eq. 'sphere_gradient') then
             select case (trim(name))
             case ('fix_tauH')
                read(value,*) fix_tauH
             case ('fix_temp')
                read(value,*) fix_temp
             case ('fix_taudust')
                read(value,*) fix_taudust
             case ('fix_vel')
                read(value,*) fix_vel
             case ('fix_box_size_cm')
                read(value,*) fix_box_size_cm
             case ('fix_dens_grad')
                read(value,*) fix_dens_grad
             end select
          end if
          !! sphere_homog_velgrad model
          if (overwrite_model .eq. 'sphere_homog_velgrad') then
             select case (trim(name))
             case ('fix_tauH')
                read(value,*) fix_tauH
             case ('fix_temp')
                read(value,*) fix_temp
             case ('fix_taudust')
                read(value,*) fix_taudust
             case ('fix_vel')
                read(value,*) fix_vel
             case ('fix_box_size_cm')
                read(value,*) fix_box_size_cm
             end select
          end if
          !! sphere_homogen_steidel
          if (overwrite_model .eq. 'sphere_homogen_steidel') then
             select case (trim(name))
             case ('fix_tauH')
                read(value,*) fix_tauH
             case ('vel_alpha')
                read(value,*) vel_alpha
             case ('fix_temp')
                read(value,*) fix_temp
             case ('fix_taudust')
                read(value,*) fix_taudust
             case ('fix_vel')
                read(value,*) fix_vel
             case ('fix_box_size_cm')
                read(value,*) fix_box_size_cm
             end select
          end if
          !! sphere_homogen_velgrad_ct_outflow_rate
          if (overwrite_model .eq. 'sphere_homogen_velgrad_ct_outflow_rate') then
             select case (trim(name))
             case ('fix_tauH')
                read(value,*) fix_tauH
             case ('vel_alpha')
                read(value,*) vel_alpha
             case ('fix_temp')
                read(value,*) fix_temp
             case ('fix_taudust')
                read(value,*) fix_taudust
             case ('fix_vel')
                read(value,*) fix_vel
             case ('fix_box_size_cm')
                read(value,*) fix_box_size_cm
             end select
          end if
          !! sphere_homogen_velgrad_rad_pressure
          if (overwrite_model .eq. 'sphere_homogen_velgrad_rad_pressure') then
             select case (trim(name))
             case ('fix_tauH')
                read(value,*) fix_tauH
             case ('fix_temp')
                read(value,*) fix_temp
             case ('fix_taudust')
                read(value,*) fix_taudust
             case ('fix_vel')
                read(value,*) fix_vel
             case ('fix_box_size_cm')
                read(value,*) fix_box_size_cm
             end select
          end if
          !! sphere_densgrad_velgrad model
          if (overwrite_model .eq. 'sphere_densgrad_velgrad') then
             select case (trim(name))
             case ('fix_dens_grad')
                read(value,*) fix_dens_grad
             case ('vel_alpha')
                read(value,*) vel_alpha
             case ('dens_alpha')
                read(value,*) dens_alpha
             case ('rho_0')
                read(value,*) rho_0
             case ('fix_temp')
                read(value,*) fix_temp
             case ('fix_taudust')
                read(value,*) fix_taudust
             case ('fix_vel')
                read(value,*) fix_vel
             case ('fix_box_size_cm')
                read(value,*) fix_box_size_cm
             end select
          end if
          !! sphere_scarlata15 model
          if (overwrite_model .eq. 'sphere_scarlata15') then
             select case (trim(name))
             case ('rsf')
                read(value,*) rsf
             case ('rw')
                read(value,*) rw
             case ('fix_temp')
                read(value,*) fix_temp
             case ('fix_taudust')
                read(value,*) fix_taudust
             case ('vel_gamma')
                read(value,*) vel_gamma
             case ('rho_rsf')
                read(value,*) rho_rsf
             case ('fix_vel')
                read(value,*) fix_vel
             case ('fix_box_size_cm')
                read(value,*) fix_box_size_cm
             end select
          end if
          !! disc_thin model
          if (overwrite_model .eq. 'disc_thin') then
             select case (trim(name))
             case ('fix_mtot_hi')
                read(value,*) fix_mtot_hi
             case ('fix_rd_disc')
                read(value,*) fix_rd_disc
             case ('fix_temp')
                read(value,*) fix_temp
             case ('fix_taudust')
                read(value,*) fix_taudust
             case ('fix_vel')
                read(value,*) fix_vel
             case ('fix_box_size_cm')
                read(value,*) fix_box_size_cm
             end select
          end if
          !! disc_thick model
          if (overwrite_model .eq. 'disc_thick') then
             select case (trim(name))
             case ('fix_mtot_hi')
                read(value,*) fix_mtot_hi
             case ('fix_rd_disc')
                read(value,*) fix_rd_disc
             case ('fix_zd_disc')
                read(value,*) fix_zd_disc
             case ('fix_temp')
                read(value,*) fix_temp
             case ('fix_taudust')
                read(value,*) fix_taudust
             case ('fix_vel')
                read(value,*) fix_vel
             case ('fix_box_size_cm')
                read(value,*) fix_box_size_cm
             end select
          end if
          !! sphere_densgrad_velfix model
          if (overwrite_model .eq. 'sphere_densgrad_velfix') then
             select case (trim(name))
!!$             case ('fix_dens_grad')
!!$                read(value,*) fix_dens_grad
!!$             case ('dens_alpha')
!!$                read(value,*) dens_alpha
!!$             case ('rho_0')
!!$                read(value,*) rho_0
             case ('fix_tauH')
                read(value,*) fix_tauH
             case ('fix_temp')
                read(value,*) fix_temp
             case ('fix_taudust')
                read(value,*) fix_taudust
             case ('fix_vel')
                read(value,*) fix_vel
             case ('fix_box_size_cm')
                read(value,*) fix_box_size_cm
             end select
          end if
       end do
    end if
    close(10)
    
  end subroutine read_overwrite_params


  subroutine print_overwrite_params(unit)
    
    implicit none

    integer(kind=4),optional,intent(in)  :: unit
    character(50)                        :: overwrite_model
    
    call read_overwrite_params(overwrite_model)

    if (present(unit)) then
       write(unit,'(a,a)')         ' overwrite_model     = ',overwrite_model
       select case (overwrite_model)
       case('uniform')
          write(unit,'(a,ES10.3)') '  fix_nhi            = ',fix_nhi
          write(unit,'(a,ES10.3)') '  fix_vth            = ',fix_vth
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
       case('sphere_homogen_velfix')
          write(unit,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(unit,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('shell_homogen_velfix')
          write(unit,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(unit,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('sphere_gradient')
          write(unit,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(unit,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
          write(unit,'(a,a)')      '  fix_dens_grad      = ',fix_dens_grad
       case('sphere_homog_velgrad')
          write(unit,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(unit,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('sphere_homogen_steidel')
          write(unit,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(unit,'(a,ES10.3)') '  vel_alpha          = ',vel_alpha
          write(unit,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('sphere_homogen_velgrad_ct_outflow_rate')
          write(unit,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(unit,'(a,ES10.3)') '  vel_alpha          = ',vel_alpha
          write(unit,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('sphere_homogen_velgrad_rad_pressure')
          write(unit,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(unit,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('sphere_densgrad_velgrad')
          write(unit,'(a,ES10.3)') '  fix_dens_grad      = ',fix_dens_grad
          write(unit,'(a,ES10.3)') '  vel_alpha          = ',vel_alpha
          write(unit,'(a,ES10.3)') '  dens_alpha         = ',dens_alpha
          write(unit,'(a,ES10.3)') '  rho_0              = ',rho_0
          write(unit,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('sphere_densgrad_velfix')
          write(unit,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(unit,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('disc_thin')
          write(unit,'(a,ES10.3)') '  fix_mtot_hi        = ',fix_mtot_hi
          write(unit,'(a,ES10.3)') '  fix_rd_disc        = ',fix_rd_disc
          write(unit,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('disc_thick')
          write(unit,'(a,ES10.3)') '  fix_mtot_hi        = ',fix_mtot_hi
          write(unit,'(a,ES10.3)') '  fix_rd_disc        = ',fix_rd_disc
          write(unit,'(a,ES10.3)') '  fix_zd_disc        = ',fix_zd_disc
          write(unit,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(unit,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('sphere_scarlata15')
          write(unit,'(a,ES10.3)') '  rsf                = ',rsf
          write(unit,'(a,ES10.3)') '  rw                 = ',rw
          write(unit,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(unit,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
          write(unit,'(a,ES10.3)') '  vel_gamma          = ',vel_gamma
          write(unit,'(a,ES10.3)') '  rho_rsf            = ',rho_rsf
          write(unit,'(a,ES10.3)') '  fix_vel            = ',fix_vel
       end select
       write(unit,'(a,ES10.3)')    '  fix_box_size_cm    = ',fix_box_size_cm
    else
       write(*,'(a,a)')         '  overwrite_model    = ',overwrite_model
       select case (overwrite_model)
       case('uniform')
          write(*,'(a,ES10.3)') '  fix_nhi            = ',fix_nhi
          write(*,'(a,ES10.3)') '  fix_vth            = ',fix_vth
          write(*,'(a,ES10.3)') '  fix_vel            = ',fix_vel
       case('sphere_homogen_velfix')
          write(*,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(*,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(*,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('shell_homogen_velfix')
          write(*,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(*,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(*,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('sphere_homog_velgrad')
          write(*,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(*,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(*,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('sphere_homogen_steidel')
          write(*,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(*,'(a,ES10.3)') '  vel_alpha          = ',vel_alpha
          write(*,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(*,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('sphere_homogen_velgrad_ct_outflow_rate')
          write(*,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(*,'(a,ES10.3)') '  vel_alpha          = ',vel_alpha
          write(*,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(*,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('sphere_homogen_velgrad_rad_pressure')
          write(*,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(*,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(*,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('sphere_densgrad_velgrad')
          write(*,'(a,ES10.3)') '  fix_dens_grad      = ',fix_dens_grad
          write(*,'(a,ES10.3)') '  vel_alpha          = ',vel_alpha
          write(*,'(a,ES10.3)') '  dens_alpha         = ',dens_alpha
          write(*,'(a,ES10.3)') '  rho_0              = ',rho_0
          write(*,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(*,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('sphere_gradient')
          write(*,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(*,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(*,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
          write(*,'(a,a)')      '  fix_dens_grad      = ',fix_dens_grad
       case('disc_thin')
          write(*,'(a,ES10.3)') '  fix_mtot_hi        = ',fix_mtot_hi
          write(*,'(a,ES10.3)') '  fix_rd_disc        = ',fix_rd_disc
          write(*,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(*,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('disc_thick')
          write(*,'(a,ES10.3)') '  fix_mtot_hi        = ',fix_mtot_hi
          write(*,'(a,ES10.3)') '  fix_rd_disc        = ',fix_rd_disc
          write(*,'(a,ES10.3)') '  fix_zd_disc        = ',fix_zd_disc
          write(*,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(*,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       case('sphere_densgrad_velfix')
          write(*,'(a,ES10.3)') '  fix_tauH           = ',fix_tauH
          write(*,'(a,ES10.3)') '  fix_temp           = ',fix_temp
          write(*,'(a,ES10.3)') '  fix_vel            = ',fix_vel
          write(*,'(a,ES10.3)') '  fix_taudust        = ',fix_taudust
       end select
       write(*,'(a,ES10.3)') '  fix_box_size_cm    = ',fix_box_size_cm
    end if
       
  end subroutine print_overwrite_params

  
end module module_idealised_models
