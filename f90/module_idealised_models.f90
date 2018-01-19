module module_idealised_models

  use module_random
  use module_constants, only : kb, mp, pi
  
  implicit none

  ! In this module, the species is unknown

  private
  real(kind=8)             :: ngas_norm           = -999.0d0         ! Gas density profile normalisation [cm-3]
  real(kind=8)             :: coldens_norm        = -999.0d0         ! Gas COLUMN density [cm-2]
  real(kind=8)             :: ngas_slope          = -999.0d0         ! slope of the density profile
  real(kind=8)             :: Vgas_norm           = -999.0d0         ! V0, Vmax, Vcirc etc [cm/s]
  real(kind=8)             :: Vgas_slope          = -999.0d0         ! slope of the velocity profile
  real(kind=8)             :: Tgas_norm           = -999.0d0         ! Gas temperature [K]
  real(kind=8)             :: ndust_norm          = -999.0d0         ! Dust density profile normalisation [cm-3]
  real(kind=8)             :: taudust_norm        = -999.0d0         ! (abs) dust opacity
  real(kind=8)             :: r_min               = -999.0d0         ! full box size [cm]
  real(kind=8)             :: r_max               = -999.0d0         ! full box size [cm]
  real(kind=8)             :: disc_rd             = -999.0d0         ! Disc radial scalelength
  real(kind=8)             :: disc_zd             = -999.0d0         ! Disc height, i.e. vertical scalelength
  real(kind=8),public      :: box_size_IM_cm      = -999.0d0         ! full box size [cm]

  character(100)           :: idealmodel_type     = 'my_model'       ! name of the idealised model
  logical                  :: MCsampling          = .false.
  logical                  :: verbose             = .false.

  
  ! public functions:
  public :: read_IdealisedModels_params, print_IdealisedModels_params, compute_idealised_gas, shell_V_rho_gradient, shellcone_V_rho_gradient, shell_chisholm, shell_V_rho_gradient_steady
  
  
contains


  subroutine compute_idealised_gas(n_dust,n_gas,b_param,v_leaf,x_leaf,leaf_level,nleaf)

    implicit none
    
    integer(kind=4),intent(in)            :: nleaf
    real(kind=8),intent(in)               :: x_leaf(nleaf,3)
    integer(kind=4),intent(in)            :: leaf_level(nleaf)
    real(kind=8),intent(inout)            :: v_leaf(3,nleaf)
    real(kind=8),intent(inout)            :: n_dust(nleaf),n_gas(nleaf),b_param(nleaf)  ! b_param is dopwidth
    real(kind=8)                          :: dx_cell
    integer(kind=4)                       :: ileaf

    dx_cell = 0.5d0**leaf_level(1) ! assumes all leaf levels are identical (regular grid)
        
    select case (idealmodel_type)
    case('shell_V_rho_gradient')
       do ileaf=1,nleaf
          call shell_V_rho_gradient(n_dust(ileaf),n_gas(ileaf),b_param(ileaf),v_leaf(1,ileaf),v_leaf(2,ileaf),v_leaf(3,ileaf),x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do
    case('shell_V_rho_gradient_steady')
       do ileaf=1,nleaf
          call shell_V_rho_gradient_steady(n_dust(ileaf),n_gas(ileaf),b_param(ileaf),v_leaf(1,ileaf),v_leaf(2,ileaf),v_leaf(3,ileaf),x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do
    case('shellcone_V_rho_gradient')
       do ileaf=1,nleaf
          call shellcone_V_rho_gradient(n_dust(ileaf),n_gas(ileaf),b_param(ileaf),v_leaf(1,ileaf),v_leaf(2,ileaf),v_leaf(3,ileaf),x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do
    case('shell_chisholm')
       do ileaf=1,nleaf
          call shell_chisholm(n_dust(ileaf),n_gas(ileaf),b_param(ileaf),v_leaf(1,ileaf),v_leaf(2,ileaf),v_leaf(3,ileaf),x_leaf(ileaf,1),x_leaf(ileaf,2),x_leaf(ileaf,3),dx_cell)
       end do
    end select
    
    return
    
  end subroutine compute_idealised_gas



  !!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ BEGINNING MODELS //////////////////////////////////////////////////////////////////////////

  !+++++++++++++++++++++++++++++++++++++++++++++++ SHELL with velocity gradient and density gradient (power-laws) +++++++++++++++++++++++++++++++++++++++++++++

  subroutine shell_V_rho_gradient(ndust_ideal,ngas_ideal,bparam_ideal,vx_ideal,vy_ideal,vz_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)
    
    implicit none

    ! Shell (or sphere if r_min=0) with "power-law" V and rho profiles 
    ! point source at center and medium transparent at R < r_min and at R > r_max

    ! Declare arguments
    real(kind=8)                          :: dx_cell
    real(kind=8),intent(inout)            :: ndust_ideal,ngas_ideal,bparam_ideal
    real(kind=8),intent(inout)            :: vx_ideal,vy_ideal,vz_ideal
    real(kind=8),intent(in)               :: xcell_ideal,ycell_ideal,zcell_ideal    
    real(kind=8)                          :: volfrac2
    real(kind=8)                          :: dist_cell,dist2,dist_cell_min,dist_cell_max
    integer(kind=4)                       :: missed_cell
    real(kind=8)                          :: n0

    
    vx_ideal    = 0.0d0
    vy_ideal    = 0.0d0
    vz_ideal    = 0.0d0
    ngas_ideal  = 0.0d0
    ndust_ideal = 0.0d0
    volfrac2    = 1.0d0
    
    missed_cell = 1 ! =1 if cell doesn't satisfy and if statements... should not happen!
    
    ! xcell, ycell and zcell are in frame with origin at bottom-left corner of box
    dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 + (zcell_ideal-0.5d0)**2  ! in frame with origin at center of box
    dist_cell = sqrt(dist2)

    if (coldens_norm .gt. 1.0d0) then 
       ! if gas norm set as column density in param file.... OK for model with n~r^-2 ONLY !!!
       n0 = coldens_norm / (r_min * box_size_IM_cm * (1.0d0 - r_min / r_max))  ! cm-3
    else
       ! I use the std density ngas_norm
       n0 = ngas_norm 
    end if
   
    if (MCsampling) then
       dist_cell_max = dist_cell + dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
       dist_cell_min = dist_cell - dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube

       if ((dist_cell_max > r_max .and. dist_cell_min < r_max) .or. (dist_cell_max > r_min .and. dist_cell_min < r_min)) then ! cell partially within shell
          volfrac2 = mc_sampling_shell(xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)
          missed_cell = 0
       end if
       if (dist_cell_max < r_max .and. dist_cell_min > r_min) then ! cell completely within shell
          volfrac2 = 1.0
          missed_cell = 0
       end if
       
       ngas_ideal  = n0 *  (r_min / dist_cell)**(ngas_slope) * volfrac2 ! ngas_slope  = +2 for P+11 fiducial model
       ndust_ideal = ndust_norm * ngas_ideal / n0

       vx_ideal = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (xcell_ideal - 0.5d0)
       vy_ideal = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (ycell_ideal - 0.5d0)
       vz_ideal = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (zcell_ideal - 0.5d0)
       
       if (dist_cell_min > r_max .or. dist_cell_max < r_min) then  ! cell completely out of shell or completely within r_min                                                       
          vx_ideal    = 0.0d0
          vy_ideal    = 0.0d0
          vz_ideal    = 0.0d0
          ngas_ideal  = 0.0d0
          ndust_ideal = 0.0d0
          missed_cell = 0
       end if
       
    else  !! Brut force: compare dist to cell center against Rmin/Rmax 
       if (dist_cell < r_max .and. dist_cell > r_min) then ! cell completely within sphere
!!$          vx_ideal    = Vgas_norm * (xcell_ideal-0.5d0) / r_max
!!$          vy_ideal    = Vgas_norm * (ycell_ideal-0.5d0) / r_max
!!$          vz_ideal    = Vgas_norm * (zcell_ideal-0.5d0) / r_max
          vx_ideal    = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (xcell_ideal - 0.5d0)
          vy_ideal    = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (ycell_ideal - 0.5d0)
          vz_ideal    = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (zcell_ideal - 0.5d0)
          
          ngas_ideal  = n0 *  (r_min / dist_cell)**(ngas_slope) ! ngas_slope  = +2 for P+11 fiducial model
          ndust_ideal = ndust_norm * ngas_ideal / n0
          missed_cell = 0
       else                                                ! cell completely out of sphere or completely within r_min               
          vx_ideal    = 0.0d0
          vy_ideal    = 0.0d0
          vz_ideal    = 0.0d0
          ngas_ideal  = 0.0d0
          ndust_ideal = 0.0d0
          missed_cell = 0
       end if
    end if
    
    bparam_ideal = sqrt(2.0d0*kb*Tgas_norm)

    if (missed_cell .eq. 1) then
       print*,'I missed a cell in module_idealised_models !'
       print*,dist_cell
       stop
    endif
    
    return
    
  end subroutine shell_V_rho_gradient


   !+++++++++++++++++++++++++++++++++++++++++++++++ SHELL with velocity gradient and density gradient (power-laws) for steady state +++++++++++++++++++++++++++++++++++++++++++++

  subroutine shell_V_rho_gradient_steady(ndust_ideal,ngas_ideal,bparam_ideal,vx_ideal,vy_ideal,vz_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)
    
    implicit none

    ! Shell (or sphere if r_min=0) with "power-law" V and rho profiles 
    ! point source at center and medium transparent at R < r_min and at R > r_max
    ! Assumes ct ejection rate: rho = const. / (R^2 V(r))
    ! Won't work for Vgas_slope=-1
    
    ! Declare arguments
    real(kind=8)                          :: dx_cell
    real(kind=8),intent(inout)            :: ndust_ideal,ngas_ideal,bparam_ideal
    real(kind=8),intent(inout)            :: vx_ideal,vy_ideal,vz_ideal
    real(kind=8),intent(in)               :: xcell_ideal,ycell_ideal,zcell_ideal    
    real(kind=8)                          :: volfrac2
    real(kind=8)                          :: dist_cell,dist2,dist_cell_min,dist_cell_max
    integer(kind=4)                       :: missed_cell
    real(kind=8)                          :: n0, coldens_dust, ndust_0

    
    vx_ideal    = 0.0d0
    vy_ideal    = 0.0d0
    vz_ideal    = 0.0d0
    ngas_ideal  = 0.0d0
    ndust_ideal = 0.0d0
    volfrac2    = 1.0d0
    
    missed_cell = 1 ! =1 if cell doesn't satisfy and if statements... should not happen!
    
    ! xcell, ycell and zcell are in frame with origin at bottom-left corner of box
    dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 + (zcell_ideal-0.5d0)**2  ! in frame with origin at center of box
    dist_cell = sqrt(dist2)

    if (coldens_norm .gt. 1.0d0) then 
       ! if gas norm set as column density in param file....
       n0 = coldens_norm * (Vgas_slope+1.0) / (r_min * box_size_IM_cm * (1.0d0 - (r_min / r_max)**(Vgas_slope+1.0)))
    else
       ! I use the std density ngas_norm
       n0 = ngas_norm 
    end if
    if (taudust_norm .gt. 0.0d0) then
       coldens_dust = taudust_norm ! assumes sigma_dust = 1 (should be taudust_norm/sigma_dust)
       ndust_0      = coldens_dust * (Vgas_slope+1.0) / (r_min * box_size_IM_cm * (1.0d0 - (r_min / r_max)**(Vgas_slope+1.0)))
    else
       ndust_0      = 0.0d0
    end if
    
    if (MCsampling) then
       dist_cell_max = dist_cell + dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
       dist_cell_min = dist_cell - dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube

       if ((dist_cell_max > r_max .and. dist_cell_min < r_max) .or. (dist_cell_max > r_min .and. dist_cell_min < r_min)) then ! cell partially within shell
          volfrac2 = mc_sampling_shell(xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)
          missed_cell = 0
       end if
       if (dist_cell_max < r_max .and. dist_cell_min > r_min) then ! cell completely within shell
          volfrac2 = 1.0
          missed_cell = 0
       end if
       
       ngas_ideal  = n0 *  (r_min / dist_cell)**(Vgas_slope+2.0) * volfrac2 
       ndust_ideal = ndust_0 / n0 * ngas_ideal

       vx_ideal = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (xcell_ideal - 0.5d0)
       vy_ideal = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (ycell_ideal - 0.5d0)
       vz_ideal = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (zcell_ideal - 0.5d0)
       
       if (dist_cell_min > r_max .or. dist_cell_max < r_min) then  ! cell completely out of shell or completely within r_min                                                       
          vx_ideal    = 0.0d0
          vy_ideal    = 0.0d0
          vz_ideal    = 0.0d0
          ngas_ideal  = 0.0d0
          ndust_ideal = 0.0d0
          missed_cell = 0
       end if
       
    else  !! Brut force: compare dist to cell center against Rmin/Rmax 
       if (dist_cell < r_max .and. dist_cell > r_min) then ! cell completely within sphere
!!$          vx_ideal    = Vgas_norm * (xcell_ideal-0.5d0) / r_max
!!$          vy_ideal    = Vgas_norm * (ycell_ideal-0.5d0) / r_max
!!$          vz_ideal    = Vgas_norm * (zcell_ideal-0.5d0) / r_max
          vx_ideal    = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (xcell_ideal - 0.5d0)
          vy_ideal    = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (ycell_ideal - 0.5d0)
          vz_ideal    = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (zcell_ideal - 0.5d0)

          ngas_ideal  = n0 *  (r_min / dist_cell)**(Vgas_slope+2.0) * volfrac2 
          ndust_ideal = ndust_0 / n0 * ngas_ideal
          missed_cell = 0
       else                                                ! cell completely out of sphere or completely within r_min               
          vx_ideal    = 0.0d0
          vy_ideal    = 0.0d0
          vz_ideal    = 0.0d0
          ngas_ideal  = 0.0d0
          ndust_ideal = 0.0d0
          missed_cell = 0
       end if
    end if
    
    bparam_ideal = sqrt(2.0d0*kb*Tgas_norm)

    if (missed_cell .eq. 1) then
       print*,'I missed a cell in module_idealised_models !'
       print*,dist_cell
       stop
    endif
    
    return
    
  end subroutine shell_V_rho_gradient_steady
  


    !+++++++++++++++++++++++++++++++++++++++++++++++ Cone with velocity gradient and density gradient (power-laws) +++++++++++++++++++++++++++++++++++++++++++++

  subroutine shellcone_V_rho_gradient(ndust_ideal,ngas_ideal,bparam_ideal,vx_ideal,vy_ideal,vz_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)
    
    implicit none

    ! Shell (or sphere if r_min=0) with "power-law" V and rho profiles 
    ! point source at center and medium transparent at R < r_min and at R > r_max

    ! Declare arguments
    real(kind=8)                          :: dx_cell
    real(kind=8),intent(inout)            :: ndust_ideal,ngas_ideal,bparam_ideal
    real(kind=8),intent(inout)            :: vx_ideal,vy_ideal,vz_ideal
    real(kind=8),intent(in)               :: xcell_ideal,ycell_ideal,zcell_ideal    
    real(kind=8)                          :: dist_cell,dist2,dist_cell_min,dist_cell_max
    integer(kind=4)                       :: missed_cell
    real(kind=8)                          :: n0,theta_coord
    real(kind=8)                          :: theta_cone

    theta_cone = 90.0d0 ! deg
    theta_cone = theta_cone * pi / 180.0d0 ! rad
    
    vx_ideal    = 0.0d0
    vy_ideal    = 0.0d0
    vz_ideal    = 0.0d0
    ngas_ideal  = 0.0d0
    ndust_ideal = 0.0d0
       
    missed_cell = 1 ! =1 if cell doesn't satisfy and if statements... should not happen!
    
    ! xcell, ycell and zcell are in frame with origin at bottom-left corner of box
    dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 + (zcell_ideal-0.5d0)**2  ! in frame with origin at center of box
    dist_cell = sqrt(dist2)
    
    if (coldens_norm .gt. 1.0d0) then 
       ! if gas norm set as column density in param file.... OK for model with n~r^-2 ONLY !!!
       n0 = coldens_norm / (r_min * box_size_IM_cm * (1.0d0 - r_min / r_max))  ! cm-3
    else
       ! I use the std density ngas_norm
       n0 = ngas_norm 
    end if
    
    theta_coord = (zcell_ideal-0.5d0) / dist_cell
    theta_coord = acos(theta_coord)
    
    !! Brut force: compare dist to cell center against Rmin/Rmax 
    if (dist_cell < r_max .and. dist_cell > r_min .and. (theta_coord .lt. theta_cone .or. theta_coord .gt. (pi-theta_cone))) then ! cell completely within sphere
!!$          vx_ideal    = Vgas_norm * (xcell_ideal-0.5d0) / r_max
!!$          vy_ideal    = Vgas_norm * (ycell_ideal-0.5d0) / r_max
!!$          vz_ideal    = Vgas_norm * (zcell_ideal-0.5d0) / r_max
       vx_ideal    = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (xcell_ideal - 0.5d0)
       vy_ideal    = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (ycell_ideal - 0.5d0)
       vz_ideal    = Vgas_norm / r_max**(Vgas_slope) * dist_cell**(Vgas_slope-1.0) * (zcell_ideal - 0.5d0)
       
       ngas_ideal  = n0 * (r_min / dist_cell)**(ngas_slope) ! ngas_slope  = +2 for P+11 fiducial model
       ndust_ideal = ndust_norm * ngas_ideal / n0
       missed_cell = 0
    else                                                ! cell completely out of sphere or completely within r_min               
       vx_ideal    = 0.0d0
       vy_ideal    = 0.0d0
       vz_ideal    = 0.0d0
       ngas_ideal  = 0.0d0
       ndust_ideal = 0.0d0
       missed_cell = 0
    end if
    
    
    bparam_ideal = sqrt(2.0d0*kb*Tgas_norm)
    
    if (missed_cell .eq. 1) then
       print*,'I missed a cell in module_idealised_models !'
       print*,dist_cell
       stop
    endif
    
    return
    
  end subroutine shellcone_V_rho_gradient


   subroutine shell_chisholm(ndust_ideal,ngas_ideal,bparam_ideal,vx_ideal,vy_ideal,vz_ideal,xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)
    
    implicit none

    ! point source at center and medium transparent at R < r_min and at R > r_max

    ! Declare arguments
    real(kind=8)                          :: dx_cell
    real(kind=8),intent(inout)            :: ndust_ideal,ngas_ideal,bparam_ideal
    real(kind=8),intent(inout)            :: vx_ideal,vy_ideal,vz_ideal
    real(kind=8),intent(in)               :: xcell_ideal,ycell_ideal,zcell_ideal    
    real(kind=8)                          :: volfrac2
    real(kind=8)                          :: dist_cell,dist2,dist_cell_min,dist_cell_max
    integer(kind=4)                       :: missed_cell
    real(kind=8)                          :: n0

    
    vx_ideal    = 0.0d0
    vy_ideal    = 0.0d0
    vz_ideal    = 0.0d0
    ngas_ideal  = 0.0d0
    ndust_ideal = 0.0d0
    volfrac2    = 1.0d0

    missed_cell = 1 ! =1 if cell doesn't satisfy and if statements... should not happen!
    
    ! xcell, ycell and zcell are in frame with origin at bottom-left corner of box
    dist2 = (xcell_ideal-0.5d0)**2 + (ycell_ideal-0.5d0)**2 + (zcell_ideal-0.5d0)**2  ! in frame with origin at center of box
    dist_cell = sqrt(dist2)

    !! ngas_norm parameter is actually a column density => convert to density n0
    ! n0 = ngas_norm / (r_min * box_size_IM_cm * (1.0d0 - r_min / r_max))  ! cm-3
    n0 = ngas_norm
    
    if (MCsampling) then
       dist_cell_max = dist_cell + dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube
       dist_cell_min = dist_cell - dx_cell * sqrt(3.0d0) / 2.0d0 ! dx_cell * sqrt(3.0) / 2. is half the longest length in a cube

       if ((dist_cell_max > r_max .and. dist_cell_min < r_max) .or. (dist_cell_max > r_min .and. dist_cell_min < r_min)) then ! cell partially within shell
          volfrac2 = mc_sampling_shell(xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)
          missed_cell = 0
       end if
       if (dist_cell_max < r_max .and. dist_cell_min > r_min) then ! cell completely within shell
          volfrac2 = 1.0
          missed_cell = 0
       end if

       ngas_ideal  = n0 *  (dist_cell / r_min)**(ngas_slope) * volfrac2 ! ngas_slope  = -5.47 for Chisholm model
       ndust_ideal = ndust_norm * ngas_ideal / n0

       vx_ideal    = Vgas_norm * (1.0d0 - (r_min / dist_cell))**(Vgas_slope) * (xcell_ideal-0.5d0) / dist_cell
       vy_ideal    = Vgas_norm * (1.0d0 - (r_min / dist_cell))**(Vgas_slope) * (ycell_ideal-0.5d0) / dist_cell
       vz_ideal    = Vgas_norm * (1.0d0 - (r_min / dist_cell))**(Vgas_slope) * (zcell_ideal-0.5d0) / dist_cell
       
       if (dist_cell_min > r_max .or. dist_cell_max < r_min) then  ! cell completely out of shell or completely within r_min                                                       
          vx_ideal    = 0.0d0
          vy_ideal    = 0.0d0
          vz_ideal    = 0.0d0
          ngas_ideal  = 0.0d0
          ndust_ideal = 0.0d0
          missed_cell = 0
       end if
       
    else  !! Brut force: compare dist to cell center against Rmin/Rmax 
       if (dist_cell < r_max .and. dist_cell > r_min) then ! cell completely within sphere
          vx_ideal    = Vgas_norm * (1.0d0 - (r_min / dist_cell))**(Vgas_slope) * (xcell_ideal-0.5d0) / dist_cell
          vy_ideal    = Vgas_norm * (1.0d0 - (r_min / dist_cell))**(Vgas_slope) * (ycell_ideal-0.5d0) / dist_cell
          vz_ideal    = Vgas_norm * (1.0d0 - (r_min / dist_cell))**(Vgas_slope) * (zcell_ideal-0.5d0) / dist_cell
          
          ngas_ideal  = n0 *  (dist_cell / r_min)**(ngas_slope) ! ngas_slope  = -5.47 for Chisholm model
          ndust_ideal = ndust_norm * ngas_ideal / n0
          missed_cell = 0
       else                                                ! cell completely out of sphere or completely within r_min               
          vx_ideal    = 0.0d0
          vy_ideal    = 0.0d0
          vz_ideal    = 0.0d0
          ngas_ideal  = 0.0d0
          ndust_ideal = 0.0d0
          missed_cell = 0
       end if
    end if
    
    bparam_ideal = sqrt(2.0d0*kb*Tgas_norm)

    if (missed_cell .eq. 1) then
       print*,'I missed a cell in module_idealised_models !'
       print*,dist_cell
       stop
    endif
    
    return
    
  end subroutine shell_chisholm
  
  !!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ END MODELS //////////////////////////////////////////////////////////////////////////
  

  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function mc_sampling_shell(xcell_ideal,ycell_ideal,zcell_ideal,dx_cell)

    !! This function returns volfrac2, i.e. weighting factor for gas density and velocity for cells partially in shell-CGM
    !! Works for shells (or spheres with r_min=0) with any rho and V profiles
    !! Velocity components are not corrected here: vx evaluated at r=dist_cell, not dist_cell_mc.... to avoid dealing with V(r) form in this function
    
    implicit none

    real(kind=8)                   :: dx_cell
    integer(kind=4)                :: nMC
    real(kind=8)                   :: volfrac
    real(kind=8)                   :: xmin,xmax,xmc,ymin,ymax,ymc,zmin,zmax,zmc
    real(kind=8)                   :: r,dist_mc,dist_mc2
    integer(kind=4)                :: localseed
    real(kind=8),intent(in)        :: xcell_ideal,ycell_ideal,zcell_ideal
    real(kind=8)                   :: mc_sampling_shell
    
    ! bounds of cell
    xmin = xcell_ideal - 0.5d0*dx_cell 
    xmax = xcell_ideal + 0.5d0*dx_cell
    ymin = ycell_ideal - 0.5d0*dx_cell 
    ymax = ycell_ideal + 0.5d0*dx_cell
    zmin = zcell_ideal - 0.5d0*dx_cell 
    zmax = zcell_ideal + 0.5d0*dx_cell
    ! use Monte Carlo to compute density of cell (i.e. fraction of its volume within the sphere)
    volfrac = 0.0d0
    nMC     = 0
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
       if (dist_mc < r_max .and. dist_mc > r_min) then ! point within shell
          volfrac = volfrac + 1.0d0
       end if
       nMC = nMC + 1
    end do
   
    mc_sampling_shell = volfrac / dble(nMC)

    return
    
  end function mc_sampling_shell

  
  
  function rho_sphere_from_tau0_T_rs_Lbox(tauH,temp,rs)
    
    implicit none 
    
    real(kind=8),intent(in)     :: tauH,temp,rs
    real(kind=8)                :: rho_sphere_from_tau0_T_rs_Lbox!,fix_box_size_cmBIS
        
    ! box_size_IM_cm is fixed by user in param file
    rho_sphere_from_tau0_T_rs_Lbox = 1.d20 * tauH * (temp/1.d4)**0.5 / (5.898d6 * rs * box_size_IM_cm) ! cm^-3
    
    return
    
  end function rho_sphere_from_tau0_T_rs_Lbox
  

  function rho_shell_from_tau0_T_rs_Lbox(tauH,temp,dr)
    
    implicit none 
    
    real(kind=8),intent(in)     :: tauH,temp,dr
    real(kind=8)                :: rho_shell_from_tau0_T_rs_Lbox!,fix_box_size_cmBIS

    rho_shell_from_tau0_T_rs_Lbox = 1.d20 * tauH * (temp/1.d4)**0.5 / (5.898d6 * dr * box_size_IM_cm) ! cm^-3

    return
    
  end function rho_shell_from_tau0_T_rs_Lbox



  
  subroutine read_IdealisedModels_params(pfile)
    
    implicit none
    
    character(*),intent(in) :: pfile
    character(1000)         :: line,name,value
    integer(kind=4)         :: err,i
    logical                 :: section_present

    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:25) == '[IdealisedModels]') then
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
          case ('idealmodel_type')
             read(value,*) idealmodel_type
          case ('MCsampling')
             read(value,*) MCsampling
          case ('ngas_norm')
             read(value,*) ngas_norm
          case ('coldens_norm')
             read(value,*) coldens_norm            
          case ('ngas_slope')
             read(value,*) ngas_slope
          case ('Vgas_norm')
             read(value,*) Vgas_norm
          case ('Vgas_slope')
             read(value,*) Vgas_slope
          case ('Tgas_norm')
             read(value,*) Tgas_norm
          case ('ndust_norm')
             read(value,*) ndust_norm
          case ('taudust_norm')
             read(value,*) taudust_norm
          case ('box_size_IM_cm')
             read(value,*) box_size_IM_cm
          case ('r_min')
             read(value,*) r_min
          case ('r_max')
             read(value,*) r_max 
          case ('verbose')
             read(value,*) verbose
          end select
       end do
    end if
    close(10)

    print*, ' idealmodel_type      = ',idealmodel_type
    print*, ' MCsampling           = ', MCsampling
    print*, ' ngas_norm            = ',ngas_norm
    print*, ' coldens_norm         = ',coldens_norm
    print*, ' ngas_slope           = ',ngas_slope
    print*, ' Vgas_norm            = ',Vgas_norm
    print*, ' Vgas_slope           = ',Vgas_slope
    print*, ' Tgas_norm            = ',Tgas_norm
    print*, ' ndust_norm           = ',ndust_norm
    print*, ' taudust_norm         = ',taudust_norm
    print*, ' box_size_IM_cm       = ',box_size_IM_cm
    print*, ' r_min                = ',r_min
    print*, ' r_max                = ',r_max
    
    
    return
    
  end subroutine read_IdealisedModels_params

  
  subroutine print_IdealisedModels_params(unit)
    
    implicit none
    
    integer(kind=4),optional,intent(in)  :: unit
    
    if (present(unit)) then
       write(unit,'(a,a)')      ' idealmodel_type       = ',idealmodel_type
       write(unit,'(a,a)')      ' MCsampling            = ', MCsampling
       write(unit,'(a,ES10.3)') '  ngas_norm            = ',ngas_norm
       write(unit,'(a,ES10.3)') '  coldens_norm         = ',coldens_norm
       write(unit,'(a,ES10.3)') '  ngas_slope           = ',ngas_slope
       write(unit,'(a,ES10.3)') '  Vgas_norm            = ',Vgas_norm
       write(unit,'(a,ES10.3)') '  Vgas_slope           = ',Vgas_slope
       write(unit,'(a,ES10.3)') '  Tgas_norm            = ',Tgas_norm
       write(unit,'(a,ES10.3)') '  ndust_norm           = ',ndust_norm
       write(unit,'(a,ES10.3)') '  taudust_norm         = ',taudust_norm
       write(unit,'(a,ES10.3)') '  box_size_IM_cm       = ',box_size_IM_cm
       write(unit,'(a,ES10.3)') '  r_min                = ',r_min
       write(unit,'(a,ES10.3)') '  r_max                = ',r_max
       write(unit,'(a,a)')      ' verbose               = ',verbose
    else
       write(*,'(a,a)')      ' idealmodel_type       = ',idealmodel_type
       write(*,'(a,a)')      ' MCsampling            = ', MCsampling
       write(*,'(a,ES10.3)') '  ngas_norm            = ',ngas_norm
       write(*,'(a,ES10.3)') '  coldens_norm         = ',coldens_norm
       write(*,'(a,ES10.3)') '  ngas_slope           = ',ngas_slope
       write(*,'(a,ES10.3)') '  Vgas_norm            = ',Vgas_norm
       write(*,'(a,ES10.3)') '  Vgas_slope           = ',Vgas_slope
       write(*,'(a,ES10.3)') '  Tgas_norm            = ',Tgas_norm
       write(*,'(a,ES10.3)') '  ndust_norm           = ',ndust_norm
       write(*,'(a,ES10.3)') '  taudust_norm         = ',taudust_norm
       write(*,'(a,ES10.3)') '  box_size_IM_cm       = ',box_size_IM_cm
       write(*,'(a,ES10.3)') '  r_min                = ',r_min
       write(*,'(a,ES10.3)') '  r_max                = ',r_max
       write(*,'(a,a)')      ' verbose               = ',verbose
    end if

  end subroutine print_IdealisedModels_params

  
end module module_idealised_models
