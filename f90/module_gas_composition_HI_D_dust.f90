module module_gas_composition

  ! mix of HI, Deuterium, and dust.
  ! - The Deuterium density is set as a fraction of H density through the parameter deut2H_nb_ratio (in module_params.f90).
  ! - The thermal velocity of D is derived from the one of H using the mass ratio (sqrt_H2Deut_mass_ratio in module_constants.f90).
  ! - The HI content is from RAMSES
  ! - The dust content is a function of metallicity and HI. 
  
  use module_HI_model
  use module_D_model
  use module_dust_model
  use module_random
  use module_ramses
  use module_constants, only : sqrt_H2Deut_mass_ratio, XH
  use module_params,    only : deut2H_nb_ratio, dust_to_metal_ratio, mH_over_mdust

  implicit none

  type gas
     ! fluid
     real(kind=8) :: v(3)      ! gas velocity [cm/s]
     ! Hydrogen 
     real(kind=8) :: nHI       ! HI numerical density [HI/cm3]
     real(kind=8) :: dopwidth  ! Doppler width [cm/s]
     ! Deuterium
     ! -> density is computed as nHI * deut2H_nb_ratio
     ! -> dopwidth is computed as dopwidth * sqrt_H2Deut_mass_ratio.
     ! dust
     real(kind=8) :: ndust     ! numerical density of dust particles [#/cm3]
  end type gas

  contains

    subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var, g)

      ! define gas contents from ramses raw data
      
      character(2000),intent(in)        :: repository 
      integer(kind=4),intent(in)        :: snapnum
      integer(kind=4),intent(in)        :: nleaf,nvar
      real(kind=8),intent(in)           :: ramses_var(nvar,nleaf)
      type(gas),dimension(:),allocatable,intent(out) :: g
      integer(kind=4)                   :: ileaf
      real(kind=8),allocatable          :: v(:,:), T(:), nhi(:), metallicity(:)
      
      ! allocate gas-element array
      allocate(g(nleaf))

      ! compute velocities in cm / s 
      allocate(v(3,nleaf))
      call ramses_get_velocity_cgs(repository,snapnum,nleaf,nvar,ramses_var,v)
      do ileaf = 1,nleaf
         g(ileaf)%v = v(:,ileaf)
      end do
      deallocate(v)
      
      ! get nHI and temperature from ramses
      allocate(T(nleaf),nhi(nleaf))
      call ramses_get_T_nhi_cgs(repository,snapnum,nleaf,nvar,ramses_var,T,nhi)
      g(:)%nHI = nhi(:)
      ! compute thermal velocity 
      ! ++++++ TURBULENT VELOCITY >>>>> parameter to add and use here
      g(:)%dopwidth = sqrt((2.0d0*kb/mp)*T) ! [ cm/s ]
      deallocate(T,nhi)

      ! get ndust
      ! Use Verhamme 2012, Eq. 3.
      allocate(metallicity(nleaf))
      call ramses_get_metallicity(repository,snapnum,nleaf,nvar,ramses_var,metallicity)
      g(:)%ndust = (dust_to_metal_ratio * mH_over_mdust / XH) * metallicity * g(:)%nHI     !! JB: To test ... 
      deallocate(metallicity)
      
      return
      
    end subroutine gas_from_ramses_leaves

    
    subroutine overwrite_gas(g,nhi,vth)
      type(gas),dimension(:),intent(inout) :: g
      real(kind=8),intent(in)              :: nhi,vth
      integer                              :: nleaf
      stop
      nleaf = size(g)
      g(:)%nhi = nhi
      g(:)%dopwidth = vth

#ifdef DEBUG
      print*,'in overwrite_gas: allocated g?',shape(g)
      print*,'in overwrite_gas: ',nhi,vth
      print*,'in overwrite_gas: ',nleaf,nhi
      print*,'in overwrite_gas: ',minval(g%nhi),maxval(g%nhi)
      print*,'in overwrite_gas: ',minval(g%dopwidth),maxval(g%dopwidth)
#endif

    end subroutine overwrite_gas



    function get_gas_velocity(cell_gas)
      type(gas),intent(in)      :: cell_gas
      real(kind=8),dimension(3) :: get_gas_velocity
      get_gas_velocity(:) = cell_gas%v(:)
      return
    end function get_gas_velocity



    function  gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs,iran)

      ! --------------------------------------------------------------------------
      ! Decide whether a scattering event occurs, and if so, on which element
      ! --------------------------------------------------------------------------
      ! INPUTS:
      ! - cell_gas : a mix of H, D, and dust
      ! - distance_to_border_cm : the maximum distance the photon may travel (before leaving the cell)
      ! - nu_cell : photon frequency in cell's frame [ Hz ]
      ! - tau_abs : optical depth at which the next scattering event will occur
      ! - iran    : random generator state of the photon 
      ! OUTPUTS:
      ! - distance_to_border_cm : comes out as the distance to scattering event (if there is an event)
      ! - tau_abs : 0 if a scatter occurs, decremented by tau_cell if photon escapes cell. 
      ! - gas_get_scatter_flag : 0 [no scatter], 1 [H scatter], 2 [D scatter], 3 [dust]
      ! --------------------------------------------------------------------------

      ! check whether scattering occurs within cell (scatter_flag > 0) or not (scatter_flag==0)
      type(gas),intent(in)                  :: cell_gas
      real(kind=8),intent(inout)            :: distance_to_border_cm
      real(kind=8),intent(in)               :: nu_cell
      real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur.
      integer,intent(inout)                 :: iran
      integer(kind=4)                       :: gas_get_scatter_flag 
      real(kind=8)                          :: tau_HI, tau_dust, tau_cell, tau_D, tirage, proba1, proba2
      
      ! compute optical depths for different components of the gas.
      tau_HI   = get_tau_HI(cell_gas%nHI, cell_gas%dopwidth, distance_to_border_cm, nu_cell)
      tau_dust = get_tau_dust(cell_gas%ndust, distance_to_border_cm)
      tau_D    = get_tau_D(cell_gas%nHI * deut2H_nb_ratio, cell_gas%dopwidth * sqrt_H2Deut_mass_ratio,distance_to_border_cm, nu_cell)
      tau_cell = tau_HI + tau_D + tau_dust
      
      if (tau_abs > tau_cell) then  ! photon is due for absorption outside the cell 
         gas_get_scatter_flag = 0   ! no scatter
         tau_abs = tau_abs - tau_cell
         if (tau_abs.lt.0.d0) then
            print*, 'tau_abs est negatif'
            stop
         endif
      else  ! the scattering happens inside the cell. 
         ! decider si HI ou dust selon rapport des tau
         proba1 = tau_HI / tau_cell         
         proba2 = proba1 + tau_D / tau_cell
         tirage = ran3(iran)
         if(tirage <= proba1)then
            gas_get_scatter_flag = 1 ! HI scatter
         else if (tirage <= proba2) then 
            gas_get_scatter_flag = 2 ! interaction with a Deuterium atom 
         else 
            gas_get_scatter_flag = 3 ! interaction with dust
         endif
         ! et transformer distance_to_border_cm en distance_to_absorption_cm ...
         distance_to_border_cm = distance_to_border_cm * (tau_abs / tau_cell)
      end if
      
      return
    end function gas_get_scatter_flag



    subroutine gas_scatter(flag,cell_gas,nu_cell,k,nu_ext,iran)

      integer, intent(in)                       :: flag
      type(gas), intent(in)                     :: cell_gas
      real(kind=8), intent(inout)               :: nu_cell, nu_ext
      real(kind=8), dimension(3), intent(inout) :: k
      integer, intent(inout)                    :: iran

      select case(flag)
      case(1)
         call scatter_HI(cell_gas%v, cell_gas%dopwidth, nu_cell, k, nu_ext, iran)
      case(2)
         call scatter_D(cell_gas%v,cell_gas%dopwidth*sqrt_H2Deut_mass_ratio, nu_cell, k, nu_ext, iran)
      case(3)
         call scatter_dust(cell_gas%v, nu_cell, k, nu_ext, iran)
      end select

    end subroutine gas_scatter



   subroutine dump_gas(unit,g)
      
      type(gas),dimension(:),intent(in) :: g
      integer,intent(in)                :: unit
      integer                           :: i,nleaf

      nleaf = size(g)
      write(unit) (g(i)%v(:), i=1,nleaf)
      write(unit) (g(i)%nHI, i=1,nleaf)
      write(unit) (g(i)%dopwidth, i=1,nleaf)
      write(unit) (g(i)%ndust, i=1,nleaf)
      write(unit) deut2H_nb_ratio 
            

    end subroutine dump_gas



    subroutine read_gas(unit,n,g)
      
      integer,intent(in)                             :: unit,n
      type(gas),dimension(:),allocatable,intent(out) :: g
      integer                                        :: i
      real(kind=8) :: D2H_ratio
      allocate(g(1:n))

      read(unit) (g(i)%v(:),i=1,n)
      read(unit) (g(i)%nHI,i=1,n)
      read(unit) (g(i)%dopwidth,i=1,n)
      read(unit) (g(i)%ndust,i=1,n)
      read(unit) D2H_ratio  ! JB: not used ... 
      if (D2H_ratio /= deut2H_nb_ratio) then
         write(*,*) 'Oh boy! D2H_ratio /= deut2H_nb_ratio'
         stop
      end if
      
    end subroutine read_gas


    subroutine gas_destructor(g)
      type(gas),dimension(:),allocatable,intent(inout) :: g
      deallocate(g)
    end subroutine gas_destructor



  end module module_gas_composition
