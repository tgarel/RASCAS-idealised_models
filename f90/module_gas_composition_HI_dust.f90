module module_gas_composition

  use module_HI_model
  use module_dust_model
  use module_ramses
  use module_random

  implicit none

  type gas
     ! fluid
     real(kind=8) :: v(3)      ! gas velocity [cm/s]
     ! Hydrogen 
     real(kind=8) :: nHI       ! HI numerical density [HI/cm3]
     real(kind=8) :: dopwidth  ! Doppler width [cm/s]
     ! dust
     real(kind=8) :: ndust     ! numerical density of dust particles [#/cm3]
  end type gas

  contains
    ! Routines list:
    ! - subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var, g)
    ! - function get_gas_velocity(cell_gas)
    ! - function gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs)
    ! - subroutine gas_scatter(flag,cell_gas,nu_cell,k,nu_ext,iran)
    ! - subroutine dump_gas(unit,g)
    ! - subroutine read_gas(unit,n,g)
    ! - subroutine gas_destructor(g)
    
    subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramses_var, g)

      ! define gas contents from ramses raw data

      character(2000),intent(in)                     :: repository 
      integer(kind=4),intent(in)                     :: snapnum
      integer(kind=4),intent(in)                     :: nleaf,nvar
      real(kind=8),dimension(:,:),intent(in)         :: ramses_var
      type(gas),allocatable,intent(out)              :: g
      integer(kind=4)                                :: ileaf

      ! make sure conversion factors are set 
      call read_conversion_scales(repository,snapnum)

      nleaf = size(ramses_var,2)
      
      ! allocate gas-element array
      allocate(g(nleaf))

      ! compute gas props. leaf by leaf
      do ileaf = 1,nleaf
         g(ileaf)%v        = get_velocity_cgs(ramses_var,ileaf)
         g(ileaf)%nHI      = get_HI_density(ramses_var,ileaf) 
         g(ileaf)%dopwidth = get_doppler_velocity_cgs(ramses_var,ileaf,vturb=10.)
         g(ileaf)%ndust    = get_dust_density(ramses_var,ileaf)
      end do

      return
      
    end subroutine gas_from_ramses_leaves



    function get_gas_velocity(cell_gas)
      type(gas),intent(in)      :: cell_gas
      real(kind=8),dimension(3) :: get_gas_velocity
      get_gas_velocity(:) = cell_gas%v(:)
      return
    end function get_gas_velocity



    function gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs)

      ! NB: also return distance to interaction (in variable distance_to_border_cm) ...
      ! LEo: and also update tau_abs if gas_get_scatter_flag = 0

      ! check whether scattering occurs within cell (scatter_flag > 0) or not (scatter_flag==0)
      type(gas),intent(in)                  :: cell_gas
      real(kind=8),intent(inout)            :: distance_to_border_cm
      real(kind=8),intent(in)               :: nu_cell
      real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur. 
      integer(kind=4)                       :: gas_get_scatter_flag 
      real(kind=8)                          :: tau_HI, tau_cell, tau_dust, one_proba
      
      ! compute optical depths for different components of the gas.
      tau_HI   = get_tau_HI(cell_gas%nHI, cell_gas%dopwidth, distance_to_border_cm, nu_cell)
      tau_dust = get_tau_dust(cell_gas%ndust, distance_to_border_cm)
      tau_cell = tau_HI + tau_dust
      
      if (tau_abs > tau_cell) then  ! photon is due for absorption outside the cell 
         gas_get_scatter_flag = 0
         tau_abs = tau_abs - tau_cell
         if (tau_abs.lt.0.d0) then
            print*, 'tau_abs est negatif'
            stop
         endif
      else  ! the scattering happens inside the cell
         ! decider si HI ou dust selon rapport des tau
         one_proba = tau_HI / tau_cell         
         tirage = ran3(iran)
         if(tirage <= one_proba)then
            gas_get_scatter_flag = 1
         else ! interaction with dust
            gas_get_scatter_flag = 3
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
         call scatter_HI(cell_gas%v, cell_gas%nHI, cell_gas%dopwidth, nu_cell, k, nu_ext, iran)
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

    end subroutine dump_gas



    subroutine read_gas(unit,n,g)
      
      integer,intent(in)                             :: unit,n
      type(gas),dimension(:),allocatable,intent(out) :: g
      integer                                        :: i

      allocate(g(1:n))

      read(unit) (g(i)%v(:),i=1,n)
      read(unit) (g(i)%nHI,i=1,n)
      read(unit) (g(i)%dopwidth,i=1,n)
      read(unit) (g(i)%ndust,i=1,n)

    end subroutine read_gas


    subroutine gas_destructor(g)
      type(gas),dimension(:),allocatable,intent(inout) :: g
      deallocate(g)
    end subroutine gas_destructor



  end module module_gas_composition
