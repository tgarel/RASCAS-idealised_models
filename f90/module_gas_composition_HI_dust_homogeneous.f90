module module_gas_composition

  ! no D and isotrope angular redistribution
  
  use module_HI_model
  use module_dust_model
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
    ! - subroutine gas_from_ramses_leaves(ramses_var, g)
    ! - subroutine overwrite_gas(g,nhi,vth,ndust)
    ! - function get_gas_velocity(cell_gas)
    ! - function gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs, iran)
    ! - subroutine gas_scatter(flag,cell_gas,nu_cell,k,nu_ext,iran)
    ! - subroutine dump_gas(unit,g)
    ! - subroutine read_gas(unit,n,g)
    ! - subroutine gas_destructor(g)

    subroutine gas_from_ramses_leaves(ramses_var, g)

      ! define gas contents from ramses raw data
      
      real(kind=8),dimension(:,:),intent(in)         :: ramses_var
      type(gas),dimension(:),allocatable,intent(out) :: g
      integer(kind=4)                                :: ileaf,nleaf

      nleaf = size(ramses_var,2)

      ! allocate gas-element array
      allocate(g(nleaf))

      ! compute gas props. leaf by leaf
      do ileaf = 1,nleaf
         g(ileaf)%v        = 0.d0
         g(ileaf)%nHI      = 1.d0
         g(ileaf)%dopwidth = 10.d0*1.e5 ! 10 km/s in cgs
         g(ileaf)%ndust    = 0.01
      end do

      return
      
    end subroutine gas_from_ramses_leaves



    subroutine overwrite_gas(g,nhi,vth,ndust)
      type(gas),dimension(:),intent(inout) :: g
      real(kind=8),intent(in)              :: nhi,vth,ndust
      integer(kind=4)                      :: nleaf

      nleaf = size(g)
      g(:)%nhi      = nhi
      g(:)%dopwidth = vth
      g(:)%ndust    = ndust
      
#ifdef DEBUG
      print*,'in overwrite_gas: allocated g?',shape(g)
      print*,'in overwrite_gas: ',nhi,vth
      print*,'in overwrite_gas: ',nleaf,nhi
      print*,'in overwrite_gas: ',minval(g%nhi),maxval(g%nhi)
      print*,'in overwrite_gas: ',minval(g%dopwidth),maxval(g%dopwidth)
      print*,'in overwrite_gas: ',minval(g%ndust),maxval(g%ndust)
#endif

    end subroutine overwrite_gas



    function get_gas_velocity(cell_gas)
      type(gas),intent(inout)      :: cell_gas
      real(kind=8),dimension(3)    :: get_gas_velocity
      get_gas_velocity(:) = cell_gas%v(:)
      return
    end function get_gas_velocity



    function gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs, iran)

      ! --------------------------------------------------------------------------
      ! Decide whether a scattering event occurs, and if so, on which element
      ! --------------------------------------------------------------------------
      ! INPUTS:
      ! - cell_gas : H and dust
      ! - distance_to_border_cm : the maximum distance the photon may travel (before leaving the cell)
      ! - nu_cell : photon frequency in cell's frame [ Hz ]
      ! - tau_abs : optical depth at which the next scattering event will occur
      ! - iran    : random generator state of the photon 
      ! OUTPUTS:
      ! - distance_to_border_cm : comes out as the distance to scattering event (if there is an event)
      ! - tau_abs : 0 if a scatter occurs, decremented by tau_cell if photon escapes cell. 
      ! - gas_get_scatter_flag : 0 [no scatter], 1 [H scatter], 3 [dust scatter]
      ! --------------------------------------------------------------------------
      

      ! check whether scattering occurs within cell (scatter_flag > 0) or not (scatter_flag==0)
      type(gas),intent(in)                  :: cell_gas
      real(kind=8),intent(inout)            :: distance_to_border_cm
      real(kind=8),intent(in)               :: nu_cell
      real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur.
      integer,intent(inout)                 :: iran      
      integer(kind=4)                       :: gas_get_scatter_flag 
      real(kind=8)                          :: tau_HI, tau_cell, tau_dust, tirage, proba

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
         proba = tau_HI / tau_cell         
         tirage = ran3(iran)
         if(tirage <= proba) then ! H scatter
            gas_get_scatter_flag = 1
         else ! interaction with dust
            gas_get_scatter_flag = 3
         endif
         ! update distance_to_border_cm to account for advance of photon up to scattering event. 
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
         call scatter_HI_isotrope(cell_gas%v, cell_gas%dopwidth, nu_cell, k, nu_ext, iran)
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

#ifdef DEBUG
      print*,'in read_gas:',n
      print*,'in read_gas: ',minval(g%nhi),maxval(g%nhi)
      print*,'in read_gas: ',minval(g%dopwidth),maxval(g%dopwidth)
      print*,'in read_gas: ',minval(g%ndust),maxval(g%ndust)
#endif

    end subroutine read_gas


    subroutine gas_destructor(g)
      type(gas),dimension(:),allocatable,intent(inout) :: g
      deallocate(g)
    end subroutine gas_destructor



  end module module_gas_composition
