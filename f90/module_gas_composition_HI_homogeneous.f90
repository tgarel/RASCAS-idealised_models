module module_gas_composition

  ! Leo: simplified model with no dust, no D and isotrope angular redistribution
  !      + for the test case: static & homogeneous sphere

  use module_HI_model

  implicit none
  
  type gas
     ! fluid
     real(kind=8) :: v(3)      ! gas velocity [cm/s]
     ! Hydrogen 
     real(kind=8) :: nHI       ! HI numerical density [HI/cm3]
     real(kind=8) :: dopwidth  ! Doppler width [cm/s]
  end type gas

  contains
    ! Routines list:
    ! - subroutine gas_from_ramses_leaves(ramses_var, g)
    ! - subroutine overwrite_gas(g,nhi,vth)
    ! - function get_gas_velocity(cell_gas)
    ! - function  gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs)
    ! - subroutine gas_scatter(flag,cell_gas,nu_cell,k,nu_ext,iran)
    ! - subroutine dump_gas(unit,g)
    ! - subroutine read_gas(unit,n,g)
    ! - subroutine gas_destructor(g)

    subroutine gas_from_ramses_leaves(ramses_var, g)

      ! define gas contents from ramses raw data

      !character(2000),intent(in)        :: repository 
      !integer(kind=4),intent(in)        :: snapnum
      !integer(kind=4),intent(in)        :: nleaf,nvar
      real(kind=8),dimension(:,:),intent(in)         :: ramses_var
      type(gas),dimension(:),allocatable,intent(out) :: g
      integer(kind=4)                   :: ileaf,nleaf

      ! make sure conversion factors are set 
      ! call read_conversion_scales(repository,snapnum)

      nleaf = size(ramses_var,2)

      ! allocate gas-element array
      allocate(g(nleaf))

      ! compute gas props. leaf by leaf
      do ileaf = 1,nleaf
         g(ileaf)%v        = 0.d0
         g(ileaf)%nHI      = 1.d0
         g(ileaf)%dopwidth = 10.d0*1.e5 ! 10 km/s in cgs
      end do

      return
      
    end subroutine gas_from_ramses_leaves



    subroutine overwrite_gas(g,nhi,vth)
      type(gas),dimension(:),intent(inout) :: g
      real(kind=8),intent(in)              :: nhi,vth
      integer                              :: nleaf

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



    function  gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs)

      ! NB: also return distance to interaction (in variable distance_to_border_cm) ...
      ! LEo: and also update tau_abs if gas_get_scatter_flag = 0

      ! check whether scattering occurs within cell (scatter_flag > 0) or not (scatter_flag==0)
      type(gas),intent(in)                  :: cell_gas
      real(kind=8),intent(inout)            :: distance_to_border_cm
      real(kind=8),intent(in)               :: nu_cell
      real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur. 
      integer(kind=4)                       :: gas_get_scatter_flag 
      real(kind=8)                          :: tau_HI, tau_cell
      
      ! compute optical depths for different components of the gas.
      tau_HI   = get_tau_HI(cell_gas%nHI, cell_gas%dopwidth, distance_to_border_cm, nu_cell)
      tau_cell = tau_HI
      ! DEBUG case
      !tau_cell = 0
      
      if (tau_abs > tau_cell) then  ! photon is due for absorption outside the cell 
         gas_get_scatter_flag = 0
         tau_abs = tau_abs - tau_cell
         if (tau_abs.lt.0.d0) then
            print*, 'tau_abs est negatif'
            stop
         endif
      else  ! the scattering happens inside the cell on a H atom (no dust)
         gas_get_scatter_flag = 1
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
         call scatter_HI_isotrope(cell_gas%v, cell_gas%nHI, cell_gas%dopwidth, nu_cell, k, nu_ext, iran)
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

    end subroutine dump_gas



    subroutine read_gas(unit,n,g)
      
      integer,intent(in)                             :: unit,n
      type(gas),dimension(:),allocatable,intent(out) :: g
      integer                                        :: i

      allocate(g(1:n))

      read(unit) (g(i)%v(:),i=1,n)
      read(unit) (g(i)%nHI,i=1,n)
      read(unit) (g(i)%dopwidth,i=1,n)

#ifdef DEBUG
      print*,'in read_gas:',n
      print*,'in read_gas: ',minval(g%nhi),maxval(g%nhi)
      print*,'in read_gas: ',minval(g%dopwidth),maxval(g%dopwidth)
#endif

    end subroutine read_gas


    subroutine gas_destructor(g)
      type(gas),dimension(:),allocatable,intent(inout) :: g
      deallocate(g)
    end subroutine gas_destructor



  end module module_gas_composition
