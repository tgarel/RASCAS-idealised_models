module module_gas_composition

  ! module deps
  ! use HI_model

  implicit none

  ! define type gas
  type gas
     ! ...
  end type gas

  contains


    subroutine gas_from_ramses_leaves(ramses_leaves, g)

      type(leaf_cell_array), intent(in) :: ramses_leaves
      type(gas),intent(out)             :: g

      print *,"There is no cause for alarm, but there probably will be"
      print *," "
      print *," this routine is not defined, better stop "
      stop

    end subroutine gas_from_ramses_leaves



    function get_gas_velocity(cell_gas)
      
      type(gas),intent(in)      :: cell_gas
      real(kind=8),dimension(3) :: get_gas_velocity
      get_gas_velocity(:) = cell_gas%v(:)
      return
    
    end function get_gas_velocity



    function gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs, iran)

      ! check whether scattering occurs within cell (scatter_flag > 0) or not (scatter_flag==0)
      type(gas),intent(in)                  :: cell_gas
      real(kind=8),intent(inout)            :: distance_to_border_cm
      real(kind=8),intent(in)               :: nu_cell
      real(kind=8),intent(inout)            :: tau_abs                ! tau at which scattering is set to occur.
      integer,intent(inout)                 :: iran 
 
      print *,"There is no cause for alarm, but there probably will be"
      print *," "
      print *," this routine is not defined, better stop "
      stop
      
      return
    end function gas_get_scatter_flag



    subroutine gas_scatter(flag,cell_gas,nu_cell,k,nu_ext,iran)

      integer, intent(in)                       :: flag
      type(gas), intent(in)                     :: cell_gas
      real(kind=8), intent(inout)               :: nu_cell, nu_ext
      real(kind=8), dimension(3), intent(inout) :: k
 
      print *,"There is no cause for alarm, but there probably will be"
      print *," "
      print *," this routine is not defined, better stop "
      stop

    end subroutine gas_scatter



    subroutine dump_gas(unit,g)
      
      type(gas),dimension(:),intent(in) :: g
      integer,intent(in)                :: unit

      print *,"There is no cause for alarm, but there probably will be"
      print *," "
      print *," this routine is not defined, better stop "
      stop

    end subroutine dump_gas



    subroutine read_gas(unit,n,g)
      
      integer,intent(in)                             :: unit,n
      type(gas),dimension(:),allocatable,intent(out) :: g
 
      print *,"There is no cause for alarm, but there probably will be"
      print *," "
      print *," this routine is not defined, better stop "
      stop

    end subroutine read_gas



    subroutine gas_destructor(g)

      type(gas),dimension(:),allocatable,intent(inout) :: g
      deallocate(g)

    end subroutine gas_destructor


  end module module_gas_composition
