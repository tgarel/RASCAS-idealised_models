module module_gas_composition

  implicit none

  type gas
     real(kind=8) :: density
     real(kind=8) :: vx
     real(kind=8) :: vy
     real(kind=8) :: vz
     real(kind=8) :: pressure
     real(kind=8) :: metallicity
  end type gas

  contains

    subroutine gas_from_ramses_leaves(ramsesvar, g)

      real(kind=8),dimension(:,:),intent(in)          :: ramsesvar
      type(gas),dimension(:),allocatable, intent(out) :: g
      integer                                         :: i,nleaf

      nleaf = size(ramsesvar,2)

      allocate(g(1:nleaf))

      do i=1,nleaf
         g(i)%density     = ramsesvar(1,i)
         g(i)%vx          = ramsesvar(2,i)
         g(i)%vy          = ramsesvar(3,i)
         g(i)%vz          = ramsesvar(4,i)
         g(i)%pressure    = ramsesvar(5,i)
         g(i)%metallicity = ramsesvar(6,i)
      end do
       
    end subroutine gas_from_ramses_leaves



    subroutine gas_destructor(g)
      type(gas),dimension(:),allocatable,intent(inout) :: g
      deallocate(g)
    end subroutine gas_destructor



    subroutine dump_gas(unit,g)
      
      type(gas),dimension(:),intent(in) :: g
      integer,intent(in)                :: unit
      integer                           :: i,nleaf

      nleaf = size(g)
      write(unit) (g(i)%density, i=1,nleaf)
      write(unit) (g(i)%vx, i=1,nleaf)
      write(unit) (g(i)%vy, i=1,nleaf)
      write(unit) (g(i)%vz, i=1,nleaf)
      write(unit) (g(i)%pressure, i=1,nleaf)
      write(unit) (g(i)%metallicity, i=1,nleaf)

    end subroutine dump_gas


    subroutine read_gas(unit,n,g)
      
      integer,intent(in)                             :: unit,n
      type(gas),dimension(:),allocatable,intent(out) :: g
      integer                                        :: i

      allocate(g(1:n))

      read(unit) (g(i)%density,i=1,n)
      read(unit) (g(i)%vx,i=1,n)
      read(unit) (g(i)%vy,i=1,n)
      read(unit) (g(i)%vz,i=1,n)
      read(unit) (g(i)%pressure,i=1,n)
      read(unit) (g(i)%metallicity,i=1,n)

    end subroutine read_gas


    ! subroutine gas_get_scatter_flag(taus)
    !   ! no need 

    ! end subroutine gas_get_taus


    !function  gas_get_scatter_flag(cell_gas, distance_to_border_cm, nu_cell, tau_abs, iran)
    !   ! no need...
      
    !end function gas_get_scatter_flag



  end module module_gas_composition
