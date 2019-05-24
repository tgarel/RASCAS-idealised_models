module module_gas_composition

  ! Pure RAMSES hydro variables
  ! - is NVAR dependant

  implicit none

  type gas
     real(kind=8) :: density
     real(kind=8) :: vx
     real(kind=8) :: vy
     real(kind=8) :: vz
     real(kind=8) :: pressure
     real(kind=8) :: metallicity
  end type gas

  real(kind=8),public :: box_size_cm   ! size of simulation box in cm. 


contains


  subroutine gas_from_ramses_leaves(repository,snapnum,nleaf,nvar,ramsesvar,g)

    ! just repack ramses raw data into gas

    use module_ramses

    character(2000),intent(in)                      :: repository 
    integer(kind=4),intent(in)                      :: snapnum
    integer(kind=4),intent(in)                      :: nleaf,nvar
    real(kind=8),dimension(nvar,nleaf),intent(in)   :: ramsesvar
    type(gas),dimension(:),allocatable, intent(out) :: g
    integer(kind=4)                                 :: i

    allocate(g(1:nleaf))

    do i=1,nleaf
       g(i)%density     = ramsesvar(1,i)
       g(i)%vx          = ramsesvar(2,i)
       g(i)%vy          = ramsesvar(3,i)
       g(i)%vz          = ramsesvar(4,i)
       g(i)%pressure    = ramsesvar(5,i)
       g(i)%metallicity = ramsesvar(6,i)
    end do

    box_size_cm = ramses_get_box_size_cm(repository,snapnum)

  end subroutine gas_from_ramses_leaves



  subroutine gas_destructor(g)
    type(gas),dimension(:),allocatable,intent(inout) :: g
    deallocate(g)
  end subroutine gas_destructor



  subroutine dump_gas(unit,g)

    type(gas),dimension(:),intent(in) :: g
    integer(kind=4),intent(in)        :: unit
    integer(kind=4)                   :: i,nleaf

    nleaf = size(g)
    write(unit) (g(i)%density, i=1,nleaf)
    write(unit) (g(i)%vx, i=1,nleaf)
    write(unit) (g(i)%vy, i=1,nleaf)
    write(unit) (g(i)%vz, i=1,nleaf)
    write(unit) (g(i)%pressure, i=1,nleaf)
    write(unit) (g(i)%metallicity, i=1,nleaf)

  end subroutine dump_gas


  subroutine read_gas(unit,n,g)

    integer(kind=4),intent(in)                     :: unit,n
    type(gas),dimension(:),allocatable,intent(out) :: g
    integer(kind=4)                                :: i

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


  ! subroutine overwrite_gas(g)
  ! no need...

  ! end subroutine overwrite_gas


  ! function get_gas_velocity(cell_gas)
    ! no need...

  ! end function get_gas_velocity


  subroutine read_gas_composition_params(pfile)

    ! ---------------------------------------------------------------------------------
    ! no parameter for this composition 
    ! subroutine here for consistency with other composition...
    ! ---------------------------------------------------------------------------------

    character(*),intent(in) :: pfile
    character(1000)         :: line
    integer(kind=4)         :: err
    logical                 :: section_present

    section_present = .false.
    open(unit=10,file=trim(pfile),status='old',form='formatted')
    ! search for section start
    do
       read (10,'(a)',iostat=err) line
       if(err/=0) exit
       if (line(1:17) == '[gas_composition]') then
          section_present = .true.
          exit
       end if
    end do

    ! print warning and stop if read section present
    if (section_present) then 
       print *,"ERROR: this gas composition do not admit parameter, better stop..."
       stop
    endif
    close(10)
    
  end subroutine read_gas_composition_params



  subroutine print_gas_composition_params(unit)

    ! ---------------------------------------------------------------------------------
    ! write parameter values to std output or to an open file if argument unit is
    ! present.
    ! ---------------------------------------------------------------------------------

    integer(kind=4),optional,intent(in) :: unit

    if (present(unit)) then 
       write(unit,'(a,a,a)') '[gas_composition]'
       write(unit,'(a)')       '# no parameter for this composition'
       write(unit,'(a)')             ' '
    else
       write(*,'(a,a,a)') '[gas_composition]'
       write(*,'(a)')       '# no parameter for this composition'
       write(*,'(a)')             ' '
    end if

  end subroutine print_gas_composition_params


end module module_gas_composition
