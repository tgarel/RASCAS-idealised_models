module module_select

  use module_gas_composition

  implicit none

  interface select_from_domain
     module procedure xleaf_select, leaflevel_select, gas_select
  end interface select_from_domain

contains

  subroutine xleaf_select(arr_in, ind_sel, arr_out)
    real(kind=8),dimension(:,:),intent(in)              :: arr_in
    integer,dimension(:),intent(in)                     :: ind_sel
    real(kind=8),dimension(:,:),allocatable,intent(out) :: arr_out
    integer ::i,n,ii

    n = size(ind_sel)
    
    allocate(arr_out(1:n,1:3))
    do i=1,n
       ii=ind_sel(i)
       arr_out(i,1:3)=arr_in(ii,1:3)
    enddo

  end subroutine xleaf_select


  subroutine leaflevel_select(arr_in, ind_sel, arr_out)
    integer,dimension(:),intent(in)              :: arr_in
    integer,dimension(:),intent(in)              :: ind_sel
    integer,dimension(:),allocatable,intent(out) :: arr_out
    integer ::i,n,ii

    n = size(ind_sel)
    
    allocate(arr_out(1:n))
    do i=1,n
       ii=ind_sel(i)
       arr_out(i)=arr_in(ii)
    enddo

  end subroutine leaflevel_select

  !!!! on beneficirait pleinement ici d'un type gas tableau(type(cell))
  !on pourrait donc selectionner ici les elements sans connaitre explicitement les composants de ce type...
  ! subroutine gas_select_oldshape(arr_in, ind_sel, arr_out)

  !   type(gas),intent(in)            :: arr_in
  !   integer,dimension(:),intent(in) :: ind_sel
  !   type(gas),intent(out)           :: arr_out
  !   integer ::i,n,ii

  !   n = size(ind_sel)
    
  !   allocate(arr_out%density(1:n))
  !   allocate(arr_out%vx(1:n))
  !   allocate(arr_out%vy(1:n))
  !   allocate(arr_out%vz(1:n))
  !   allocate(arr_out%pressure(1:n))
  !   allocate(arr_out%metallicity(1:n))
  !   do i=1,n
  !      ii=ind_sel(i)
  !      arr_out%density(i)=arr_in%density(ii)
  !      arr_out%vx(i)=arr_in%vx(ii)
  !      arr_out%vy(i)=arr_in%vy(ii)
  !      arr_out%vz(i)=arr_in%vz(ii)
  !      arr_out%pressure(i)=arr_in%pressure(ii)
  !      arr_out%metallicity(i)=arr_in%metallicity(ii)
  !   enddo

  ! end subroutine gas_select_oldshape

  
  subroutine gas_select(arr_in, ind_sel, arr_out)

    type(gas),intent(in),dimension(:)            :: arr_in
    integer,dimension(:),intent(in) :: ind_sel
    type(gas),intent(out),dimension(:),allocatable           :: arr_out
    integer ::i,n,ii


    n = size(ind_sel)
    
    allocate(arr_out(1:n))
    do i=1,n
       ii=ind_sel(i)
       arr_out(i) = arr_in(ii)
    enddo

  end subroutine gas_select



end module module_select
    
