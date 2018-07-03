module module_mock

  private

  ! direction of observation
  real(kind=8) :: kobs(3),kobs_perp_1(3),kobs_perp_2(3)

  ! center of observation
  real(kind=8) :: mock_center(3)
  
  ! detector (image for now)
  integer(kind=4) :: npix
  real(kind=8)    :: image(npix,npix)
  
  
  

  public mock_init
  
contains

  subroutine mock_init()

    implicit none
        
    image = 0.0d0
    ! define direction of observation (normalise vector)
    kobs = (/0.31,0.2,0.42/)
    kobs = kobs / sqrt(kobs(1)*kobs(1)+kobs(2)*kobs(2)+kobs(3)*kobs(3))
    ! define basis for sky plane
    if (kobs(0) < 1.) then
       ! kobs_perp_1 is cross prod. of kobs and x 
       kobs_perp_1(1) = 0.
       kobs_perp_1(2) = kobs(3)
       kobs_perp_1(3) = -kobs(2)
       kobs_perp_1 = kobs_perp_1 / sqrt(kobs_perp_1(1)*kobs_perp_1(1)+kobs_perp_1(2)*kobs_perp_1(2)+kobs_perp_1(3)*kobs_perp_1(3))
       ! kobs_perp_2 = kobs x kobs_perp_1
       kobs_perp_2(1) = kobs(2)*kobs_perp_1(3) - kobs(3)*kobs_perp_1(2)
       kobs_perp_2(2) = -kobs(1)*kobs_perp_1(3)
       kobs_perp_2(3) = kobs(1)*kobs_perp_1(2)
       kobs_perp_2 = kobs_perp_2 / sqrt(kobs_perp_2(1)*kobs_perp_2(1)+kobs_perp_2(2)*kobs_perp_2(2)+kobs_perp_2(3)*kobs_perp_2(3))       
    else
       ! kobs is aligned with x -> use kobs_perp_1 = y, kobs_perp_2 = z
       kobs_perp_1 = (/0.0d0,1.0d0,0.0d0/)
       kobs_perp_2 = (/0.0d0,0.0d0,1.0d0/)
    end if

    ! define mock center (needed to rotate coords ...)
    
       
    return
  end subroutine mock_init
  

  function projected_pos(pos)
    implicit none
    real(kind=8),intent(in) :: pos(3)
    real(kind=8),intent(out) :: projected_pos(2)
    projected_pos(1) = (pos(1)-mock_center(1))*kobs_perp_1(1) + (pos(2)-mock_center(2))*kobs_perp_1(2) + (pos(3)-mock_center(3))*kobs_perp_1(3)
    projected_pos(2) = (pos(1)-mock_center(1))*kobs_perp_2(1) + (pos(2)-mock_center(2))*kobs_perp_2(2) + (pos(3)-mock_center(3))*kobs_perp_2(3)
    return
  end function projected_pos
    
  
  subroutine peel_to_map
    implicit none

    ! get exp(-tau)
    ! if exp(-tau) > ridiculous then 
    !    get ix,iy
    !    increment map

    
    return
  end subroutine peel_to_map
  
end module module_mock
