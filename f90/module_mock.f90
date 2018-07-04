module module_mock

  public

  ! direction of observation
  real(kind=8) :: kobs(3),kobs_perp_1(3),kobs_perp_2(3)

  ! center of observation
  real(kind=8) :: mock_center(3)

  ! extent of observation
  real(kind=8) :: mock_image_side ! box units
  
  ! detector (image for now)
  integer(kind=4),parameter :: npix = 1000
  real(kind=8)    :: image(npix,npix)


contains

  subroutine mock_init()

    implicit none
        
    image = 0.0d0
    ! define direction of observation (normalise vector)
    kobs = (/0.31,0.2,0.42/)
    kobs = kobs / sqrt(kobs(1)*kobs(1)+kobs(2)*kobs(2)+kobs(3)*kobs(3))
    ! define basis for sky plane
    if (kobs(1) < 1.) then
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
  
  
  subroutine projected_pos(pos,pp)
    implicit none
    real(kind=8),intent(in) :: pos(3)
    real(kind=8),intent(inout) :: pp(2)
    
    pp(1) = (pos(1)-mock_center(1))*kobs_perp_1(1) + (pos(2)-mock_center(2))*kobs_perp_1(2) + (pos(3)-mock_center(3))*kobs_perp_1(3)
    pp(2) = (pos(1)-mock_center(1))*kobs_perp_2(1) + (pos(2)-mock_center(2))*kobs_perp_2(2) + (pos(3)-mock_center(3))*kobs_perp_2(3)

    return
  end subroutine projected_pos
    
  
  subroutine peel_to_map(peel_pos,peel_contribution)
    implicit none

    real(kind=8),intent(in) :: peel_pos(3),peel_contribution
    real(kind=8)    :: pp(2) 
    integer(kind=4) :: ix,iy

    call projected_pos(peel_pos,pp)
    ix = int((pp(1) + 0.5d0 * mock_image_side) /mock_image_side)
    iy = int((pp(2) + 0.5d0 * mock_image_side) /mock_image_side)
    if (ix>0 .and. ix<=npix .and. iy>0 .and. iy<=npix) image(ix,iy) = image(ix,iy) + peel_contribution
    
    return
  end subroutine peel_to_map
  
end module module_mock
