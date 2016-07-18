module module_utils
  
  ! general-purpose functions : 
  ! 
  ! - voigt_fit
  ! - isotropic_direction
  ! - anisotropic_direction1 
  
  use module_constants, only : pi, sqrtpi, twopi
  use module_random, only : ran3
  
  public

contains

  function voigt_fit(x,a)

    ! returns ... what exactly?
    ! REF to where the fit is taken from ?
    ! comment on accuracy ?
    
    implicit none
    
    real(kind=8),intent(in) :: x,a
    real(kind=8)            :: voigt_fit 
    real(kind=8)            :: q,z,x2
    
    x2 = x**2
    z  = (x2 - 0.855d0) / (x2 + 3.42d0)
    if (z > 0) then 
       q = z * (1.0d0 + 21.0d0/x2) * a / pi / (x2 + 1.0d0)
       q = q * (((5.674d0*z - 9.207d0)*z + 4.421d0)*z + 0.1117)
    else
       q = 0.0d0 
    end if
    voigt_fit = q + exp(-x2) / sqrtpi
    
    return
    
  end function voigt_fit

  subroutine isotropic_direction(k,iran)

    ! return k vector with a random direction
    ! send back updated state of random number generator (iran) 
    
    implicit none

    real(kind=8),intent(out)      :: k(3)
    integer(kind=4),intent(inout) :: iran
    real(kind=8)                  :: cos_theta,sin_theta,phi
    
    phi   = twopi*ran3(iran)
    cos_theta = 1.0d0 - 2.0d0 * ran3(ira)  ! in [-1,1]
    sin_theta = sqrt(1.0d0 - cos_theta**2) ! in [0,1]
    k(1) = sin_theta * cos(phi)   !x
    k(2) = sin_theta * sin(phi)   !y
    k(3) = cos_theta              !z
    
  end subroutine isotropic_direction
  
end module module_utils
