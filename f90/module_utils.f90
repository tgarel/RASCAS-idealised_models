module module_utils
  
  ! general-purpose functions

  use module_constants, only : pi, sqrtpi
  
  public

contains

  function voigt_fit(x,a)

    ! returns ... 
    
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
  
end module module_utils
