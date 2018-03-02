module module_voigt
  
  use module_constants, only : pi, sqrtpi
  
  private
  
  character(20) :: method = 'tasitsiomi'  ! could be 'colt', 'tasitsiomi' or 'humlicek_w4'
  
  public :: voigt_fit

contains

  function voigt_fit(x,a)

    ! returns H(x,a) = a/pi * integral(exp(-z**2) dz / (a**2+(z-x)**2))
    
    implicit none
    
    real(kind=8),intent(in) :: x,a
    real(kind=8)            :: voigt_fit 

    select case(trim(method))
    case('tasitsiomi')
       voigt_fit = tasitsiomi_fit(x,a)
    end select

    return

  end function voigt_fit

  function tasitsiomi_fit(x,a)
    ! Fit from Tasitsiomi 2006
    implicit none
    real(kind=8),intent(in) :: x,a
    real(kind=8)            :: tasitsiomi_fit 
    real(kind=8)            :: q,z,x2
    x2 = x**2
    z  = (x2 - 0.855d0) / (x2 + 3.42d0)
    if (z > 0.0d0) then 
       q = z * (1.0d0 + 21.0d0/x2) * a / pi / (x2 + 1.0d0)
       q = q * (((5.674d0*z - 9.207d0)*z + 4.421d0)*z + 0.1117)
    else
       q = 0.0d0 
    end if
    tasitsiomi_fit = sqrtpi*q + exp(-x2)
    return
   end function tasitsiomi_fit
   
 end module module_voigt
