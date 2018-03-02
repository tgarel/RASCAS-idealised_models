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
    case('colt')
       voigt_fit = colt_approx(x,a)
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


  function colt_approx(x,a)
    ! from Smith+15, appendix A1
    implicit none
    real(kind=8), intent(in) :: x,a
    real(kind=8)             :: A0,A1,A2,A3,A4,A5,A6
    real(kind=8)             :: B0,B1,B2,B3,B4,B5,B6,B7,B8
    real(kind=8)             :: z, colt_approx
    !--- coefficients for the rational function (Table A1, Smith+15) ---
    A0 = 15.75328153963877d0
    A1 = 286.9341762324778d0
    A2 = 19.05706700907019d0
    A3 = 28.22644017233441d0
    A4 = 9.526399802414186d0
    A5 = 35.29217026286130d0
    A6 = 0.8681020834678775d0
    B0 = 0.0003300469163682737d0
    B1 = 0.5403095364583999d0
    B2 = 2.676724102580895d0
    B3 = 12.82026082606220d0
    B4 = 3.21166435627278d0
    B5 = 32.032981933420d0
    B6 = 9.0328158696d0
    B7 = 23.7489999060d0
    B8 = 1.82106170570d0
    !------------------------------
    z = x*x
    if(z <= 3.0d0)then
       colt_approx = exp(-z) * (1.0d0-a*(A0+A1/(z-A2+A3/(z-A4+A5/(z-A6)))))
    else
       if(z < 25.0d0) then
          colt_approx = exp(-z) + a * (B0+B1/(z-B2+B3/(z+B4+B5/(z-B6+B7/(z-B8)))))
       else
          colt_approx = a/sqrtpi/(z-1.5d0-1.5d0/(z-3.5d0-5.0d0/(z-5.5d0)))
       end if
    end if
    return
  end function colt_approx

end module module_voigt
