module module_constants
  
  public

  real(kind=8),parameter :: clight  = 2.9979250d+10       ![cm/s] light speed
  real(kind=8),parameter :: cmtoA   = 1.d8                ! from cm to A
  real(kind=8),parameter :: evtoerg = 1.6d-12             ! from eV to erg
  real(kind=8),parameter :: kb      = 1.3806200d-16       ![erg/K] Boltzman constant
  real(kind=8),parameter :: mp      = 1.66d-24            ![g] proton mass
  real(kind=8),parameter :: mpc     = 3.08d24             ![cm] from Mpc to cm
  real(kind=8),parameter :: me      = 9.11d-28            ![g] electron mass
  real(kind=8),parameter :: pi      = 3.1415926535898     ! pi
  real(kind=8),parameter :: sqrtpi  = 1.77245387556703    ! sqrt(pi)
  real(kind=8),parameter :: grtoev  = 1.78d-33*clight**2  ! from gr to eV
  real(kind=8),parameter :: e_ch    = 4.8067e-10          ![esu] electron charge
  real(kind=8),parameter :: planck  = 6.626196e-27        ![erg s] Planck's constant
  real(kind=8),parameter :: sqrt_H2Deut_mass_ratio = 0.7071067811865d0  ! == sqrt(mp/mdeut) = 1/sqrt(2).

end module module_constants
