pro lbox_from_tau0_temp

  tau0          = 1.d7
  temp          = 1.d2          ; K

  radius_sphere = 0.00100000d0
  nh            = 1.d0 ; cm^-3
  
  fix_box_size_cm = 1.d20 * tau0 * (temp/1.d4)^0.5 / (5.898d6 * nh * radius_sphere) ; cm
  vth             = 12.85*(temp / 1.d4)^0.5 * 1.d5 ; cm/s

  print,'radius_sphere = ',radius_sphere
  print,'fix_nhi = ',nh,' cm^-3'
  print,'fix_temp = ',temp,' K'
  print,'vth = ',vth,' cm/s'
  print,'fix_box_size_cm = ',fix_box_size_cm,' cm'

end
