pro lbox_from_taudust_temp_slab


  tau0          = 1.d6
  temp          = 1.d2          ; K

  taudust       = 0.01
  grain_radius  = 1.0d-5        ; cm
  albedo        = 0.46
  
  thickness     = 0.3  ; total thickness: use half of it to compute tau0, following Neufeld
  nh            = 1.d0 ; cm^-3

  thickness_over_two = thickness / 2.
  
  fix_box_size_cm = 1.d20 * tau0 * (temp/1.d4)^0.5 / (5.898d6 * nh * thickness_over_two) ; cm
  vth             = 12.85*(temp / 1.d4)^0.5 * 1.d5 ; cm/s

  print,'Slab thickness = ',thickness
  print,'fix_nhi = ',nh,' cm^-3'
  print,'fix_temp = ',temp,' K'
  print,'vth = ',vth,' cm/s'
  print,'fix_box_size_cm = ',fix_box_size_cm,' cm'

  ;; Wrong
  ; ndust         = taudust / (!pi * grain_radius*grain_radius) / thickness_over_two / fix_box_size_cm / (1.-albedo)

  ;; Here, taudust is TauD (Taud = TauA+TauS)
  ndust         = taudust / (!pi * grain_radius*grain_radius) / thickness_over_two / fix_box_size_cm * (1.-albedo)

  print,'fix_ndust = ',ndust,' cm^-3'

end
