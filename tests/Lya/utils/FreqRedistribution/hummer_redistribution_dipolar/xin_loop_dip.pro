pro xin_loop_dip

  xin_arr = findgen(10)
 ; xin_arr = [8.]
  for i=0,n_elements(xin_arr)-1 do begin
     xinput = xin_arr[i]
     hummer_integ_dipolar,xinput
  endfor




end
