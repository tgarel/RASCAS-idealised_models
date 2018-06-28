pro read_photons,myfile,vth,x_out,tau0,status,xlast,nb_abs,time,ids,k_ext

  openu,11,myfile
  nphotons = 0L
  r=0L
  readu,11,r,nphotons,r
  print,nphotons
  
  ids    = lonarr(nphotons)
  status = lonarr(nphotons)
  xlast  = dblarr(nphotons,3)
  nu_ext = dblarr(nphotons)
  k_ext  = dblarr(nphotons,3)
  nb_abs = lonarr(nphotons)
  time   = dblarr(nphotons)

  r=0L
  readu,11,r
  xx = -1L
  for i=0L,nphotons-1L do begin
     readu,11,xx & ids[i]=xx
  endfor  
  r=0L
  readu,11,r,r
  xx = -1L
  for i=0L,nphotons-1L do begin
     readu,11,xx & status[i]=xx
  endfor
  r=0L
  readu,11,r,r
  xx = dblarr(3) 
  for i=0L,nphotons-1L do begin
     readu,11,xx & xlast[i,*] = xx
  endfor
  r=0L
  readu,11,r,r
  xx = -1.d0
  for i=0L,nphotons-1L do begin
     readu,11,xx & nu_ext[i] = xx
  endfor

  r=0L
  readu,11,r,r
  xx = dblarr(3) 
  for i=0L,nphotons-1L do begin
     readu,11,xx & k_ext[i,*] = xx
     print,xx
  endfor

  r=0L
  readu,11,r,r
  xx = -1L
  for i=0L,nphotons-1L do begin
     readu,11,xx & nb_abs[i]=xx
  endfor
  
  r=0L
  readu,11,r,r
  xx = -1.d0
  for i=0L,nphotons-1L do begin
     readu,11,xx & time[i] = xx
  endfor

  close,11
  
  lambda_0_cm      = 1215.67 * 1.d-8
  delta_nu_doppler = vth / lambda_0_cm
  clight           = 2.9979250d10 
  nu_0             = clight / lambda_0_cm 
  x_out            = (nu_ext - nu_0) / delta_nu_doppler

end
